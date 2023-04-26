# voronoi polygon neighbourhood function ----------------------------------


neighbours_define_voronoi <- function(real_point_map = route_map, #sf map of strata
                              species = "",
                              plot_dir = "",
                              plot_file = "_route_maps",
                              save_plot_data = FALSE,
                              strata_map = NULL,
                              strat_indicator = "routeF",
                              lab_size = 3,
                              concavity = 1){
  
  require(spdep)
  require(sf)
  require(tidyverse)
  require(concaveman)
  
  if(is.null(strata_map)){stop("Strata map is missing and is required to constrain neighbours \n at edges of species distribution")}

  
    # function to prep spatial data for stan model ----------------------------
  mungeCARdata4stan = function(adjBUGS,numBUGS) {
    N = length(numBUGS);
    nn = numBUGS;
    N_edges = length(adjBUGS) / 2;
    node1 = vector(mode="numeric", length=N_edges);
    node2 = vector(mode="numeric", length=N_edges);
    iAdj = 0;
    iEdge = 0;
    for (i in 1:N) {
      for (j in 1:nn[i]) {
        iAdj = iAdj + 1;
        if (i < adjBUGS[iAdj]) {
          iEdge = iEdge + 1;
          node1[iEdge] = i;
          node2[iEdge] = adjBUGS[iAdj];
        }
      }
    }
    return (list("N"=N,"N_edges"=N_edges,"node1"=node1,"node2"=node2));
  }
  
  # neighbourhood define ----------------------------------------------------
  
  real_point_map <- real_point_map %>% rename_with(.,
                                                     ~ gsub(pattern = strat_indicator, 
                                                            replacement = "strat_lab",
                                                            .x, fixed = TRUE)) %>% 
    group_by(strat_lab) %>% 
    summarise() %>% 
    arrange(strat_lab)
  
  
  centres <- suppressWarnings(st_centroid(real_point_map))
  
  coords = st_coordinates(centres)
  
  check_concave <- TRUE # this and the associated while loop are required
  # because in some rare situations the concaveman function seems to fail
  # it will exclude one point near the edge. Will submit an issue to the 
  # package gitHub with reprex, hopefully will get solved.
  while(check_concave){
  #concave hull of route locations for clipping later in the function
  cov_hull_fill <- concaveman::concaveman(centres,
                                          concavity = concavity) %>% 
    st_buffer(.,50000) %>% #buffer by 50km to ensure all of route is included 
    mutate(concave = TRUE) 

  miss <- real_point_map %>%
    st_join(.,cov_hull_fill,
            join = st_intersects,
            left_join = TRUE) %>% 
    filter(is.na(concave))
check_concave <- ifelse(nrow(miss) > 0,TRUE,FALSE)
concavity <- concavity+0.5
  }
  
  # tmpp <- ggplot(data = cov_hull_fill)+
  #   geom_sf()+
  #   geom_sf(data = real_point_map)+
  #   geom_sf_text(data = real_point_map,
  #                aes(label = strat_lab))
  # tmpp
  
  
    # Voronoi polygons from centres -----------------------------------
    box <- st_as_sfc(st_bbox(centres))

    xb <- range(st_coordinates(box)[,"X"])
    yb <- range(st_coordinates(box)[,"Y"])
    
    v <- st_cast(st_voronoi(st_union(centres), envelope = box))
    v <- st_join(st_sf(v),centres,join = st_contains) # join to label the polygons with route-indices
    
    # clip voronoi polygon with intersection of sstrata map and concave hull
    # simplify strata map (strata_map)
      add_fill <- strata_map %>%
        st_buffer(.,10000) %>% #add 10 km buffer to ensure all routes are covered and simplify boundary
        ungroup() %>% 
        summarise() %>% 
        st_cast(.,to = "POLYGON")
     
      # intersect strata map and concave hull of route locations 
      add_fill2 <- st_intersection(add_fill,cov_hull_fill) %>% 
        ungroup() %>% 
        summarise() %>% 
        st_cast(.,to = "POLYGON")
      
    #clip the voronoi polygons with the strata+concavehull intersection
      vintj <- suppressWarnings(st_intersection(v,add_fill2)) %>% 
        group_by(strat_lab) %>% 
        summarise() %>% 
        arrange(.,strat_lab)

    #create neighbours with spdep::poly2nb
    nb_db <- spdep::poly2nb(vintj,row.names = vintj$strat_lab,queen = FALSE)#polygon to neighbour definition
    #create neighbour matrix
    nb_mat <- spdep::nb2mat(nb_db, style = "B",
                           zero.policy = TRUE) #binary adjacency matrix
    nb_info = spdep::nb2WB(nb_db) #summaries of neighbour info to check all points are connected
    
    # if some points have no neighbours, then fill based on nearest neighbours according to their centroids
    if(min(nb_info$num) == 0){
      nn_fill <- TRUE
      message("Some strata have no neighbours, filling by 2 nearest neighbours by centroids")
      
      
      nn = knearneigh(centres, k=2)
      
      w_rep = which(nb_info$num == 0)
      
      for(i in w_rep){
        wm <- nn[[1]][i,c(1,2)]
        
        for(jjt in c(1,2)){
          wwm <- wm[jjt]
          
          nb_db[[i]] <- as.integer(unique(c(nb_db[[i]],wwm)))
          if(nb_db[[i]][1] == 0){nb_db[[i]] <- nb_db[[i]][-1]}
          nb_db[[wwm]] <- as.integer(unique(c(nb_db[[wwm]],i)))
          if(nb_db[[wwm]][1] == 0){nb_db[[wwm]] <- nb_db[[wwm]][-1]}
        }
        
      }
    }
    
    
    #identify the number of "islands" in the neighbour graph (i.e., fully connected graph has 1 island)
    n_islands <- n.comp.nb(nb_db)$nc
    #distance matrix among points
    distnc <- st_distance(centres)
    
    #if the graph is not fully conencted, then force a full connection between the islands of points
    while(n_islands > 1){
      message(paste(n_islands-1,"groups of nodes are isolated, linking by distance between centroids"))

      isls <- n.comp.nb(nb_db)


      ww1 <- which(isls$comp.id == 1)
      tmp <- distnc[ww1,-c(ww1)]

      clstn <- apply(tmp,1,min)
      clst <- as.numeric(clstn) #minimum values in each row (for each of the sites in ww1)
      wwcl <- (names(which.min(clstn))) #which row includes the minumum values (which site is closest)
      ww2 <- which(as.numeric(distnc[wwcl,]) == clstn[wwcl] |
                     (as.numeric(distnc[wwcl,]) > clstn[wwcl] &
                        as.numeric(distnc[wwcl,]) < clstn[wwcl]*1.1))

      if(any(ww1 %in% ww2)){ww2 <- ww2[-which(ww2 %in% ww1)]}
      #ww2 are the strata that should be linked to the isolated group
      for(i in ww2){

        nb_db[[i]] <- unique(c(nb_db[[i]],as.integer(wwcl)))
        nb_db[[as.integer(wwcl)]] <- unique(c(nb_db[[as.integer(wwcl)]],i))


      }

      n_islands <- n.comp.nb(nb_db)$nc #check that n_islands now == 1

    }


    ## final adjacency matrix and neighbours
    nb_mat <- spdep::nb2mat(nb_db, style = "B",
                            zero.policy = TRUE) #binary adjacency matrix
    
    
    nb_info = spdep::nb2WB(nb_db)
    
    
    # this shouldn't happen, but here as a check
    if(min(nb_info$num) == 0){stop("ERROR some strata have no neighbours")}
    
    
    
    
    #plot the neighbourhood relationships 
 
  
    species_dirname <- gsub(pattern = " ",
                            replacement = "_",
                            x = species)
    
    
    plot_file_name = paste0(plot_dir,species_dirname,plot_file,".pdf")
    
    
    
    
    nb_l <- nb2listw(nb_db)
    nt = length(attributes(nb_l$neighbours)$region.id)
    DA = data.frame(
      from = rep(1:nt,sapply(nb_l$neighbours,length)),
      to = unlist(nb_l$neighbours)
    )
    DA = cbind(DA,coords[DA$from,c("X","Y")],coords[DA$to,c("X","Y")])
    colnames(DA)[3:6] = c("long","lat","long_to","lat_to")
    
    ggp <- ggplot(data = centres)+ 
      geom_sf(aes(col = strat_lab,alpha = 0.5)) 

    
    if(!is.null(strata_map)){
      
      ggp <- ggp +
        geom_sf(data = strata_map,alpha = 0,colour = grey(0.9))
    }
  
    ggp <- ggp + 
      geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),
                   inherit.aes = FALSE,linewidth=0.3,alpha=0.4) +
      geom_sf(data = vintj,alpha = 0,colour = grey(0.95))+ 
      geom_sf(data = real_point_map,alpha = 0,colour = grey(0.85))+
      geom_sf_text(aes(label = strat_lab,
                       colour = strat_lab),size = lab_size,alpha = 0.7)+
      labs(title = species)+
        theme_minimal() +
        coord_sf(xlim = xb,ylim = yb)+
        theme(legend.position = "none")
    
  
    
  
    if(save_plot_data){
      #print the map
      pdf(file = plot_file_name,
          width = 11,
          height = 8.5)
      print(ggp)
      dev.off()
      
      #save the data in a .RData file
      save_file_name = paste0(plot_dir,species_dirname,plot_file,"_data.RData")
      
      save(list = c("centres",
                    "real_point_map",
                    "vintj",
                    "nb_db",
                    "coords",
                    "nb_mat",
                    "DA",
                    "ggp",
                    "distnc"),
           file = save_file_name)
    }
  
  
  
  ### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
  car_stan <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                numBUGS = nb_info$num)
  

  return(list(N = car_stan$N,
              N_edges = car_stan$N_edges,
              node1 = car_stan$node1,
              node2 = car_stan$node2,
              adj_matrix = nb_mat,
              map = ggp,
              distance_matrix = distnc))
} ### end of function




# Distance_matrix ---------------------------------------------------------

dist_matrix <- function(
    points_sf = route_starts, #simple feature points
    strat_indicator = "route"
    ){
  require(spdep)
  require(sf)
  require(tidyverse)
  
  real_map <- points_sf %>% rename_with(.,
                                        ~ gsub(pattern = strat_indicator,
                                               replacement = "site_lab",
                                               .x, fixed = TRUE)) %>% 
    group_by(site_lab) %>% 
    summarise() %>% 
    arrange(site_lab)
  
  #export distance matrix in 1000km
  mat <- as.matrix(st_distance(real_map))/1000000
  attr(mat,"units") <- "1000km"
  return(mat)
  
}


