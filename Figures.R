#Figures

library(tidyverse)
library(patchwork)


# Figure_1 ----------------------------------------------------------------

# demonstration of neighbours for BHVI

# load saved neighbours data from neighbours_define_voronoi(..., save_plot_data = TRUE)

load("data/Lazuli_Bunting_route_maps_data.RData")

box <- st_as_sfc(st_bbox(strata_map))

xb <- range(st_coordinates(box)[,"X"])
yb <- range(st_coordinates(box)[,"Y"])


ggp1 <- ggplot()+ 
  geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+ 
  #geom_sf(data = cov_hull_clip,alpha = 0,colour = grey(0.8))+
  #geom_sf(data = v,alpha = 0,colour = grey(0.9))+
  geom_sf(data = full_clip,alpha = 0,colour = grey(0.1))+
  geom_sf(data = vintj,alpha = 0,colour = grey(0.5))+
  geom_sf(data = centres,alpha = 1,size = 0.6)+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none")
ggp1

ggp2 <- ggplot()+ 
  geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+ 
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),
               inherit.aes = FALSE,linewidth=0.3,colour = grey(0.6)) +
  geom_sf(data = centres,alpha = 1,size = 1.3)+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none")
ggp2


pdf("Figures/Figure_1.pdf",
    width = 7.5,
    height = 10)
print(ggp1 + ggp2 + plot_layout(nrow = 2))
dev.off()






# Figure 2 - 4 model comparison -------------------------------------------

firstYear <- 2006
lastYear <- 2021

ppy <- function(x){
  p <- exp(x)-1
  return(p*100)
}

species <- "Lazuli Bunting" 
  

species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  out_base_temp <- paste0(species_f,"_GP_",firstYear,"_",lastYear)
  

  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
  
  
  load(sp_data_file)
  
  
  outboth <- NULL
  for(spp1 in c("BYM","iCAR","nonspatial","GP")){
    
    spp <- paste0("_",spp1,"_")
    
    
    out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
    
    
    
    
    
    #stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
    
    summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))
    
    abundance <- summ %>% 
      filter(grepl(pattern = "alpha[",
                   variable,
                   fixed = TRUE),
             !grepl(pattern = "gp",
                    variable,
                    fixed = TRUE)) %>% 
      mutate(routeF = row_number(),
             across(mean:q95,exp)) %>% 
      rename_with(.,
                  ~ paste0("abundance_",.x),
                  .cols = mean:ess_tail)
    
    
    slope <- summ %>% 
      filter(grepl(pattern = "beta[",
                   variable,
                   fixed = TRUE),
             !grepl(pattern = "gp",
                    variable,
                    fixed = TRUE)) %>% 
      mutate(routeF = row_number(),
             across(mean:q95, ppy )) %>% 
      rename_with(.,
                  ~ paste0("trend_",.x),
                  .cols = mean:ess_tail) %>% 
      select(-variable)
    
    both <- inner_join(slope,abundance,
                       by = "routeF") %>% 
      mutate(model = spp1)
    
    outboth <- bind_rows(outboth,both) %>% 
      mutate(species = species)
    
    
    
  } #end models loop
  
  strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
  bb = st_bbox(strata_bounds)
  xlms = as.numeric(c(bb$xmin,bb$xmax))
  ylms = as.numeric(c(bb$ymin,bb$ymax))
  
  
  plot_map <- route_map %>% 
    left_join(.,outboth,
              by = "routeF") %>% 
    mutate(model = factor(model,
                          levels = c("iCAR","GP","BYM","nonspatial"),
                          ordered = TRUE))
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  plot_map$Tplot <- cut(plot_map$trend_mean,breaks = c(-Inf, breaks, Inf),labels = labls)
  
  
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  
  
  map <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = Tplot,
                size = abundance_mean))+
    scale_size_continuous(range = c(0.05,2),
                          name = "Mean Count")+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0(lgnd_head,firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    labs(title = paste(species, "trends"))+
    facet_wrap(vars(model))
  
  pdf(paste0("Figures/Figure_2.pdf"),
      height = 7,
      width = 7.5)
  print(map)
  dev.off()
  
  map_se <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = trend_sd,
                size = abundance_sd))+
    scale_size_continuous(range = c(0.05,2),
                          name = "SE of Mean Count",
                          trans = "reverse")+
    scale_colour_viridis_c(aesthetics = c("colour"),
                           guide = guide_legend(reverse=TRUE),
                           name = paste0("SE of Trend",firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    labs(title = paste(species, "Standard error"))+
    facet_wrap(vars(model))
  
 # version with teh SE of hte maps
  pdf(paste0("Figures/Figure_S2.pdf"),
      height = 10.5,
      width = 7.5)
  print(map / map_se)
  dev.off()
  
