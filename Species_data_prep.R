## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
source("functions/neighbours_define_voronoi.R") ## function to define neighbourhood relationships
source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
## above source() call over-writes the bbsBayes prepare_data() function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output



output_dir <- "output"

# 
# if(!dir.exists("Figures")){
#   dir.create("Figures")
# }
# if(!dir.exists("data/maps/")){
#   if(!dir.exists("data/")){
#     dir.create("data/")
#   }
#   dir.create("data/maps/")
# }
# if(!dir.exists(output_dir)){
#   dir.create(output_dir)
# }
# if(!dir.exists("trends")){
#   dir.create("trends")
# }


# load and stratify CASW data ---------------------------------------------
strat = "bbs_usgs"
model = "slope"



## this list should include all of the species that we're interested in for the grasslands project
species_list = c("Dickcissel",
                 "Chestnut-collared Longspur",
                 "Thick-billed Longspur",
                 "Eastern Meadowlark",
                 "Western Meadowlark",
                 "Savannah Sparrow")




strat_data <- bbsBayes::stratify(by = strat)



spp <- "_base_"

spans <- data.frame(ly = c(2021), #last year of the time-span
                    fy = c(2001)) # first year of the time-span

# SPECIES LOOP ------------------------------------------------------------

# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
#species <- species_list[2]


  
  
  
  
  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)



#for(ii in 1:nrow(spans)){
ii <- 1
firstYear <- spans[ii,"fy"]
  lastYear <- spans[ii,"ly"]
  
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")



### this is the alternate prepare data function - modified from bbsBayes
jags_data = prepare_data(strat_data = strat_data,
                              species_to_run = species,
                              model = model,
                              #n_knots = 10,
                              min_year = firstYear,
                              max_year = lastYear,
                              min_n_routes = 1,
                         min_max_route_years = 1)# spatial neighbourhood define --------------------------------------------

# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours,
## NOT directly used in the model
strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)


#create list of routes and locations to ID routes that are not inside of original strata (some off-shore islands)
route_map1 = unique(data.frame(route = jags_data$route,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))

#create spatial object of above
route_map1 = st_as_sf(route_map1,coords = c("Longitude","Latitude"))
st_crs(route_map1) <- 4269 #NAD83 commonly used by US federal agencies

#reconcile the projections of routes and base bbs strata
route_map1 = st_transform(route_map1,crs = st_crs(strata_map))

#drops the routes geographically outside of the strata (some offshore islands) 
# and adds the strat indicator variable to link to model output
strata_map_buf <- strata_map %>% 
  filter(ST_12 %in% route_map1$strat) %>% 
  summarise() %>% 
  st_buffer(.,10000) #drops any routes with start-points > 10 km outside of strata boundaries
realized_routes <- route_map1 %>% 
  st_join(.,strata_map_buf,
          join = st_within,
          left = FALSE) 



# reorganizes data after routes were dropped outside of strata
new_data <- data.frame(strat_name = jags_data$strat_name,
                       strat = jags_data$strat,
                       route = jags_data$route,
                       strat = jags_data$strat_name,
                       Latitude = jags_data$Latitude,
                       Longitude = jags_data$Longitude,
                       count = jags_data$count,
                       year = jags_data$year,
                       firstyr = jags_data$firstyr,
                       ObsN = jags_data$ObsN,
                       r_year = jags_data$r_year) %>% 
  filter(route %in% realized_routes$route)




strata_list <- data.frame(ST_12 = unique(new_data$strat_name),
                          strat = unique(new_data$strat))


realized_strata_map <- strata_map %>%
  filter(ST_12 %in% strata_list$ST_12)

# Spatial neighbours set up --------------------

new_data$routeF <- as.integer(factor((new_data$route))) #main route-level integer index

#create a data frame of each unique route in the species-specific dataset
route_map = unique(data.frame(route = new_data$route,
                              routeF = new_data$routeF,
                              strat = new_data$strat_name,
                              Latitude = new_data$Latitude,
                              Longitude = new_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
# these duplicates happen when a "new" route is established (i.e., a route is re-named) 
# because some large proportion
# of the end of a route is changed, but the start-point remains the same
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
while(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
  
}
dups = which(duplicated(route_map[,c("Latitude","Longitude")])) 
if(length(dups) > 0){stop(paste(spec,"At least one duplicate route remains"))}



#create spatial object from route_map dataframe
route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies


#reproject teh routes spatial object ot match the strata-map in equal area projection
route_map = st_transform(route_map,crs = st_crs(strata_map))

## custom function that returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir

## this function generates a Voronoi tesselation with each route start-point
## at teh centre of the voronoi polygons.
## it then clips the tesselated surface using an intersection of the strata map 
## and a concave polygon surrounding the route start locations
## this clipping ensures that routes at the edges of the species range are not able
## to be considered neighbours with other very-distant edge-routes and retains
## the complex shape of many species distributions (e.g., species with boreal and 
## mountainous distributions, such as Pileated Woodpecker)
## it also forces a completely connected network of routes. So if some portions of the 
## species range are separated by gaps (e.g., intervening BBS strata without any routes)
## it forces a neighbourhood connection between the closest pair of routes that would
## bridge the gap (and any other routes within 10% of the distance of the closest pair)
car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                  species = species,
                                  strat_indicator = "routeF",
                                  strata_map = realized_strata_map)


#print(car_stan_dat$map)



# set up cross-validation folds -------------------------------------------

# new_data$folds <- cv_folds(new_data,
#                   fold_groups = "routeF")
# 
# check_folds <- new_data %>% 
#   group_by(routeF,ObsN,folds) %>% 
#   summarise(n = n())


stan_data <- list()
stan_data[["count"]] <- new_data$count
stan_data[["ncounts"]] <- length(new_data$count)
stan_data[["strat"]] <- new_data$strat
stan_data[["route"]] <- new_data$routeF
stan_data[["year"]] <- new_data$year
stan_data[["firstyr"]] <- new_data$firstyr
stan_data[["fixedyear"]] <- jags_data$fixedyear


stan_data[["nyears"]] <- max(new_data$year)
stan_data[["observer"]] <- as.integer(factor((new_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] <- car_stan_dat$N_edges
stan_data[["node1"]] <- car_stan_dat$node1
stan_data[["node2"]] <- car_stan_dat$node2
stan_data[["nroutes"]] <- max(stan_data$route)


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

dist_matrix_km <- dist_matrix(route_map,
                              strat_indicator = "routeF")
save(list = c("stan_data",
              "new_data",
              "route_map",
              "realized_strata_map",
              "firstYear",
              "car_stan_dat",
              "dist_matrix_km"),
     file = sp_data_file)


} #end species loop




