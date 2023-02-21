## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/BBS_iCAR_route_trends")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
source("functions/neighbours_define.R") ## function to define neighbourhood relationships
source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
## above source() call over-writes the bbsBayes prepare_data() function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output



output_dir <- "output"


if(!dir.exists("Figures")){
  dir.create("Figures")
}
if(!dir.exists("data/maps/")){
  if(!dir.exists("data/")){
    dir.create("data/")
  }
  dir.create("data/maps/")
}
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}
if(!dir.exists("trends")){
  dir.create("trends")
}


# load and stratify CASW data ---------------------------------------------
#species = "Pacific Wren"
strat = "bbs_usgs"
model = "slope"


firstYear = 1995
lastYear = 2021

## this list should include all of the species that we're interested in for the grasslands project
species_list = c("Chestnut-collared Longspur",
                 "Thick-billed Longspur",
                 "Eastern Meadowlark",
                 "Western Meadowlark",
                 "Savannah Sparrow")


# these two can be set to FALSE if it makes sense just to fit the models and store the model output from Stan
# e.g., if it's easier to run the models and then to produce maps and trends after the fact
produce_trends <- TRUE
produce_maps <- TRUE # this is only relevant if produce_trends == TRUE


strat_data <- bbsBayes::stratify(by = strat)

# I've got this running as a simple species loop
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp <- "_BYM_"

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)


# SPECIES LOOP ------------------------------------------------------------




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")



### this is the alternate prepare data function
jags_data = prepare_data(strat_data = strat_data,
                              species_to_run = species,
                              model = model,
                              #n_knots = 10,
                              min_year = firstYear,
                              max_year = lastYear,
                              min_n_routes = 1)# spatial neighbourhood define --------------------------------------------

# strata map of one of the bbsBayes base maps
# helps group and set boundaries for the route-level neighbours,
## NOT directly used in the model
strata_map  <- get_basemap(strata_type = strat,
                           transform_laea = TRUE,
                           append_area_weights = FALSE)

strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))

#drops the strata that are not included for this species
# and adds the strat indicator variable to link to model output
realized_strata_map <- filter(strata_map,ST_12 %in% unique(jags_data$strat_name)) %>% 
  inner_join(.,strata_list, by = "ST_12")

# Spatial boundaries set up --------------------

jags_data[["routeF"]] <- as.integer(factor((jags_data$route)))

#create a data frame of each unique route in the species-specific dataset
route_map = unique(data.frame(route = jags_data$route,
                              routeF = jags_data$routeF,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
# adhoc way of separating different routes with the same starting coordinates
# this shifts the starting coordinates of teh duplicates by ~1.5km to the North East 
# ensures that the duplicates have a unique spatial location, but remain very close to
# their original location and retain reasonable neighbourhood relationships
# these duplicates happen when a "new" route is established because some large proportion
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

#reconcile the projections of routes and base bbs strata
route_map = st_transform(route_map,crs = st_crs(realized_strata_map))


## custom function that returns the adjacency data necessary for the stan model
## also exports maps and saved data objects to plot_dir

car_stan_dat <- neighbours_define(real_strata_map = route_map,
                                  #strat_link_fill = 100000,
                                  plot_neighbours = TRUE,
                                  species = species,
                                  plot_dir = "data/maps/",
                                  plot_file = paste0("_",firstYear,"_",lastYear,"_route_maps"),
                                  save_plot_data = TRUE,
                                  voronoi = TRUE,
                                  strat_indicator = "routeF",
                                  add_map = realized_strata_map)






stan_data = jags_data[c("ncounts",
                        #"nstrata",
                        #"nobservers",
                        "count",
                        #"strat",
                        #"obser",
                        "year",
                        "firstyr",
                        "fixedyear")]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["observer"]] <- as.integer(factor((jags_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2
stan_data[["route"]] = jags_data$routeF
stan_data[["nroutes"]] = max(jags_data$routeF)


if(car_stan_dat$N != stan_data[["nroutes"]]){stop("Some routes are missing from adjacency matrix")}

save(list = c("stan_data",
              "jags_data",
              "route_map",
              "realized_strata_map",
              "firstYear"),
     file = sp_data_file)




}



