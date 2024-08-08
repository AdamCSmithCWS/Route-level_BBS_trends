

#setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")





library(bbsBayes2)#updated version of bbsBayes package
library(bbsBayes) # archived earlier version of bbsBayes
library(tidyverse)
library(sf)
library(cmdstanr)
library(patchwork)

output_dir <- "output"
# species <- "Horned Grebe" 
# species <- "Red-winged Blackbird" 
source("functions/neighbours_define_voronoi.R") ## function to define neighbourhood relationships
source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
source("functions/GAM_basis_function_mgcv.R")
strat = "bbs_usgs"
model = "slope"


crs_use <- readRDS("functions/custom_crs_for_maps.rds")

base_strata_map <- bbsBayes2::load_map("bbs_usgs")%>% 
  st_transform(.,crs_use)

state_prov <- bbsBayes2::load_map("prov_state") %>% 
  st_transform(.,crs_use)



strat_data <- bbsBayes::stratify(by = strat)


# species_list <- c("American Coot",
#                   "Pied-billed Grebe",
#                   "Northern Shoveler",
#                   "Black Tern",
#                   "Eared Grebe",
#                   "Horned Grebe",
#                   "Sora")

species_list <- c("Horned Grebe")



for(species in species_list){


species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp <- "_moisture_"

exp_t <- function(x){
  y <- (exp(x)-1)*100
}

firstYear <- 1975
lastYear <- 2017


out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

sp_data_file <- paste0("data_open/",species_f,"_",firstYear,"_",lastYear,"_covariate_stan_data.RData")



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
                           append_area_weights = FALSE) %>% 
  filter(COUNTRY == "CA")


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


# Link covariates and counts ----------------------------------------------

# Load covariate data -------------------------------------------------------

clim <- read.csv("data_open/HOGR_climate_covariates.csv")

clim2 <- read.csv("data_open/HOGR_ponds_buffer.csv")
### include the ponds variable
clim <- clim %>% 
  inner_join(.,clim2,by = c("Route" = "ProvRoute",
                            "Year")) %>% 
  filter(Year >= firstYear,
         Year <= lastYear) %>% 
  mutate(route =  paste(trunc(Route/1000),as.integer(str_sub(Route,nchar(Route)-2)),
                        sep = "-")) %>%
  mutate(prec_scale = prec/sd(prec),
         temp_scale = maxt/sd(maxt),
         pond_scale = log(Ponds+1)) %>% 
  group_by(route) %>% 
  mutate(prec_c = prec_scale-mean(prec_scale),
         temp_c = temp_scale-mean(temp_scale),
         pond_c = pond_scale-mean(pond_scale,na.rm = TRUE)) %>% 
  rename(r_year = Year) %>% 
  select(-Route) %>% 
  drop_na()

test1 <- clim %>% 
  group_by(route) %>% 
  summarise(m_prec = mean(prec_c),
            m_temp = mean(temp_c),
            m_pond = mean(pond_c),
            sd_pred = sd(prec_c),
            sd_temp = sd(temp_c),
            sd_pond = sd(pond_c)) 

## drop all BBS data with no climate info
new_data <- new_data %>% 
  inner_join(.,clim,
             by = c("route","r_year"))

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
## mountainous distributions, or forest species that span either
## side of the great plains such as Pileated Woodpecker)
## it also forces a completely connected network of routes. So if some portions of the 
## species range are separated by gaps (e.g., intervening BBS strata without any routes)
## it forces a neighbourhood connection between the closest pair of routes that would
## bridge the gap (and any other routes within 10% of the distance of the closest pair)
car_stan_dat <- neighbours_define_voronoi(real_point_map = route_map,
                                          species = species,
                                          strat_indicator = "routeF",
                                          strata_map = realized_strata_map,
                                          concavity = 1,
                                          save_plot_data = TRUE,
                                          plot_dir = "data")#concavity argument from concaveman()

## a relative measure of concavity. 1 results in a relatively detailed shape, Infinity results in a convex hull.


pdf(paste0("data/maps/route_map_",firstYear,"-",lastYear,"_",species_f,".pdf"))
print(car_stan_dat$map)
dev.off()


# # generate GAM predictor basis
# gam_prec <- gam_basis(new_data$prec_c,
#                       sm_name = "prec")
# gam_temp <- gam_basis(new_data$temp_c,
#                       sm_name = "temp")


stan_data <- list()
stan_data[["count"]] <- new_data$count
stan_data[["ncounts"]] <- length(new_data$count)
stan_data[["strat"]] <- new_data$strat
stan_data[["route"]] <- new_data$routeF
stan_data[["year"]] <- new_data$year
stan_data[["firstyr"]] <- new_data$firstyr
stan_data[["pond"]] <- new_data$pond_c
# stan_data[["temp"]] <- new_data$temp_c
stan_data[["fixedyear"]] <- jags_data$fixedyear


stan_data[["nyears"]] <- max(new_data$year)
stan_data[["observer"]] <- as.integer(factor((new_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)



stan_data[["N_edges"]] <- car_stan_dat$N_edges
stan_data[["node1"]] <- car_stan_dat$node1
stan_data[["node2"]] <- car_stan_dat$node2
stan_data[["nroutes"]] <- max(stan_data$route)
stan_data[["sd_alpha_prior"]] <- 2

# stan_data[["temp_pred_basis"]] <- gam_temp$temp_basis
# stan_data[["nknots_temp"]] <- gam_temp$nknots_temp
# stan_data[["temp_vis_basis"]] <- gam_temp$temp_basispred


# stan_data[["prec_pred_basis"]] <- gam_prec$prec_basis
# stan_data[["nknots_prec"]] <- gam_prec$nknots_prec
# stan_data[["prec_vis_basis"]] <- gam_prec$prec_basispred


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


mod.file = paste0("models/slope",spp,"route_NB.stan")

slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

re_fit <- FALSE # set to TRUE if rerunning Stan model

if(re_fit){
stanfit <- slope_model$sample(
  data=stan_data,
  refresh=400,
  iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 4)

summ <- stanfit$summary()
print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))

}else{
  stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
  summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))
}


summ %>% arrange(-rhat)

RHO <- summ %>% 
  filter(variable == "RHO")
BETA <- summ %>% 
  filter(variable == "BETA") %>% 
  mutate(across(2:7,~exp_t(.x)))


load(sp_data_file)

mod.file = paste0("models/slope_iCAR_route_NB.stan")

slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stan_data[["pond"]] <- NULL

if(refit){
stanfit_simple <- slope_model$sample(
  data=stan_data,
  refresh=400,
  iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 4)

summ_simple <- stanfit_simple$summary()
print(paste(species, stanfit_simple$time()[["total"]]))

saveRDS(stanfit_simple,
        paste0(output_dir,"/",out_base,"_simple_stanfit.rds"))

saveRDS(summ_simple,
        paste0(output_dir,"/",out_base,"_simple_summ_fit.rds"))
}else{
  stanfit_simple <- readRDS(paste0(output_dir,"/",out_base,"_simple_stanfit.rds"))
  summ_simple <- readRDS(paste0(output_dir,"/",out_base,"_simple_summ_fit.rds"))
}



BETA_simple <- summ_simple %>% 
  filter(variable == "BETA") %>% 
  mutate(across(2:7,~exp_t(.x)))





pond_eff <- summ %>% 
  filter(variable == "RHO") %>% 
  select(mean,q5,q95) %>% 
  unlist() %>% 
  signif(.,2)
pond_eff <- paste0(pond_eff["mean"]," ","[",pond_eff["q5"],
                     " : ",pond_eff["q95"],"]")



BETA_eff <- paste0(signif(BETA["mean"],2)," %/year (moisture removed) : ",signif(BETA_simple["mean"],2)," %/year")



exp_t <- function(x){
  y <- (exp(x)-1)*100
}


# plot trends -------------------------------------------------------------



strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))

betas1 <- summ %>% 
  filter(grepl("beta[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "After removing pond effect") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

rho1 <- summ %>% 
  filter(grepl("rho[",variable,fixed = TRUE)) %>% 
  mutate(routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Pond effect") %>% 
  select(routeF,median,sd,q5) %>% 
  rename(ponds = median,
         ponds_se = sd,
         ponds_lci = q5)
betas1 <- betas1 %>% 
  inner_join(.,rho1,
             by = "routeF")

betas2 <- summ_simple %>% 
  filter(grepl("beta[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Full with pond effect") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

betas <- bind_rows(betas1,betas2)

plot_map <- route_map %>% 
  left_join(.,betas,
            by = "routeF",
            multiple = "all") 

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
plot_map$Tplot <- cut(plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)


map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                 "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
names(map_palette) <- labls



map <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map,
          aes(colour = Tplot,
              size = trend_se))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  guides(size = "none")+
  labs(title = paste(species,firstYear,"-",lastYear,BETA_eff))+
  facet_wrap(vars(parameter))

plot_map2 <- plot_map %>% 
  filter(!is.na(ponds))
map_ponds <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map2,
          aes(colour = ponds))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("Pond effect"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  labs(title = paste(firstYear,"-",lastYear))
# 
map_ponds_lc <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map2,
          aes(colour = ponds_lci))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("LCL of Pond effect"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  labs(title = paste(firstYear,"-",lastYear))

# 
# 
# map_se <- ggplot()+
#   geom_sf(data = base_strata_map,
#           fill = NA,
#           colour = grey(0.75))+
# geom_sf(data = state_prov,
#         fill = NA,
#         colour = grey(0.5))+
#   geom_sf(data = plot_map,
#           aes(colour = trend_se,
#               size = abundance_se))+
#   scale_size_continuous(range = c(0.05,2),
#                         name = "SE of Mean Count",
#                         trans = "reverse")+
#   scale_colour_viridis_c(aesthetics = c("colour"),
#                          guide = guide_legend(reverse=TRUE),
#                          name = paste0("SE of Trend"))+
#   coord_sf(xlim = xlms,ylim = ylms)+
#   theme_bw()+
#   guides(size = "none")+
#   labs(title = paste(firstYear,"-",lastYear))+
#   facet_wrap(vars(parameter))
# 
# 



map_save <- map / map_ponds + plot_layout(guides = "collect")
saveRDS(map_save,paste0("Figures/saved_map_",species_f,"_covariate_",firstYear,".rds"))

pdf(paste0("Figures/Figure_S14_",species_f,"_",firstYear,".pdf"),
    height = 8.5,
    width = 8.5)


print(map/ (map_ponds + map_ponds_lc))


dev.off()

}


