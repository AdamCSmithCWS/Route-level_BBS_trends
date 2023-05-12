## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)

source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output



output_dir <- "F:/iCAR_route_2021/output"



# load and stratify CASW data ---------------------------------------------
strat = "bbs_usgs"
model = "slope"


species_list <- readRDS("data/species_to_include_4_model_comparison.rds")


firstYear <- 2006
lastYear <- 2021
base_year <- lastYear - floor((lastYear-firstYear)/2) 



# SPECIES LOOP ------------------------------------------------------------

# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options

for(species in species_list){
#species <- species_list[2]

species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)




  spp <- "_GP_"
  
out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")

load(sp_data_file)


units(dist_matrix_km) <- NULL



stan_data[["distances"]] = dist_matrix_km

stan_data[["N_edges"]] <- NULL
stan_data[["node1"]] <- NULL
stan_data[["node2"]] <- NULL





mod.file = "models/slope_GP_route_NB.stan"



slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data=stan_data,
  refresh=100,
  chains=4, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 4,
  adapt_delta = 0.8,
  max_treedepth = 10)

print(species)
 print(stanfit$time())

 vars <- stanfit$metadata()
 vars <- vars$variables
 vars_sum <- vars[!grepl("_LK",vars)]
 
 summ <- stanfit$summary(variables = vars_sum)
 
 
saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))

}


