## Fitting the BYM model to 1995 - 2021 BBS data on Azure




# load and stratify BBS data ---------------------------------------------

#setwd("C:/GitHub/BBS_iCAR_route_trends")
library(bbsBayes)
library(tidyverse)
#library(cmdstanr)
#source("functions/neighbours_define.R") ## function to define neighbourhood relationships
#source("functions/prepare-data-alt.R") ## small alteration of the bbsBayes function
## above source() call over-writes the bbsBayes prepare_data() function
# source("functions/get_basemap_function.R") ## loads one of the bbsBayes strata maps
# source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output



output_dir <- ""



strat = "bbs_usgs"
model = "slope"

species_list <- c("Western Meadowlark",
                  "Savannah Sparrow")
firstYear = 1995
lastYear = 2021

# strat_data <- bbsBayes::stratify(by = strat)
stan_data_list <- vector(mode = "list",length = length(species_list))
names(stan_data_list) <- gsub(gsub(species_list,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

for(species in species_list){
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  spp <- "_BYM_"
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  load(paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData"))
  
  stan_data_list[[species_f]] <- stan_data
}

mod.file = "models/slope_BYM_route_NB.stan"

mod_string <- read_file(mod.file)














# parallel function -------------------------------------------------------

run_model_azure <- function(species){
  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp <- "_BYM_"

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

stan_data <- stan_data_list[[species_f]]


mod_file <- cmdstanr::write_stan_file(code = mod_string)

slope_model <- cmdstanr::cmdstan_model(mod_file, stanc_options = list("Oexperimental"))


# save(list = c("stan_data",
#               "jags_data",
#               "route_map",
#               "realized_strata_map",
#               "firstYear"),
#      file = sp_data_file)






  
  init_def <- function(){ list(alpha_raw = rnorm(stan_data$nroutes,0,0.1),
                               sdalpha = runif(1,0.01,0.1),
                               ALPHA = 0,
                               BETA = 0,
                               eta = 0,
                               obs_raw = rnorm(stan_data$nobservers,0,0.1),
                               sdnoise = 0.2,
                               sdobs = 0.1,
                               sdbeta_rand = runif(1,0.01,0.1),
                               beta_raw_rand = rnorm(stan_data$nroutes,0,0.01),
                               sdbeta_space = runif(1,0.01,0.1),
                               beta_raw_space = rnorm(stan_data$nroutes,0,0.01))} 




stanfit <- slope_model$sample(
  data=stan_data,
  refresh=200,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  # output_dir = output_dir,
  # output_basename = out_base,
  init = init_def
  )

tmp <- stanfit$draws()
tmp <- stanfit$sampler_diagnostics()

ret <- list(stanfit = stanfit)

# save(list = c("stanfit","summ"),
#      file = paste0(output_dir,"/",out_base,"_stanfit.RData"))
return(ret)
}





# AZURE Run ---------------------------------------------------------------


# Load the doAzureParallel library 
library(doAzureParallel) 

setCredentials("credentials.json")

# generate your cluster in the cloud; this takes a few minutes
cluster <- makeCluster("cluster.json") 
#this .json file currently includes instructions to use a rocker Bayesian image


# Register your parallel backend 
registerDoAzureParallel(cluster) 

# Check that the nodes are running 
getDoParWorkers() 














opt <- list(wait = TRUE, #you can set this to false if you want to free-up your local machine, but you then have to frequently check back with the cluster to download the results
            enableCloudCombine = TRUE) #this combines the results into a single R-object list before downloading


t1 = Sys.time()


results <- foreach(i = species_list,
                   .packages = c("cmdstanr"),
                   .errorhandling = "pass", #this option passes any errors back to the console
                   .options.azure = opt) %dopar% {
                     
                     # This code is executed, in parallel, across your cluster.
                     run_model_azure(species = i)
                   }

output_dir <- "output"

for(species in species_list){
  i <- which(species_list == species)
  stanfit <- results[[i]][["stanfit"]]
  
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  spp <- "_BYM_"
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
# stanfit$save_object(dir = output_dir, basename = out_base, timestamp = FALSE, random = FALSE)
stanfit$save_object(file = paste0(output_dir,"/",out_base,"_stanfit.RDS"))
# 
# save(list = c("summ"),
#      file = paste0(output_dir,"/",out_base,"_summ.RData"))
  
}
#results

t2 = Sys.time()
t2-t1
### when finished, be sure to run this stopCluster function, so you avoid paying for idle virtual machines
stopCluster(cluster)






# Post cluster summaries --------------------------------------------------
library(tidyverse)
library(cmdstanr)

species_list <- c("Chestnut-collared Longspur",
                  "Thick-billed Longspur")
firstYear = 1995
lastYear = 2021

source("Functions/posterior_summary_functions.R")

output_dir <- "output"

for(species in species_list){

  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  spp <- "_BYM_"
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  # stanfit$save_object(dir = output_dir, basename = out_base, timestamp = FALSE, random = FALSE)
  stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.RDS"))
  
  
  slopes <- posterior_samples(fit = stanfit,
                              parm = "beta",
                              dims = "ste") 
  intercepts <- posterior_samples(fit = stanfit,
                              parm = "alpha",
                              dims = "ste")

  slope_sum <- posterior_sums(slopes,dims = "ste")
  
}












