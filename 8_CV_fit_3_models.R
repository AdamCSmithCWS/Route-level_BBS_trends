## 1-step ahead, cross-validation of three route-level trend models for the BBS
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(sf)
source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020
output_dir <- "output"



strat = "bbs_usgs"
model = "slope"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")


spp <- "_cv_"

firstYear <- 2006
lastYear <- 2021
base_year <- lastYear - floor((lastYear-firstYear)/2) 


rerun <- TRUE

if(rerun){
  # SPECIES LOOP ------------------------------------------------------------
  
  conv_rerun <- readRDS("data/convergence_fail_rerun.rds") %>% 
    filter(model != "GP") %>% 
    rename(speciesl = species) %>% 
    arrange(speciesl)
  
  sp_list_rerun <- unique(conv_rerun$speciesl)
}



for(species in sp_list_rerun){

  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  

# CROSS-VALIDATION loop through the annual re-fitting --------------------------------------
  
  if(rerun){
    wr <- conv_rerun %>% 
      filter(speciesl == species)
    modsToRun <- as.character(wr$model)
    
  }else{
    modsToRun <- c("BYM","iCAR","nonspatial")
  }
  

  for(sppn in modsToRun){
    

  spp <- paste0("_",sppn,"_")
  
  if(!rerun){
  if(file.exists(paste0("output/",species_f,spp,"_pred_save.rds"))){
    next
  }
  }
  
  out_base_1 <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  
  
  
  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_CV_data.RData")
  
  
 load(sp_data_file)
  

  predictions_save <- NULL
  

for(ynext in (base_year+1):lastYear){
  if(ynext == 2020){next} #there are no BBS data in 2020 to predict
  out_base <- paste0(species_f,spp,firstYear,"_",ynext,"_CV")
  
  sp_file <- paste0(output_dir,"/",out_base,".RData")
  
  # setting up the fitting data ------------------------------------------
  
  
  obs_df_fit <- full_data %>% 
    filter(r_year <= ynext-1) %>% 
    mutate(observer = as.integer(factor(ObsN)))
  
  stan_data <- list(count = obs_df_fit$count,
                    year = obs_df_fit$year,
                    route = obs_df_fit$routeF,
                    firstyr = obs_df_fit$firstyr,
                    observer = obs_df_fit$observer,
                    nobservers = max(obs_df_fit$observer),
                    nyears = max(obs_df_fit$year),
                    nroutes = nrow(route_map),
                    ncounts = length(obs_df_fit$count),
                    fixedyear = floor(max(obs_df_fit$year)/2))
  obs_df <- obs_df_fit %>% 
    select(observer,ObsN) %>% 
    distinct()

  if(spp != "_nonspatial_"){
    
  stan_data[["N_edges"]] = car_stan_dat$N_edges
  stan_data[["node1"]] = car_stan_dat$node1
  stan_data[["node2"]] = car_stan_dat$node2
  }  
  
  # setting up the prediction data ------------------------------------------
  
  # slightly more informative prior on sdalpha intercept sd
  # this seems to help for low-abundance species where the 
  # half-normal sd = 2 prior includes a lot of prior mass at extreme values
  # e.g., 100 fold changes in mean abundance, when the counts range from 0 - 2
  # alternate prior sets the sd of the intercept variance to the observed sd
  # of mean counts among all routes. This prior helps convergence for some
  # species (some of the birds of prey and waterbirds that have very low abundance)
  if(rerun){
    tmp <- data.frame(route = stan_data$route,count = stan_data$count) %>% 
      group_by(route) %>% 
      summarise(mean_count = mean(count))
    sd_alpha_prior <- sd(tmp$mean_count)
  }else{
    sd_alpha_prior <- 2 
  }
  
  stan_data[["sd_alpha_prior"]] <- sd_alpha_prior
  
  
  obs_df_predict <- full_data %>% 
    filter(r_year == ynext) %>% 
    left_join(.,obs_df,
              by = "ObsN") %>% 
    mutate(observer = ifelse(!is.na(observer),observer,0))
  
  stan_data[["route_pred"]] <- obs_df_predict$routeF
  stan_data[["count_pred"]] <- obs_df_predict$count
  stan_data[["firstyr_pred"]] <- obs_df_predict$firstyr
  stan_data[["observer_pred"]] <- obs_df_predict$observer
  stan_data[["ncounts_pred"]] <- length(obs_df_predict$count)
  
  
  mod.file = paste0("models/slope",spp,"route_NB_CV.stan")
  
  ## compile model
  slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))
  
  

  slope_stanfit <- slope_model$sample(
    data=stan_data,
    refresh=400,
    chains=3, iter_sampling=1000,
    iter_warmup=1000,
    parallel_chains = 3,
    #pars = parms,
    adapt_delta = 0.8,
    max_treedepth = 10)
  


  

  
  log_lik_samples_full <- posterior_samples(fit = slope_stanfit,
                                            parm = "log_lik",
                                            dims = "i") 
  
  log_lik_samples <- log_lik_samples_full %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(log_lik_samples) <- paste0("log_lik_",names(log_lik_samples))
  
  
  E_pred_samples_full <- posterior_samples(fit = slope_stanfit,
                                           parm = "E_pred",
                                           dims = "i") 
  
  E_pred_samples <- E_pred_samples_full %>% 
    posterior_sums(.,quantiles = NULL,dims = "i") 
  names(E_pred_samples) <- paste0("E_pred_",names(E_pred_samples))
  
  
  obs_df_predict_out <- bind_cols(obs_df_predict,log_lik_samples)
  obs_df_predict_out <- bind_cols(obs_df_predict_out,E_pred_samples)
  obs_df_predict_out$species <- species
  obs_df_predict_out$model <- sppn
  obs_df_predict_out$base <- out_base
  
  
  
  predictions_save <- bind_rows(predictions_save,obs_df_predict_out)
  

  print(paste("Finished",sppn,ynext))
  
  
  saveRDS(predictions_save,file = paste0("output/",species_f,spp,"_pred_save.rds"))
  
  
}
  
  
 




}
  
}#end species loop
