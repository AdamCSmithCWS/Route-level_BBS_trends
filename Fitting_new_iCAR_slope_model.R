## Fitting an new parameterisation of the iCAR model that takes advantage of the
## Stan data type for som-to-zero constrained values
## 
## script currently written to fit the model then save the Stan output to a directory
## 
library(tidyverse)
library(cmdstanr)


output_dir <- "output"
## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")


firstYear <- 2006
lastYear <- 2021


# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
  #species <- species_list[2]
  
  
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
 
spp1  <- "iCAR"
    spp <- paste0("_",spp1,"_")
    
    
    out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
    
    
    
    
    sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
    
    
    load(sp_data_file)
    
    

    # slightly more informative prior on sdalpha intercept 
    # based on an aproximation of the sd of mean log-scale observed counts
    # 
      tmp <- data.frame(route = stan_data$route,count = stan_data$count) %>% 
        group_by(route) %>% 
        summarise(mean_count = mean(log(count+1))) # plus 1 to ensure approximate sd values are positive 
      sd_alpha_prior <- sd(tmp$mean_count) #approximates the observed mean, log-scale variation among routes

    
    stan_data[["sd_alpha_prior"]] <- sd_alpha_prior
    
    mod.file = paste0("models/slope",spp,"route_NB_New.stan")
    
    
    slope_model <- cmdstan_model(mod.file, stanc_options = list("O1"))
    
    stanfit <- slope_model$sample(
      data=stan_data,
      refresh=400,
      chains=4, iter_sampling=2000,
      iter_warmup=2000,
      parallel_chains = 4,
      #pars = parms,
      adapt_delta = 0.8,
      max_treedepth = 10,
      show_exceptions = TRUE)# Change to FALSE if you want to stop Stan spewing the exceptions during the initialisation
                              # but when developing a new model, change to TRUE because these messages are informative.
    
    summ <- stanfit$summary()
  
    stanfit$save_object(paste0(output_dir,"/",out_base,"_stanfit.rds"))
    # saveRDS(stanfit,
    #         paste0(output_dir,"/",out_base,"_stanfit.rds"))
    
    saveRDS(summ,
            paste0(output_dir,"/",out_base,"_summ_fit.rds"))
    
    print(paste(species, round(stanfit$time()[["total"]]/60),"minutes"))
    
}#end species loop

