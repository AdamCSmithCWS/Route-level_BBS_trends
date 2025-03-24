## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(cmdstanr)




output_dir <- "F:/iCAR_route_2021/output"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")


firstYear <- 2006
lastYear <- 2021

skip_bym <- FALSE
rerun <- FALSE

if(rerun){
# SPECIES LOOP ------------------------------------------------------------

conv_rerun <- readRDS("data/convergence_fail_rerun.rds") %>% 
  filter(model != "GP") %>% 
  rename(speciesl = species) %>% 
  arrange(speciesl)

species_list <- unique(conv_rerun$speciesl)
}
# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
#species <- species_list[2]

  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)


if(rerun){
  wr <- conv_rerun %>% 
    filter(speciesl == species)
  modsToRun <- as.character(wr$model)
  
}else{
modsToRun <- c("BYM","iCAR","nonspatial")
}

for(spp1 in modsToRun){

  if(skip_bym & spp1 == "BYM"){next}
  
  spp <- paste0("_",spp1,"_")
  

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)


if(spp1 == "nonspatial"){
  
  stan_data[["N_edges"]] <- NULL
  stan_data[["node1"]] <- NULL
  stan_data[["node2"]] <- NULL
  
}

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
  
   mod.file = paste0("models/slope",spp,"route_NB.stan")


slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

stanfit <- slope_model$sample(
  data=stan_data,
  refresh=200,
  chains=4, iter_sampling=2000,
  iter_warmup=2000,
  parallel_chains = 4,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 10)

summ <- stanfit$summary()
print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))


} #end models loop

}#end species loop

