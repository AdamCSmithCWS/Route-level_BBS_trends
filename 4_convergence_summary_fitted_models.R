## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(cmdstanr)
library(sf)
library(patchwork)


output_dir <- "output"
output_dir <- "F:/iCAR_route_2021/output"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
species_list <- c(species_list,"Blue-headed Vireo")

species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")

species_list <- unique(c(species_list,species_list_broad))
firstYear <- 2006
lastYear <- 2021


conv_save <- NULL
# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
#species <- species_list[2]

  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
 # out_base_temp <- paste0(species_f,"_GP_",firstYear,"_",lastYear)
  
 # if(!file.exists(paste0(output_dir,"/",out_base_temp,"_summ_fit.rds"))){next}
    
#     
# 
# sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
# 
# 
# load(sp_data_file)



for(spp1 in c("BYM","iCAR","nonspatial","GP")){

  spp <- paste0("_",spp1,"_")
  

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)



if(!file.exists(paste0(output_dir,"/",out_base,"_summ_fit.rds"))){next}

#stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
 
summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))

summ <- summ %>% 
  mutate(species = species,
         model = spp1)

conv_save <- bind_rows(conv_save,summ)

} #end models loop


}#end species loop



conv_sum <- conv_save %>% 
  filter(!grepl(variable,pattern = "_LK")) %>% 
  mutate(variable_type = str_extract(pattern = "([[:alpha:]]|_){1,}",
                                     string = variable))

#saveRDS(conv_sum,"data/convergence_summary_May31.rds")

conv_sum <- readRDS("data/convergence_summary_May31.rds")
# tmp <- conv_sum %>% 
#   filter(model == "GP",
#          species == "Lazuli Bunting")



sdobs <- conv_sum %>% 
  filter(grepl("sdobs",variable_type)) %>% 
  select(species,model,mean) %>% 
  rename(sdobs = mean)
sdroute <- conv_sum %>% 
  filter(grepl("sdalpha",variable_type)) %>% 
  select(species,model,mean) %>% 
  rename(sdroute = mean)

sd_comp <- sdobs %>% 
  left_join(.,sdroute,
            by = c("species","model")) %>% 
  mutate(sddiff = sdroute - sdobs)

hists <- ggplot(data = sd_comp) +
  geom_freqpoly(aes(sddiff,colour = model))
hists


n_fail_rhat <- function(x,th = 1.02){
  length(which(x > th))
}
n_fail_ess <- function(x,th = 100){
  length(which(x < th))
}

conv_summary <- conv_sum %>% 
  group_by(variable_type,species,model) %>% 
  summarise(mean_ess = mean(ess_bulk),
            min_ess = min(ess_bulk),
            mean_rhat = mean(rhat,na.rm = T),
            max_rhat = max(rhat,na.rm = T),
            n_fail_rhat = n_fail_rhat(rhat),
            n_fail_ess = n_fail_ess(ess_bulk))

conv_summary2 <- conv_sum %>% 
  group_by(species,model) %>% 
  summarise(mean_ess = mean(ess_bulk),
            min_ess = min(ess_bulk),
            mean_rhat = mean(rhat,na.rm = T),
            max_rhat = max(rhat,na.rm = T),
            n_fail_rhat = n_fail_rhat(rhat,1.02),
            n_fail_ess = n_fail_ess(ess_bulk,100))


conv_fail <- conv_summary2 %>% 
  filter((n_fail_ess > 0 ) |
         n_fail_rhat > 0)



tmp <- conv_sum %>% 
  filter(grepl("gp_sq",variable_type))

tmp2 <- conv_sum %>% 
  filter(grepl("gp_eta_beta",variable_type))
