
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(sf)
library(patchwork)

source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020
output_dir <- "output"



strat = "bbs_usgs"
model = "slope"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include.rds")




firstYear <- 2006
lastYear <- 2021
base_year <- lastYear - floor((lastYear-firstYear)/2) 



# Collect the summaries of each species -----------------------------------


cv_sum = NULL
n_obs <- NULL

for(species in species_list){
  #species <- species_list[2]
  

  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  
  

 

  out_base <- paste0(species_f,"_cv_",firstYear,"_",lastYear)
  
  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_CV_data.RData")
  
  
  load(sp_data_file)
  
  for(sppn in c("iCAR","BYM","nonspatial","GP")){
    

    output_dir <- "output"
    spp <- paste0("_",sppn,"_")
    predictions_save <- readRDS(paste0("output/",species_f,spp,"_pred_save.rds")) 
    
    cv_sum = bind_rows(cv_sum,predictions_save)
    
    
    
    
  }
  n_obst <- full_data %>% 
    group_by(route,ObsN) %>% 
    summarise(n_years = n(),
              .groups = "keep") %>% 
    group_by(route) %>% 
    summarise(sum_n_years = sum(n_years),
              max_nyears = max(n_years),
              mean_nyears = mean(n_years),
              n_obs = n(),
              .groups = "keep") %>% 
    mutate(species = species)
  
  n_obs <- bind_rows(n_obs,n_obst)
  
}




# explore and plot comparison ---------------------------------------------




# point wise differences among models -------------------------------------
cv_sum$Year <- factor(cv_sum$r_year)

simpl_sum <- cv_sum %>% 
  group_by(species,model) %>% 
  summarise(mean = mean(log_lik_mean),
            meanp = mean(E_pred_mean),
            se = sd(log_lik_mean)/sqrt(n()),
            lci = mean - se*1.96,
            uci = mean + se*1.96)

sum_plot <- ggplot(simpl_sum,
                   aes(x = species,y = mean, colour = model))+
  geom_pointrange(aes(ymin = lci,ymax = uci),
                  position = position_dodge(width = 0.8))+
  coord_flip()+
  scale_colour_viridis_d()
sum_plot



diffs <- cv_sum %>% 
  select(species,count,Year,r_year,route,observer,model,log_lik_mean) %>% 
  pivot_wider(.,names_from = model,
              values_from = c(log_lik_mean))  %>% 
  mutate(iCAR_BYM = iCAR - BYM,
         iCAR_GP = iCAR - GP,
         iCAR_nonspatial = iCAR-nonspatial,
         BYM_nonspatial = BYM-nonspatial,
         GP_nonspatial = GP - nonspatial) %>% 
  left_join(.,n_obs,by = c("route","species"))

cv_sum <- cv_sum %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(model = factor(model,levels = c("nonspatial","BYM","iCAR","GP"),
                        ordered = FALSE),
         I = paste(species,E_pred_i,sep = "-"))

# save(list = c("diffs","cv_sum"),
#      file = "data/cv_summary_25_data.RData")


lpos = function(x){
  p = length(which(x > 0))/length(x)
}
mndiffs = diffs %>% 
  group_by(Year,species) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_GP = median(iCAR_GP),
            m_iCAR_nonspatial = median(iCAR_nonspatial),
            m_BYM_nonspatial = median(BYM_nonspatial),
            m_GP_nonspatial = median(GP_nonspatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial))
mndiffs

y_diffs_bym <- ggplot(data = mndiffs,
                  aes(y = m_iCAR_BYM,x = species,
                      colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))
y_diffs_icar <- ggplot(data = mndiffs,
                      aes(y = m_iCAR_nonspatial,x = species,
                          colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))
y_diffs_gp <- ggplot(data = mndiffs,
                       aes(y = m_iCAR_GP,x = species,
                           colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp + y_diffs_icar + y_diffs_bym + plot_layout(guides = "collect")


mndiffs = diffs %>% 
  group_by(species) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_GP = median(iCAR_GP),
            m_iCAR_nonspatial = median(iCAR_nonspatial),
            m_BYM_nonspatial = median(BYM_nonspatial),
            m_GP_nonspatial = median(GP_nonspatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial))
mndiffs










