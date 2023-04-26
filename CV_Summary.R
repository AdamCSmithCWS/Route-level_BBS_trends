
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
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")




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
    
    if(!file.exists(paste0("output/",species_f,spp,"_pred_save.rds"))){next}
    predictions_save <- readRDS(paste0("output/",species_f,spp,"_pred_save.rds")) 
    
    cv_sum = bind_rows(cv_sum,predictions_save)
    
    
    
    
  }
 
  min_pos <- function(x){
    y = x[which(x > 0)]
    return(min(y,na.rm = TRUE))
  }
  
  mean_pos <- function(x){
    y = x[which(x > 0)]
    return(mean(y,na.rm = TRUE))
  }
  
  adj_mat <- car_stan_dat$adj_matrix
  neig_dist_mat <- dist_matrix_km
  for(i in 1:nrow(neig_dist_mat)){
    neig_dist_mat[i,] <- as.numeric(neig_dist_mat[i,])*as.numeric(adj_mat[i,])
  }

  dist_min <- apply(dist_matrix_km,
                    MARGIN = 2,
                    FUN = min_pos)
  
  
  dist_mean <- apply(neig_dist_mat,
                    MARGIN = 2,
                    FUN = mean_pos)
  dist_df <- data.frame(routeF = 1:length(dist_min),
                        min_distance = as.numeric(dist_min),
                        mean_distance = as.numeric(dist_mean))
  
  n_obst <- full_data %>% 
    group_by(route,routeF,ObsN) %>% 
    summarise(n_years = n(),
              .groups = "keep") %>% 
    group_by(route,routeF) %>% 
    summarise(sum_n_years = sum(n_years),
              max_nyears = max(n_years),
              mean_nyears = mean(n_years),
              n_obs = n(),
              .groups = "drop") %>% 
    mutate(species = species) %>% 
    inner_join(dist_df,
               by = "routeF") %>% 
    select(-routeF)
  
  n_obs <- bind_rows(n_obs,n_obst)
  
}




# explore and plot comparison ---------------------------------------------




# point wise differences among models -------------------------------------
cv_sum$Year <- factor(cv_sum$r_year)

simpl_sum <- cv_sum %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  group_by(species,model) %>% 
  summarise(mean = mean(log_lik_mean),
            meanp = mean(E_pred_mean),
            se = sd(log_lik_mean)/sqrt(n()),
            lci = mean - se*1.96,
            uci = mean + se*1.96,
            mean_dist = mean(mean_distance)) %>% 
  mutate(species = fct_reorder(species,mean_dist) )

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
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(dist_cat = ifelse(min_distance > 0.1,"Isolated","Close"),
         dist_cat_mean = ifelse(mean_distance > 0.2,"Isolated","Close"))

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
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            mean_dist = mean(mean_distance)) %>% 
  mutate(species = fct_reorder(species,mean_dist) )
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

y_diffs_gp2 <- ggplot(data = mndiffs,
                     aes(y = m_GP_nonspatial,x = species,
                         colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp2 + y_diffs_icar + y_diffs_gp + y_diffs_bym + plot_layout(guides = "collect")


mndiffs_sp = diffs %>% 
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
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            mean_dist = mean(mean_distance)) %>% 
  mutate(species = fct_reorder(species,mean_dist) )
mndiffs_sp



y_diffs_bym <- ggplot(data = mndiffs_sp,
                      aes(y = m_iCAR_BYM,x = species))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_icar <- ggplot(data = mndiffs_sp,
                       aes(y = m_iCAR_nonspatial,x = species))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp <- ggplot(data = mndiffs_sp,
                     aes(y = m_iCAR_GP,x = species))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp2 <- ggplot(data = mndiffs_sp,
                      aes(y = m_GP_nonspatial,x = species))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp2 + y_diffs_icar + y_diffs_gp + y_diffs_bym + plot_layout(guides = "collect")





# By route distances ------------------------------------------------------

mndiffs = diffs %>% 
  group_by(species,dist_cat_mean) %>% 
#  group_by(species,route,min_distance,sum_n_years,n_obs,max_nyears) %>% 
  summarise(m_iCAR_BYM = median(iCAR_BYM),
            m_iCAR_GP = median(iCAR_GP),
            m_iCAR_nonspatial = median(iCAR_nonspatial),
            m_BYM_nonspatial = median(BYM_nonspatial),
            m_GP_nonspatial = median(GP_nonspatial),
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            n_routes = n(),
            mean_dist = mean(mean_distance)) %>% 
  mutate(species = factor(species,ordered = TRUE,
                          levels = levels(mndiffs_sp$species)),
            .groups = "keep") 
mndiffs

# y_diffs_bym <- ggplot(data = mndiffs,
#                       aes(y = m_iCAR_BYM,x = min_distance,
#                           colour = max_nyears))+
#   geom_point()+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   scale_x_continuous(trans = "log")+
#   theme(legend.position = "none")+
#   facet_wrap(vars(species),
#              scales = "free")
# y_diffs_bym


y_diffs_gp <- ggplot(data = mndiffs,
                      aes(y = m_iCAR_GP,x = species,
                          colour = dist_cat_mean))+
  geom_point()+
  scale_colour_viridis_d(direction = -1,
                         begin = 0.1,end = 0.8)+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip()#ylim = c(-0.1,0.1))

y_diffs_gp


y_diffs_icar <- ggplot(data = mndiffs,
                     aes(y = m_iCAR_nonspatial,x = species,
                         colour = dist_cat))+
  geom_point()+
  scale_colour_viridis_d(direction = -1,
                         begin = 0.1,end = 0.8)+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip()#ylim = c(-0.1,0.1))

y_diffs_icar



y_diffs_icar_gp <- ggplot(data = mndiffs,
                      aes(y = m_iCAR_GP,colour = mean_distance,
                          x = max_nyears))+
  geom_point(position = position_dodge(width = 0.25),
             alpha = 0.2)+
  scale_colour_viridis_c(direction = -1)+
  geom_hline(yintercept = 0)+
  #theme(legend.position = "none")+
  facet_wrap(vars(species),
             scales = "free")

y_diffs_icar_gp



y_diffs_icar <- ggplot(data = mndiffs,
                       aes(y = m_iCAR_GP,x = min_distance,
                           colour = max_nyears))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  theme(legend.position = "none")+
  facet_wrap(vars(species),
             scales = "free")

y_diffs_icar






y_diffs_gp <- ggplot(data = mndiffs,
                     aes(y = m_iCAR_GP,x = species,
                         colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp2 <- ggplot(data = mndiffs,
                      aes(y = m_GP_nonspatial,x = species,
                          colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip(ylim = c(-0.1,0.1))

y_diffs_gp2 + y_diffs_icar + y_diffs_gp + y_diffs_bym + plot_layout(guides = "collect")






