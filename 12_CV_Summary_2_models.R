
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(sf)
library(patchwork)
library(ggdist)

source("functions/posterior_summary_functions.R") ## functions similar to tidybayes that work on cmdstanr output
## changes captured in a commit on Nov 20, 2020
output_dir <- "output"



strat = "bbs_usgs"
model = "slope"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")

species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")

sp_list_full <- c(species_list,species_list_broad)


firstYear <- 2006
lastYear <- 2021
base_year <- lastYear - floor((lastYear-firstYear)/2) 



resummarise <- FALSE

# Collect the summaries of each species -----------------------------------

if(resummarise){

cv_sum = NULL
n_obs <- NULL

for(species in sp_list_full){
  #species <- species_list[2]
  

  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  
  if(!file.exists(paste0("output/",species_f,"_nonspatial__pred_save.rds"))){next}
  

 

  out_base <- paste0(species_f,"_cv_",firstYear,"_",lastYear)
  
  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_CV_data.RData")
  
  
  load(sp_data_file)
  
  for(sppn in c("iCAR","nonspatial")){
    

    output_dir <- "output"
    spp <- paste0("_",sppn,"_")
    
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



saveRDS(cv_sum,"output/saved_cross_validation_summary_2_models.rds")
saveRDS(n_obs,"output/saved_nobs_cross_validation_summary_2_models.rds")


}else{
  cv_sum <- readRDS("output/saved_cross_validation_summary_2_models.rds")
    n_obs <- readRDS("output/saved_nobs_cross_validation_summary_2_models.rds")
  }
#end if resummarise
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
  mutate(iCAR_nonspatial = iCAR-nonspatial) %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(dist_cat = ifelse(min_distance > 0.1,"Isolated","Close"),
         dist_cat_mean = ifelse(mean_distance > 0.1,"Isolated","Close"))

cv_sum <- cv_sum %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(model = factor(model,levels = c("nonspatial","BYM","iCAR","GP"),
                        ordered = FALSE),
         I = paste(species,E_pred_i,sep = "-"))




lpos = function(x){
  p = length(which(x > 0))/length(x)
}
mndiffs_y = diffs %>% 
  #filter(species %in% c("Dickcissel","White-crowned Sparrow","Sharp-shinned Hawk")) %>% 
  group_by(Year,species) %>% 
  summarise(n_routes = length(unique(route)),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
           nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
           mean_dist = mean(mean_distance),
           p_isolated = sum(ifelse(dist_cat == "Isolated",TRUE,FALSE))/n_routes,
           p_isolated_mean = sum(ifelse(dist_cat_mean == "Isolated",TRUE,FALSE))/n_routes) %>% 
  mutate(species = fct_reorder(species,mean_dist) )
mndiffs_y


y_diffs_icar <- ggplot(data = mndiffs_y,
                      aes(y = z_iCAR_nonspatial,x = species,
                          colour = as.integer(Year)))+
  geom_point()+
  scale_colour_viridis_c()+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip()#ylim = c(-0.1,0.1))

y_diffs_icar



mndiffs_sp = diffs %>% 
  group_by(species) %>% 
  summarise(n_routes = length(unique(route)),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            mean_dist = mean(mean_distance),
            p_isolated = sum(ifelse(dist_cat == "Isolated",TRUE,FALSE))/n(),
            p_isolated_mean = sum(ifelse(dist_cat_mean == "Isolated",TRUE,FALSE))/n()) %>% 
  mutate(species = fct_reorder(species,z_iCAR_nonspatial))


mndiffs_sp

# 
# 
# y_diffs_icar <- ggplot(data = mndiffs_sp,
#                        aes(y = m_iCAR_nonspatial,x = species))+
#   geom_point()+
#   geom_errorbar(aes(ymin = lci_iCAR_nonspatial,
#                     ymax = uci_iCAR_nonspatial),
#                 width = 0,
#                 alpha = 0.2)+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   theme_classic()+
#   ylab("Mean difference in pointwise lpd (iCAR - Non-spatial)")+
#   xlab("")+
#   #coord_cartesian()+
#   coord_flip(ylim = c(-0.1,0.1))
# 
# pdf("Figures/CV_comparison_iCAR_nonspatial.pdf",
#     height = 10,
#     width = 7.5)
# y_diffs_icar
# dev.off()

mndiffs_sp <- mndiffs_sp %>% 
  mutate(nbet_d = cut(nbet_iCAR_nonspatial,breaks = seq(0,10,1)/10),
         base = 0)

y_diffs_icar <- ggplot(data = mndiffs_sp,
                       aes(y = z_iCAR_nonspatial,x = species))+
  #geom_point(aes(colour = nbet_d))+
  geom_errorbar(aes(colour = nbet_d,ymin = base,ymax = z_iCAR_nonspatial),
            linewidth = 0.5,
            width = 0)+
  scale_color_viridis_d(direction = -1,
                        begin = 0.2,end = 0.9)+
  guides(colour = guide_legend(title = "Proportion supporting \n spatial model"))+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = c(-2,2),
             alpha = 0.5)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 4))+
  ylab("Z-score difference in pointwise lpd (iCAR - Non-spatial)")+
  xlab("")+
  #coord_cartesian()+
  coord_flip()#ylim = c(-0.1,0.1))

pdf("Figures/CV_comparison_iCAR_nonspatial.pdf",
    height = 10,
    width = 7.5)
y_diffs_icar
dev.off()





# Raincloud plot ----------------------------------------------------------


bns <- seq(-1,22,by = 1)


y_diffs_hist <- ggplot(data = mndiffs_sp,
                       aes(y = z_iCAR_nonspatial,
                           colour = nbet_iCAR_nonspatial))+
  ggdist::stat_halfeye(.width = 0,
                       adjust = 0.5,
                       justification = 0,
                       point_color = NA)+
  # geom_boxplot(width = 0.12,
  #              alpha = 0.5)+
  ggdist::stat_dots(side = "left",
                    justification = 1.1,
                    binwidth = 0.2)+
  theme_bw()+
  ylab("Z-score difference in pointwise lpd (iCAR - Non-spatial)")+
  xlab("")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  geom_hline(yintercept = 2, alpha = 0.4)+
  geom_hline(yintercept = 0, alpha = 1)+
  coord_flip(xlim = c(-3,1))

y_diffs_hist


# Spatial pattern in route level differences ------------------------------





mndiffs_route <- diffs %>% 
  group_by(species,route) %>% 
  summarise(m_iCAR_nonspatial = mean(iCAR_nonspatial),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            mean_dist = mean(mean_distance))



pdf("Figures/iCAR_vs_NonSpatial_preference_map.pdf")
for(sp in species_list_broad){
  #species <- species_list[2]
  
  
  species_f <- gsub(gsub(sp,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_CV_data.RData")
  
  
  load(sp_data_file)
  
  cv_comp <- mndiffs_route %>% 
    filter(species == sp)
  
  
  cv_map <- route_map %>% 
    inner_join(.,cv_comp,
               by = "route")
  
  bks <- c(0,0.33,0.67,1)
  
  icar_ns <- ggplot(data = cv_map)+
    geom_sf(aes(colour = nbet_iCAR_nonspatial))+
    colorspace::scale_color_continuous_diverging(n_interp = 3,
                                      breaks = bks,
                                      rev = TRUE,
                                      mid = 0.5,
                                      name = "Proportion")+
    labs(title = paste(sp,"Support for iCAR over Non-spatial by route"),
         subtitle = "Porportion of predictions supporting iCAR over Non-spatial \n values > 0.5 (blue) support first model grey values indicate similar support")+
    theme_bw()
  
  print(icar_ns)
  
  
  
  
}
dev.off()







