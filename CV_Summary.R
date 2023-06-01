
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

species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")

species_list <- c(species_list,"Blue-headed Vireo")
species_list <- c(species_list,"Lazuli Bunting")
species_list <- c(species_list,"Scissor-tailed Flycatcher")


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
  mutate(dist_cat = ifelse(min_distance > 0.2,"Isolated","Close"),
         dist_cat_mean = ifelse(mean_distance > 0.1,"Isolated","Close"))

cv_sum <- cv_sum %>% 
  left_join(.,n_obs,by = c("route","species")) %>% 
  mutate(model = factor(model,levels = c("nonspatial","BYM","iCAR","GP"),
                        ordered = FALSE),
         I = paste(species,E_pred_i,sep = "-"))

# save(list = c("diffs","cv_sum"),
#      file = "data/cv_summary_25_data.RData")


# 
lpos = function(x){
  p = length(which(x > 0))/length(x)
}
# annual cv summary
# mndiffs = diffs %>%
#   group_by(Year,species) %>%
#   summarise(m_iCAR_BYM = mean(iCAR_BYM),
#             m_iCAR_GP = mean(iCAR_GP),
#             m_iCAR_nonspatial = mean(iCAR_nonspatial),
#             m_BYM_nonspatial = mean(BYM_nonspatial),
#             m_GP_nonspatial = mean(GP_nonspatial),
#             nbet_iCAR_BYM = lpos(iCAR_BYM),
#             nbet_iCAR_GP = lpos(iCAR_GP),
#             nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
#             nbet_BYM_nonspatial = lpos(BYM_nonspatial),
#             nbet_GP_nonspatial = lpos(GP_nonspatial),
#             mean_dist = mean(mean_distance)) %>%
#   mutate(species = fct_reorder(species,mean_dist) )
# mndiffs
# 
# y_diffs_bym <- ggplot(data = mndiffs,
#                   aes(y = m_iCAR_BYM,x = species,
#                       colour = as.integer(Year)))+
#   geom_point()+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian()+
#   coord_flip(ylim = c(-0.1,0.1))
# 
# y_diffs_icar <- ggplot(data = mndiffs,
#                       aes(y = m_iCAR_nonspatial,x = species,
#                           colour = as.integer(Year)))+
#   geom_point()+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian()+
#   coord_flip(ylim = c(-0.1,0.1))
# 
# y_diffs_gp <- ggplot(data = mndiffs,
#                        aes(y = m_iCAR_GP,x = species,
#                            colour = as.integer(Year)))+
#   geom_point()+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian()+
#   coord_flip(ylim = c(-0.25,0.25))
# 
# y_diffs_gp2 <- ggplot(data = mndiffs,
#                      aes(y = m_GP_nonspatial,x = species,
#                          colour = as.integer(Year)))+
#   geom_point()+
#   scale_colour_viridis_c()+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian()+
#   coord_flip(ylim = c(-0.1,0.1))
# 
# y_diffs_gp2 + y_diffs_icar + y_diffs_gp + y_diffs_bym + plot_layout(guides = "collect")


mndiffs_sp <- diffs %>% 
  group_by(species) %>% 
  summarise(m_iCAR_BYM = mean(iCAR_BYM),
            m_iCAR_GP = mean(iCAR_GP),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            m_BYM_nonspatial = mean(BYM_nonspatial),
            m_GP_nonspatial = mean(GP_nonspatial),
            se_iCAR_BYM = sd(iCAR_BYM)/sqrt(n()),
            se_iCAR_GP = sd(iCAR_GP)/sqrt(n()),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            se_BYM_nonspatial = sd(BYM_nonspatial)/sqrt(n()),
            se_GP_nonspatial = sd(GP_nonspatial)/sqrt(n()),
            z_iCAR_BYM = mean(iCAR_BYM)/se_iCAR_BYM,
            z_iCAR_GP = mean(iCAR_GP)/se_iCAR_GP,
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            z_BYM_nonspatial = mean(BYM_nonspatial)/se_BYM_nonspatial,
            z_GP_nonspatial = mean(GP_nonspatial)/se_GP_nonspatial,
            lci_iCAR_BYM = m_iCAR_BYM - se_iCAR_BYM*1.96,
            lci_iCAR_GP = m_iCAR_GP - se_iCAR_GP*1.96,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            lci_BYM_nonspatial = m_BYM_nonspatial - se_BYM_nonspatial*1.96,
            lci_GP_nonspatial = m_GP_nonspatial - se_GP_nonspatial*1.96,
            uci_iCAR_BYM = m_iCAR_BYM + se_iCAR_BYM*1.96,
            uci_iCAR_GP = m_iCAR_GP + se_iCAR_GP*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
            uci_BYM_nonspatial = m_BYM_nonspatial + se_BYM_nonspatial*1.96,
            uci_GP_nonspatial = m_GP_nonspatial + se_GP_nonspatial*1.96,
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            max_dist = max(mean_distance),
            n_routes = length(unique(route)),
            p_isolated = sum(ifelse(dist_cat == "Isolated",TRUE,FALSE))/n_routes,
            p_isolated_mean = sum(ifelse(dist_cat_mean == "Isolated",TRUE,FALSE))/n_routes) %>% 
  mutate(species = fct_reorder(species,p_isolated) )
mndiffs_sp

z_diffs_plot <- mndiffs_sp %>% 
  select(species,
         starts_with("z_")) %>% 
  pivot_longer(names_to = "model_comparison",
               values_to = "z",
               cols = -species,
               names_prefix = "z_")

mn_diffs_plot <- mndiffs_sp %>% 
  select(species,
         starts_with("m_")) %>% 
  pivot_longer(names_to = "model_comparison",
               values_to = "mean",
               cols = -species,
               names_prefix = "m_")

lci_diffs_plot <- mndiffs_sp %>% 
  select(species,
         starts_with("lci_")) %>% 
  pivot_longer(names_to = "model_comparison",
               values_to = "lci",
               names_prefix = "lci_",
               cols = -species)


uci_diffs_plot <- mndiffs_sp %>% 
  select(species,
         starts_with("uci_")) %>% 
  pivot_longer(names_to = "model_comparison",
               values_to = "uci",
               names_prefix = "uci_",
               cols = -species)


diffs_plot <- inner_join(mn_diffs_plot,
                         lci_diffs_plot) %>% 
  inner_join(.,uci_diffs_plot) %>% 
  inner_join(.,z_diffs_plot) %>% 
  #filter(species != "White Ibis") %>% 
  mutate(model_comparison = factor(model_comparison,ordered = TRUE,
                                   levels = c("iCAR_nonspatial",
                                              "BYM_nonspatial",
                                              "GP_nonspatial",
                                              "iCAR_GP",
                                              "iCAR_BYM")))


saveRDS(diffs_plot,"data_cv_summary_4models_plotting_data.rds")

# y_diffs_all <- ggplot(data = diffs_plot,
#                       aes(y = mean,x = species))+
#   geom_point()+
#   geom_errorbar(aes(ymin = lci,
#                     ymax = uci),
#                 alpha = 0.3,
#                 width = 0)+
#   geom_hline(yintercept = 0)+
#   #coord_cartesian()+
#   ylab("Mean difference in pointwise lpd")+
#   xlab("")+
#   coord_flip(ylim = c(-0.5,0.5))+
#   facet_wrap(vars(model_comparison),
#              nrow = 1,
#              scales = "free_x")
# 
# y_diffs_all


diffs_plot_sel <- diffs_plot %>% 
  filter(grepl("nonspatial",model_comparison))

z_diffs_all <- ggplot(data = diffs_plot_sel,
                      aes(y = z,x = species,
                          colour = model_comparison))+
  geom_point()+
  geom_hline(yintercept = 0)+
  scale_colour_viridis_d(direction = -1,
                         end = 0.8)+
  ylab("Z-score difference in pointwise lpd")+
  xlab("")+
  geom_hline(yintercept = c(-2,2),
             alpha = 0.5)+
  theme_bw()+
  coord_flip()#ylim = c(-0.5,0.5))

pdf("Figures/Spatial_vs_non_3models.pdf",
    height = 10,
    width = 7.5)
z_diffs_all
dev.off()



diffs_plot_sel2 <- diffs_plot %>% 
  filter(!grepl("nonspatial",model_comparison))

z_diffs_spat <- ggplot(data = diffs_plot_sel2,
                      aes(y = z,x = species,
                          colour = model_comparison))+
  geom_point()+
  geom_hline(yintercept = 0)+
  scale_colour_viridis_d(direction = -1,
                         end = 0.8)+
  ylab("Z-score difference in pointwise lpd")+
  xlab("")+
  geom_hline(yintercept = c(-2,2),
             alpha = 0.5)+
  theme_bw()+
  coord_flip()#ylim = c(-0.5,0.5))

pdf("Figures/Spatial_options_cv.pdf",
    height = 10,
    width = 7.5)
z_diffs_spat
dev.off()








# By route distances ------------------------------------------------------


mndiffs_sp <- diffs %>% 
  group_by(species) %>% 
  summarise(m_iCAR_BYM = mean(iCAR_BYM),
            m_iCAR_GP = mean(iCAR_GP),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            m_BYM_nonspatial = mean(BYM_nonspatial),
            m_GP_nonspatial = mean(GP_nonspatial),
            se_iCAR_BYM = sd(iCAR_BYM)/sqrt(n()),
            se_iCAR_GP = sd(iCAR_GP)/sqrt(n()),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            se_BYM_nonspatial = sd(BYM_nonspatial)/sqrt(n()),
            se_GP_nonspatial = sd(GP_nonspatial)/sqrt(n()),
            z_iCAR_BYM = mean(iCAR_BYM)/se_iCAR_BYM,
            z_iCAR_GP = mean(iCAR_GP)/se_iCAR_GP,
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            z_BYM_nonspatial = mean(BYM_nonspatial)/se_BYM_nonspatial,
            z_GP_nonspatial = mean(GP_nonspatial)/se_GP_nonspatial,
            lci_iCAR_BYM = m_iCAR_BYM - se_iCAR_BYM*1.96,
            lci_iCAR_GP = m_iCAR_GP - se_iCAR_GP*1.96,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            lci_BYM_nonspatial = m_BYM_nonspatial - se_BYM_nonspatial*1.96,
            lci_GP_nonspatial = m_GP_nonspatial - se_GP_nonspatial*1.96,
            uci_iCAR_BYM = m_iCAR_BYM + se_iCAR_BYM*1.96,
            uci_iCAR_GP = m_iCAR_GP + se_iCAR_GP*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
            uci_BYM_nonspatial = m_BYM_nonspatial + se_BYM_nonspatial*1.96,
            uci_GP_nonspatial = m_GP_nonspatial + se_GP_nonspatial*1.96,
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            max_dist = max(mean_distance),
            n_routes = length(unique(route)),
            p_isolated = sum(ifelse(dist_cat == "Isolated",TRUE,FALSE))/n_routes,
            p_isolated_mean = sum(ifelse(dist_cat_mean == "Isolated",TRUE,FALSE))/n_routes) %>% 
  mutate(species = fct_reorder(species,z_iCAR_GP) )

mndiffs = diffs %>% 
  group_by(species,dist_cat) %>% 
#  group_by(species,route,min_distance,sum_n_years,n_obs,max_nyears) %>% 
  summarise(m_iCAR_BYM = mean(iCAR_BYM),
            m_iCAR_GP = mean(iCAR_GP),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            m_BYM_nonspatial = mean(BYM_nonspatial),
            m_GP_nonspatial = mean(GP_nonspatial),
            se_iCAR_BYM = sd(iCAR_BYM)/sqrt(n()),
            se_iCAR_GP = sd(iCAR_GP)/sqrt(n()),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            se_BYM_nonspatial = sd(BYM_nonspatial)/sqrt(n()),
            se_GP_nonspatial = sd(GP_nonspatial)/sqrt(n()),
            z_iCAR_BYM = mean(iCAR_BYM)/se_iCAR_BYM,
            z_iCAR_GP = mean(iCAR_GP)/se_iCAR_GP,
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            z_BYM_nonspatial = mean(BYM_nonspatial)/se_BYM_nonspatial,
            z_GP_nonspatial = mean(GP_nonspatial)/se_GP_nonspatial,
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            n_routes = length(unique(route)),
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
                          colour = dist_cat))+
  geom_point()+
  scale_colour_viridis_d(direction = -1,
                         begin = 0.1,end = 0.8)+
  geom_hline(yintercept = 0)+
  #coord_cartesian()+
  coord_flip()#ylim = c(-0.1,0.1))

pdf("Figures/CV_GP_iCAR_isolated.pdf",
    height = 10,
    width = 7.5)
y_diffs_gp
dev.off()









# Spatial pattern in route level differences ------------------------------





mndiffs_route <- diffs %>% 
  group_by(species,route) %>% 
  summarise(m_iCAR_BYM = mean(iCAR_BYM),
            m_iCAR_GP = mean(iCAR_GP),
            m_iCAR_nonspatial = mean(iCAR_nonspatial),
            m_BYM_nonspatial = mean(BYM_nonspatial),
            m_GP_nonspatial = mean(GP_nonspatial),
            se_iCAR_BYM = sd(iCAR_BYM)/sqrt(n()),
            se_iCAR_GP = sd(iCAR_GP)/sqrt(n()),
            se_iCAR_nonspatial = sd(iCAR_nonspatial)/sqrt(n()),
            se_BYM_nonspatial = sd(BYM_nonspatial)/sqrt(n()),
            se_GP_nonspatial = sd(GP_nonspatial)/sqrt(n()),
            z_iCAR_BYM = mean(iCAR_BYM)/se_iCAR_BYM,
            z_iCAR_GP = mean(iCAR_GP)/se_iCAR_GP,
            z_iCAR_nonspatial = mean(iCAR_nonspatial)/se_iCAR_nonspatial,
            z_BYM_nonspatial = mean(BYM_nonspatial)/se_BYM_nonspatial,
            z_GP_nonspatial = mean(GP_nonspatial)/se_GP_nonspatial,
            lci_iCAR_BYM = m_iCAR_BYM - se_iCAR_BYM*1.96,
            lci_iCAR_GP = m_iCAR_GP - se_iCAR_GP*1.96,
            lci_iCAR_nonspatial = m_iCAR_nonspatial - se_iCAR_nonspatial*1.96,
            lci_BYM_nonspatial = m_BYM_nonspatial - se_BYM_nonspatial*1.96,
            lci_GP_nonspatial = m_GP_nonspatial - se_GP_nonspatial*1.96,
            uci_iCAR_BYM = m_iCAR_BYM + se_iCAR_BYM*1.96,
            uci_iCAR_GP = m_iCAR_GP + se_iCAR_GP*1.96,
            uci_iCAR_nonspatial = m_iCAR_nonspatial + se_iCAR_nonspatial*1.96,
            uci_BYM_nonspatial = m_BYM_nonspatial + se_BYM_nonspatial*1.96,
            uci_GP_nonspatial = m_GP_nonspatial + se_GP_nonspatial*1.96,
            nbet_iCAR_BYM = lpos(iCAR_BYM),
            nbet_iCAR_GP = lpos(iCAR_GP),
            nbet_iCAR_nonspatial = lpos(iCAR_nonspatial),
            nbet_BYM_nonspatial = lpos(BYM_nonspatial),
            nbet_GP_nonspatial = lpos(GP_nonspatial),
            max_dist = max(mean_distance))



pdf("Figures/iCAR_vs_GP_preference_map.pdf")
for(sp in species_list){
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

  icar_gp <- ggplot(data = cv_map)+
    geom_sf(aes(colour = nbet_iCAR_GP))+
    colorspace::scale_colour_continuous_diverging(n_interp = 3,
                                      breaks = bks,
                                      rev = TRUE,
                                      mid = 0.5,
                                      name = "Proportion")+
    labs(title = paste(sp,"Support for iCAR over GP by route"),
         subtitle = "Porportion of predictions supporting iCAR over GP \n values > 0.5 (blue) support first model grey values indicate similar support")+
     theme_bw()

  print(icar_gp)
  

  
  
  }
dev.off()

