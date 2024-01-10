#Figures

library(tidyverse)
library(patchwork)
library(sf)
 species_sort <- readRDS("data/species_sort.rds")
# Figure_1 ----------------------------------------------------------------

# demonstration of neighbours for BHVI

# load saved neighbours data from neighbours_define_voronoi(..., save_plot_data = TRUE)

load("data/Baird's_Sparrow_route_maps_data.RData")

strata_map2 <- bbsBayes2::load_map("bbs_usgs")

strata_map2 <- st_transform(strata_map2,st_crs(strata_map))

box <- st_as_sfc(st_bbox(strata_map))

xb <- range(st_coordinates(box)[,"X"])
yb <- range(st_coordinates(box)[,"Y"])


ggp1 <- ggplot()+ 
  geom_sf(data = strata_map2,alpha = 0,colour = grey(0.9))+ 
  geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+ 
  #geom_sf(data = cov_hull_clip,alpha = 0,colour = grey(0.8))+
  #geom_sf(data = v,alpha = 0,colour = grey(0.9))+
  geom_sf(data = full_clip,alpha = 0,colour = grey(0.1))+
  geom_sf(data = vintj,alpha = 0,colour = grey(0.5))+
  geom_sf(data = centres,alpha = 1,size = 0.6)+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none",
        text = element_text(family = "serif",
                          size = 11),
      panel.grid = element_line(colour = grey(0.95)),
      plot.caption = element_text(hjust = 0))
ggp1

ggp2 <- ggplot()+ 
  geom_sf(data = strata_map2,alpha = 0,colour = grey(0.9))+ 
  geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+ 
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),
               inherit.aes = FALSE,linewidth=0.3,colour = grey(0.6)) +
  geom_sf(data = centres,alpha = 1,size = 1.3)+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none",
        text = element_text(family = "serif",
                            size = 11),
        panel.grid = element_line(colour = grey(0.95)),
        plot.caption = element_text(hjust = 0))
ggp2


pdf("Figures/Figure_1.pdf",
    width = 7.5,
    height = 10)
print(ggp1 + ggp2 + plot_layout(nrow = 2))
dev.off()






# Alternate GP posterior plot ---------------------------------------------

output_dir <- "output"
species <- "Baird's Sparrow" 
firstYear <- 2006
lastYear <- 2021

species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

out_base_temp <- paste0(species_f,"_GP_",firstYear,"_",lastYear)


sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)

spp <- paste0("_iCAR_")


out_base <- paste0(species_f,spp,firstYear,"_",lastYear)


summ_CAR <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds")) %>% 
  filter(variable %in% c("sdalpha","sdbeta_space"))

summ_CAR
  


spp <- paste0("_GP_")
  
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  
  
  summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds")) %>% 
    filter(variable %in% c("gp_sq_rho_beta",
                           "gp_sq_rho_alpha",
                           "gp_sq_alpha_beta",
                           "gp_sq_alpha_alpha"))
  
  summ
  

  
  stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
  
  
  dist_sim <- expand_grid(dist = seq(0,3,by = 0.02),
                          .draw = 1:4000)
  
  gp_sim_function <- function(alpha2,rho2,dist){
    cv <- alpha2 * exp(-1*rho2 * dist^2)
    return(cv)
  }
  gp_params <- stanfit$draws(variables = c("gp_sq_rho_beta",
                                          "gp_sq_rho_alpha",
                                          "gp_sq_alpha_beta",
                                          "gp_sq_alpha_alpha"),
                            format = "draws_df") %>% 
    sample_n(.,50) %>% 
    left_join(.,dist_sim,by = ".draw") %>% 
    mutate(pred_cov_beta = gp_sim_function(gp_sq_alpha_beta,gp_sq_rho_beta,dist),
           pred_cov_alpha = gp_sim_function(gp_sq_alpha_alpha,gp_sq_rho_alpha,dist),
           km = dist*1000)
    
  col_1 <- RColorBrewer::brewer.pal(4,"Paired")
  
  species_latin <- bbsBayes2::search_species(species)
  species_latin <- paste(species_latin[1,"genus"],species_latin[1,"species"])
  
  capt_tmp <- paste("Figure S1. Sample of posterior draws for the distance-based decay of covariance of
                    abundance and trend among BBS routes for Baird's Sparrow (",species_latin,") estimated using an
                    isotropic spatial Gaussian Process model. The modeled estimates show that the
                    covariance among routes in relative abundance of the species decreases very quickly
                    with increasing distance (green lines), while the covariance among routes in trends
                    decreases very slowly with distance (blue lines).")
  
  spag <- ggplot(data = gp_params)+
    geom_line(aes(x = km,y = pred_cov_beta,
                  group = .draw),
              alpha = 0.9,
              colour = col_1[2])+
    geom_line(aes(x = km,y = pred_cov_alpha,
                  group = .draw),
              alpha = 0.9,
              colour = col_1[3])+
    theme_bw()+
    labs(caption = capt_tmp)+
    theme(plot.caption = element_text(hjust = 0),
          text = element_text(family = "serif",
                              size = 11))+
    ylab("Covariance between pairs of routes")+
    xlab("Distance between routes (km)")
    
  pdf("Figures/Figure_S1.pdf",
      width = 8.5,
      height = 11)
print(spag)
  dev.off()
  
# Figure 2 and 3 - 4 model comparison -------------------------------------------
output_dir <- "F:/iCAR_route_2021/output"
base_strata_map <- bbsBayes2::load_map("bbs_usgs")

  
firstYear <- 2006
lastYear <- 2021

ppy <- function(x){
  p <- exp(x)-1
  return(p*100)
}

species <- "Baird's Sparrow" 
  

species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  out_base_temp <- paste0(species_f,"_GP_",firstYear,"_",lastYear)
  

  sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
  
  
  load(sp_data_file)
  
  
  outboth <- NULL
  for(spp1 in c("BYM","iCAR","nonspatial","GP")){
    
    spp <- paste0("_",spp1,"_")
    
    
    out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
    
    
    
    
    
    #stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
    
    summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))
    
    abundance <- summ %>% 
      filter(grepl(pattern = "alpha[",
                   variable,
                   fixed = TRUE),
             !grepl(pattern = "gp",
                    variable,
                    fixed = TRUE)) %>% 
      mutate(routeF = row_number(),
             across(mean:q95,exp)) %>% 
      rename_with(.,
                  ~ paste0("abundance_",.x),
                  .cols = mean:ess_tail)
    
    
    slope <- summ %>% 
      filter(grepl(pattern = "beta[",
                   variable,
                   fixed = TRUE),
             !grepl(pattern = "gp",
                    variable,
                    fixed = TRUE)) %>% 
      mutate(routeF = row_number(),
             across(mean:q95, ppy )) %>% 
      rename_with(.,
                  ~ paste0("trend_",.x),
                  .cols = mean:ess_tail) %>% 
      select(-variable)
    
    both <- inner_join(slope,abundance,
                       by = "routeF") %>% 
      mutate(model = spp1)
    
    outboth <- bind_rows(outboth,both) %>% 
      mutate(species = species)
    
    both_wide <- both %>% 
      select(routeF,trend_mean,trend_sd,abundance_mean,abundance_sd) %>% 
      rename_with(.cols = -routeF,
                  .fn = ~paste(spp1,.x,sep = "_"))
    
    if(spp1 == "BYM"){
      both_wide_out <- both_wide
    }else{
    both_wide_out <- inner_join(both_wide_out,both_wide,
                                by = "routeF")
    }
    
    
  } #end models loop
  
  strata_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
  bb = st_bbox(strata_bounds)
  xdif <- as.numeric(bb$xmax - bb$xmin)*1.1
  ydif <- as.numeric(bb$ymax - bb$ymin)*1.1
  
  xlms = as.numeric(c(ifelse(bb$xmin > 0,bb$xmin*0.9,bb$xmin*1.1),
                      bb$xmin + xdif))
  ylms = as.numeric(c(ifelse(bb$ymin > 0,bb$ymin*0.9,bb$ymin*1.1),
                      bb$ymin + ydif))
  
  
  
  plot_map <- route_map %>% 
    left_join(.,outboth,
              by = "routeF") %>% 
    mutate(model = ifelse(model == "nonspatial","Non-spatial",model),
           model = factor(model,
                          levels = c("iCAR","GP","BYM","Non-spatial"),
                          ordered = TRUE),
           abundance_cv = abundance_sd/abundance_mean)
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  lgnd_head <- "Mean Trend\n"
  trend_title <- "Mean Trend"
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %/year")
  plot_map$Tplot <- cut(plot_map$trend_mean,breaks = c(-Inf, breaks, Inf),labels = labls)
  
  
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  
  
  
  map <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = Tplot,
                size = abundance_mean))+
    scale_size_continuous(range = c(0.05,2),
                          name = "Mean Count")+
    scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                        guide = guide_legend(reverse=TRUE),
                        name = paste0(lgnd_head,firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    theme(text = element_text(family = "serif",
                              size = 11),
          panel.grid = element_line(colour = grey(0.95)),
          plot.caption = element_text(hjust = 0))+
    #labs(title = paste(species, "trends"))+
    facet_wrap(vars(model))
  
  pdf(paste0("Figures/Figure_2.pdf"),
      height = 7,
      width = 7.5)
  print(map)
  dev.off()
  
 #  map_se <- ggplot()+
 #    geom_sf(data = base_strata_map,
 #            fill = NA,
 #            colour = grey(0.75))+
 #    geom_sf(data = plot_map,
 #            aes(colour = trend_sd,
 #                size = abundance_sd))+
 #    scale_size_continuous(range = c(0.05,2),
 #                          name = "SE of Mean Count",
 #                          trans = "reverse")+
 #    scale_colour_viridis_c(aesthetics = c("colour"),
 #                           guide = guide_legend(reverse=TRUE),
 #                           name = paste0("SE of Trend",firstYear,"-",lastYear))+
 #    coord_sf(xlim = xlms,ylim = ylms)+
 #    theme_bw()+
 #    #labs(title = paste(species, "Standard error"))+
 #    facet_wrap(vars(model))
 #  
 # # version with teh SE of hte maps
 #  pdf(paste0("Figures/Figure_S2.pdf"),
 #      height = 10.5,
 #      width = 7.5)
 #  print(map / map_se)
 #  dev.off()
 #  

  
  lgnd_head <- "Mean Abundance\n"
  map_abund <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = abundance_mean,
                size = abundance_sd/abundance_mean),
            alpha = 0.7)+
    scale_size_continuous(range = c(0.05,2),
                          name = "CV Abundance",
                          trans = "reverse")+
    scale_colour_viridis_c(guide = guide_colourbar(),
                        name = paste0(lgnd_head,firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    theme(text = element_text(family = "serif",
                              size = 11),
          panel.grid = element_line(colour = grey(0.95)),
          plot.caption = element_text(hjust = 0))+
    #labs(title = paste(species, "abundance"))+
    facet_wrap(vars(model))
  
  pdf(paste0("Figures/Figure_3.pdf"),
      height = 7,
      width = 7.5)
  print(map_abund)
  dev.off()
  
  
# 
#   ab_iCAR <- ggplot(data = both_wide_out)+
#     geom_point(aes(x = nonspatial_abundance_mean,
#                    y = iCAR_abundance_mean))+
#     geom_abline(slope = 1,intercept = 0)+
#     scale_y_continuous(limits = c(0.05,20),
#                        trans = "log10")+
#     scale_x_continuous(trans = "log10")+
#     theme_bw()
#   
#   ab_BYM <- ggplot(data = both_wide_out)+
#     geom_point(aes(x = nonspatial_abundance_mean,
#                    y = BYM_abundance_mean))+
#     geom_abline(slope = 1,intercept = 0)+
#     scale_y_continuous(limits = c(0.05,20),
#                        trans = "log10")+
#     scale_x_continuous(trans = "log10")+
#     theme_bw()
#   
#   ab_GP <- ggplot(data = both_wide_out)+
#     geom_point(aes(x = nonspatial_abundance_mean,
#                    y = GP_abundance_mean))+
#     geom_abline(slope = 1,intercept = 0)+
#     scale_y_continuous(limits = c(0.05,20),
#                        trans = "log10")+
#     scale_x_continuous(trans = "log10")+
#     theme_bw()
# 
# ab_all <- ab_iCAR + ab_BYM + ab_GP 
# 
# ab_all 
  
  ## showing that sdobs is lower with the increased smoothing of the iCAR
#   conv_sum <- readRDS("data/convergence_summary_May31.rds")
#   sdobs_BASP <- conv_sum %>% filter(species == "Baird's Sparrow", variable == "sdobs")
  # sdobs_BASP %>% select(variable,mean,sd,rhat,model)
  # # A tibble: 4 × 5
  # variable  mean    sd  rhat model     
  # <chr>    <dbl> <dbl> <dbl> <chr>     
  # 1 sdobs     1.09 0.130  1.01 BYM       
  # 2 sdobs     1.08 0.125  1.00 iCAR      
  # 3 sdobs     1.30 0.142  1.00 nonspatial
  # 4 sdobs     1.22 0.135  1.01 GP 
  
  
  
  

# Figure 4 ----------------------------------------------------------------

  ## CV summary for all 4 models

  diffs_plot <- readRDS("data_cv_summary_4models_plotting_data.rds")
 
  BASP <- diffs_plot %>% 
    filter(species == "Baird's Sparrow")
  
  ### specific values of CV summary included in text
  ### describing the CV results for Baird's sparrow
  BASP
  
  
  CATO <- diffs_plot %>% 
    filter(species == "Canyon Towhee")
  
  ### specific values of CV summary included in text
  ### describing the CV results for Canyon Towhee
  CATO
  WEBL <- diffs_plot %>% 
    filter(species == "Western Bluebird")
  
  ### specific values of CV summary included in text
  ### describing the CV results for Western Bluebird
  WEBL
  
  
  cv_comparisons <- diffs_plot %>% 
    select(model_comparison) %>% 
    distinct() %>% 
    mutate(model_comp_plot = str_replace(model_comparison,
                                         "_"," - "),
           model_comp_plot = str_replace(model_comp_plot,
                                         "nonspatial","Non-spatial"),
           model_comp_plot = factor(model_comp_plot,
                                    levels = rev(c("iCAR - Non-spatial",
                                               "GP - Non-spatial",
                                               "BYM - Non-spatial",
                                               "iCAR - GP",
                                               "iCAR - BYM")),
                                    ordered = TRUE))
  
  diffs_plot_sel <- diffs_plot %>% 
    filter(grepl("nonspatial",model_comparison)) %>% 
    inner_join(.,cv_comparisons,
               by = "model_comparison") %>% 
    mutate(species = factor(species,
                            levels = species_sort$english,
                            ordered = TRUE))
  
  
  capt_tmp <- paste("Figure S4. Leave Future Out (LFO) cross-validation results for 71 small-range species
                    from the BBS database, comparing three explicitly spatial models (iCAR, BYM, and GP)
                    to an otherwise identical non-spatial model. The values represent z-score summaries
                    of the difference between each of the spatial models and the non-spatial model. For all
                    of the comparisons, positive values indicate that the spatial model had higher out-of-
                    sample predictive accuracy (higher log point-wise predictive density, lppd) than the 
                    non-spatial model.")
  z_diffs_all <- ggplot(data = diffs_plot_sel,
                        aes(y = z,x = species,
                            colour = model_comp_plot))+
    geom_point()+
    geom_hline(yintercept = 0)+
    scale_colour_viridis_d(direction = -1,
                           end = 0.8)+
    ylab("Z-score difference in pointwise lpd")+
    xlab("")+
    geom_hline(yintercept = c(-2,2),
               alpha = 0.5)+
    theme_bw()+
    labs(caption = capt_tmp)+
    theme(plot.caption = element_text(hjust = 0),
          text = element_text(family = "serif",
                              size = 11))+
    coord_flip()#ylim = c(-0.5,0.5))
  
  ### species-level detail from Figure 4 comparing 3 spatial models with nonspatial
  pdf("Figures/Figure_S4.pdf",
      height = 11,
      width = 8.5)
  z_diffs_all
  dev.off()
  

  ## raincloud plot of z-score differences from nonspatial
  
  
  

  diffs_plot_sel <- diffs_plot %>% 
   # filter(grepl("nonspatial",model_comparison)| model_comparison == "iCAR_GP") %>% 
    inner_join(.,cv_comparisons,
               by = "model_comparison")
  
  
  y_diffs_hist <- ggplot(data = diffs_plot_sel,
                         aes(y = z,
                             x = model_comp_plot,
                             #colour = model_comp_plot,
                             group = model_comp_plot))+
    # ggdist::stat_halfeye(.width = 0,
    #                      #adjust = 0.3,
    #                      justification = 0,
    #                      point_color = NA,
    #                      alpha = 0.3)+
    # geom_boxplot(width = 0.12,
    #              alpha = 0.5)+
    ggdist::stat_dots(side = "right",
                      justification = 0,
                      binwidth = 0.2)+
    theme_bw()+
    ylab("Z-score difference in pointwise lppd")+
    xlab("")+
    scale_x_discrete(guide = guide_axis(position = "right"))+
    # theme(axis.text.y = element_blank(),
    #       axis.ticks.y = element_blank())+
    geom_hline(yintercept = c(-2,2), alpha = 0.4)+
    geom_hline(yintercept = 0, alpha = 1)+
    # guides(colour = guide_legend(title = "Model comparison",
    #                              reverse = TRUE),
    #        fill = guide_legend(title = "Model comparison",
    #                            reverse = TRUE))+
    # scale_colour_viridis_d(begin = 0.2,end = 0.9,
    #                        aesthetics = c("fill","colour"))+
    coord_flip()
  
  pdf("Figures/Figure_4.pdf",
      height = 4,
      width = 7.5)
  y_diffs_hist
  dev.off()
  
  
  
  

# Figure 5 ----------------------------------------------------------------

  #read in the saved cross validation summary from file "CV_summary_2_models.R"
  cv_sum <- readRDS("output/saved_cross_validation_summary_2_models.rds")
  n_obs <- readRDS("output/saved_nobs_cross_validation_summary_2_models.rds")
  
  
  cv_sum$Year <- factor(cv_sum$r_year)
  
  
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
  
  
  
  # diffs_plot_all <- mndiffs_sp %>% 
  #   # filter(grepl("nonspatial",model_comparison)| model_comparison == "iCAR_GP") %>% 
  #   inner_join(.,cv_comparisons,
  #              by = "model_comparison")
  # 
  # 
  difs_sp <- ggplot(data = mndiffs_sp,
                         aes(y = z_iCAR_nonspatial))+
    ggdist::stat_dots(side = "right",
                      justification = 0,
                      binwidth = 0.3)+
    theme_bw()+
    ylab("Z-score difference in pointwise lppd iCAR - Non-spatial")+
    xlab("")+
    scale_x_continuous(limits = c(0,1), expand = expansion(add = 0))+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())+
    geom_hline(yintercept = c(-2,2), alpha = 0.4)+
    geom_hline(yintercept = 0, alpha = 1)+
    coord_flip()
  
  pdf("Figures/Figure_5.pdf",
      height = 3,
      width = 7.5)
  difs_sp
  dev.off()
  
  
  

# Figure 6 ----------------------------------------------------------------

  ### iCAR vs nonspatial comparison maps for 4 species (2x4 maps)
  
  species_sel <- c("American Robin","Hairy Woodpecker",
                   "Common Yellowthroat","Bald Eagle")

  output_dir <- "F:/iCAR_route_2021/output"
  base_strata_map <- bbsBayes2::load_map("bbs_usgs")
  
  
  firstYear <- 2006
  lastYear <- 2021
  
  ppy <- function(x){
    p <- exp(x)-1
    return(p*100)
  }
  
  plot_map_out <- NULL
  
  for(species in species_sel){
   
    species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
    
    out_base_temp <- paste0(species_f,"_nonspatial_",firstYear,"_",lastYear)
    
    if(!file.exists(paste0(output_dir,"/",out_base_temp,"_summ_fit.rds"))){next}
    
    
    
    sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
    
    
    load(sp_data_file)
    
    
    outboth <- NULL
    for(spp1 in c("iCAR","nonspatial")){
      
      spp <- paste0("_",spp1,"_")
      
      
      out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
      
      
      
      
      
      #stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
      summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))
      
      abundance <- summ %>% 
        filter(grepl(pattern = "alpha[",
                     variable,
                     fixed = TRUE),
               !grepl(pattern = "gp",
                      variable,
                      fixed = TRUE)) %>% 
        mutate(routeF = row_number(),
               across(mean:q95,exp)) %>% 
        rename_with(.,
                    ~ paste0("abundance_",.x),
                    .cols = mean:ess_tail)
      
      
      slope <- summ %>% 
        filter(grepl(pattern = "beta[",
                     variable,
                     fixed = TRUE),
               !grepl(pattern = "gp",
                      variable,
                      fixed = TRUE)) %>% 
        mutate(routeF = row_number(),
               across(mean:q95, ppy )) %>% 
        rename_with(.,
                    ~ paste0("trend_",.x),
                    .cols = mean:ess_tail) %>% 
        select(-variable)
      
      both <- inner_join(slope,abundance,
                         by = "routeF") %>% 
        mutate(model = spp1)
      
      outboth <- bind_rows(outboth,both)
      
    } #end models loop
    

    
    
    plot_map <- route_map %>% 
      left_join(.,outboth,
                by = "routeF") %>% 
      mutate(model = ifelse(model == "nonspatial","Non-spatial",model),
             model = factor(model,
                            levels = c("iCAR","GP","BYM","Non-spatial"),
                            ordered = TRUE),
             sp = species) %>% 
      select(-c(routeF,variable))
    
    plot_map_out <- bind_rows(plot_map_out,plot_map)
    
  }#end species loop
    
    breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
    lgnd_head <- "Mean Trend\n"
    trend_title <- "Mean Trend"
    labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
    labls = paste0(labls, " %/year")
    plot_map_out$Tplot <- cut(plot_map_out$trend_mean,breaks = c(-Inf, breaks, Inf),labels = labls)
    
    
    map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                     "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
    names(map_palette) <- labls
    

    
    map <- ggplot()+
      geom_sf(data = base_strata_map,
              fill = NA,
              colour = grey(0.75))+
      geom_sf(data = plot_map_out,
              aes(colour = Tplot),
              size = 0.01,
              shape = 20)+
      # scale_size_continuous(range = c(0.05,2),
 #                           name = "Mean Count")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE,
                                               override.aes = list(size = 3)),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      guides(size = "none")+
      theme_bw()+
      theme(panel.grid = element_line(colour = grey(0.95)))+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    
    pdf(paste0("Figures/Figure_6.pdf"),
        height = 10,
        width = 7.5)
    print(map)
    dev.off()
    
    

# Figure S5 ---------------------------------------------------------------

    
    capt_tmp <- paste0("Figure S5. Map of standard error of route-level trend estimates for four broad-ranging
                       species from an iCAR spatial model and an otherwise identical non-spatial model.")
    
    map_se <- ggplot()+
      geom_sf(data = base_strata_map,
              fill = NA,
              colour = grey(0.75))+
      geom_sf(data = plot_map_out,
              aes(colour = trend_sd,
                  size = abundance_sd))+
      scale_size_continuous(range = c(0.05,1),
                            name = "SE of Mean Count")+
      scale_colour_viridis_c(aesthetics = c("colour"),
                             guide = guide_colourbar(reverse=TRUE),
                             name = paste0("SE of Trend \n",firstYear,"-",lastYear))+
      guides(size = "none")+
      theme_bw()+
      labs(caption = capt_tmp)+
      theme(plot.caption = element_text(hjust = 0),
            text = element_text(family = "serif",
                                size = 11),
            panel.grid = element_line(colour = grey(0.95)))+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    # 
    pdf(paste0("Figures/Figure_S5.pdf"),
        height = 10,
        width = 7.5)
    print(map_se)
    dev.off()
    
    
    


# Figure 7 ----------------------------------------------------------------

  
  ### add a GP vs iCAR map comparison for two extreme cv difference species
  ### Canyon Towhee - GP favoured
  ### Western Bluebird - iCAR favoured
  
  
    
    species_sel <- c("Canyon Towhee","Western Bluebird")
    
    output_dir <- "F:/iCAR_route_2021/output"
    base_strata_map <- bbsBayes2::load_map("bbs_usgs")
    
    species_latin1 <- bbsBayes2::search_species(species_sel[1])
    species_latin1 <- paste(species_latin1[1,"genus"],species_latin1[1,"species"])
    
    species_latin2 <- bbsBayes2::search_species(species_sel[2])
    species_latin2 <- paste(species_latin2[1,"genus"],species_latin2[1,"species"])
    
    
    firstYear <- 2006
    lastYear <- 2021
    
    ppy <- function(x){
      p <- exp(x)-1
      return(p*100)
    }
    
    plot_map_out <- NULL
    
    for(species in species_sel){
      
      species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
      
      out_base_temp <- paste0(species_f,"_nonspatial_",firstYear,"_",lastYear)
      
      if(!file.exists(paste0(output_dir,"/",out_base_temp,"_summ_fit.rds"))){next}
      
      
      
      sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")
      
      
      load(sp_data_file)
      
      
      outboth <- NULL
      for(spp1 in c("iCAR","GP")){
        
        spp <- paste0("_",spp1,"_")
        
        
        out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
        
        
        
        
        
        #stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.rds"))
        summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))
        
        abundance <- summ %>% 
          filter(grepl(pattern = "alpha[",
                       variable,
                       fixed = TRUE),
                 !grepl(pattern = "gp",
                        variable,
                        fixed = TRUE)) %>% 
          mutate(routeF = row_number(),
                 across(mean:q95,exp)) %>% 
          rename_with(.,
                      ~ paste0("abundance_",.x),
                      .cols = mean:ess_tail)
        
        
        slope <- summ %>% 
          filter(grepl(pattern = "beta[",
                       variable,
                       fixed = TRUE),
                 !grepl(pattern = "gp",
                        variable,
                        fixed = TRUE)) %>% 
          mutate(routeF = row_number(),
                 across(mean:q95, ppy )) %>% 
          rename_with(.,
                      ~ paste0("trend_",.x),
                      .cols = mean:ess_tail) %>% 
          select(-variable)
        
        both <- inner_join(slope,abundance,
                           by = "routeF") %>% 
          mutate(model = spp1)
        
        outboth <- bind_rows(outboth,both)
        
      } #end models loop
      
      plot_map <- route_map %>% 
        left_join(.,outboth,
                  by = "routeF") %>% 
        mutate(model = ifelse(model == "nonspatial","Non-spatial",model),
               model = factor(model,
                              levels = c("iCAR","GP","BYM","Non-spatial"),
                              ordered = TRUE),
               sp = species) %>% 
        select(-c(routeF,variable))
      
      plot_map_out <- bind_rows(plot_map_out,plot_map)
      
      
    }#end species loop
    
    
   
    
      strata_bounds <- st_union(plot_map_out) #union to provide a simple border of the realised strata
      bb = st_bbox(strata_bounds)
      xdif <- as.numeric(bb$xmax - bb$xmin)*1.1
      ydif <- as.numeric(bb$ymax - bb$ymin)*1.1
      
      xlms = as.numeric(c(ifelse(bb$xmin > 0,bb$xmin*0.9,bb$xmin*1.1),
                          bb$xmin + xdif))
      ylms = as.numeric(c(ifelse(bb$ymin > 0,bb$ymin*0.9,bb$ymin*1.1),
                          bb$ymin + ydif))
      
      
  
    
    
    breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
    lgnd_head <- "Mean Trend\n"
    trend_title <- "Mean Trend"
    labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
    labls = paste0(labls, " %/year")
    plot_map_out$Tplot <- cut(plot_map_out$trend_mean,breaks = c(-Inf, breaks, Inf),labels = labls)
    
    
    map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                     "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
    names(map_palette) <- labls
    
    
    
    map <- ggplot()+
      geom_sf(data = base_strata_map,
              fill = NA,
              colour = grey(0.75))+
      geom_sf(data = plot_map_out,
              aes(colour = Tplot,
                  size = abundance_mean))+
      scale_size_continuous(range = c(0.05,2),
                            name = "Mean Count")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms,
               expand = TRUE)+
      #guides(size = "none")+
      theme_bw()+
      theme(panel.grid = element_line(colour = grey(0.95)))+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    
    pdf(paste0("Figures/Figure_7.pdf"),
        height = 6,
        width = 7.5)
    print(map)
    dev.off()
    
    
    

# Figure S6 ---------------------------------------------------------------

    
    capt_tmp <- paste0("Figure S6. Map of standard error of route-level trend estimates for two species from two spatial models.
                       Interestingly, the standard errors of the GP model's estimates are smaller than those of the iCAR model
                       for both species. However, this higher estimated precision does not reflect higher accuracy
                       because the out-of-sample predictive accuracy suggests that the best model varies between these two species.
                       For Canyon Towhee (",species_latin1,") the GP model has higher accuracy and 
                       for Western Bluebird (",species_latin2,") the iCAR model has higher accuracy.")
    

    
    map_se <- ggplot()+
      geom_sf(data = base_strata_map,
              fill = NA,
              colour = grey(0.75))+
      geom_sf(data = plot_map_out,
              aes(colour = trend_sd,
                  size = abundance_sd))+
      scale_size_continuous(range = c(0.05,3),
                            name = "SE of Mean Count")+
      scale_colour_viridis_c(aesthetics = c("colour"),
                             guide = guide_colorbar(reverse=TRUE),
                             name = paste0("SE of Trend \n",firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      guides(size = "none")+
      theme_bw()+
      labs(caption = capt_tmp)+
      theme(plot.caption = element_text(hjust = 0),
            text = element_text(family = "serif",
                                size = 11),
            panel.grid = element_line(colour = grey(0.95)))+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    

    # 
    pdf(paste0("Figures/Figure_S6.pdf"),
        height = 6,
        width = 7.5)
    print(map_se)
    dev.off()
    
    

   
    
  

# Species tables ----------------------------------------------------------

    
    
    sps_all <- bbsBayes2::load_bbs_data()$species
    
    species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
    
    species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")
    
    sp_s <- data.frame(Species = species_list,
                       Models_fit = c("Non-spatial, iCAR, BYM, and GP"))
    sp_l <- data.frame(Species = species_list_broad,
                       Models_fit = c("Non-spatial, iCAR"))
    
    # sp_all <- bind_rows(sp_s,sp_l) %>% 
    #   rowwise() %>% 
    #   mutate(latin_name = latin_function(species),
    #          aou = aou_function(species))
      
    sp_all <- bind_rows(sp_s,sp_l) %>% 
      left_join(.,sps_all,
                by = c("Species" = "english")) %>% 
      mutate(latin_name = paste(genus,species),
             Models_fit = factor(Models_fit,
                                  levels = c("Non-spatial, iCAR, BYM, and GP",
                                             "Non-spatial, iCAR"),ordered = TRUE)) %>% 
      select(Species,latin_name,seq,Models_fit) %>% 
      distinct() %>% 
      arrange(Models_fit,seq) %>% 
      select(-seq)

        
   write.csv(sp_all,"figures/species_list.csv")
  
    
   

# Figure S8 supplement ----------------------------------------------------

load("data/cv_summary_4models_data.RData") #four model cross-validation results
   
   lpos = function(x){
     p = length(which(x > 0))/length(x)
   }
   

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
                             levels = levels(species_sort$english)),
            .groups = "keep") 

   
   mndiffs_sel_plot <- mndiffs %>% 
     filter(dist_cat == "Isolated")
   mn_diffs_plot <- mndiffs %>% 
     filter(species %in% mndiffs_sel_plot$species)
   
   
   capt_tmp <- paste("Figure S8. Mean point-wise differences in log posterior predictive density (lppd, iCAR-GP)
                     between iCAR and GP spatial models for BBS routes that are Isolated from other routes
                     (i.e., greater than 200 km from the nearest neighbour, in dark purple) and other 
                     routes that are more central to the species range and therefor have closer (and usually
                     more) neighbors. For most of the species here, the iCAR based simplification of spatial 
                     relationships has higher predictive accuracy for the trends on Isolated routes than the 
                     GP model (dark points to the right of the vertal line at 0) that uses precise distance 
                     information and therefore necessarily treats these isolated routes as having lower covariance
                     in trends and abundances than the iCAR model that treats them as immediate neighbours of 
                     the nearest routes that are 200 km away.")
   
   
   y_diffs_gp <- ggplot(data = mn_diffs_plot,
                        aes(y = m_iCAR_GP,x = species,
                            colour = dist_cat))+
     geom_point()+
     scale_colour_viridis_d(direction = -1,
                            begin = 0.1,end = 0.8)+
     geom_hline(yintercept = 0)+
     theme_bw()+
     xlab("")+
     ylab("Mean difference (iCAR-GP) \n positive values support iCAR over GP")+
     guides(guide_legend(title = "Distance to \n nearest neighbor"))+
     labs(caption = capt_tmp)+
     theme(plot.caption = element_text(hjust = 0),
           text = element_text(family = "serif",
                               size = 11),
           panel.grid = element_line(colour = grey(0.95)))+
     coord_flip()
   
   pdf("Figures/Figure_S8.pdf",
       height = 11,
       width = 8.5)
   y_diffs_gp
   dev.off()
   
   
   
  