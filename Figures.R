#Figures

library(tidyverse)
library(patchwork)
library(sf)

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

species <- "Baird's Sparrow" 


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
    ylab("Covariance between pairs of routes")+
    xlab("Distance between routes (km)")
    
  pdf("Figures/Figure_S1.pdf",
      width = 4,
      height = 3)
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
               by = "model_comparison")
  
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
    coord_flip()#ylim = c(-0.5,0.5))
  
  ### species-level detail from Figure 4 comparing 3 spatial models with nonspatial
  pdf("Figures/Figure_S4.pdf",
      height = 10,
      width = 7.5)
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
    ylab("Z-score difference in pointwise lpd")+
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
  
  
  
  diffs_plot_all <- mndiffs_sp %>% 
    # filter(grepl("nonspatial",model_comparison)| model_comparison == "iCAR_GP") %>% 
    inner_join(.,cv_comparisons,
               by = "model_comparison")
  
  
  difs_sp <- ggplot(data = mndiffs_sp,
                         aes(y = z_iCAR_nonspatial))+
    ggdist::stat_dots(side = "right",
                      justification = 0,
                      binwidth = 0.3)+
    theme_bw()+
    ylab("Z-score difference in pointwise lpd iCAR - Non-spatial")+
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
              aes(colour = Tplot,
                  size = abundance_mean))+
      scale_size_continuous(range = c(0.05,2),
                            name = "Mean Count")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      guides(size = "none")+
      theme_bw()+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    
    pdf(paste0("Figures/Figure_6.pdf"),
        height = 10,
        width = 7.5)
    print(map)
    dev.off()
    
    map_se <- ggplot()+
      geom_sf(data = base_strata_map,
              fill = NA,
              colour = grey(0.75))+
      geom_sf(data = plot_map_out,
              aes(colour = trend_sd,
                  size = abundance_sd))+
      scale_size_continuous(range = c(0.05,2),
                            name = "SE of Mean Count")+
      scale_colour_viridis_c(aesthetics = c("colour"),
                             guide = guide_legend(reverse=TRUE),
                             name = paste0("SE of Trend",firstYear,"-",lastYear))+
      guides(size = "none")+
      theme_bw()+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    # 
    pdf(paste0("Figures/Figure_S6.pdf"),
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
      scale_size_continuous(range = c(0.05,3),
                            name = "Mean Count")+
      scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                          guide = guide_legend(reverse=TRUE),
                          name = paste0(lgnd_head,firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms,
               expand = TRUE)+
      guides(size = "none")+
      theme_bw()+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    
    pdf(paste0("Figures/Figure_7.pdf"),
        height = 6,
        width = 7.5)
    print(map)
    dev.off()
    
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
                             guide = guide_legend(reverse=TRUE),
                             name = paste0("SE of Trend",firstYear,"-",lastYear))+
      coord_sf(xlim = xlms,ylim = ylms)+
      guides(size = "none")+
      theme_bw()+
      facet_grid(cols = vars(model),
                 rows = vars(sp))
    # 
    pdf(paste0("Figures/Figure_S7.pdf"),
        height = 6,
        width = 7.5)
    print(map_se)
    dev.off()
    
    

   
    
  
  
  
    
  