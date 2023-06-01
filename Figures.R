#Figures

library(tidyverse)
library(patchwork)
library(sf)

# Figure_1 ----------------------------------------------------------------

# demonstration of neighbours for BHVI

# load saved neighbours data from neighbours_define_voronoi(..., save_plot_data = TRUE)

load("data/Baird's_Sparrow_route_maps_data.RData")

box <- st_as_sfc(st_bbox(strata_map))

xb <- range(st_coordinates(box)[,"X"])
yb <- range(st_coordinates(box)[,"Y"])


ggp1 <- ggplot()+ 
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
  theme(legend.position = "none")
ggp1

ggp2 <- ggplot()+ 
  geom_sf(data = strata_map,alpha = 0,colour = grey(0.8))+ 
  geom_segment(data=DA,aes(x = long, y = lat,xend=long_to,yend=lat_to),
               inherit.aes = FALSE,linewidth=0.3,colour = grey(0.6)) +
  geom_sf(data = centres,alpha = 1,size = 1.3)+
  theme_minimal() +
  xlab("")+
  ylab("")+
  coord_sf(xlim = xb,ylim = yb)+
  theme(legend.position = "none")
ggp2


pdf("Figures/Figure_1.pdf",
    width = 7.5,
    height = 10)
print(ggp1 + ggp2 + plot_layout(nrow = 2))
dev.off()






# Figure 2 - 4 model comparison -------------------------------------------
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
  xlms = as.numeric(c(bb$xmin,bb$xmax))*1.05
  ylms = as.numeric(c(bb$ymin,bb$ymax))*1.05
  
  
  plot_map <- route_map %>% 
    left_join(.,outboth,
              by = "routeF") %>% 
    mutate(model = factor(model,
                          levels = c("iCAR","GP","BYM","nonspatial"),
                          ordered = TRUE))
  
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
    #labs(title = paste(species, "trends"))+
    facet_wrap(vars(model))
  
  pdf(paste0("Figures/Figure_2.pdf"),
      height = 7,
      width = 7.5)
  print(map)
  dev.off()
  
  map_se <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = trend_sd,
                size = abundance_sd))+
    scale_size_continuous(range = c(0.05,2),
                          name = "SE of Mean Count",
                          trans = "reverse")+
    scale_colour_viridis_c(aesthetics = c("colour"),
                           guide = guide_legend(reverse=TRUE),
                           name = paste0("SE of Trend",firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
    #labs(title = paste(species, "Standard error"))+
    facet_wrap(vars(model))
  
 # version with teh SE of hte maps
  pdf(paste0("Figures/Figure_S2.pdf"),
      height = 10.5,
      width = 7.5)
  print(map / map_se)
  dev.off()
  

  
  lgnd_head <- "Mean Abundance\n"
  map_abund <- ggplot()+
    geom_sf(data = base_strata_map,
            fill = NA,
            colour = grey(0.75))+
    geom_sf(data = plot_map,
            aes(colour = abundance_mean,
                size = abundance_sd/abundance_mean))+
    scale_size_continuous(range = c(0.05,2),
                          name = "CV Abundance",
                          trans = "reverse")+
    scale_colour_viridis_c(guide = guide_legend(reverse=TRUE),
                        name = paste0(lgnd_head,firstYear,"-",lastYear))+
    coord_sf(xlim = xlms,ylim = ylms)+
    theme_bw()+
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
  
  pdf("Figures/Figure_S4.pdf",
      height = 10,
      width = 7.5)
  z_diffs_all
  dev.off()
  

  ## raincloud plot of z-score differences from nonspatial
  
  
  # Raincloud plot ----------------------------------------------------------
  
  

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

  
  
  
  
  
  
  
  
  
  
  
  
    
  