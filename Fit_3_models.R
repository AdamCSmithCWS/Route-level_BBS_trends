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

skip_bym <- TRUE
# SPECIES LOOP ------------------------------------------------------------

# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in rev(species_list_broad)){
#species <- species_list[2]

  
species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

for(spp1 in c("BYM","iCAR","nonspatial")){

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

# 
# 
# # extract trends and abundances -------------------------------------------
# if(produce_trends){
# routes_df <- data.frame(routeF = new_data$routeF,
#                         route = new_data$route,
#                         year = new_data$r_year,
#                         obs = new_data$ObsN,
#                         stratum = new_data$strat_name,
#                         Latitude = new_data$Latitude,
#                         Longitude = new_data$Longitude) # %>% 
#   group_by(route,routeF,obs,stratum,Longitude,Latitude) %>%
#   summarise(n_yr_obs = n(),
#             .groups = "drop") %>% 
#   group_by(route,routeF,stratum,Longitude,Latitude) %>% 
#   summarise(n_obs = n(),
#             mean_y_obs = mean(n_yr_obs),
#             max_y_obs = max(n_yr_obs),
#             .groups = "drop") %>% 
#   distinct()
# 
# # function to turn a log-scale slope parameter into a trend estimate in %/year
# tr_f <- function(x){
#   t <- (exp(x)-1)*100
# }
# 
# 
# # route-level trends and abundances ------------------------------------------
# 
# 
# abunds <- posterior_samples(fit = stanfit,
#                                   parm = "alpha",
#                                   dims = "routeF") %>% 
#   posterior_sums(.,
#                  dims = "routeF")%>% 
#   left_join(.,routes_df,by = "routeF") %>% 
#   mutate(abund = exp(mean),
#          abund_lci = exp(lci),
#          abund_uci = tr_f(uci),
#          abund_se = tr_f(sd),
#          abund_Wci = abund_uci-abund_lci)%>% 
#   select(route,abund,abund_lci,abund_uci,abund_se,abund_Wci)
# 
# trends <- posterior_samples(fit = stanfit,
#                              parm = "beta",
#                              dims = "routeF") %>%
#   posterior_sums(.,
#                  dims = "routeF")%>% 
#   left_join(.,routes_df,by = "routeF") %>% 
#   mutate(trend = tr_f(mean),
#          trend_lci = tr_f(lci),
#          trend_uci = tr_f(uci),
#          trend_se = tr_f(sd),
#          trend_Wci = trend_uci-trend_lci,
#          log_scale_slope = mean,
#          log_scale_slope_lci = lci,
#          log_scale_slope_uci = uci,
#          log_scale_slope_se = sd,
#          log_scale_slope_Wci = log_scale_slope_uci-log_scale_slope_lci) %>% 
#   select(route,stratum,Latitude,Longitude,trend,trend_lci,trend_uci,trend_se,trend_Wci,
#          log_scale_slope,log_scale_slope_lci,log_scale_slope_uci,log_scale_slope_se,log_scale_slope_Wci)
# 
# est_table <- inner_join(trends,
#                         abunds,
#                         by = "route") %>% 
#   mutate(species = species,
#          trend_time = paste(firstYear,lastYear,sep = "-"))
# 
# 
# #writes a csv file for each species
# write.csv(est_table,
#           paste0("trends/","trends_",out_base,".csv"),
#           row.names = FALSE)
# #also appends species results to an all trends file
# # if(file.exists(paste0("trends/","All_trends_",firstYear,"-",lastYear,".csv"))){
# # write.table(est_table,
# #           paste0("trends/","All_trends_",firstYear,"-",lastYear,".csv"),
# #           sep = ",",
# #           row.names = FALSE,
# #           append = TRUE,
# #           col.names = FALSE)
# #   }else{
# #             write.table(est_table,
# #                         paste0("trends/","All_trends_",firstYear,"-",lastYear,".csv"),
# #                         sep = ",",
# #                         row.names = FALSE,
# #                         append = FALSE,
# #                         col.names = TRUE)  
# #           }
# # if it makes sense to output the pdf maps of species trends
# if(produce_maps){
#   
#   other_routes <- strat_data$route_strat %>% 
#     select(strat_name,Latitude,Longitude,rt.uni,Year) %>%
#     filter(Year %in% c(firstYear:lastYear)) %>% 
#     distinct() %>% 
#     rename(route = rt.uni) %>% 
#     filter(!route %in% route_map$route) %>% 
#     st_as_sf(.,coords = c("Longitude","Latitude"))
#   
#   st_crs(other_routes) <- 4269 #NAD83 commonly used by US federal agencies
#   
#   
#   
# route_bounds <- st_union(route_map) #union to provide a simple border of the realised strata
# bb = st_bbox(route_bounds)
# xlms = as.numeric(c(bb$xmin,bb$xmax))
# ylms = as.numeric(c(bb$ymin,bb$ymax))
# 
# 
# # merging route map with the estimates for this species -------------------
# trend_plot_map <- route_map %>% 
#   left_join(.,est_table,by = "route") 
# 
# 
# breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
# lgnd_head <- "Mean Trend\n"
# trend_title <- "Mean Trend"
# labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
# labls = paste0(labls, " %/year")
# trend_plot_map$Tplot <- cut(trend_plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)
# 
# 
# map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
#                  "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
# names(map_palette) <- labls
# 
# 
# tmap = ggplot(trend_plot_map)+
#   #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(data = other_routes,colour = grey(0.7),inherit.aes = FALSE,size = 0.15,
#           show.legend = FALSE)+
#   geom_sf(aes(colour = Tplot,size = abund))+
#   scale_size_continuous(range = c(0.4,2),
#                         name = "Mean Count")+
#   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
#                       guide = guide_legend(reverse=TRUE),
#                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
#   theme_minimal()+
#   coord_sf(xlim = xlms,ylim = ylms)
# 
# 
# pdf(file = paste0("Figures/",out_base,"_trend_map.pdf"),
#     width = 10,
#     height = 8)
# print(tmap)
# dev.off()
# 
# }#end if produce_maps 
# 
# }# end if produce_trends
# 
# }# end of species loop


