
#setwd("C:/GitHub/iCAR_route_2021")
#setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(cmdstanr)
library(sf)
library(patchwork)


#output_dir <- "output"
output_dir <- "e:/iCAR_route_2021/output"

crs_use <- readRDS("functions/custom_crs_for_maps.rds")

base_strata_map <- bbsBayes2::load_map("bbs_usgs")%>% 
  st_transform(.,crs_use)

state_prov <- bbsBayes2::load_map("prov_state") %>% 
  st_transform(.,crs_use)




species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
#species_list <- c(species_list,"Blue-headed Vireo")

species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")


firstYear <- 2006
lastYear <- 2021

ppy <- function(x){
  p <- exp(x)-1
  return(p*100)
}
# SPECIES LOOP ------------------------------------------------------------

#outboth_save <- NULL
re_summarize <- FALSE # set to true if running for the first time
estimates_compile <- NULL


pdf(paste0("Figures/Figure_S2.pdf"),
    height = 11,
    width = 8)
# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list){
#species <- species_list[2]
  species_latin <- bbsBayes2::search_species(species)
  species_latin <- paste(species_latin[1,"genus"],species_latin[1,"species"])
  
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  out_base_temp <- paste0(species_f,"_GP_",firstYear,"_",lastYear)
  
  if(!file.exists(paste0(output_dir,"/",out_base_temp,"_summ_fit.rds"))){
    print(paste("missing fitted model for",species))
    next}
    
    

sp_data_file <- paste0("Data/",species_f,"_",firstYear,"_",lastYear,"_stan_data.RData")


load(sp_data_file)

if(re_summarize){
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

saveRDS(outboth,
        paste0("output/predictions_4_models_",species_f,".RDS"))
# outboth_save <- bind_rows(outboth_save,outboth)


} #end models loop
}else{#end if re_summarize

  outboth <- readRDS(paste0("output/predictions_4_models_",species_f,".RDS"))
  
}
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
         abundance_cv = (abundance_q95-abundance_q5)/(abundance_median*4))

plot_map_df <- plot_map %>%
  sf::st_transform(.,st_crs(4269)) %>% 
  dplyr::mutate(longitude = sf::st_coordinates(.)[,1],
                latitude = sf::st_coordinates(.)[,2],
                species = species) %>% 
  sf::st_drop_geometry()
plot_map_df$species <- species

estimates_compile <- bind_rows(estimates_compile,plot_map_df)

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
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map,
          aes(colour = Tplot,
              size = abundance_mean))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  labs(title = paste0(species," (",species_latin,")"),
       subtitle = paste("Trends"))+
  theme_bw()+
  theme(text = element_text(family = "serif",
                            size = 11,
                            hjust = 0.5),
        panel.grid = element_line(colour = grey(0.95)))+
  coord_sf(xlim = xlms,ylim = ylms)+
  facet_wrap(vars(model))
  
if(species == species_list[1]){
sp_caption <- paste0("Figure S2. Comparison of the predictions for ",species," (",species_latin,") \n from three spatially explicit models (iCAR, BYM, GP) and a non-spatial model of route-level \n abundance (size of dots) and trends (colours) on individual survey routes from the \n North American Breeding Bird Survey. Each point represents the starting location (first 3-minute point count) \n of the 50-point count 40 km long roadside survey transect")
}else{
  sp_caption <- paste0("Figure S2 (continued). Comparison of the predictions for ",species," (",species_latin,") \n from three spatially explicit models (iCAR, BYM, GP) and a non-spatial model of route-level \n abundance (size of dots) and trends (colours) on individual survey routes from the \n North American Breeding Bird Survey. Each point represents the starting location (first 3-minute point count) \n of the 50-point count 40 km long roadside survey transect")
}
map_se <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map,
          aes(colour = trend_sd,
              size = abundance_cv))+
  scale_size_continuous(range = c(0.05,2),
                        name = "CV of Mean Count")+
  scale_colour_viridis_c(aesthetics = c("colour"),
                         guide = guide_colorbar(reverse=TRUE),
                         name = paste0("SE of Trend ",firstYear,"-",lastYear))+
  labs(subtitle = paste("Standard error"))+
  labs(caption = sp_caption) +
  theme_bw()+
  theme(text = element_text(family = "serif",
                            size = 11),
        panel.grid = element_line(colour = grey(0.95)),
        plot.caption = element_text(hjust = 0))+
  coord_sf(xlim = xlms,ylim = ylms)+
  facet_wrap(vars(model))

fullmap <- map / map_se
print(fullmap)





}#end species loop

dev.off()


estimates_out <- estimates_compile %>% 
  select(species,model,route,latitude,longitude,
         trend_median,trend_sd,trend_q5,trend_q95,
         abundance_median,abundance_q5,abundance_q95,abundance_cv)

write_csv(estimates_out,
          "data_open/All_estimates_from_species_w_four_models.csv")



# exporting the icar trends for all species and routes --------------------


icar2m <- read_csv("data_open/All_estimates_from_species_w_two_models.csv") %>% 
  filter(model == "iCAR")

icar4m <- read_csv("data_open/All_estimates_from_species_w_four_models.csv") %>% 
  filter(model == "iCAR")


bbs_all_data <- bbsBayes2::load_bbs_data(release = 2022)
bbs_routes_strata <- bbs_all_data$routes %>%
  select(country,route_name,state_num,route,
         bcr,st_abrev) %>% 
  distinct() %>% 
  mutate(route_num = route,
         route = paste0(state_num,"-",route_num))
  


icar_out <- bind_rows(icar2m,icar4m) %>% 
  left_join(.,bbs_routes_strata,
            by = "route") %>% 
  select(-c(abundance_median:abundance_cv))

write_csv(icar_out,"data_open/all_route_level_bbs_trends.csv")




