
library(bbsBayes2)
library(tidyverse)
library(sf)
library(cmdstanr)
library(patchwork)

output_dir <- "output"
species <- "Rufous Hummingbird" 

crs_use <- readRDS("functions/custom_crs_for_maps.rds")

base_strata_map <- bbsBayes2::load_map("bbs_usgs")%>% 
  st_transform(.,crs_use)

state_prov <- bbsBayes2::load_map("prov_state") %>% 
  st_transform(.,crs_use)



species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)

spp <- "_habitat_"

exp_t <- function(x){
  y <- (exp(x)-1)*100
}

firstYear <- 2006
lastYear <- 2021


out_base <- paste0(species_f,spp,firstYear,"_",lastYear)




sp_data_file <- paste0("Data_open/",species_f,"_",firstYear,"_",lastYear,"_covariate_stan_data.RData")


load(sp_data_file)

mod.file = paste0("models/slope",spp,"route_NB.stan")

# trend habitat effects are not changed, but the intercept effect is
# removes the optional spatial components for intercepts 
stan_data[["fit_spatial"]] <- 0 # this sets an option in the model
# to estimate the residual intercept component using a simple random
# effect, instead of a spatial one. This allows the model to estimate
# variation in abundance that is not predicted by local habitat suitability
# but does not fit an inherently spatial residual structure
# setting this fit_spatial value to 1 uses the iCAR structure to model
# a spatially explicit residual term



slope_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))



stanfit <- slope_model$sample(
  data=stan_data,
  refresh=400,
  iter_sampling=2000,
  iter_warmup=2000,
  parallel_chains = 4)

summ <- stanfit$summary()
print(paste(species, stanfit$time()[["total"]]))

saveRDS(stanfit,
        paste0(output_dir,"/",out_base,"_stanfit.rds"))

saveRDS(summ,
        paste0(output_dir,"/",out_base,"_summ_fit.rds"))

summ %>% arrange(-rhat)





# graphing ----------------------------------------------------------------


firstYear <- 2006
lastYear <- 2021

out_base <- paste0(species_f,spp,firstYear,"_",lastYear)

sp_data_file <- paste0("data_open/",species_f,"_",firstYear,"_",lastYear,"_covariate_stan_data.RData")

load(sp_data_file)

summ <- readRDS(paste0(output_dir,"/",out_base,"_summ_fit.rds"))

hab_slope <- stan_data$route_habitat_slope

hab_mean <- stan_data$route_habitat

hab_data <- data.frame(routeF = 1:length(hab_slope),
                       habitat_slope = hab_slope,
                       habitat_mean = hab_mean)


# mn0 <- new_data %>% 
#   group_by(routeF) %>% 
#   summarise(mn = mean(count),
#             mx = max(count),
#             ny = n(),
#             fy = min(year),
#             ly = max(year),
#             sp = max(year)-min(year))

route_map_2006 <- route_map 

exp_t <- function(x){
  y <- (exp(x)-1)*100
}


# plot trends -------------------------------------------------------------




strata_bounds <- st_buffer(st_union(route_map),
                           dist = 20000)#union to provide a simple border of the realised strata
bb = st_bbox(strata_bounds)
xlms = as.numeric(c(bb$xmin,bb$xmax))
ylms = as.numeric(c(bb$ymin,bb$ymax))

betas1 <- summ %>% 
  filter(grepl("beta[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Full with Habitat-Change") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

alpha1 <- summ %>% 
  filter(grepl("alpha[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Full with Habitat") %>% 
  select(routeF,median,sd) %>% 
  rename(abundance = median,
         abundance_se = sd)
betas1 <- betas1 %>% 
  inner_join(.,alpha1)

alpha2 <- summ %>% 
  filter(grepl("alpha_resid[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Residual") %>% 
  select(routeF,median,sd) %>% 
  rename(abundance = median,
         abundance_se = sd)


betas2 <- summ %>% 
  filter(grepl("beta_resid[",variable,fixed = TRUE)) %>% 
  mutate(across(2:7,~exp_t(.x)),
         routeF = as.integer(str_extract(variable,"[[:digit:]]{1,}")),
         parameter = "Residual") %>% 
  select(routeF,mean,sd,parameter) %>% 
  rename(trend = mean,
         trend_se = sd)

betas2 <- betas2 %>% 
  inner_join(.,alpha2,by = "routeF") 


betas <- bind_rows(betas1,betas2)

plot_map <- route_map_2006 %>% 
  left_join(.,betas,
            by = "routeF",
            multiple = "all") 

breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
lgnd_head <- "Mean Trend\n"
trend_title <- "Mean Trend"
labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
labls = paste0(labls, " %/year")
plot_map$Tplot <- cut(plot_map$trend,breaks = c(-Inf, breaks, Inf),labels = labls)


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
              size = abundance))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  theme(panel.grid = element_line(colour = grey(0.95)))+
  facet_wrap(vars(parameter))


map_abund <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map,
          aes(colour = abundance))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_legend(reverse=TRUE),
                         name = paste0("Relative Abundance"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  theme(panel.grid = element_line(colour = grey(0.95)))+
  labs(title = "Relative abundance")+
  facet_wrap(vars(parameter))

plot_hab_map <- route_map_2006 %>% 
  left_join(.,hab_data,
            by = "routeF",
            multiple = "all") 

capt_tmp <- paste0("Figure S7. Map of route-level habitat covariates for Rufous Hummingbird from 2006-2021.
                   The left plot shows the relative distribution of mean annual habitat amount. The right plot
                   shows the distribution of the changes in habitat between 2006-2021. These maps demonstrate
                   the general east-west pattern in both habitat amount and habitat change, where habitat has 
                   decreased in western portion of the species' range and increased in the east.")
map_hab <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_hab_map,
          aes(colour = habitat_mean))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9,
                         guide = guide_colourbar(reverse=FALSE),
                         name = paste0("Relative Habitat Amount"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  theme(text = element_text(family = "serif",size = 11),
        panel.grid = element_line(colour = grey(0.95)))+
  labs(title = "Relative habitat amount")


map_hab_slope <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_hab_map,
          aes(colour = habitat_slope))+
  colorspace::scale_colour_continuous_diverging(name = paste0("Change in Habitat"),
                                                rev = TRUE,
                                                palette = "Blue-Red 2")+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  theme(text = element_text(family = "serif",
                            size = 11),
        panel.grid = element_line(colour = grey(0.95)))+
  labs(title = "Observed change in habitat")

map_hab_slope

#map

# pdf(paste0("Figures/Four_trends_model_comparison_",species_f,".pdf"),
#     height = 8,
#     width = 8)
# print(map)
# dev.off()

map_se <- ggplot()+
  geom_sf(data = base_strata_map,
          fill = NA,
          colour = grey(0.75))+
  geom_sf(data = state_prov,
          fill = NA,
          colour = grey(0.5))+
  geom_sf(data = plot_map,
          aes(colour = trend_se,
              size = abundance_se))+
  scale_size_continuous(range = c(0.05,2),
                        name = "SE of Mean Count",
                        trans = "reverse")+
  scale_colour_viridis_c(aesthetics = c("colour"),
                         guide = guide_colourbar(reverse=TRUE),
                         name = paste0("SE of Trend"))+
  coord_sf(xlim = xlms,ylim = ylms)+
  theme_bw()+
  theme(panel.grid = element_line(colour = grey(0.95)))+
  guides(size = "none")+
  labs(title  = "Standard error")+
  facet_wrap(vars(parameter))




#print(map2 / map_se2)

# pdf(paste0("Figures/Figure_supplement_1_Trend_map_w_habitat_and_withing_",species_f,".pdf"),
#     height = 10.5,
#     width = 7.5)
# 
# 
# print(map / map_se + plot_layout(guides = "collect"))
# 
# 
# dev.off()

pdf(paste0("Figures/Figure_11.pdf"),
    height = 5,
    width = 7)


print(map)


dev.off()

map <- map+
  labs(title = "Trend")
map_save <- map / map_abund + plot_layout(guides = "collect")
saveRDS(map_save,paste0("Figures/saved_map_RUHU_covariate_",firstYear,".rds"))



maphaball <- map_hab + map_hab_slope + plot_layout(guides = "collect") +
  plot_annotation(caption = capt_tmp,
                  theme = theme(plot.caption = element_text(hjust = 0)))


pdf(paste0("Figures/Figure_S7.pdf"),
    height = 10,
    width = 7)
maphaball


dev.off()

  
#   



  