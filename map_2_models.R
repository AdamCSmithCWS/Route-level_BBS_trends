## Fitting the BYM model to 1995 - 2021 BBS data
## script currently written to fit the model then save the Stan output to a directory
## 
#setwd("C:/GitHub/iCAR_route_2021")
setwd("C:/Users/SmithAC/Documents/GitHub/iCAR_route_2021")
library(tidyverse)
library(cmdstanr)
library(sf)
library(patchwork)


output_dir <- "output"
output_dir <- "F:/iCAR_route_2021/output"

## this list should include all of the species that we're interested in for the grasslands project
species_list <- readRDS("data/species_to_include_4_model_comparison.rds")
species_list_broad <- readRDS("data/species_to_include_2_model_comparison.rds")


firstYear <- 2006
lastYear <- 2021

ppy <- function(x){
  p <- exp(x)-1
  return(p*100)
}
# SPECIES LOOP ------------------------------------------------------------

pdf(paste0("Figures/Two_model_comparison_all_species.pdf"),
    height = 10.5,
    width = 7.5)

# I've got this running as a species loop with a time-spans loop nested within
# it would be more efficient to run it in parallel using the foreach and parallel packages, but I can't seem to get Stan to work using these parallel options
for(species in species_list_broad){
#species <- species_list[2]

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


# tmap = ggplot(trend_plot_map)+
#   #geom_sf(data = realized_strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(data = strata_map,colour = gray(0.8),fill = NA)+
#   geom_sf(aes(colour = Tplot,size = abund))+
#   scale_size_continuous(range = c(0.05,2),
#                         name = "Mean Count")+
#   scale_colour_manual(values = map_palette, aesthetics = c("colour"),
#                       guide = guide_legend(reverse=TRUE),
#                       name = paste0(lgnd_head,firstYear,"-",lastYear))+
#   coord_sf(xlim = xlms,ylim = ylms)+
#   facet_wrap(vars(version),nrow = 3)

map <- ggplot()+
  geom_sf(data = plot_map,
          aes(colour = Tplot,
              size = abundance_mean))+
  scale_size_continuous(range = c(0.05,2),
                        name = "Mean Count")+
  scale_colour_manual(values = map_palette, aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0(lgnd_head,firstYear,"-",lastYear))+
  labs(title = paste(species, "trends"))+
  facet_wrap(vars(model))
  
# pdf(paste0("Figures/Two_model_comparison_",species_f,".pdf"),
#     height = 5,
#     width = 8)
# print(map)
# dev.off()

map_se <- ggplot()+
  geom_sf(data = plot_map,
          aes(colour = trend_sd,
              size = abundance_sd))+
  scale_size_continuous(range = c(0.05,2),
                        name = "SE of Mean Count")+
  scale_colour_viridis_c(aesthetics = c("colour"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("SE of Trend",firstYear,"-",lastYear))+
  labs(title = paste(species, "Standard error"))+
  facet_wrap(vars(model))
# 
# pdf(paste0("Figures/Two_SE_model_comparison_",species_f,".pdf"),
#     height = 5,
#     width = 8)
# print(map_se)
# dev.off()




# pdf(paste0("Figures/Two_model_comparison_",species_f,".pdf"),
#     height = 5,
#     width = 8)
print(map / map_se)
# dev.off()



}#end species loop


dev.off()



