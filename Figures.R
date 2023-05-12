#Figures

library(tidyverse)
library(patchwork)


# Figure_1 ----------------------------------------------------------------

# demonstration of neighbours for BHVI

# load saved neighbours data from neighbours_define_voronoi(..., save_plot_data = TRUE)

load("data/Blue-headed_Vireo_route_maps_data.RData")

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



