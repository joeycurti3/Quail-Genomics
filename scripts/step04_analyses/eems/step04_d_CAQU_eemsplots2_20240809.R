######################################################
###### R script to visualize results from EEMS #######
######################################################

# Author: Joey Curti 
# Date: SUN AUG 11 2024
# Description: Make maps of migration surface using eems outputs
# References:
# https://devonderaad.github.io/aph.rad/eems/run.eems.html
# https://github.com/dipetkov/reemsplots2

## Clean Workspace and set working directory

setwd("~/Downloads/")
rm(list = ls())

## Dependencies

library("devtools")
install_github("dipetkov/reemsplots2")
library(reemsplots2)
library(patchwork)

## Read in data

# Quail points 
points <- read.csv("~/Downloads/20240729_CAQU_Locations_Numbers.csv",sep = ",", header = T, stringsAsFactors = F)
points_samo <- points[1:21,]

# Polygon of eems output to clip roads to
coords <- read.csv("~/Downloads/GCA_023055505.1_bCalCai1.0.p_PassSNPs_unrelated_SAMO_chr1.outer",stringsAsFactors = F, header = F, sep="\t")
coords <- coords[-120,] # Weird last row

polygon <- coords %>%
  st_as_sf(coords = c("V1", "V2"),crs="+proj=longlat +datum=WGS84") %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# Road data - Los Angeles County 
roads_la <- st_read("~/Downloads/tl_2023_06037_roads/")
roads_la_sf <- st_transform(roads_la,crs="+proj=longlat +datum=WGS84")
roads_la_crop <- st_intersection(roads_la_sf,polygon)
road_la_crop_PrimaryAndSecondary <- roads_la_crop[roads_la_crop$MTFCC %in% c("S1100","S1200"),]

# Road data - Ventura County
roads_vent <- st_read("~/Downloads/tl_2023_06111_roads/")
roads_vent_sf <- st_transform(roads_vent,crs="+proj=longlat +datum=WGS84")
roads_vent_crop <- st_intersection(roads_vent_sf,polygon)
roads_vent_crop_PrimaryAndSecondary <- roads_vent_crop[roads_vent_crop$MTFCC %in% c("S1100","S1200"),]

# Combine road data 
roads_sf_crop_PrimaryAndSecondary <- do.call(rbind,list(road_la_crop_PrimaryAndSecondary,roads_vent_crop_PrimaryAndSecondary))

## Plot it - nDemes = 100

plots_nDemes100 <- make_eems_plots(c("~/Downloads/nIndiv56-nSites25199223-nDemes100-chain1","~/Downloads/nIndiv56-nSites25199223-nDemes100-chain2","~/Downloads/nIndiv56-nSites25199223-nDemes100-chain3"),
                         longlat = TRUE,
                         add_demes = FALSE,
                         add_grid = TRUE,
                         add_outline = TRUE,
                         col_outline = "darkgrey")  


gg_nDemes100 <- plots_nDemes100$mrates01 + 
  labs(title = "nDemes=100") +
  theme(plot.title = element_text(hjust=.9, vjust=-8),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill=NA, colour = "black"),
        panel.background = element_rect(fill=NA))

gg_nDemes100_roads <- gg_nDemes100 + geom_sf(data = roads_sf_crop_PrimaryAndSecondary,colour="grey40", linewidth=.6,inherit.aes = FALSE) + # inherit.aes=false gets around 'object 'y' not found' error
  geom_point(data=points_samo,aes(longitude,latitude),colour="black", fill="black",size=5,pch=19) +
  geom_text(data=points_samo,aes(longitude,latitude),label=points_samo$sample_num, colour="white",size=3, fontface="bold") +
  theme(axis.text.x = element_text(colour = "white"))
  
## Plot it - nDemes = 200

plots_nDemes200 <- make_eems_plots(c("~/Downloads/nIndiv56-nSites25199223-nDemes200-chain1","~/Downloads/nIndiv56-nSites25199223-nDemes200-chain2","~/Downloads/nIndiv56-nSites25199223-nDemes200-chain3"),
                                   longlat = TRUE,
                                   add_demes = FALSE,
                                   add_grid = TRUE,
                                   add_outline = TRUE,
                                   col_outline = "darkgrey") 

gg_nDemes200 <- plots_nDemes200$mrates01 + 
  labs(title = "nDemes=200") +
  theme(plot.title = element_text(hjust=.9, vjust=-8),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill=NA, colour = "black"),
        panel.background = element_rect(fill=NA))

gg_nDemes200_roads <- gg_nDemes200 + geom_sf(data = roads_sf_crop_PrimaryAndSecondary,colour="grey40", linewidth=.6,inherit.aes = FALSE) + # inherit.aes=false gets around 'object 'y' not found' error
  geom_point(data=points_samo,aes(longitude,latitude),colour="black", fill="black",size=5,pch=19) +
  geom_text(data=points_samo,aes(longitude,latitude),label=points_samo$sample_num, colour="white",size=3, fontface="bold") +
  theme(legend.position = "none", axis.text.x = element_text(colour = "white"))

## Plot it - nDemes = 300

plots_nDemes300 <- make_eems_plots(c("~/Downloads/nIndiv56-nSites25199223-nDemes300-chain1","~/Downloads/nIndiv56-nSites25199223-nDemes300-chain2","~/Downloads/nIndiv56-nSites25199223-nDemes300-chain3"),
                                   longlat = TRUE,
                                   add_demes = FALSE,
                                   add_grid = TRUE,
                                   add_outline = TRUE,
                                   col_outline = "darkgrey") 

gg_nDemes300 <- plots_nDemes300$mrates01 + 
  labs(title = "nDemes=300") +
  theme(plot.title = element_text(hjust=.9, vjust=-8),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill=NA, colour = "black"),
        panel.background = element_rect(fill=NA))

gg_nDemes300_roads <- gg_nDemes300 + geom_sf(data = roads_sf_crop_PrimaryAndSecondary,colour="grey40", linewidth=.6,inherit.aes = FALSE) + # inherit.aes=false gets around 'object 'y' not found' error
  geom_point(data=points_samo,aes(longitude,latitude),colour="black", fill="black",size=5,pch=19) +
  geom_text(data=points_samo,aes(longitude,latitude),label=points_samo$sample_num, colour="white",size=3, fontface="bold") +
  theme(legend.position = "none")

# Arrange them all and save

gg_final <- gg_nDemes100_roads / gg_nDemes200_roads / gg_nDemes300_roads


ggsave(filename="20240811_CAQU_eems_plots.png", plot = gg_final, height = 12, width = 6,device="png",dpi=600,bg="white")
