#####################################
## Create inputs for Circuitscape ##
####################################

# Author: Joseph Curti (jcurti3@g.ucla.edu)
# Date: TUE JUN 10 2025
# Version: V1

## References: 

## Set up workspace

rm(list = ls())
setwd("~/Desktop/MYYU Landscape Genomics/")

library(terra)
library(sf)
library(ggplot2)

## Import HSM and spatial points 

# Load resistance raster - factor of 3
suitability_raster_agr3 <- rast("RF_pred_occ_all_mean_no_isl_agg_3.tiff")
crs(suitability_raster_agr3)
plot(suitability_raster_agr3)

# invert so have a map of resistance, not suitability
resistance_raster_agr3 <- 1 - suitability_raster_agr3
crs(suitability_raster_agr3) == crs(resistance_raster_agr3) #TRUE
plot(resistance_raster_agr3)

# output new raster for Circuitscape
writeRaster(resistance_raster_agr3, "resistance_surface_agr3.tif", overwrite = TRUE)

# Load spatial points
points <- read.csv("MYYU_2025_CA.csv", sep=",", stringsAsFactors = F) 
points_sf <- st_as_sf(points, coords = c("longitude", "latitude"), crs="EPSG:4326")
points_trans <- st_transform(points_sf, crs(resistance_raster))
crs(points_trans) == crs(resistance_raster) #TRUE

# Output points with row numbers 
points <- st_coordinates(points_trans$geometry)
ids <- seq(1:45)
write.table(cbind(ids,points), "points_trans.txt",sep="\t", col.names = F, row.names = F)

## Plotting

# plot raster and points to make sure they look alright
df_resistance_agr3 <- as.data.frame(resistance_raster_agr3, xy=T) # need to convert raster to df

pp1 <- ggplot() +
  geom_raster(data = df_resistance_agr3, aes(x = x, y = y, fill = `RF_pred_occ_all_mean_no_isl_agg_3`)) +
  scale_fill_gradient(name = "Resistance",low = "white", high = "black") +
  geom_sf(data=points_trans, color="gold") +
  theme_void()

ggsave(filename = "myyu_resistance_map_agr3.tif", plot=pp1, device = "tif",dpi = 300)

