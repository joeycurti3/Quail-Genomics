##############################################################
## Spatial Analyses and Reciprical Causal Modeling for CAQU ##
##############################################################

# Author: Zac MacDonald (zmacdonald@ioes.ucla.edu)
# Description: create spatial layers and environmental/resistance distances for input into Reciprocal
# Version: V1
# Date: Oct 14 2024

## References
# https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15604?casa_token=2DUf0mKA5toAAAAA%3A7W8mACi3cktR-G_rJf8HCnOIIlHJyKIwIUJiBhaOw-xrn3KLVwvXU0PeFrLiAJ-WrH-QrxSuaFVk

## Clean Workspace

rm(list = ls())

## Dependencies 

### R version 4.1.2
### List of required packages
required.pkg <- c("terra", # replaces raster
                  "rgdal",
                  "raster", # most raster-based functions
                  "sp", # most spatial points functions
                  "data.table",# basic data manipulation
                  "spdplyr",
                  "spatialEco", # heatload index
                  "sf", # many spatial things
                  "dismo", # implements MaxEnt
                  "gdistance", # estimate resistance distance
                  "vegan", # estimate correlation coefficients
                  "ENMeval", # various MaxEnt extensions
                  "maptools",
                  "rgeos",
                  "prettymapr"
                  )

pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
if (length(pkgs.not.installed) > 0) install.packages(pkgs.not.installed, dependencies=TRUE)

# Load the required libraries
lapply(required.pkg, library, character.only = TRUE) 

# # I couldn't get first function to work, so use devtools
# # devtools::install_local("spdplyr_0.4.0.tar.gz") 
# library("spdplyr")
# # devtools::install_github("danddr/USE")
# library(USE)

### Set wd
setwd(<insert working directory>)

################### Load specimens and gbif occurrences ###################
### Import Joey's specimen metadata
californica_wgs <- fread("<insert path to metadata file>") 
head(californica_wgs)
californica_wgs$ID <- californica_wgs$individuals
length(unique(californica_wgs$latitude))
unique(californica_wgs$locality)
# only include individuals for which Joey generated genetic distance
genetic_dist_nei_ind <- fread("<insert path to output file of step04_c_CAQU_IndGeneticDist_20240820.R>")
head(genetic_dist_nei_ind)
ind_list <- data.frame("ID"=genetic_dist_nei_ind$V1)
nrow(ind_list)
californica_wgs_46_inds <- merge(californica_wgs,ind_list, by="ID")
nrow(californica_wgs_46_inds)
### Make it spatial, baby
californica_wgs.sp <- vect(californica_wgs_46_inds, geom=c("longitude", "latitude"), crs="+init=epsg:4326", keepgeom=FALSE)
plot(californica_wgs.sp)
### reproject to epsg:3310 and crop to study area
californica_wgs.sp <- terra::project(californica_wgs.sp, CRS("+init=epsg:3310"))
plot(californica_wgs.sp)

### Import iNat locations
californica_iNat <- fread("<insert path to iNat data>")
californica_iNat <- californica_iNat %>% filter(californica_iNat$positional_accuracy <= 100)
# californica_iNat <- californica_iNat %>% filter(californica_iNat$quality_grade == "research")
californica_iNat
### Make it spatial, baby
californica_iNat.sp <- vect(californica_iNat, geom=c("longitude","latitude"), crs="+init=epsg:4326", keepgeom=FALSE)
plot(californica_iNat.sp)
### reproject to epsg:3310 and crop to study area
californica_iNat.sp <- terra::project(californica_iNat.sp, CRS("+init=epsg:3310"))
plot(californica_iNat.sp)
plot(californica_wgs.sp, add=T, col="red", cex=1)

### Define study extent (MCP of WGS samples buffered by 50km)
Euclidean_dist <- terra::distance(californica_wgs.sp, unit="m")
range(Euclidean_dist) # max separation of samples is ~43 km
study_extent_californica <- terra::buffer(convHull(californica_wgs.sp), 25000)
californica_iNat.sp <- terra::intersect(californica_iNat.sp, study_extent_californica)
plot(study_extent_californica)
plot(californica_iNat.sp, add=T, col="black", cex=1)
plot(californica_wgs.sp, add=T, col="red", cex=0.5)
writeVector(x = study_extent_californica, file = "<insert study extent buffer>", overwrite=T)


################### Load spatial data ###################
### Get predictor layers
# Note: predictor layers are from a study in review, MacDonald et al. 2025
predictor_list <- c("TRI", 
                    "heatload", 
                    "landcover", 
                    "hydro_rivers_lakes_raster_100m_buffer",
                    "dist_water_hydrosheds",
                    "EVI_median",
                    "night_light",
                    "ClimateNA_Normal_1991_2020_SpatRaster_uncorrelated")
predictors_spatrast <- rast(paste("spatial_data_layers/ClimateNA/ClimateNA/", predictor_list, ".tif", sep=""))
names(predictors_spatrast)
layer_names <- c("TRI",
                 "heatload",
                 "landcover",
                 "rivers_lakes",
                 "dist_water",
                 "EVI_median",
                 "night_light",
                 "continentality",
                 "temp_winter",
                 "temp_summer",
                 "precip_winter",
                 "precip_summer")
names(predictors_spatrast) <- layer_names
### Crop down the predictor layers:
predictors_spatrast_californica <- mask(crop(predictors_spatrast, study_extent_californica), study_extent_californica)
plot(predictors_spatrast_californica$TRI)
plot(californica_iNat.sp, add=T, col="red", cex=1)
plot(californica_wgs.sp, add=T)
plot(study_extent_californica, add=T)
writeRaster(predictors_spatrast_californica, filename="Callipepla_californica/ClimateNA_Normal_1991_2020_Callipepla_californica.tif", overwrite=TRUE)

### Add NLCD 2016 Percent Developed Imperviousness (CONUS): https://www.mrlc.gov/data/nlcd-2016-percent-developed-imperviousness-conus
imperv <- rast("~/nlcd_2016_impervious_l48_20210604.img")
# crop data before reprojecting to save time (make super big buffer to avoid issues)
imperv_crop <- terra::crop(imperv, terra::project(terra::buffer(convHull(study_extent_californica), 10000), imperv), snap='out')
plot(imperv_crop) # looks good
# Now that it is cropped, we can reproject, mask, and write the raster to hard drive, and have another beer
imperv_californica <-  terra::crop(terra::project(imperv_crop, predictors_spatrast_californica$TRI), predictors_spatrast_californica$TRI, mask=TRUE)
plot(imperv_californica) # looks good
names(imperv_californica) <- "imperviousness"
writeRaster(imperv_californica, filename="~/imperv_californica.tif", overwrite=TRUE)
# Add imperviousness to spatraster
predictors_spatrast_californica <- c(predictors_spatrast_californica,imperv_californica)
global(predictors_spatrast_californica, fun="isNA")
writeRaster(predictors_spatrast_californica, filename="~/ClimateNA_Normal_1991_2020_Callipepla_californica.tif", overwrite=TRUE)


### Make occ file
# merge locations of sequenced inds and gbif records
occ.sp <- rbind(californica_wgs.sp, californica_iNat.sp)
plot(occ.sp)
# extract values
occ <- terra::extract(as.numeric(predictors_spatrast), occ.sp)
# apply isNA for all columns (exclude any NA)
occ$anyNA <- rowSums(occ)
occ.sp$hasNA = occ$anyNA
# get rid of na values
occ_noNA.sp <- occ.sp[complete.cases(occ.sp$hasNA),]
# thin occ data (keep one occurrence point per cell)
# prep occ_species for cellfromxy
occ_noNA_df <- as.data.frame(occ_noNA.sp, geom = "XY")
occ_noNA_df <- occ_noNA_df[,c("x","y")]
cells <- cellFromXY(predictors_spatrast[[1]], occ_noNA_df)
dups <- duplicated(cells)
occ_noNA_thin <- occ_noNA_df[!dups, ]
cat(nrow(occ_noNA_df) - nrow(occ_noNA_thin), "records are removed")
nrow(occ_noNA_thin)
colnames(occ_noNA_thin) <- c("easting", "northing")
# turn back into spatvec points (might not need this though)
occ_noNA_thin.sp <- vect(occ_noNA_thin, geom = c("easting", "northing"))
crs(occ_noNA_thin.sp) <- "epsg:3310"
plot(occ_noNA_thin.sp)

### MaxEnt via dismo requires raster objects
predictors <- raster::stack(x="~/ClimateNA_Normal_1991_2020_Callipepla_californica.tif")
names(predictors)
predictors$rivers_lakes <- ratify(predictors$rivers_lakes)
predictors$landcover <- ratify(predictors$landcover)
is.factor(predictors)
### Write the final compilation of layers to HD for Joey
writeRaster(predictors, filename="~/layers_for_Joey/.tif", bylayer=TRUE, overwrite=TRUE, suffix=names(predictors))


### Running MaxEnt via dismo package (example http://spatialecology.weebly.com/r-code--data/category/sdm-maxent)
ENM_ClimateNA_californica <- maxent(predictors, occ_noNA_thin, nbg = 10000, factors=c("landcover", "rivers_lakes"), removeDuplicates = T, path="Callipepla_californica/maxent_models/ClimateNA_100m_riverlakes")
### Plot showing importance of each variable
plot(ENM_ClimateNA_californica)
### Predict habitat suitability across extent
maxent_pred_ENM_ClimateNA_californica <- predict(ENM_ClimateNA_californica, predictors)
maxent_pred <- maxent_pred_ENM_ClimateNA_californica
plot(maxent_pred)
plot(californica_iNat.sp, add=T, col="red", cex=1)
plot(californica_wgs.sp, add=T)
plot(study_extent_californica, add=T)
writeRaster(maxent_pred, filename="~/maxent_pred_ClimateNA_californica.tif", overwrite=TRUE)
writeRaster(maxent_pred, filename="~/maxent_pred_ClimateNA_californica.tif", overwrite=TRUE)
writeRaster(maxent_pred, filename="~/layers_for_Joey/maxent_pred_ClimateNA_californica.tif", overwrite=TRUE)

### Make habitat suitability figure
pdf(file = "~/figures/maxent_pred_ENM_ClimateNA_californica.pdf", width = 6, height = 4, family = "serif")
par(oma=c(0,0,0,0))
par(mar = c(0.55, 0.55, 1.5, 2)) 
# bottom, left, top and right margins
plot(maxent_pred, main="", cex.main=2,
     xlab="", ylab="", cex.lab=1, bty="n", box=TRUE, legend=FALSE, axis = FALSE, col.axis = "white", col.lab = "white", tck = 0)
title(main=expression(italic("Callipepla californica")* " habitat suitability"), adj = 0, line = 0.5, cex.main = 1)
plot(maxent_pred, legend.only=T, legend.args=list(text='  index', cex=1, line=1), axis.args=list(cex.axis=0.8))
plot(californica_iNat.sp, pch=16, col="blue", cex=0.5, add=T)
plot(californica_wgs.sp, pch=16, col="red", cex=0.75, add=T)
legend('bottom', bg="white", expression("sequenced indivudals", "iNaturalist records", "AUC = 0.88"),
       pch=c(16,16, 0), col=c("red", "blue", "white"), bty="n", cex=c(0.75))
addnortharrow(pos = "bottomright", padin = c(0.1, 0.1), scale = 0.5,
              lwd = 1, border = "black", cols = c("white", "black"),
              text.col = "black")
dev.off()

### Variable contribution plot
pdf(file = "~/maxent_pred_var_cont.pdf", width = 6, height = 6, family = "serif")
plot(ENM_ClimateNA_californica)
dev.off()

### Extent indicator plot
state_border <- vect("<insert CA boundary shapefile>")
crs(state_border)
state_border <- terra::project(state_border, CRS("+init=epsg:3310"))

pdf(file = "~/extent_indicator.pdf", width = 6, height = 6, family = "serif")
plot(state_border,axis = FALSE, col.axis = "white", col.lab = "white") 
plot(maxent_pred, add=T, legend=FALSE, axis = FALSE, col.axis = "white", col.lab = "white")
dev.off()

################### Resistance, LCP, and Euclidean distances ###################
### Import MaxEnt predicted surface
maxent_pred <- raster(x="~/maxent_pred_ClimateNA_californica.tif") # reading in as RasterLayer instead of SpatRast makes a better plot (legends go in wrong place otherwise, wtf)
crs(maxent_pred) <- "+init=epsg:3310"
# transform into a transition objects
tr1 <- transition(maxent_pred, transitionFunction=mean, directions=8)
# correct for the fact that diagonal neighbors are more remote from each other than orthogonal neighbors
tr1C <- geoCorrection(tr1, type="c")
plot(raster(tr1C), main="raster(tr1C)") # looks good

### Calculate Resistance, LCP, and Euclidean distances for 1) individuals then 3) locations
### Individuals
# need SpatialPoints, won't take SpatVec
occurrences <- rbind(geom(californica_wgs.sp, df=T))
occurrences <- occurrences[,3:4]
location.coords <- as.matrix(occurrences)
coordinates(occurrences) <- c("x", "y")
californica.SpatialPoints <- SpatialPoints(coordinates(occurrences), proj4string=CRS("+init=epsg:3310"))
plot(californica.SpatialPoints)
# run commuteDistance
RD1 <- commuteDistance(tr1C, californica.SpatialPoints)
RD1_matrix <- as.matrix(RD1)
ind_names <- californica_wgs.sp$ID
rownames(RD1_matrix) <- ind_names
colnames(RD1_matrix) <- ind_names
dim(RD1_matrix)
# View(RD1_matrix)
write.table(RD1_matrix, file="~/hs_resistance_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)

### Calculate pairwise least-cost distances between sequenced individuals
LCD1 <- costDistance(tr1C, californica.SpatialPoints)
# # visualize
# for (i in 1:nrow(location.coords)){
#   for (n in 2:nrow(location.coords)){
#     plot(shortestPath(tr1C, location.coords[i,], location.coords[n,], output="SpatialLines"), add=TRUE)
#   }}
# get matrix
LCD1_matrix <- as.matrix(LCD1)
rownames(LCD1_matrix) <- ind_names
colnames(LCD1_matrix) <- ind_names
dim(LCD1_matrix)
# View(LCD1_matrix)
write.table(LCD1_matrix, file="~/hs_least_cost_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)

### Calculate Euclidean distances between locations of sequenced individuals
euclidean_dist_matrix <- as.matrix(terra::distance(californica_wgs.sp, unit="m"))
rownames(euclidean_dist_matrix) <- ind_names
colnames(euclidean_dist_matrix) <- ind_names
# View(euclidean_dist_matrix)
dim(euclidean_dist_matrix)
write.table(euclidean_dist_matrix, file="~/euclidean_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)

### Run quick Mantel test between resistance, least-cost, and Euclidean distance
vegan::mantel(RD1_matrix, LCD1_matrix, method="pearson", permutations=999)
vegan::mantel(RD1_matrix, euclidean_dist_matrix, method="pearson", permutations=999)
vegan::mantel(LCD1_matrix, euclidean_dist_matrix, method="pearson", permutations=999)


### Extract suitability and predictor variable values for each sequenced individual
suitability_values_ind <- terra::extract(rast(maxent_pred), californica_wgs.sp)
colnames(suitability_values_ind) <- c("ID", "habitat_suitability")
predictor_values_ind <- terra::extract(rast(predictors), californica_wgs.sp)
Callipepla_californica_extracted_values_ind <- cbind(californica_wgs.sp$ID, suitability_values_ind, predictor_values_ind)
head(Callipepla_californica_extracted_values_ind)
Callipepla_californica_extracted_values_ind <- Callipepla_californica_extracted_values_ind[ , !(names(Callipepla_californica_extracted_values_ind) %in% "ID")]
names(Callipepla_californica_extracted_values_ind)[names(Callipepla_californica_extracted_values_ind) == 'californica_wgs.sp$ID'] <- 'ID'
write.csv(Callipepla_californica_extracted_values_ind, file="~/Callipepla_californica_extracted_values_ind.csv", row.names = FALSE)

### Calculate environmental distances
# sequenced individuals
Callipepla_californica_extracted_values_ind <- fread("~/Callipepla_californica_extracted_values_ind.csv")
colnames(Callipepla_californica_extracted_values_ind)
continentality_matrix  <- as.matrix(dist(Callipepla_californica_extracted_values_ind$continentality, diag = TRUE, upper = TRUE, p = 2))
temp_winter_matrix  <- as.matrix(dist(Callipepla_californica_extracted_values_ind$temp_winter, diag = TRUE, upper = TRUE, p = 2))
temp_summer_matrix <- as.matrix(dist(Callipepla_californica_extracted_values_ind$temp_summer, diag = TRUE, upper = TRUE, p = 2))
precip_winter_matrix <- as.matrix(dist(Callipepla_californica_extracted_values_ind$precip_winter, diag = TRUE, upper = TRUE, p = 2))
precip_summer_matrix <- as.matrix(dist(Callipepla_californica_extracted_values_ind$precip_summer, diag = TRUE, upper = TRUE, p = 2))
# name rows/columns
rownames(continentality_matrix) <- colnames(continentality_matrix) <-  Callipepla_californica_extracted_values_ind$ID
rownames(temp_winter_matrix) <- colnames(temp_winter_matrix) <-  Callipepla_californica_extracted_values_ind$ID
rownames(temp_summer_matrix) <- colnames(temp_summer_matrix) <-  Callipepla_californica_extracted_values_ind$ID
rownames(precip_winter_matrix) <- colnames(precip_winter_matrix) <-  Callipepla_californica_extracted_values_ind$ID
rownames(precip_summer_matrix) <- colnames(precip_summer_matrix) <-  Callipepla_californica_extracted_values_ind$ID
# write matrices
write.table(continentality_matrix, file="~/continentality_dist_ind.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(temp_winter_matrix, file="~/temp_winter_dist_ind.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(temp_summer_matrix, file="~/temp_summer_dist_ind.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(precip_summer_matrix, file="~/precip_summer_dist_ind.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(precip_winter_matrix, file="~/precip_winter_dist_ind.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)


################### Road rasters ###################
### Rasterize road shapefiles that Joey provided
# for some reason can't combine "lines" with "union" so do LA and Ventura independently then combine rasters after
# LA
roads_LA <- vect("~/tl_2023_06037_roads.shp")
roads_LA <- terra::project(roads_LA, CRS("+init=epsg:3310"))
# plot(roads_LA)
roads_LA <- subset(roads_LA, roads_LA$MTFCC == "S1100" | roads_LA$MTFCC == "S1200")
values(roads_LA)<-1
roads_LA_rast <- terra::mask(terra::crop(rasterize(roads_LA, predictors_spatrast_californica$TRI, "value", background=0, touches=T), predictors_spatrast_californica$TRI), predictors_spatrast_californica$TRI)
plot(roads_LA_rast)
# Ventura
roads_Ventura <- vect("~/tl_2023_06111_roads.shp")
roads_Ventura <- terra::project(roads_Ventura, CRS("+init=epsg:3310"))
roads_Ventura <- subset(roads_Ventura, roads_Ventura$MTFCC == "S1100" | roads_Ventura$MTFCC == "S1200")
values(roads_Ventura)<-1
roads_Ventura_rast <- terra::mask(terra::crop(rasterize(roads_Ventura, predictors_spatrast_californica$TRI, "value", background=0, touches=T), predictors_spatrast_californica$TRI), predictors_spatrast_californica$TRI)
plot(roads_Ventura_rast)
# combine
roads_all<- mosaic(roads_LA_rast, roads_Ventura_rast, fun = "sum")
roads_all <- subst(roads_all, 2, 1, others=NULL, raw=T)
plot(roads_all)
writeRaster(roads_all, filename="~/roads_new.tif", overwrite=TRUE)

################### Resistance and LCP based on roads ###################
# import road raster
roads_all <- rast("~/roads_new.tif")
plot(roads_all)
roads_all_conductance <- subst(roads_all, 0, 10, raw=T) # make landscape 10x more conductive than roads
plot(roads_all_conductance)
roads_all <- raster(roads_all_conductance) # reading in as RasterLayer instead of SpatRast makes a better plot (legends go in wrong place otherwise, wtf)
plot(roads_all)
# transform into a transition objects
tr1 <- transition(roads_all, transitionFunction=mean, directions=8)
# correct for the fact that diagonal neighbors are more remote from each other than orthogonal neighbors
tr1C <- geoCorrection(tr1, type="c")
plot(raster(tr1C), main="raster(tr1C)") # looks good
writeRaster(raster(tr1C), filename="~/road_conductance_layer.asc", overwrite=TRUE)

# Make a figure showing this transition object
pdf(file = "~/road_conductance_surface.pdf", width = 6, height = 4, family = "serif")
par(oma=c(0,0,0,0))
par(mar = c(0.55, 0.55, 1.5, 2)) 
# bottom, left, top and right margins
plot(raster(tr1C), main="", cex.main=2,
     xlab="", ylab="", cex.lab=1, bty="n", box=TRUE, legend=FALSE, axis = FALSE, col.axis = "white", col.lab = "white", tck = 0)
title(main="Conductance based on major roads", adj = 0, line = 0.5, cex.main = 1)
plot(raster(tr1C), legend.only=T, legend.args=list(text='  index', cex=1, line=1), axis.args=list(cex.axis=0.8))
plot(californica_iNat.sp, pch=16, col="blue", cex=0.5, add=T)
plot(californica_wgs.sp, pch=16, col="red", cex=0.75, add=T)
legend('bottom', bg="white", expression("sequenced indivudals", "iNaturalist records", "AUC = 0.88"),
       pch=c(16,16, 0), col=c("red", "blue", "white"), bty="n", cex=c(0.75))
addnortharrow(pos = "bottomright", padin = c(0.1, 0.1), scale = 0.5,
              lwd = 1, border = "black", cols = c("white", "black"),
              text.col = "black")
dev.off()


### Calculate the resistance distance between sequenced individuals
occurrences <- rbind(geom(californica_wgs.sp, df=T))
occurrences <- occurrences[,3:4]
location.coords <- as.matrix(occurrences)
coordinates(occurrences) <- c("x", "y")
californica.SpatialPoints <- SpatialPoints(coordinates(occurrences), proj4string=CRS("+init=epsg:3310"))
plot(californica.SpatialPoints)
# run commuteDistance
RD1 <- commuteDistance(tr1C, californica.SpatialPoints)
RD1_matrix <- as.matrix(RD1)
ind_names <- californica_wgs.sp$ID
rownames(RD1_matrix) <- ind_names
colnames(RD1_matrix) <- ind_names
dim(RD1_matrix)
# View(RD1_matrix)
write.table(RD1_matrix, file="~/roads_resistance_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)

### Calculate pairwise least-cost distances between sequenced individuals
LCD1 <- costDistance(tr1C, californica.SpatialPoints)
# # visualize
# for (i in 1:nrow(location.coords)){
#   for (n in 2:nrow(location.coords)){
#     plot(shortestPath(tr1C, location.coords[i,], location.coords[n,], output="SpatialLines"), add=TRUE)
#   }}
# get matrix
LCD1_matrix <- as.matrix(LCD1)
ind_names <- californica_wgs.sp$ID
rownames(LCD1_matrix) <- ind_names
colnames(LCD1_matrix) <- ind_names
dim(LCD1_matrix)
# View(LCD1_matrix)
write.table(LCD1_matrix, file="~/roads_least_cost_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)


################### Traffic Volume Layer  ###################
### Import MaxEnt predicted surface
roads_ryan <- rast("~/InterpolatedSoCalTVw1.tif")
crs(roads_ryan) <- "+init=epsg:3857"
# Reproject, mask, and write the raster to hard drive
roads_ryan_californica <-  terra::crop(terra::project(roads_ryan, predictors_spatrast_californica$TRI), predictors_spatrast_californica$TRI, mask=TRUE)
plot(roads_ryan_californica) # looks good
writeRaster(roads_ryan_californica, filename="~/roads_ryan_californica.tif", overwrite=TRUE)

### Resistance and LCP distance
# transform into a transition objects
tr1 <- transition(raster(roads_ryan_californica), transitionFunction=mean, directions=8)
# correct for the fact that diagonal neighbors are more remote from each other than orthogonal neighbors
tr1C <- geoCorrection(tr1, type="c")
plot(raster(tr1C), main="raster(tr1C)") # looks good

### Calculate Resistance, LCP, and Euclidean distances for 1) individuals then 3) locations
### Individuals
# need SpatialPoints, won't take SpatVec
occurrences <- rbind(geom(californica_wgs.sp, df=T))
occurrences <- occurrences[,3:4]
location.coords <- as.matrix(occurrences)
coordinates(occurrences) <- c("x", "y")
californica.SpatialPoints <- SpatialPoints(coordinates(occurrences), proj4string=CRS("+init=epsg:3310"))
plot(californica.SpatialPoints)
# run commuteDistance
RD1 <- commuteDistance(tr1C, californica.SpatialPoints)
RD1_matrix <- as.matrix(RD1)
ind_names <- californica_wgs.sp$ID
rownames(RD1_matrix) <- ind_names
colnames(RD1_matrix) <- ind_names
dim(RD1_matrix)
# View(RD1_matrix)
write.table(RD1_matrix, file="~/roads_traffic_resistance_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)

### Calculate pairwise least-cost distances between sequenced individuals
LCD1 <- costDistance(tr1C, californica.SpatialPoints)
# # visualize
# for (i in 1:nrow(location.coords)){
#   for (n in 2:nrow(location.coords)){
#     plot(shortestPath(tr1C, location.coords[i,], location.coords[n,], output="SpatialLines"), add=TRUE)
#   }}
# get matrix
LCD1_matrix <- as.matrix(LCD1)
rownames(LCD1_matrix) <- ind_names
colnames(LCD1_matrix) <- ind_names
dim(LCD1_matrix)
# View(LCD1_matrix)
write.table(LCD1_matrix, file="~/roads_traffic_least_cost_dist_ind.txt", 
            sep = " ", dec = ".", row.names = TRUE, col.names = TRUE, quote = F)


################### Reciprocal causal modelling ###################
### Assess correlation between resistance distance, least-cost distance (LCP), and Euclidean distance
# import geographic distance matrices, and make sure values are numeric
# note that these all import a little differently (R messes up row/column names), so check each one
ind_names <- californica_wgs.sp$ID
# loc_names <- californica_wgs_loc.sp$Pop
files <- list.files(path = paste("<insert path to output distance matrices from above>", sep = ""), pattern = ".txt")
files
for(i in 1:length(files)){  
  temp <- fread(paste("Callipepla_californica/distance_matrices/", files[i], sep = ""), sep = " ", dec = ".", header = FALSE)
  temp <- temp[,-1] # removing row names
  # rownames(temp) <-  ind_names
  # colnames(temp) <- ind_names
  name <- substr(files[i],1,nchar(files[i])-4)
  assign(name, temp)  
}
# View(hs_least_cost_dist_ind)

### Import genetic distances
# Nei's
genetic_dist_nei_ind <- fread("<insert path to output file of step04_c_CAQU_IndGeneticDist_20240820.R>")
genetic_dist_nei_ind <- genetic_dist_nei_ind[,-1] # note that the rownames are included as a variable; have to delete
head(genetic_dist_nei_ind)
# Euclidean
genetic_dist_euclidean_ind <- fread("<insert path to output file of step04_c_CAQU_IndGeneticDist_20240820.R>")
head(genetic_dist_euclidean_ind) # note that the rownames are included as a variable; have to delete
genetic_dist_euclidean_ind <- genetic_dist_euclidean_ind[,-1]
# Dps
genetic_dist_dps_ind <- fread("<insert path to output file of step04_c_CAQU_IndGeneticDist_20240820.R>")
head(genetic_dist_dps_ind) # note that the rownames are included as a variable; have to delete
genetic_dist_dps_ind <- genetic_dist_dps_ind[,-1]
# assess correlations
vegan::mantel(genetic_dist_nei_ind, genetic_dist_euclidean_ind, method="pearson", permutations=999) # redundant distances
vegan::mantel(genetic_dist_nei_ind, genetic_dist_dps_ind, method="pearson", permutations=999) 
vegan::mantel(genetic_dist_dps_ind, euclidean_dist_ind, method="pearson", permutations=999) 

# ### Import linearized FST (made by Joey)
# FST_loc <- fread("Callipepla_californica/linearizedFST_CAQU_20240820.csv") 
# FST_loc <- FST_loc[,-1]
# dim(FST_loc)
# vegan::mantel(FST_loc, euclidean_dist_loc, method="pearson", permutations=999) # no IBD here

RCM_data <- list(genetic_dist_euclidean_ind,
                 euclidean_dist_ind,
                 hs_least_cost_dist_ind,
                 hs_resistance_dist_ind, 
                 roads_least_cost_dist_ind,
                 roads_resistance_dist_ind,
                 roads_traffic_least_cost_dist_ind,
                 roads_traffic_resistance_dist_ind,
                 continentality_dist_ind,
                 temp_summer_dist_ind, 
                 temp_winter_dist_ind,
                 precip_summer_dist_ind,
                 precip_winter_dist_ind)

names(RCM_data) <- c("genetic_dist",
                      "euclidean_dist",
                      "hs_least_cost_dist",
                      "hs_resistance_dist",
                      "roads_least_cost_dist",
                      "roads_resistance_dist",
                      "traffic_least_cost_dist",
                      "traffic_resistance_dist",
                      "continentality_dist",
                      "temp_summer_dist", 
                      "temp_winter_dist",
                      "precip_summer_dist",
                      "precip_winter_dist")

test <- RCM_data[[1]] #worked
require(vegan)
x <- mantel.partial(RCM_data[[1]], RCM_data[[2]], RCM_data[[3]], method = "pearson", permutations = 999)
x$statistic
x$signif


### now try the loop
variables = length(RCM_data)-1 # dependent variable (genetic distance) not included!
mantel_output <- matrix(ncol=variables, nrow=variables)
rownames(mantel_output) <- c("euclidean_dist",
                             "hs_least_cost_dist",
                             "hs_resistance_dist",
                             "roads_least_cost_dist",
                             "roads_resistance_dist",
                             "traffic_least_cost_dist",
                             "traffic_resistance_dist",
                             "continentality_dist",
                             "temp_summer_dist", 
                             "temp_winter_dist",
                             "precip_summer_dist",
                             "precip_winter_dist")
colnames(mantel_output) <- rownames(mantel_output)

# print(variables+1) can we chunk this in below?
for (i in 2:length(RCM_data)){
  for (n in 2:length(RCM_data)){
    x <- mantel.partial(RCM_data[[1]], RCM_data[[i]], RCM_data[[n]], method = "pearson", permutations = 999)
    mantel_output[i-1,n-1] <- x$statistic
  }}

# View(mantel_output) # need to change diagonal values to zero 
for (i in 1:12){
  mantel_output[i,i] <- 0
}
mantel_output

### Now need to generate matrix of difference between first mantel and reciprocal mantel (upper vs lower triangle)
heatmap_value_table <- matrix(ncol=variables, nrow=variables)
rownames(heatmap_value_table) <- c("euclidean_dist",
                                   "habitat_least_cost_dist",
                                   "habitat_resistance_dist",
                                   "roads_least_cost_dist",
                                   "roads_resistance_dist",
                                   "traffic_least_cost_dist",
                                   "traffic_resistance_dist",
                                   "continentality_dist",
                                   "temp_summer_dist", 
                                   "temp_winter_dist",
                                   "precip_summer_dist",
                                   "precip_winter_dist")
colnames(heatmap_value_table) <- rownames(heatmap_value_table)

for (i in 1:12){
  for (n in 1:12){
    x <- mantel_output[i,n] - mantel_output[n,i]
    heatmap_value_table[n,i] <- x
  }}

# View(heatmap_value_table) 
levelplot(heatmap_value_table)
RCM_relative_support <- cbind.data.frame(rownames(heatmap_value_table), as.numeric(colSums(heatmap_value_table)))
RCM_relative_support
colnames(RCM_relative_support) <- c("distance_measure", "relative_support")
RCM_relative_support <- RCM_relative_support[order(RCM_relative_support$relative_support, decreasing=TRUE),]
RCM_relative_support
write.csv(RCM_relative_support, file = "~/RCM_relative_support_ind.csv", quote = F, row.names = F, col.names = T)

pdf(file = "~/RCM_heatmap_ind.pdf", width = 10, height = 10, family = "serif")
library("pals")
p1 <- coolwarm(255)
levelplot(heatmap_value_table, 
          col.regions=p1,
          colorkey=list(labels=list(cex=1.25)),
          scales=list(x=list(at=seq(12:1), rot=45, cex=1.4),
                      y=list(at=seq(12:1), cex=1.4),
                      labels=c("Euclidean dist.",
                               "Habitat LCP dist.",
                               "Habitat resistance dist.", 
                               "Roads LCP dist.",
                               "Roads resistance dist.",
                               "Traffic LCP dist.",
                               "Traffic resistance dist.",
                               "Continentality dist.",
                               "Summer temp. dist.", 
                               "Winter temp. dist.",
                               "Summer precip. dist.",
                               "Winter precip. dist."),
                      tck = c(1,0),
                      col = "black"),
          par.settings = list(axis.line = list(col = "transparent")),
          # xlab=list('Alternative variable', cex=1.5), ylab=list('Focal variable', cex=1.5),
          xlab=NULL, ylab=NULL, 
          ylim=c(12.75,0.5))#, xlim=c(12.5,0.5))

dev.off()
