library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v07.RData")

mydatasighting <- read.csv("mydatasighting.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)

# Analyse data for 2011:
id <- mydata$yearnew == 2011
table(id)
mydata <- mydata[id,]
dim(mydata)
names(mydata)


# 354 rows

# Create a new variable called "caller". If home phone number is missing, use the work number. If work number is missing we use the caller name:
# mydata$caller <- ifelse(mydata$homephone == " ", as.character(mydata$workphone), as.character(mydata$homephone))
# If either home or work phone number is missing, use the callername:
# mydata$caller <- ifelse(mydata$caller == " ", as.character(mydata$callername), as.character(mydata$caller))
# sort(table(mydata$caller))

# Drop those records without a valid caller identifier (n = 16):
# id <- mydata$caller != " "
# mydata <- mydata[id,]

# Create a unique 1 to n caller identifier:
# tcaller <- unique(mydata$caller)
# tid <- 1:length(tcaller)
# lookup <- data.frame(id = tid, caller = tcaller)

# Add the unique caller identifier to the mydata data frame:
# mydata$uniq <- lookup$id[match(mydata$caller, lookup$caller)]
# sort(table(mydata$uniq))

# Drop the multiple reporters:
# id <- mydata$uniq != 182 & mydata$uniq != 19 & mydata$uniq != 64 & mydata$uniq != 65 & mydata$uniq != 9 & mydata$uniq != 62 & mydata$uniq != 23 & mydata$uniq != 18 & mydata$uniq != 88 
# mydata <- mydata[id,]


# Bring in MGA56 study area boundary map:
unzip("vector\\AU_Qld_study_area.zip")
studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 
plot(studyarea.shp)
id <- mydata$lat > -30 # remove record with wrong coordinate
mydata <- mydata[id,]
head(mydata)
# mydata in lat-lon:
coords <- SpatialPoints(mydata[, c("lng", "lat")])
mydata.ll <- SpatialPointsDataFrame(coords, mydata)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))
plot(mydata.mga)
# Add the MGA coordinates to mydata:
mydata$x <- coordinates(mydata.mga)[,1]
mydata$y <- coordinates(mydata.mga)[,2]

windows(); plot(mydata.mga, cex = 0.3, col = 'blue', pch = 15)

# Create a 3 km by 3 km grid:
mydata.r <- raster(studyarea.shp)
res(mydata.r) <- 1000
mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset
#writeRaster(mydata.r,"testing.asc")
 #windows(); plot(mydata.r)
# create the disance matrix
#####Go to 109 for grid based sampling. START Distaance based selection using the function
# Function to estimate the distance to the nearest point of the same vector
source("Lib_DistEstimatesToItself.r")
# start distance based selection method
mydata.n <- mydata.mga[c(3,4)]
mydata.n = as.data.frame(mydata.n)
myInPtDF = SpatialPointsDataFrame(mydata.n[c("x","y")], mydata.n)
windows();plot(myInPtDF, axes = T)
myInPtDF$disttoitself = Lib_DistEstimatesToItself(myInPtDF$x, myInPtDF$y)
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 500)
windows();plot(myInPtDF2, axes = T)
acsel = myInPtDF[,1:2]
acsel$all = 1
myFD2 = myInPtDF2[,1:2]
myFD2$all = 0
myFD = rbind(acsel,myFD2)
windows();plot(myFD, axes = T,col=c("red","blue")[myFD$all+1])
myFD3<-as.data.frame(myFD)
#########END of distance based selection#####
# Sample points from within each grid cell:
acsel <- gridSample(mydata.mga, mydata.r, n = 1)
mydata.p <- rasterToPolygons(mydata.r)

windows(); plot(mydata.p, border = 'gray')
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(mydata.mga)
# Selected points in red:
points(acsel, cex = 0.2, col = 'red', pch = 15)

dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

windows(); plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(x = mydata$x, y = mydata$y,cex = 0.2, col = 'blue', pch = 15)

# kdate as a date:
# dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
# dat$kmonth <- format(dat$kdate, format = "%m")

# =================================================================================================================================
#1.total phosphorous:% Read in explanatory variable raster maps - 

tpo.r <- raster("raster\\AU_Qld_soil_total_phosphorous_000_005cm-wgs84.tif")
projection(tpo.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
tpo.r <- projectRaster(from = tpo.r, crs = CRS("+init=epsg:28356"))
windows(); plot(tpo.r, main="tpo")
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#chnaged reolution 34*34 to 100 *102
r <- raster(ncol = 100, nrow = 100)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
tpo.crop <- crop(x = tpo.r, y = r, snap = "near")
windows(); plot(tpo.crop)
# Resample:okay. lets keep this resolution for all rasters. why do we need to change this match r which is our grid 
tpo.r <- resample(x =tpo.crop , y = r)
writeRaster(tpo.r,"raster_crop\\tpo.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(tpo.r, main="tpo")
points(acsel, cex = 0.4, col = 'red', pch = 15)

#new dem file in tif formt("raster\\AU_Qld_dem-MGA65.tif")
#new rainfall in tif formt ("raster\\AU_Qld_rainfall-MGA65.tif") #created from ascii to raster tool in archgis.
#new temp file in tif formt <- raster("raster\\AU_Qld_temp-MGA65.tif")
##new hpop file("AU_Qld_hpop-wgs-84")
##new "AU_Qld_FPC-MGA65.tif"

#2.
awc.r<- raster("raster\\AU_Qld_soil_available_water_capacity_005_015cm-wgs84.tif")
projection(awc.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
awc.r <- projectRaster(from = awc.r, crs = CRS("+init=epsg:28356"))
windows(); plot(awc.r)
awc.crop <- crop(x = awc.r, y = r, snap = "near")
windows(); plot(awc.crop)
# Resample:
awc.r <- resample(x = awc.crop, y = r)
projection(awc.r) <- CRS("+init=epsg:28356")
writeRaster(awc.r,"raster_crop\\sbd.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(awc.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
points(mydata.mga)
#rhohat for thhis not calculated
-------------------------------------------------------------------------------------------------------------------------------
#3. Foliage projective cover (fpc) Read in explanatory variable raster maps -
fpc.r<- raster("raster\\AU_Qld_fcp_2014-wgs-84.tif")
projection(fpc.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
fpc.r <- projectRaster(from = fpc.r, crs = CRS("+init=epsg:28356"))
windows(); plot(fpc.r)
fpc.crop <- crop(x = fpc.r, y = r, snap = "near")
windows(); plot(fpc.crop)
# Resample:
fpc.r <- resample(x = fpc.crop, y = r)
projection(fpc.r) <- CRS("+init=epsg:28356")
writeRaster(fpc.r,"raster_crop\\fpc.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(fpc.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
points(mydata.mga)

# ---------------------------------------------------------------------------------------------------------------------------------
#2.Soil bulk density:g/cm3

sbd.r <- raster("raster\\AU_Qld_soil_bulk_density_000_005cm-wgs84.tif")
projection(sbd.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
sbd.r <- projectRaster(from = sbd.r, crs = CRS("+init=epsg:28356"))
windows(); plot(sbd.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
sbd.crop <- crop(x = sbd.r, y = r, snap = "near")
windows(); plot(sbd.crop)

# Resample:
sbd.r <- resample(x = sbd.crop, y = r)
projection(sbd.r) <- CRS("+init=epsg:28356")
writeRaster(sbd.r,"raster_crop\\sbd.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(sbd.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
points(mydata.mga)


# ---------------------------------------------------------------------------------------------------------------------------------
#3. Soil clay:%
fig.clay <- function() {
clay.r <- raster("raster\\AU_Qld_soil_clay_000_005cm-wgs84.tif")
projection(clay.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
clay.r <- projectRaster(from = clay.r, crs = CRS("+init=epsg:28356"))
windows(); plot(clay.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
#mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
clay.crop <- crop(x = clay.r, y = r, snap = "near")
windows(); plot(clay.crop)

# Resample:
clay.r <- resample(x = clay.crop, y = r)
projection(clay.r) <- CRS("+init=epsg:28356")
writeRaster(clay.r,"raster_crop\\clay.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(clay.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
}
fig.clay()

#fig.clay()
#pdf("clay.pdf", width=6, height=8)
#fig.clay()
#dev.off()
# ---------------------------------------------------------------------------------------------------------------------------------
#4. Soil water capacity:
fig.water <- function() {
  
water.r <- raster("raster\\AU_Qld_soil_available_water_capacity_000_005cm-wgs84.tif")
projection(water.r) <- CRS("+init=epsg:4326")

# Reproject to MGA-56:
water.r <- projectRaster(from = water.r, crs = CRS("+init=epsg:28356"))
windows(); plot(water.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
water.crop <- crop(x = water.r, y = r, snap = "near")
windows(); plot(water.crop)

# Resample:
water.r <- resample(x = water.crop, y = r)
projection(water.r) <- CRS("+init=epsg:28356")
writeRaster(water.r,"raster_crop\\water.TIF",overwrite=TRUE)
points(acsel, cex = 0.4, col = 'red', pch = 15)
# Plot to check:
windows(); plot(water.r)
plot(mydata.p, add = TRUE)
}
fig.water()
# ---------------------------------------------------------------------------------------------------------------------------------
#5. Soil nitrogen:%
fig.nitro <- function() {
nitro.r <- raster("raster\\AU_Qld_soil_nitrogen_000_005cm-wgs84.tif")
projection(nitro.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
nitro.r <- projectRaster(from = nitro.r, crs = CRS("+init=epsg:28356"))
windows(); plot(nitro.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
nitro.crop <- crop(x = nitro.r, y = r, snap = "near")
windows(); plot(nitro.crop)

# Resample:
nitro.r <- resample(x = nitro.crop, y = r)
projection(nitro.r) <- CRS("+init=epsg:28356")
writeRaster(nitro.r,"raster_crop\\nitro.TIF",overwrite=TRUE)

# Plot to check:
windows(); plot(nitro.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
}
fig.nitro()
# ---------------------------------------------------------------------------------------------------------------------------------
#6. Road network (expressed as the number of metres of road per square km grid):
fig.roaden <- function() {
unzip("vector\\AU_Qld_road_length.zip")

roaden.shp <- readShapePoly("AU_Qld_road_length-MGA65.shp")

proj4string(roaden.shp) <- CRS("+init=epsg:28356") 

# Vector to raster conversion:

roaden.r <- rasterize(x = roaden.shp, y = mydata.r, field = "Sum_lengh")
windows(); plot(roaden.r2)
windows(); plot(roaden.r)
plot(mydata.p, add = TRUE)
#r <- raster(ncol = 100, nrow = 100)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
roaden.crop <- crop(x = roaden.r, y = r, snap = "near")
windows(); plot(roaden.crop)

# Resample:
roaden.r <- resample(x = roaden.crop, y = r)
projection(roaden.r) <- CRS("+init=epsg:28356")
writeRaster(roaden.r,"raster_crop\\roaden.TIF",overwrite=TRUE)
# Plot to check:
windows(); plot(roaden.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
#plot(mydata.p, add = TRUE)
writeRaster(roaden.r,"roaden.TIF", overwrite=TRUE)
}
fig.roaden()
# ---------------------------------------------------------------------------------------------------------------------------------
#7. Human population density
fig.hpop <- function() {
hpop.r <- raster("raster\\AU_Qld_hpop-wgs84.tif")
projection(hpop.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
hpop.r <- projectRaster(from = hpop.r, crs = CRS("+init=epsg:28356"))
windows(); plot(hpop.r)

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
hpop.crop <- crop(x = hpop.r, y = r, snap = "near")
windows(); plot(hpop.crop)s

# Resample:
hpop.r <- resample(x = hpop.crop, y = r)
projection(hpop.r) <- CRS("+init=epsg:28356")
writeRaster(hpop.r,"raster_crop\\hpop.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(hpop.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
}
fig.hpop()
# ---------------------------------------------------------------------------------------------------------------------------------
#8. Rainfall
fig.rainfall <- function() {
rain.r <- raster("raster\\AU_Qld_rainfallmeanannual-MGA56.TIF")
# Reproject to MGA-56:
windows(); plot(rain.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
rain.crop <- crop(x = rain.r, y = r, snap = "near")
windows(); plot(rain.crop)

# Resample:
rain.r <- resample(x = rain.crop, y = r)
projection(rain.r) <- CRS("+init=epsg:28356")
writeRaster(rain.r,"raster_crop\\rain.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(rain.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
}
fig.rainfall()

# ---------------------------------------------------------------------------------------------------------------------------------
#9. Temperature:http://www.bom.gov.au/jsp/ncc/climate_averages/temperature/index.jsp
fig.temp <- function() {
temp.r <- raster("raster\\AU_Qld_temp-MGA65.tif")

projection(temp.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
temp.r <- projectRaster(from = temp.r, crs = CRS("+init=epsg:28356"))

windows(); plot(temp.r)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 52, nrow = 56)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
temp.crop <- crop(x = temp.r, y = r, snap = "near")
windows(); plot(temp.crop)

# Resample:
temp.r <- resample(x = temp.crop, y = r)
projection(temp.r) <- CRS("+init=epsg:28356")
#contour(temp.r)
#contourplot(temp.r)
# Plot to check:
windows(); plot(temp.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
writeRaster(temp.r,"raster_crop\\temp.TIF",overwrite=TRUE)
}
fig.temp()
# ---------------------------------------------------------------------------------------------------------------------------------
#10. Elevation
fig.elev <- function() {
elev.r <- raster("raster\\AU_Qld_dem-MGA65.TIF")
projection(elev.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
elev.r<- projectRaster(from = elev.r, crs = CRS("+init=epsg:28356"))
windows(); plot(elev.r, axes = TRUE)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)
# Crop the raster to the shape file spatial extent:
elev.crop <- crop(x = elev.r, y = r, snap = "near")
windows(); plot(elev.crop)
# Resample:
elev.r <- resample(x = elev.crop, y = r)
projection(elev.r) <- CRS("+init=epsg:28356")
writeRaster(elev.r,"raster_crop\\elev.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(elev.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
}
fig.elev()
---------------------------------------------------------------------------------
library(rasterVis)  #https://pakillo.github.io/R-GIS-tutorial/#plotraster
#https://geoscripting-wur.github.io/AdvancedRasterAnalysis/
#load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_CAR_model_v01neditedRavi.RData")


histogram(aspect.r)
density(elev.r)
density(stack.n)
persp(slope.r)
#contour(stack.r)
contourplot(roaden.r)
#levelplot(topo.r)
#plot3D(aspect.r)
#library(rgl)
#plot3D(elev.r,col=rainbow)#
#elev.r.Moran <- MoranLocal(elev.r)
#plot(elev.r.Moran)
#ccr<- terrain(elev.r, opt = c("slope", "aspect"), unit = "degrees")
#slope.r<- terrain(elev.r, opt = c("slope"), unit = "degrees")
#aspect.r<- terrain(elev.r, opt = c("aspect"), unit = "degrees")
#writeRaster(slope.r,"raster\\slope.TIF", overwrite=TRUE)
#writeRaster(aspect.r,"raster\\aspect.TIF", overwrite=TRUE)
#windows (); plot(slope.r)
#windows (); plot(aspect.r)
#contourplot(aspect.r)
#head(elev.r) 
#points(acsel, cex = 0.4, col = 'red', pch = 15)
#plot(mydata.p, add = TRUE)
#-------------------------------------------------------------------------------------------
#The Topographic Position Index (TPI)#Terrain Attributes twi,tpi,slope,aspect
#compares the elevation of each cell in a
#DEM to the mean elevation of a specified neighborhood around that cell.
#Positive TPI values represent locations that are higher than the average of their
#surroundings, as defined by the neighborhood (ridges). Negative TPI values
#represent locations that are lower than their surroundings (valleys). TPI values
#near zero are either flat areas (where the slope is near zero) or areas of constant
#slope (where the slope of the point is significantly greater than zero).
#SRTM_TopographicPositionIndex_141
#11.Topographic Position Index (TPI)
fig.topo <- function()
topo.r <- raster("raster\\AU_Qld_TopographicPositionIndex-wgs-84.TIF")
projection(topo.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
topo.r <- projectRaster(from = topo.r, crs = CRS("+init=epsg:28356"))
windows(); plot(topo.r, axes = TRUE)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)
# Crop the raster to the shape file spatial extent:
topo.crop <- crop(x = topo.r, y = r, snap = "near")
windows(); plot(topo.crop)
# Resample:
topo.r <- resample(x = topo.crop, y = r)
writeRaster(topo.r,"raster_crop\\topo.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(topo.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)



fig.topo()
------------------
#12.slope:percentage change in the elevation calcuated by difference of elevation of two points diveded by lateral distance.
fig.slope <- function() {
  
#slope.r <- raster("raster\\slope.TIF") 
#slope.r <- raster("raster\\Slope percent_61.TIF")
slope.r <- raster("raster\\AU_Qld_slope_reliefclass81-wgs-84.TIF")
projection(slope.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
slope.r <- projectRaster(from = slope.r, crs = CRS("+init=epsg:28356"))
windows(); plot(slope.r, axes = TRUE)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)
# Crop the raster to the shape file spatial extent:
slope.crop <- crop(x = slope.r, y = r, snap = "near")
windows(); plot(slope.crop)
# Resample:
slope.r <- resample(x = slope.crop, y = r)
writeRaster(slope.r,"raster_crop\\slope.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(slope.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)

}
fig.slope()
#------------------------------------------------------------------------------------------
#13. 
fig.aspect <- function() {
aspect.r <- raster("raster\\AU_Qld_aspect91-wgs-84.TIF")
projection(aspect.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
aspect.r <- projectRaster(from = aspect.r, crs = CRS("+init=epsg:28356"))
windows(); plot(aspect.r, axes = TRUE)
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
  #r <- raster(ncol = 33, nrow = 33)
  #extent(r) <- extent(mydata.r)
  # Crop the raster to the shape file spatial extent:
  aspect.crop <- crop(x = aspect.r, y = r, snap = "near")
  windows(); plot(aspect.crop)
  # Resample:
  aspect.r <- resample(x = aspect.crop, y = r)
  writeRaster(aspect.r,"raster_crop\\aspect.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(aspect.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)
  
}
fig.aspect()
-------------------------------------------------------------
#14.Topographic_Wetness_Index  is a steady state wetness index. It is commonly used to quantify 
#topographic control on hydrological processes  
 fig.twi <- function() {
  
  twi.r <- raster("raster\\Topographic_Wetness_Index-wgs84.tif")
  projection(twi.r) <- CRS("+init=epsg:4326")
  # Reproject to MGA-56:
  twi.r <- projectRaster(from = twi.r, crs = CRS("+init=epsg:28356"))
  windows(); plot(twi.r, axes = TRUE)
  # Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
  mydata.r
  #r <- raster(ncol = 33, nrow = 33)
  #extent(r) <- extent(mydata.r)
  # Crop the raster to the shape file spatial extent:
  twi.crop <- crop(x = twi.r, y = r, snap = "near")
  windows(); plot(twi.crop)
  # Resample:
  twi.r <- resample(x = twi.crop, y = r)
  writeRaster( twi.r,"raster_crop\\ twi.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(twi.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)
}
fig.aspect()  


# ---------------------------------------------------------------------------------------------------------------------------------
#15. Habitat

unzip("vector\\AU_Qld_koala_habitat_suitability.zip")
habit.shp <- readShapePoly("AU_Qld_koala_habitat_suitability_MGA56.shp")
proj4string(habit.shp) <- CRS("+init=epsg:28356") 

# Vector to raster conversion:
habit.r <- rasterize(x = habit.shp, y = mydata.r, field = "suitabilit")

windows(); plot(habit.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
mydata.r
#r <- raster(ncol = 33, nrow = 33)
#extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
habit.crop <- crop(x = habit.r, y = r, snap = "near")
windows(); plot(habit.crop)

# Resample:
habit.r <- resample(x = habit.crop, y = r)
projection(habit.r) <- CRS("+init=epsg:28356")

writeRaster(habit.r,"raster_crop\\habit.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(habit.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)

# =================================================================================================================================
# Make a ppp object for spatstat:

# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:
source("owin2sp_source.r")

dat.w <- as(as(studyarea.shp, "SpatialPolygons"), "owin")
dat.spw <- owin2SP(dat.w)

# Set projection of the owin as GDA94 / SA Lambert:
proj4string(dat.spw)
proj4string(dat.spw) <- CRS("+init=epsg:28356") 
dat.spw <- raster::intersect(x = dstudyarea.shp, y = dat.spw)

# Convert the sp object into a window:
dat.w <- as.owin(dat.spw)

# set the dimensions for the grid and the bandwidth:
dimyx = c(200,200)
bandwidth <- 2000

# Select only those koalas within the owin:
# id <- inside.owin(x = mydata$x, y = mydata$y, w = dat.w)
# mydata <- mydata[id,]

windows(); plot(studyarea.shp, axes = TRUE)
points(x = acsel[,1], y = acsel[,2])

# Make a ppp object:
dat.ppp <- ppp(x = acsel[,1], y = acsel[,2], window = dat.w)

windows(); plot(dat.ppp, axes = TRUE)

# Work out density for dat.ppp (using Diggle's edge correction):
dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)

# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre:
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))

windows(); plot(dat.den)

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 0.01, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
points(acsel, cex = 0.4, col = 'red', pch = 15)

## rhohat - total phosphorous:Computes a smoothing estimate of the intensity of a point process, as a function of a (continuous) spatial covariate.
#Foliage projective cover
# Convert tpo.r into a spatstat image object:
fpc.im <- as.im(fpc.r)
windows(); plot(fpc.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
fpc.rho <- rhohat(object = dat.ppp, covariate = fpc.im)

windows(); plot(fpc.rho, xlab = "FPC (units)", main = "")
savePlot(filename = "fpc",type = c( "png"),device = dev.cur())


# =================================================================================================================================
twi.im <- as.im(twi.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
twi.rho <- rhohat(object = dat.ppp, covariate = twi.im)

windows(); plot(twi.rho, xlab = "twi (units)", main = "")
savePlot(filename = "twi",type = c( "png"),device = dev.cur())
-------
# Convert tpo.r into a spatstat image object:
aspect.im <- as.im(aspect.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
aspect.rho <- rhohat(object = dat.ppp, covariate = aspect.im)

windows(); plot(aspect.rho, xlab = "aspect (units)", main = "")
savePlot(filename = "aspect",type = c( "png"),device = dev.cur())

-----
# rhohat - total phosphorous:Computes a smoothing estimate of the intensity of a point process, as a function of a (continuous) spatial covariate.

# Convert tpo.r into a spatstat image object:
tpo.im <- as.im(tpo.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
tpo.rho <- rhohat(object = dat.ppp, covariate = tpo.im)

windows(); plot(tpo.rho, xlab = "Total phosphorous (units)", main = "")
savePlot(filename = "tpo",type = c( "png"),device = dev.cur())
------------------------------------------------------------------------------------------------------------------------------
# rhohat - total nitrogen:
  
# Convert tpo.r into a spatstat image object:
nitro.im <- as.im(nitro.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total Nitrogen?
nitro.rho <- rhohat(object = dat.ppp, covariate = nitro.im)

windows(); plot(nitro.rho, xlab = "Total soil nitrogen (units)", main = "")
savePlot(filename = "nitro",type = c( "png"),device = dev.cur())

# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat - soil bulk density:

# Convert tpo.r into a spatstat image object:
sbd.im <- as.im(sbd.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations a soil bulk density?
sbd.rho <- rhohat(object = dat.ppp, covariate = sbd.im)

windows(); plot(sbd.rho, xlab = "Soil bulk density (units)", main = "")
savePlot(filename = "sbd",type = c( "png"),device = dev.cur())


-----------------------------------------------------------------------------------------------------------------------------------
# rhohat - clay:

# Convert tpo.r into a spatstat image object:
clay.im <- as.im(clay.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and clay?
clay.rho <- rhohat(object = dat.ppp, covariate = clay.im)

windows(); plot(clay.rho, xlab = "clay (units)", main = "")
savePlot(filename = "clay",type = c( "png"),device = dev.cur())



------------------------------------------------------------------------------------------------------------------------------
# rhohat - elevation:

# Convert .ph.r into a spatstat image object:
elev.im <- as.im(elev.r)

windows(); plot(elev.im, axes = TRUE)

# What is the nature of the association between koala sight locations and elevation?
elev.rho <- rhohat(object = dat.ppp, covariate = elev.im)

windows(); plot(elev.rho, xlim = c(0, 600), xlab = "Elevation (metres)", main = "")
savePlot(filename = "elev",type = c( "png"),device = dev.cur())

# Most koalas sighted at elevations of less than 100 metres. Plot raster map to show areas less than 100 m:
elev100.r <- elev.r < 100

windows(); plot(elev100.r, axes = TRUE)
points(x = mydata$x, y = mydata$y)


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - human population density:

# Convert .ph.r into a spatstat image object:
hpop.im <- as.im(hpop.r)

windows(); plot(hpop.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)

windows(); plot(hpop.rho, xlim = c(0, 3000), ylim = c(0, 0.2e-06), xlab = "Human population density", main = "")
savePlot(filename = "hpop",type = c( "png"),device = dev.cur())



# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - annual rainfall (mm):

# Convert .ph.r into a spatstat image object:chnaged rainfall.r to rain.r
rainfall.im <- as.im(rain.r)

windows(); plot(rainfall.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
rainfall.rho <- rhohat(object = dat.ppp, covariate = rainfall.im)

windows(); plot(rainfall.rho, xlab = "Annual rainfall (mm)", main = "")
savePlot(filename = "rainfall",type = c( "png"),device = dev.cur())

# Strong association with rainfall between 1200 and 1400 mm per year:
crainfall.r <- rainfall.r >1200 & rainfall.r < 1400
windows(); plot(crainfall.r)

crainfall.im <- as.im(crainfall.r)

# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - road density: This was catogorsed now in to groups.

roaden.im <- as.im(roaden.r)

windows(); plot(roaden.im, axes = TRUE)

roaden.rho <- rhohat(object = dat.ppp, covariate = roaden.im)

windows(); plot(roaden.rho, xlab = "Road density (metres per square km)", main = "")
savePlot(filename = "roaden",type = c( "png"),device = dev.cur())


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - habitat:

unzip("vector\\AU_Qld_koala_habitat_suitability.zip")
habit.shp <- readShapePoly("AU_Qld_koala_habitat_suitability_MGA56.shp")
proj4string(habit.shp) <- CRS("+init=epsg:28356") 
windows(); plot(habit.shp)
table(habit.shp$suitabilit)

# Take habitat suitability area 1 only:
id <- habit.shp$suitabilit == 1; habit1.shp <- habit.shp[id,]
habit1.w <- as(as(habit1.shp, "SpatialPolygons"), "owin")
habit1.im <- distmap(habit1.w)
# Take habitat suitability area 2 only:
id <- habit.shp$suitabilit == 2; habit2.shp <- habit.shp[id,]
habit2.w <- as(as(habit2.shp, "SpatialPolygons"), "owin")
habit2.im <- distmap(habit2.w)

# Take habitat suitability area 3 only:
id <- habit.shp$suitabilit == 3; habit3.shp <- habit.shp[id,]
habit3.w <- as(as(habit3.shp, "SpatialPolygons"), "owin")
habit3.im <- distmap(habit3.w)
windows(): plot(habit)

# Turn habit1.im into a raster and re-sample to match the dimensions of the other rasters:
habit1.r <- raster(habit1.im)
#r <- raster(ncol = 52, nrow = 56); extent(r) <- extent(mydata.r)
habit1.crop <- crop(x = habit1.r, y = r, snap = "near")
habit1.r <- resample(x = habit1.crop, y = r)
projection(habit1.r) <- CRS("+init=epsg:28356")
writeRaster(habit1.r,"raster_crop\\habit1.TIF", overwrite=TRUE)
habit01.im <- as.im(habit1.r)
windows(); plot(habit01.im)
plot(mydata.p, add = TRUE)

# Turn habit2.im into a raster and re-sample to match the dimensions of the other rasters:
habit2.r <- raster(habit2.im)
#r <- raster(ncol = 52, nrow = 56); extent(r) <- extent(mydata.r)
habit2.crop <- crop(x = habit2.r, y = r, snap = "near")
habit2.r <- resample(x = habit2.crop, y = r)
projection(habit2.r) <- CRS("+init=epsg:28356")
writeRaster(habit2.r,"raster_crop\\habit2.TIF", overwrite=TRUE)
habit02.im <- as.im(habit2.r)
windows(); plot(habit02.im)
plot(mydata.p, add = TRUE)

# Turn habit3.im into a raster and re-sample to match the dimensions of the other rasters:
habit3.r <- raster(habit3.im)
#r <- raster(ncol = 52, nrow = 56); extent(r) <- extent(mydata.r)
habit3.crop <- crop(x = habit3.r, y = r, snap = "near")
habit3.r <- resample(x = habit3.crop, y = r)
projection(habit3.r) <- CRS("+init=epsg:28356")
writeRaster(habit3.r,"raster_crop\\habit3.TIF", overwrite=TRUE)
habit03.im <- as.im(habit3.r)
windows(); plot(habit03.im)
plot(mydata.p, add = TRUE)
---------------#rhohat suitability
  #habitat suitability"suitable catogory1"
  habit01.im <- as.im(habit1.r)
windows(); plot(habit01.im, axes = TRUE)
habit01.im.rho <- rhohat(object = dat.ppp, covariate = habit01.im)
windows(); plot(habit01.im.rho, xlab = "distance to habitat suitability low suitable", main = "")
savePlot(filename = "habit01",type = c( "png"),device = dev.cur())

--------------------------------------
  #habitat suitability"suitable catogory2"
  habit02.im <- as.im(habit2.r)
windows(); plot(habit02.im, axes = TRUE)
habit02.im.rho <- rhohat(object = dat.ppp, covariate = habit02.im)
windows(); plot(habit02.im.rho, xlab = "distance to habitat suitability /suitable", main = "")
savePlot(filename = "habit02",type = c( "png"),device = dev.cur())

savePlot(filename = "habit02",type = c( "png"),device = dev.cur())
-----------------------------------------------------------------------
  #habitat suitability"suitable catogory3"
  habit03.im <- as.im(habit3.r)
windows(); plot(habit03.im, axes = TRUE)
habit03.im.rho <- rhohat(object = dat.ppp, covariate = habit03.im)
windows(); plot(habit03.im.rho, xlab = "distance to habitat suitability /high suitable", main = "")
savePlot(filename = "habit03",type = c( "png"),device = dev.cur())

# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - water:

water.im <- as.im(water.r)
windows(); plot(water.im, axes = TRUE)

water.rho <- rhohat(object = dat.ppp, covariate = water.im)
windows(); plot(water.rho, xlab = "Available water capacity (units)", main = "")
savePlot(filename = "water",type = c( "png"),device = dev.cur())


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - temp:

temp.im <- as.im(temp.r)
windows(); plot(temp.im, axes = TRUE)

temp.rho <- rhohat(object = dat.ppp, covariate = temp.im)
windows(); plot(temp.rho, xlab = "Average annual temperature (degrees Celcius)", main = "")
savePlot(filename = "temp",type = c( "png"),device = dev.cur())

save.image("AU_Qld_koala_habitat_v08editedRavi.RData")


# =================================================================================================================================
# Poisson point process model:

# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409
plot(dat.ppp)

# Backward stepwise variable selection:
dat.ppm01 <- ppm(dat.ppp, trend = ~ tpo + sbd + water + nitro + roaden + hpop + temp + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
# Values of the covariates 'roaden', 'temp' were NA or undefined at 0.5% (5 out of 1005) of the quadrature points. 
# See http://gis.stackexchange.com/questions/161893/creating-polygonal-windows-for-spatial-analysis-in-r for a work around.

# Extract the quadrature scheme from dat.ppm01:
Qz <- quad.ppm(dat.ppm01, drop = TRUE)
# This step will drop points within the quadrature scheme that had NA-values.

# Run the model again, but this time, use the corrected quadrature scheme:
dat.ppm02 <- ppm(Qz, trend = ~ tpo + sbd + water + nitro + roaden + hpop + temp + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm02)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm02)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm02) # 3413.045
step(dat.ppm02)
# Drop tpo:temp +
dat.ppm03 <- ppm(Qz, trend = ~ sbd +tpo+ water + nitro + roaden +  hpop + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm03)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm03)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm03) # 3411.272


# Drop habitat03:tpo
dat.ppm04 <- ppm(Qz, trend = ~ sbd + roaden+water + nitro + habit03+ hpop + elev + habit01 + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm04)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm04)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm04) # 3460.872


# Drop habit1:
dat.ppm05 <- ppm(Qz, trend = ~ sbd + roaden+water + nitro + habit03+ hpop + elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm05)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm05)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm05) # 3459.164


# Drop roadden:
dat.ppm06 <- ppm(Qz, trend = ~ sbd + water + nitro + habit03+ hpop + elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm06)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm06)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm06) # 3457.573


# Drop hpop:
dat.ppm07 <- ppm(Qz, trend = ~ sbd + water + nitro + habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm07)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm07)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm07) # 3455.898


# Drop nitro:
dat.ppm08 <- ppm(Qz, trend = ~ sbd + water +  habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm08)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm08)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm08) # 3455.163


# Drop water
dat.ppm09 <- ppm(Qz, trend = ~ sbd + habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm09)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm09)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm09) # 3453.953


# Drop habit03
dat.ppm10 <- ppm(Qz, trend = ~ sbd +  elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm10)$coefs.SE.CI
AIC(dat.ppm10) # 3411.58

windows(); plot(dat.ppm10, how = "image", se = TRUE)

# The lurking variable plot gives you an idea of mis-specification of spatial trend:
windows(); diagnose.ppm(dat.ppm10, type = "raw")
savePlot(filename = "datppm10",type = c( "png"),device = dev.cur())

# Now add a Strauss spatial interaction term:
dat.ppm11 <- ppm(Qz, trend = ~ sbd +  elev  + habit02, interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm11)$coefs.SE.CI
AIC(dat.ppm11) # 2459.733

windows(); plot(dat.ppm11, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm11, type = "raw")


# Drop +  elev:
dat.ppm12 <- ppm(Qz, trend = ~ sbd   + habit02, interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm12)$coefs.SE.CI
AIC(dat.ppm12) # 2459.799

windows(); plot(dat.ppm12, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm12, type = "raw")


# Drop + habit02
dat.ppm13 <- ppm(Qz, trend = ~ sbd   , interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm13)$coefs.SE.CI
AIC(dat.ppm13) # 2459.806

windows(); plot(dat.ppm13, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm13, type = "raw")


# QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
# dat.qq11 <- qqplot.ppm(fit = dat.ppm11, nsim = 10)

# windows(); plot(dat.qq11, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# =================================================================================================================================
# Model with spatial interaction term:
dat.pred13 <- predict(dat.ppm13); pred <- dat.pred13

# Express predicted koala intensity in terms of the number of koalas per hectare:
pred$v <- pred$v * 1E04
summary(as.vector(pred$v))
windows(); hist(as.vector(pred$v))



# Predictions from the model:
breaks <- seq(from = 0, to = 0.00013, length = 6)
col <- brewer.pal(n = 5, name = "Reds")

windows(); plot(x = xlim, y = ylim, main="Predictions from the model", type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "blue")

plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Predictions from the model",type = c( "png"),device = dev.cur())

# Actual kernel smoothed data:
breaks <- seq(from = 0, to = 0.0015, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n",main="Actual kernel smoothed data", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Actual kernel smoothed data",type = c( "png"),device = dev.cur())


# =================================================================================================================================


# Save the workspace because these analyses take forever to run:
save.image("AU_Qld_koala_habitat_v08.RData")

xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 100, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.sden$xcol, y = dat.sden$yrow, z = t(dat.sden$v), col = col, breaks = breaks, add = TRUE)
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


# =================================================================================================================================
