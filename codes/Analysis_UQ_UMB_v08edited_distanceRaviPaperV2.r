library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v08edited_distance_Ravi.RData")

mydatasighting <- read.csv("sightingdata\\mydatasighting_cleaned.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)
mydata <- subset(mydata, yearnew >2010) #select 1997-1999/ then yearnew>=2000|yearnew<2003

head(mydata)
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))

#analyse data for 1997-1999

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
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 200)
windows();plot(myInPtDF2, axes = T)
acsel= myInPtDF2[,1:2]
windows();plot(acsel,axes = T)
#acsel$all = 1
#myFD2 = myInPtDF2[,1:2]
#myFD2$all = 0
#myFD = rbind(acsel,myFD2)
#windows();plot(myFD, axes = T,col=c("red","blue")[myFD$all+1])
#myFD3<-as.data.frame(myFD)
#########END of distance based selection#####
# Sample points from within each grid cell:
acsel <- gridSample(acsel, mydata.r, n = 1)
mydata.p <- rasterToPolygons(mydata.r)
windows(); plot(mydata.p, border = 'gray')
# Selected points in red:
points(acsel, cex = 0.5, col = 'red', pch = 15)

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
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 100
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 0.009, length = 4)
col <- brewer.pal(n = 3, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
points(acsel, cex = 0.4, col = 'red', pch = 15)

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
projection(tpo.r ) <- CRS("+init=epsg:28356")
#writeRaster(tpo.r,"raster_crop\\tpo.TIF", overwrite=TRUE)
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

windows(); plot(awc.r)
awc.crop <- crop(x = awc.r, y = r, snap = "near")
windows(); plot(awc.crop)
# Resample:
awc.r <- resample(x = awc.crop, y = r)
#writeRaster(awc.r,"raster_crop\\sbd.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(awc.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)

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
#writeRaster(fpc.r,"raster_crop\\fpc.TIF", overwrite=TRUE)
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

#writeRaster(sbd.r,"raster_crop\\sbd.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(sbd.r)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
points(mydata.mga)

# ---------------------------------------------------------------------------------------------------------------------------------
#3. Soil clay:%

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

  #writeRaster(clay.r,"raster_crop\\clay.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(clay.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)

#fig.clay()
#pdf("clay.pdf", width=6, height=8)
#fig.clay()
#dev.off()
# ---------------------------------------------------------------------------------------------------------------------------------
#4. Soil water capacity:

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
  #writeRaster(water.r,"raster_crop\\water.TIF",overwrite=TRUE)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  # Plot to check:
  windows(); plot(water.r)
  plot(mydata.p, add = TRUE)

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
  #writeRaster(nitro.r,"raster_crop\\nitro.TIF",overwrite=TRUE)
  
  # Plot to check:
  windows(); plot(nitro.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)

# ---------------------------------------------------------------------------------------------------------------------------------
#6. Road network (expressed as the number of metres of road per square km grid):

  unzip("vector\\AU_Qld_road_length.zip")
  roaden.shp <- readShapePoly("AU_Qld_road_length-MGA65.shp")
  proj4string(roaden.shp) <- CRS("+init=epsg:28356") 
    # Vector to raster conversion:
    roaden.r <- rasterize(x = roaden.shp, y = mydata.r, field = "Sum_lengh")
    windows(); plot(roaden.r)
  plot(mydata.p, add = TRUE)
  #r <- raster(ncol = 100, nrow = 100)
  extent(r) <- extent(mydata.r)
    # Crop the raster to the shape file spatial extent:
  roaden.crop <- crop(x = roaden.r, y = r, snap = "near")
  windows(); plot(roaden.crop)
    # Resample:
  roaden.r <- resample(x = roaden.crop, y = r)
  
 # writeRaster(roaden.r,"raster_crop\\roaden.TIF",overwrite=TRUE)
  # Plot to check:
  windows(); plot(roaden.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  #plot(mydata.p, add = TRUE)
  writeRaster(roaden.r,"roaden.TIF", overwrite=TRUE)

# ---------------------------------------------------------------------------------------------------------------------------------
#7. Human population density

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
  windows(); plot(hpop.crop)
  
  # Resample:
  hpop.r <- resample(x = hpop.crop, y = r)
  projection(hpop.r) <- CRS("+init=epsg:28356")
  writeRaster(hpop.r,"raster_crop\\hpop.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(hpop.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)

# ---------------------------------------------------------------------------------------------------------------------------------
#8. Rainfall

  rain.r <- raster("raster\\AU_Qld_rainfall-MGA65.tif")
  projection(rain.r) <- CRS("+init=epsg:4326")
  # Reproject to MGA-56:
  rain.r <- projectRaster(from = rain.r, crs = CRS("+init=epsg:28356"))
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
  #writeRaster(rain.r,"raster_crop\\rain.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(rain.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
#9. Temperature:http://www.bom.gov.au/jsp/ncc/climate_averages/temperature/index.jsp

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
  
  #contour(temp.r)
  #contourplot(temp.r)
  # Plot to check:
  windows(); plot(temp.r)
  points(acsel, cex = 0.4, col = 'red', pch = 15)
  plot(mydata.p, add = TRUE)
  writeRaster(temp.r,"raster_crop\\temp.TIF",overwrite=TRUE)

# ---------------------------------------------------------------------------------------------------------------------------------
#10. Elevation

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
-----------------------------------------------------------------------------------------------------------------
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
      #writeRaster(aspect.r,"raster_crop\\aspect.TIF", overwrite=TRUE)
      # Plot to check:
      windows(); plot(aspect.r)
      points(acsel, cex = 0.4, col = 'red', pch = 15)
      plot(mydata.p, add = TRUE)
      
      
      -------------------------------------------------------------
        #14.Topographic_Wetness_Index  is a steady state wetness index. It is commonly used to quantify 
        #topographic control on hydrological processes  
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
      #writeRaster( twi.r,"raster_crop\\ twi.TIF", overwrite=TRUE)
      # Plot to check:
      windows(); plot(twi.r)
      points(acsel, cex = 0.4, col = 'red', pch = 15)
      plot(mydata.p, add = TRUE)
  -----------------------------------------------------------------------------------------------------------------
        #11.Topographic Position Index (TPI)
        fig.topo <- function() {
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
          #writeRaster(topo.r,"raster_crop\\topo.TIF", overwrite=TRUE)
          # Plot to check:
          windows(); plot(topo.r)
          points(acsel, cex = 0.4, col = 'red', pch = 15)
          plot(mydata.p, add = TRUE)
          
        }
      
      fig.topo()
      graphics.off()
      ------------------
        #12.slope:percentage change in the elevation calcuated by difference of elevation of two points diveded by lateral distance.
        
        
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
      #writeRaster(slope.r,"raster_crop\\slope.TIF", overwrite=TRUE)
      # Plot to check:
      windows(); plot(slope.r)
      points(acsel, cex = 0.4, col = 'red', pch = 15)
      plot(mydata.p, add = TRUE)
      -------------------------------------------------------------------------------------------------------------------
temp_maximum.r <- raster("raster\\AU_Qld_temp_maximum1990-2000-MGA56.TIF")  
temp_maximum.crop <- crop(x = temp_maximum.r, y = r, snap = "near")
windows(); plot(temp_maximum.crop)
writeRaster(temp_maximum.crop,"raster_crop\\temp_maximum.TIF", overwrite=TRUE)
-----------------------------------------------------------------------------------------------------------------
  temp_minimum.r <- raster("raster\\AU_Qld_temp_minimum1996_2005-MGA56.TIF")  
temp_minimum.crop <- crop(x = temp_minimum.r, y = r, snap = "near")
windows(); plot(temp_minimum.crop)
writeRaster(temp_minimum.crop,"raster_crop\\temp_minimum.TIF", overwrite=TRUE)

----------------------------------------------------------------------------------------------------------------
  rain96_2005.r <- raster("raster\\AU_Qld_rain_96_20051-MGA56.TIF")  
rain96_2005.crop <- crop(x = rain96_2005.r, y = r, snap = "near")
windows(); plot(rain96_2005.crop)
writeRaster(rain96_2005.r.crop,"raster_crop\\rain96_2005.r.TIF", overwrite=TRUE)
-----------------------------------------------------------------------------------------------------------------------
pto100_200.r <- raster("raster\\PTO_100_200_EV_N_P_AU_NAT_C_161.TIF")  
pto100_200.crop <- crop(x = pto100_200.r, y = r, snap = "near")
windows(); plot(pto100_200.crop)
writeRaster(pto100_200.crop,"raster_crop\\pto100_200.TIF", overwrite=TRUE)

------------------------------------------------------------------------------------------------------------------
pto060_100.r <- raster("raster\\PTO_060_100_EV_N_P_AU_NAT_C_131.TIF")  
pto060_100.crop <- crop(x = pto060_100.r, y = r, snap = "near")
windows(); plot(pto060_100.crop)
writeRaster(pto060_100.crop,"raster_crop\\pto060_100.TIF", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------
  pto030_060.r <- raster("raster\\PTO_030_060_EV_N_P_AU_NAT_C_101.TIF")  
pto030_060.crop <- crop(x = pto030_060.r, y = r, snap = "near")
windows(); plot(pto030_060.crop)
writeRaster(pto030_060.crop,"raster_crop\\pto030_060.TIF", overwrite=TRUE)
-------------------------------------------------------------------------------------------------------------------------------  
pto015_030.r <- raster("raster\\PTO_015_030_EV_N_P_AU_NAT_C_71.TIF")  
pto015_030.crop <- crop(x = pto015_030.r, y = r, snap = "near")
windows(); plot(pto015_030.crop)
writeRaster(pto015_030.crop,"raster_crop\\pto015_030.TIF", overwrite=TRUE)  
--------------------------------------------------------------------------------------------------------------------------------
pto005_015.r <- raster("raster\\PTO_005_015_EV_N_P_AU_NAT_C_41.TIF")  
pto005_015.crop <- crop(x = pto005_015.r, y = r, snap = "near")
windows(); plot(pto005_015.crop)
writeRaster(pto005_015.crop,"raster_crop\\pto005_015.TIF", overwrite=TRUE)  
---------------------------------------------------------------------------------------------------------------------------------  
#nitrogen
nto100_200.r <- raster("raster\\NTO_100_200_EV_N_P_AU_NAT_C_161.TIF")  
nto100_200.crop <- crop(x = nto100_200.r, y = r, snap = "near")
windows(); plot(nto100_200.crop)
writeRaster(nto100_200.crop,"raster_crop\\nto100_200.TIF", overwrite=TRUE)   
----------------------------------------------------------------------------------------------------------------------------------------------------
  nto060_100.r <- raster("raster\\NTO_060_100_EV_N_P_AU_NAT_C_131.TIF")  
nto060_100.crop <- crop(x = nto060_100.r, y = r, snap = "near")
windows(); plot(nto060_100.crop)
writeRaster(nto060_100.crop,"raster_crop\\nto060_100.TIF", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------------------------------------------
nto030_060.r <- raster("raster\\NTO_030_060_EV_N_P_AU_NAT_C_101.TIF")  
nto030_060.crop <- crop(x = nto030_060.r, y = r, snap = "near")
windows(); plot(nto030_060.crop)
writeRaster(nto030_060.crop,"raster_crop\\nto030_060.TIF", overwrite=TRUE)
-------------------------------------------------------------------------------------------------------------------------------------------------
  nto015_030.r <- raster("raster\\NTO_015_030_EV_N_P_AU_NAT_C_71.TIF")  
projection(nto015_030.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
nto015_030.r<- projectRaster(from = nto015_030.r, crs = CRS("+init=epsg:28356"))
writeRaster(nto015_030.r,"raster\\NTO_015_030_EV_N_P_AU_NAT_C_71.TIF", overwrite=TRUE)
nto015_030.crop <- crop(x = nto015_030.r, y = r, snap = "near")
windows(); plot(nto015_030.crop)
writeRaster(nto015_030.crop,"raster_crop\\nto030_060.TIF", overwrite=TRUE)  
-------------------------------------------------------------------------------------------------------------------------------------------
  nto005_015.r <- raster("raster\\NTO_005_015_EV_N_P_AU_NAT_C_41.TIF")  # file is missing
projection(nto005_015.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
nto005_015.r<- projectRaster(from = nto005_015.r, crs = CRS("+init=epsg:28356"))
writeRaster(nto005_015.r,"raster\\NTO_005_015_EV_N_P_AU_NAT_C_41.TIF", overwrite=TRUE)
nto005_015.crop <- crop(x = nto005_015.r, y = r, snap = "near")
windows(); plot(nto005_015.crop)
writeRaster(nto005_015.crop,"raster_crop\\nto005_015.TIF", overwrite=TRUE)   
----------------------------------------------------------------------------------------------------------
clay100_200.r <- raster("raster\\CLY_100_200_EV_N_P_AU_NAT_C_161.tif")  
clay100_200.crop <- crop(x = clay100_200.r, y = r, snap = "near")
windows(); plot(clay100_200.crop)
writeRaster(clay100_200.crop,"raster_crop\\clay100_200.TIF", overwrite=TRUE)   
-----------------------------------------------------------------------------
clay060_100.r <- raster("raster\\CLY_060_100_EV_N_P_AU_NAT_C_131.tif")  
clay060_100.crop <- crop(x = clay060_100.r, y = r, snap = "near")
windows(); plot(clay060_100.crop)
writeRaster(clay060_100.crop,"raster_crop\\clay060_100.TIF", overwrite=TRUE)
-------------------------------------------------------------------------------------------------------
  clay030_060.r <- raster("raster\\CLY_030_060_EV_N_P_AU_NAT_C_101.tif")  
projection(clay030_060.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
clay030_060.r<- projectRaster(from = clay030_060.r, crs = CRS("+init=epsg:28356"))
writeRaster(clay030_060.r,"raster\\CLY_030_060_EV_N_P_AU_NAT_C_101.TIF", overwrite=TRUE)
clay030_060.crop <- crop(x = clay030_060.r, y = r, snap = "near")
windows(); plot(clay030_060.crop)
writeRaster(clay030_060.crop,"raster_crop\\clay030_060.TIF", overwrite=TRUE)
---------------------------------------------------------------------------------------------
  clay015_030.r <- raster("raster\\CLY_015_030_EV_N_P_AU_NAT_C_71.tif")  
projection(clay015_030.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
clay015_030.r<- projectRaster(from = clay015_030.r, crs = CRS("+init=epsg:28356"))
writeRaster(clay015_030.r,"raster\\CLY_030_060_EV_N_P_AU_NAT_C_101.TIF", overwrite=TRUE)
clay015_030.crop <- crop(x = clay015_030.r, y = r, snap = "near")
windows(); plot(clay015_030.crop)
writeRaster(clay015_030.crop,"raster_crop\\clay015_030.TIF", overwrite=TRUE)
--------------------------------------------------------------------------------------------- 
clay005_015.r <- raster("raster\\CLY_005_015_EV_N_P_AU_NAT_C_41.tif")  
clay005_015.crop <- crop(x = clay005_015.r, y = r, snap = "near")
windows(); plot(clay005_015.crop)
writeRaster(clay005_015.crop,"raster_crop\\clay005_015.TIF", overwrite=TRUE)
---------------------------------------------------------------------------------------------
BDW_100_200.r <- raster("raster\\BDW_100_200_EV_N_P_AU_NAT_C_161.tif")  
BDW_100_200.crop <- crop(x = BDW_100_200.r, y = r, snap = "near")
windows(); plot(BDW_100_200.crop)
writeRaster(BDW_100_200.crop,"raster_crop\\BDW_100_200.tif", overwrite=TRUE)
-------------------------------------------------------------------------------------------------------------------------
BDW_060_100.r <- raster("raster\\BDW_060_100_EV_N_P_AU_NAT_C_131.tif")  
BDW_060_100.crop <- crop(x = BDW_060_100.r, y = r, snap = "near")
windows(); plot(BDW_060_100.crop)
writeRaster(BDW_060_100.crop,"raster_crop\\BDW_060_100.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------------------------
BDW_030_060.r <- raster("raster\\BDW_030_060_EV_N_P_AU_NAT_C_101.tif")  
BDW_030_060.crop <- crop(x = BDW_030_060.r, y = r, snap = "near")
windows(); plot(BDW_030_060.crop)
writeRaster(BDW_030_060.crop,"raster_crop\\BDW_030_060.tif", overwrite=TRUE) 
--------------------------------------------------------------------------------------------------------------------------
BDW_005_015.r <- raster("raster\\BDW_005_015_EV_N_P_AU_NAT_C_41.tif")  
BDW_005_015.crop <- crop(x = BDW_005_015.r, y = r, snap = "near")
windows(); plot(BDW_005_015.crop)
writeRaster(BDW_005_015.crop,"raster_crop\\BDW_005_015.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------
awc_005_015.r <- raster("raster\\AU_Qld_soil_available_water_capacity_100_200cm-wgs84.tif")  
awc_005_015.crop <- crop(x = awc_005_015.r, y = r, snap = "near")
windows(); plot(awc_005_015.crop)
writeRaster(awc_005_015.crop,"raster_crop\\awc_005_015.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------
  awc_060_100.r <- raster("raster\\AU_Qld_soil_available_water_capacity_060_100cm-wgs84.tif")  

awc_060_100.crop <- crop(x = awc_060_100.r, y = r, snap = "near")
windows(); plot(awc_060_100.crop)
writeRaster(awc_060_100.crop,"raster_crop\\awc_060_100.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------
  awc_030_060.r <- raster("raster\\AU_Qld_soil_available_water_capacity_030_060cm-wgs84.tif")  
awc_030_060.crop <- crop(x = awc_030_060.r, y = r, snap = "near")
windows(); plot(awc_030_060.crop)
writeRaster(awc_030_060.crop,"raster_crop\\awc_030_060.tif", overwrite=TRUE)
------------------------------------------------------------------------------------------
awc_100_200.r <- raster("raster\\AU_Qld_soil_available_water_capacity_100_200cm-wgs84.tif")  
awc_100_200.crop <- crop(x = awc_100_200.r, y = r, snap = "near")
windows(); plot(awc_005_015.crop)
writeRaster(awc_100_200.crop,"raster_crop\\awc_100_200.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------------------
  awc_005_015.r <- raster("raster\\AU_Qld_soil_available_water_capacity_005_015cm-wgs84.tif")  
awc_005_015.crop <- crop(x = awc_005_015.r, y = r, snap = "near")
windows(); plot(awc_005_015.crop)
writeRaster(awc_005_015.crop,"raster_crop\\awc_005_015.tif", overwrite=TRUE)
---------------------------------------------------------------------------------------------------------------
  
  library(rasterVis)  #https://pakillo.github.io/R-GIS-tutorial/#plotraster
#https://geoscripting-wur.github.io/AdvancedRasterAnalysis/
#load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_CAR_model_v01neditedRavi.RData")
#histogram(aspect.r)
#density(elev.r)
#density(stack.n)
#persp(slope.r)
#contour(stack.r)
#contourplot(roaden.r)
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

# ---------------------------------------------------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------------------------------
  #lot density number of land plots in 1k grid.
  lot_density.r <- raster("raster\\lot_density.tif")

lot_density.crop <- crop(x = lot_density.r, y = r, snap = "near")
windows(); plot(lot_density.crop)
# Plot to check:
windows(); plot(lot_density.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
---------------------------------------------------------------------------------------------------------------
#Road legth in 1k2rhohat to be done
  unclassified.r <- raster("roads_raster\\unclassified.tif")

unclassified.crop <- crop(x = unclassified.r, y = r, snap = "near")
windows(); plot(unclassified.crop)
# Resample:
unclassified.r <- resample(x = unclassified.crop, y = r)
# Plot to check:
windows(); plot(unclassified.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
------------------------------------------------------------------------------------------------------------------
  track.r <- raster("roads_raster\\track.tif")

track.crop <- crop(x = track.r, y = r, snap = "near")
windows(); plot(track.crop)

# Plot to check:
windows(); plot(track.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
----------------------------------------------------------------------------------------------------------------
  trunck.r <- raster("roads_raster\\trunck.tif")
trunck.crop <- crop(x = trunck.r, y = r, snap = "near")
windows(); plot(trunck.crop)
# Plot to check:
windows(); plot(trunck.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
-----------------------------------------------------------------------------------------------------------------
  motorwaylink.r <- raster("roads_raster\\motorwaylink.tif")
motorwaylink.crop <- crop(x = motorwaylink.r, y = r, snap = "near")
windows(); plot(motorwaylink.crop)
# Plot to check:
windows(); plot(motorwaylink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
-----------------------------------------------------------------------------------------------------------------
  motorway.r <- raster("roads_raster\\motorway.tif")
motorway.crop <- crop(x = motorway.r, y = r, snap = "near")
windows(); plot(motorway.crop)

# Plot to check:
windows(); plot(motorway.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
----------------------------------------------------------------------------------------------------------------
  tertiary_liink.r <- raster("roads_raster\\tertiary_liink.tif")
tertiary_liink.crop <- crop(x = tertiary_liink.r, y = r, snap = "near")
windows(); plot(tertiary_liink.crop)
# Plot to check:
windows(); plot(tertiary_liink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
----------------------------------------------------------------------------------------------------------------- 
  tertiary.r <- raster("roads_raster\\tertiary.tif")
tertiary.crop <- crop(x = tertiary.r, y = r, snap = "near")
windows(); plot(tertiary.crop)
# Plot to check:
windows(); plot(tertiary.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
  
--------------------------------------------------------------------------------------------------------------------
  secondry_link.r <- raster("roads_raster\\secondry_link.tif")
secondry_link.crop <- crop(x = secondry_link.r, y = r, snap = "near")
windows(); plot(secondry_link.crop)
writeRaster( secondry_link.crop,"raster_crop\\ RDsecondry_link.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(secondry_link.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)
--------------------------------------------------------------------------------------------------------------------
  secondry.r <- raster("roads_raster\\secondry.tif")
secondry.crop <- crop(x = secondry.r, y = r, snap = "near")
windows(); plot(secondry.crop)
# Plot to check:
windows(); plot(secondry.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
--------------------------------------------------------------------------------------------------------------------
  roads_other.r <- raster("roads_raster\\roads_other.tif")
roads_other.crop <- crop(x = roads_other.r, y = r, snap = "near")
windows(); plot(roads_other.crop)
# Plot to check:
windows(); plot(roads_other.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
--------------------------------------------------------------------------------------------------------------------
  residential.r <- raster("roads_raster\\residential.tif")
residential.crop <- crop(x = residential.r, y = r, snap = "near")
windows(); plot(residential.crop)
# Plot to check:
windows(); plot(residential.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
  
------------------------------------------------------------------------------------------------------------------- 
  primary.r <- raster("roads_raster\\primary.tif")
primary.crop <- crop(x = primary.r, y = r, snap = "near")
windows(); plot(primary.crop)
# Plot to check:
windows(); plot(primary.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
-------------------------------------------------------------------------------------------------------------------
  pedestrian.r <- raster("roads_raster\\pedestrian.tif")
pedestrian.crop <- crop(x = pedestrian.r, y = r, snap = "near")
windows(); plot(pedestrian.crop)
# Plot to check:
windows(); plot(pedestrian.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
-------------------------------------------------------------------------------------------------------------------
  path.r <- raster("roads_raster\\path.tif")
path.crop <- crop(x = path.r, y = r, snap = "near")
windows(); plot(path.crop)
# Plot to check:
windows(); plot(path.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
-----------------------------------------------------------------------------
  cycle.r <- raster("roads_raster\\cycle.tif")
cycle.crop <- crop(x = cycle.r, y = r, snap = "near")
windows(); plot(cycle.crop)
# Plot to check:
windows(); plot(cycle.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
-------------------------------------------------------------------------------------------------------------------  
  bridle.r <- raster("roads_raster\\bridle.tif")
bridle.crop <- crop(x = bridle.r, y = r, snap = "near")
windows(); plot(bridle.crop)
# Plot to check:
windows(); plot(bridle.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE) 
-------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------
# now add distance to roads
  distance_pedestrian.r <- raster("raster\\distance_pedestrian.tif")
distance_pedestrian.crop <- crop(x = distance_pedestrian.r, y = r, snap = "near")
windows(); plot(distance_pedestrian.crop)
# Plot to check:
windows(); plot(distance_pedestrian.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
-------------------------------------------------------------------------------------------------------------------  
  distance_service.r <- raster("raster\\distance_service.tif")
distance_service.crop <- crop(x = distance_service.r, y = r, snap = "near")
windows(); plot(distance_service.crop)
# Plot to check:
windows(); plot(distance_service.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
-------------------------------------------------------------------------------------------------------------------  
  distance_bridleway.r <- raster("raster\\distance_bridleway.tif")
distance_bridleway.crop <- crop(x = distance_bridleway.r, y = r, snap = "near")
windows(); plot(distance_bridleway.crop)
# Plot to check:
windows(); plot(distance_bridleway.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
-------------------------------------------------------------------------------------------------------------------
  distance_cycleway.r <- raster("raster\\distance_cycleway.tif")
distance_cycleway.crop <- crop(x = distance_cycleway.r, y = r, snap = "near")
windows(); plot(distance_cycleway.crop)
# Plot to check:
windows(); plot(distance_cycleway.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_path.r <- raster("raster\\distance_path.tif")
distance_path.crop <- crop(x = distance_path.r, y = r, snap = "near")
windows(); plot(distance_path.crop)
# Plot to check:
windows(); plot(distance_path.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------ 
  distance_steps.r <- raster("raster\\distance_steps.tif")
distance_steps.crop <- crop(x = distance_steps.r, y = r, snap = "near")
windows(); plot(distance_steps.crop)
# Plot to check:
windows(); plot(distance_steps.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_residentil.r <- raster("raster\\distance_residentil.tif")
distance_residentil.crop <- crop(x = distance_residentil.r, y = r, snap = "near")
windows(); plot(distance_residentil.crop)
# Plot to check:
windows(); plot(distance_residentil.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_unclassified.r <- raster("raster\\distance_unclassified.tif")
distance_unclassified.crop <- crop(x = distance_unclassified.r, y = r, snap = "near")
windows(); plot(distance_unclassified.crop)
# Plot to check:
windows(); plot(distance_unclassified.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_footway.r <- raster("raster\\distance_footway.tif")
distance_footway.crop <- crop(x = distance_footway.r, y = r, snap = "near")
windows(); plot(distance_footway.crop)
writeRaster( distance_footway.crop,"raster_crop\\ distance_footway.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(distance_footway.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_tertiaryandlink.r <- raster("raster\\distance_tertiaryandlink.tif")
distance_tertiaryandlink.crop <- crop(x = distance_tertiaryandlink.r, y = r, snap = "near")
windows(); plot(distance_tertiaryandlink.crop)
# Plot to check:
windows(); plot(distance_tertiaryandlink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_secondaryandlink.r <- raster("raster\\distance_secondaryandlink.tif")
distance_secondaryandlink.crop <- crop(x = distance_secondaryandlink.r, y = r, snap = "near")
windows(); plot(distance_secondaryandlink.crop)
# Plot to check:
windows(); plot(distance_secondaryandlink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_primaryandlink.r <- raster("raster\\distance_primaryandlink.tif")
distance_primaryandlink.crop <- crop(x = distance_primaryandlink.r, y = r, snap = "near")
windows(); plot(distance_primaryandlink.crop)
# Plot to check:
windows(); plot(distance_primaryandlink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_motorwayandlink.r <- raster("raster\\distance_motorwayandlink.tif")
distance_motorwayandlink.crop <- crop(x = distance_motorwayandlink.r, y = r, snap = "near")
windows(); plot(distance_motorwayandlink.crop)
# Plot to check:
windows(); plot(distance_motorwayandlink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_trunkandlink.r <- raster("raster\\distance_trunkandlink.tif")
distance_trunkandlink.crop <- crop(x = distance_trunkandlink.r, y = r, snap = "near")
windows(); plot(distance_trunkandlink.crop)
# Plot to check:
windows(); plot(distance_trunkandlink.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  distance_track.r <- raster("raster\\distance_track.tif")
distance_track.crop <- crop(x = distance_track.r, y = r, snap = "near")
windows(); plot(distance_track.crop)
# Plot to check:
windows(); plot(distance_track.crop)
points(acsel, cex = 0.4, col = 'red', pch = 15)
plot(mydata.p, add = TRUE)  
------------------------------------------------------------------------------------------------------------------
  

------------------------------------------------------------------------------------------------------------------------------------
## rhohat - :Computes a smoothing estimate of the intensity of a point process, as a function of a (continuous) spatial covariate.
#Foliage projective cover
# Convert tpo.r into a spatstat image object:
--------------------------------------------------------------------------------------------------------------------------------
 lot_density.im <- as.im(lot_density.r)
windows(); plot(lot_density.im, axes = TRUE)

# What is the nature of the association between koala sight locations and lot densities?
lot_density.rho <- rhohat(object = dat.ppp, covariate = lot_density.im)
windows(); plot(lot_density.rho, xlab = "lot_density (units)", main = "")
savePlot(filename = "lot_density",type = c( "png"),device = dev.cur())
----------------------------------------------------------------------------------------------------------------------
  
  suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif")
suitable_3.im <- as.im(suitable_3.r)
windows(); plot(suitable_3.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_3.rho <- rhohat(object = dat.ppp, covariate = suitable_3.im)
windows(); plot(suitable_3.rho,xlab = "suitable_3(units)", main = "")
savePlot(filename = "suitable_3(units)",type = c( "png"),device = dev.cur())
------------------------------------------------------------------------------------------------
  suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif")
suitable_2.im <- as.im(suitable_2.r)
windows(); plot(suitable_2.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_2.rho <- rhohat(object = dat.ppp, covariate = suitable_2.im)
windows(); plot(suitable_2.rho,xlab = "suitable_2(units)", main = "")
savePlot(filename = "suitable_2(units)",type = c( "png"),device = dev.cur())
------------------------------------------------------------------------------------------------
  suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif")
suitable_1.im <- as.im(suitable_1.r)
windows(); plot(suitable_1.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_1.rho <- rhohat(object = dat.ppp, covariate = suitable_1.im)
windows(); plot(suitable_1.rho,xlab = "suitable_1(units)", main = "")
savePlot(filename = "suitable_1(units)",type = c( "png"),device = dev.cur())
------------------------------------------------------------------------------------------------
 suitable_0.r<- raster("habit_vector_raster\\suitable_0.tif")
suitable_0.im <- as.im(suitable_0.r)
windows(); plot(suitable_0.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_0.rho <- rhohat(object = dat.ppp, covariate = suitable_0.im)
windows(); plot(suitable_0.rho,xlab = "suitable_0(units)", main = "")
savePlot(filename = "suitable_0(units)",type = c( "png"),device = dev.cur())
----------------------------------------------------------------------------------------
#next step is to add road length density and road distance rasters and rohoat. files already preapred. in folders
#roads_raster= road length, 
-------------------------------------------------------------------------------------------------------------------------
 #rhohat for road densities. and distace to roads follows.
  unclassified.im <- as.im(unclassified.r )
windows(); plot(unclassified.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
unclassified.rho <- rhohat(object = dat.ppp, covariate = unclassified.im)
windows(); plot(unclassified.rho, xlab = "unclassified road density (units)", main = "")
savePlot(filename = "unclassifiedRD",type = c( "png"),device = dev.cur())
---------------------------------------------------------------------------------------------------------------------------
  track.im <- as.im(track.r  )
windows(); plot(track.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
track.rho <- rhohat(object = dat.ppp, covariate = track.im)
windows(); plot(track.rho, xlab = "track density (units)", main = "")
savePlot(filename = "trackRD",type = c( "png"),device = dev.cur())
---------------------------------------------------------------------------------------------------------------------------  
  trunck.im <- as.im(trunck.r  )
windows(); plot(trunck.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
trunck.rho <- rhohat(object = dat.ppp, covariate = trunck.im)
windows(); plot(trunck.rho, xlab = "trunck density (units)", main = "")
savePlot(filename = "trunckRD",type = c( "png"),device = dev.cur())
--------------------------------------------------------------------------------------------------------------------- 
  motorwaylink.im <- as.im(motorwaylink.r )
windows(); plot(motorwaylink.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
motorwaylink.rho <- rhohat(object = dat.ppp, covariate = motorwaylink.im)
windows(); plot(motorwaylink.rho, xlab = "motorwaylink density (units)", main = "")
savePlot(filename = "motorwaylinkRD",type = c( "png"),device = dev.cur())
---------------------------------------------------------------------------------------------------------------------  
  motorway.im <- as.im(motorway.r )
windows(); plot(motorway.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
motorway.rho <- rhohat(object = dat.ppp, covariate = motorway.im)
windows(); plot(motorway.rho, xlab = "motorway density (units)", main = "")
savePlot(filename = "motorwayRD",type = c( "png"),device = dev.cur())
--------------------------------------------------------------------------------------------------------------------- 
tertiary_liink.im <- as.im(tertiary_liink.r )
windows(); plot(tertiary_liink.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
tertiary_liink.rho <- rhohat(object = dat.ppp, covariate = tertiary_liink.im)
windows(); plot(tertiary_liink.rho, xlab = "tertiary_liink density (units)", main = "")
savePlot(filename = "tertiary_liinkRD",type = c( "png"),device = dev.cur())
--------------------------------------------------------------------------------------------------------------------- 
  tertiary.im <- as.im(tertiary.r )
windows(); plot(tertiary.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
tertiary.rho <- rhohat(object = dat.ppp, covariate = tertiary.im)
windows(); plot(tertiary.rho, xlab = "tertiary density (units)", main = "")
savePlot(filename = "tertiaryRD",type = c( "png"),device = dev.cur())
---------------------------------------------------------------------------------------------------------------------
  secondry_link.im <- as.im(secondry_link.r )
windows(); plot(secondry_link.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
secondry_link.rho <- rhohat(object = dat.ppp, covariate = secondry_link.im)
windows(); plot(secondry_link.rho, xlab = "secondry_link density (units)", main = "")
savePlot(filename = "secondry_linkRD",type = c( "png"),device = dev.cur())
---------------------------------------------------------------------------------------------------------------------
  secondry.im <- as.im(secondry.r )
windows(); plot(secondry.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
secondry.rho <- rhohat(object = dat.ppp, covariate = secondry.im)
windows(); plot(secondry.rho, xlab = "secondry density (units)", main = "")
savePlot(filename = "secondryRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  roads_other.im <- as.im(roads_other.r )
windows(); plot(roads_other.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
roads_other.rho <- rhohat(object = dat.ppp, covariate = roads_other.im)
windows(); plot(roads_other.rho, xlab = "roads_other density (units)", main = "")
savePlot(filename = "roads_otherRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  residential.im <- as.im(residential.r )
windows(); plot(residential.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
residential.rho <- rhohat(object = dat.ppp, covariate = residential.im)
windows(); plot(residential.rho, xlab = "residential density (units)", main = "")
savePlot(filename = "residentialRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  primary.im <- as.im(primary.r )
windows(); plot(primary.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
primary.rho <- rhohat(object = dat.ppp, covariate = primary.im)
windows(); plot(primary.rho, xlab = "primary density (units)", main = "")
savePlot(filename = "primaryRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  pedestrian.im <- as.im(pedestrian.r )
windows(); plot(pedestrian.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
pedestrian.rho <- rhohat(object = dat.ppp, covariate = pedestrian.im)
windows(); plot(pedestrian.rho, xlab = "pedestrian density (units)", main = "")
savePlot(filename = "pedestrianRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  path.im <- as.im(path.r )
windows(); plot(path.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
path.rho <- rhohat(object = dat.ppp, covariate = path.im)
windows(); plot(path.rho, xlab = "path density (units)", main = "")
savePlot(filename = "pathRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  cycle.im <- as.im(cycle.r )
windows(); plot(cycle.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
cycle.rho <- rhohat(object = dat.ppp, covariate = cycle.im)
windows(); plot(cycle.rho, xlab = "cycle density (units)", main = "")
savePlot(filename = "cycleRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  bridle.im <- as.im(bridle.r )
windows(); plot(bridle.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
bridle.rho <- rhohat(object = dat.ppp, covariate = bridle.im)
windows(); plot(bridle.rho, xlab = "bridle density (units)", main = "")
savePlot(filename = "bridleRD",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
#new set of varissables
  
 
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
  fig.aspect.rho <- function() {
  # Convert tpo.r into a spatstat image object:
  aspect.im <- as.im(aspect.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
aspect.rho <- rhohat(object = dat.ppp, covariate = aspect.im)

windows(); plot(aspect.rho, xlab = "aspect (units)", main = "")
savePlot(filename = "aspect",type = c( "png"),device = dev.cur())
  }
fig.aspect.rho () 
graphics.off()
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

windows(); plot(hpop.rho, xlim = c(0, 3000), ylim = c(0, 0.5e-06), xlab = "Human population density", main = "")
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
# rhohat analyses - road density:

roaden.im <- as.im(roaden.r)

windows(); plot(roaden.im, axes = TRUE)

roaden.rho <- rhohat(object = dat.ppp, covariate = roaden.im)

windows(); plot(roaden.rho, xlab = "Road density (metres per square km)", main = "")
savePlot(filename = "roaden",type = c( "png"),device = dev.cur())


# ---------------------------------------------------------------------------------------------------------------------------------

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
--------------------------------------------------------------------------------------------
  awc.im <- as.im(awc.r)
windows(); plot(awc.im, axes = TRUE)

awc.rho <- rhohat(object = dat.ppp, covariate = awc.im)
windows(); plot(awc.rho, xlab = "available water capacity", main = "")
savePlot(filename = "awc",type = c( "png"),device = dev.cur())

-------------------------------------------------------------------------------------------------------------------------------------------
  topo.im <- as.im(topo.r)
windows(); plot(topo.im, axes = TRUE)

topo.rho <- rhohat(object = dat.ppp, covariate = topo.im)
windows(); plot(topo.rho, xlab = "topo ()", main = "")
savePlot(filename = "topo",type = c( "png"),device = dev.cur())
-----------------------------------------------------------------------------------------------------------------------------------
  slope.im <- as.im(slope.r)
windows(); plot(slope.im, axes = TRUE)

slope.rho <- rhohat(object = dat.ppp, covariate = slope.im)
windows(); plot(slope.rho, xlab = "slope()", main = "")
savePlot(filename = "slope",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------- 
  temp_maximum.im <- as.im(temp_maximum.r)
windows(); plot(temp_maximum.im, axes = TRUE)

temp_maximum.rho <- rhohat(object = dat.ppp, covariate = temp_maximum.im)
windows(); plot(temp_maximum.rho, xlab = "temp_maximum()", main = "")
savePlot(filename = "temp_maximum",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------
  temp_minimum.im <- as.im(temp_minimum.r)
windows(); plot(temp_minimum.im, axes = TRUE)

temp_minimum.rho <- rhohat(object = dat.ppp, covariate = temp_minimum.im)
windows(); plot(temp_minimum.rho, xlab = "temp_minimum()", main = "")
savePlot(filename = "temp_minimum",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------
  rain96_2005.im <- as.im(rain96_2005.r)
windows(); plot(rain96_2005.im, axes = TRUE)

rain96_2005.rho <- rhohat(object = dat.ppp, covariate = rain96_2005.im)
windows(); plot(rain96_2005.rho, xlab = "rain96_2005()", main = "")
savePlot(filename = "rain96_2005",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------  
  pto100_200.im <- as.im(pto100_200.r)
windows(); plot(pto100_200.im, axes = TRUE)

pto100_200.rho <- rhohat(object = dat.ppp, covariate = pto100_200.im)
windows(); plot(pto100_200.rho, xlab = "pto100_200()", main = "")
savePlot(filename = "pto100_200",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------
  pto060_100.im <- as.im(pto060_100.r)
windows(); plot(pto060_100.im, axes = TRUE)

pto060_100.rho <- rhohat(object = dat.ppp, covariate = pto060_100.im)
windows(); plot(pto060_100.rho, xlab = "pto060_100()", main = "")
savePlot(filename = "pto060_100",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------
  pto030_060.im <- as.im(pto030_060.r)
windows(); plot(pto030_060.im, axes = TRUE)

pto030_060.rho <- rhohat(object = dat.ppp, covariate = pto030_060.im)
windows(); plot(pto030_060.rho, xlab = "pto030_060()", main = "")
savePlot(filename = "pto030_060",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------------------------------------------
  pto015_030.im <- as.im(pto015_030.r)
windows(); plot(pto015_030.im, axes = TRUE)

pto015_030.rho <- rhohat(object = dat.ppp, covariate = pto015_030.im)
windows(); plot(pto015_030.rho, xlab = "pto015_030()", main = "")
savePlot(filename = "pto015_030",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  pto005_015.im <- as.im(pto005_015.r)
windows(); plot(pto005_015.im, axes = TRUE)

pto005_015.rho <- rhohat(object = dat.ppp, covariate = pto005_015.im)
windows(); plot(pto005_015.rho, xlab = "pto005_015()", main = "")
savePlot(filename = "pto005_015",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  nto100_200.im <- as.im(nto100_200.r)
windows(); plot(nto100_200.im, axes = TRUE)

nto100_200.rho <- rhohat(object = dat.ppp, covariate = nto100_200.im)
windows(); plot(nto100_200.rho, xlab = "nto100_200()", main = "")
savePlot(filename = "nto100_200",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  nto060_100.im <- as.im(nto060_100.r)
windows(); plot(nto060_100.im, axes = TRUE)

nto060_100.rho <- rhohat(object = dat.ppp, covariate = nto060_100.im)
windows(); plot(nto060_100.rho, xlab = "nto060_100()", main = "")
savePlot(filename = "nto060_100",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  nto030_060.im <- as.im(nto030_060.r)
windows(); plot(nto030_060.im, axes = TRUE)

nto030_060.rho <- rhohat(object = dat.ppp, covariate = nto030_060.im)
windows(); plot(nto030_060.rho, xlab = "nto030_060()", main = "")
savePlot(filename = "nto030_060",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  nto015_030.im <- as.im(nto015_030.r)
windows(); plot(nto015_030.im, axes = TRUE)

nto015_030.rho <- rhohat(object = dat.ppp, covariate = nto015_030.im)
windows(); plot(nto015_030.rho, xlab = "nto015_030()", main = "")
savePlot(filename = "nto015_030",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  nto005_015.im <- as.im(nto005_015.r) # not avaibale
windows(); plot(nto005_015.im, axes = TRUE)

nto005_015.rho <- rhohat(object = dat.ppp, covariate = nto005_015.im)
windows(); plot(nto005_015.rho, xlab = "nto005_015()", main = "")
savePlot(filename = "nto005_015",type = c( "png"),device = dev.cur()) 
--------------------------------------------------------------------------------------------------------------------------------- 
  clay100_200.im <- as.im(clay100_200.r) 
windows(); plot(clay100_200.im, axes = TRUE)

clay100_200.rho <- rhohat(object = dat.ppp, covariate = clay100_200.im)
windows(); plot(clay100_200.rho, xlab = "clay100_200()", main = "")
savePlot(filename = "clay100_200",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  clay060_100.im <- as.im(clay060_100.r) 
windows(); plot(clay060_100.im, axes = TRUE)

clay060_100.rho <- rhohat(object = dat.ppp, covariate = clay060_100.im)
windows(); plot(clay060_100.rho, xlab = "clay060_100()", main = "")
savePlot(filename = "clay060_100",type = c( "png"),device = dev.cur()) 



---------------------------------------------------------------------------------------------------------------------------------
  clay030_060.im <- as.im(clay030_060.r) 
windows(); plot(clay030_060.im, axes = TRUE)

clay030_060.rho <- rhohat(object = dat.ppp, covariate = clay030_060.im)
windows(); plot(clay030_060.rho, xlab = "clay030_060()", main = "")
savePlot(filename = "clay030_060",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
clay015_030.im <- as.im(clay015_030.r) 
windows(); plot(clay015_030.im, axes = TRUE)

clay015_030.rho <- rhohat(object = dat.ppp, covariate = clay015_030.im)
windows(); plot(clay015_030.rho, xlab = "clay015_030()", main = "")
savePlot(filename = "clay015_030",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  clay005_015.im <- as.im(clay005_015.r) 
windows(); plot(clay005_015.im, axes = TRUE)

clay005_015.rho <- rhohat(object = dat.ppp, covariate = clay005_015.im)
windows(); plot(clay005_015.rho, xlab = "clay005_015()", main = "")
savePlot(filename = "clay005_015",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  BDW_100_200.im <- as.im(BDW_100_200.r) 
windows(); plot(BDW_100_200.im, axes = TRUE)

BDW_100_200.rho <- rhohat(object = dat.ppp, covariate = BDW_100_200.im)
windows(); plot(BDW_100_200.rho, xlab = "BDW_100_200()", main = "")
savePlot(filename = "BDW_100_200",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  BDW_060_100.im <- as.im(BDW_060_100.r) 
windows(); plot(BDW_060_100.im, axes = TRUE)

BDW_060_100.rho <- rhohat(object = dat.ppp, covariate = BDW_060_100.im)
windows(); plot(BDW_060_100.rho, xlab = "BDW_060_100()", main = "")
savePlot(filename = "BDW_060_100",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  BDW_030_060.im <- as.im(BDW_030_060.r) 
windows(); plot(BDW_030_060.im, axes = TRUE)

BDW_030_060.rho <- rhohat(object = dat.ppp, covariate = BDW_030_060.im)
windows(); plot(BDW_030_060.rho, xlab = "BDW_030_060()", main = "")
savePlot(filename = "BDW_030_060",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  BDW_005_015.im <- as.im(BDW_005_015.r) 
windows(); plot(BDW_005_015.im, axes = TRUE)

BDW_005_015.rho <- rhohat(object = dat.ppp, covariate = BDW_005_015.im)
windows(); plot(BDW_005_015.rho, xlab = "BDW_005_015()", main = "")
savePlot(filename = "BDW_005_015",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
  awc_005_015.im <- as.im(awc_005_015.r) 
windows(); plot(awc_005_015.im, axes = TRUE)

awc_005_015.rho <- rhohat(object = dat.ppp, covariate = awc_005_015.im)
windows(); plot(awc_005_015.rho, xlab = "awc_005_015()", main = "")
savePlot(filename = "awc_005_015",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------
 awc_015_030.im <- as.im(awc_015_030.r) # not tested as there is no relationship.
windows(); plot(awc_015_030.im, axes = TRUE)
awc_015_030.rho <- rhohat(object = dat.ppp, covariate =awc_015_030.im)
windows(); plot(awc_015_030.rho, xlab = "awc_015_030()", main = "")
savePlot(filename = "awc_015_030",type = c( "png"),device = dev.cur()) 
------------------------------------------------------------------------- 
awc_030_060.im <- as.im(awc_030_060.r) 
windows(); plot(awc_030_060.im, axes = TRUE)

awc_030_060.rho <- rhohat(object = dat.ppp, covariate = awc_030_060.im)
windows(); plot(awc_030_060.rho, xlab = "awc_030_060()", main = "")
savePlot(filename = "awc_030_060",type = c( "png"),device = dev.cur()) 
--------------------------------------------------------------------------------------------------------------------------------

    awc_060_100.im <- as.im(awc_060_100.r) 
windows(); plot(awc_060_100.im, axes = TRUE)

awc_060_100.rho <- rhohat(object = dat.ppp, covariate = awc_060_100.im)
windows(); plot(awc_060_100.rho, xlab = "awc_060_100()", main = "")
savePlot(filename = "awc_060_100",type = c( "png"),device = dev.cur()) 

---------------------------------------------------------------------------------------------------------------------------------
  awc_100_200.im <- as.im(awc_100_200.r) 
windows(); plot(awc_100_200.im, axes = TRUE)

awc_100_200.rho <- rhohat(object = dat.ppp, covariate = awc_100_200.im)
windows(); plot(awc_100_200.rho, xlab = "awc_100_200()", main = "")
savePlot(filename = "awc_100_200",type = c( "png"),device = dev.cur()) 
---------------------------------------------------------------------------------------------------------------------------------  
  
save.image("AU_Qld_koala_habitat_v08edited_distance_Ravi.RData")


# =================================================================================================================================
# Poisson point process model:

# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, 
                                                         tpo = tpo.im, sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                                         hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409#19912.25
windows();plot(dat.ppp)

#save.image("AU_Qld_koala_habitat_v08edited_distance16_03.RData")
# Backward stepwise variable selection:
dat.ppm01 <- ppm(dat.ppp, trend = ~  clay+nitro+ fpc+lot+ twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit1+ suit2+ suit3+ hpop,
                 covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                   sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                   hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
# Values of the covariates 'roaden', 'temp' were NA or undefined at 0.5% (5 out of 1005) of the quadrature points. 
# See http://gis.stackexchange.com/questions/161893/creating-polygonal-windows-for-spatial-analysis-in-r for a work around.

# Extract the quadrature scheme from dat.ppm01:
Qz <- quad.ppm(dat.ppm01, drop = TRUE)
# This step will drop points within the quadrature scheme that had NA-values.

# Run the model again, but this time, use the corrected quadrature scheme:
dat.ppm02 <- ppm(Qz,trend = ~  clay+nitro+ fpc+lot+ twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit1+ suit2+ suit3+ hpop,
                 covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                   sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                   hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm02)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm02)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm02) # 3413.045#19580.54
-----------------------------------------------------------------------------------------------------------------
#this is the automating of the model selction. hwoever when use the most significant varibales it fails to give same results when run independetly.
  step(dat.ppm02) # automate 
-----------------------------------------------------------------------------------------------------------------------------------------
  #drop:suit1
dat.ppm03 <- ppm(Qz,trend = ~  clay+nitro+ fpc+lot+ twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit2+ suit3+ hpop,
                 covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                   sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                   hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm03)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm03)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm03)#19578.55
----------------------------------------------------------------------------------------------------------------------
  #drop:nitro
  dat.ppm04 <- ppm(Qz,trend = ~  clay+fpc+lot+ twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit2+ suit3+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm04)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm04)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm04)#19576.72
------------------------------------------------------------------------------------------------------------------------
  #drop:suit3
  dat.ppm05 <- ppm(Qz,trend = ~  clay+fpc+lot+ twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit2+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm05)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm05)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm05)#119575.59
------------------------------------------------------------------------------------------------------------------------
  #drop:lot
  dat.ppm06 <- ppm(Qz,trend = ~  clay+fpc+twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ suit2+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm06)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm06)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm06)#19574.52
------------------------------------------------------------------------------------------------------------------------
  #drop:suit2
  dat.ppm07 <- ppm(Qz,trend = ~  clay+fpc+twi+ tpo+ sbd+ water+residential+ secondary+ thertiary+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm07)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm07)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm07)#19573.68
----------------------------------------------------------------------------------------------
  #drop:secondary
  dat.ppm08 <- ppm(Qz,trend = ~  clay+fpc+twi+ tpo+ sbd+ water+residential+ thertiary+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm08)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm08)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm08)#19572.95
----------------------------------------------------------------------------------------------
  #drop:secondary
  dat.ppm09 <- ppm(Qz,trend = ~  clay+fpc+twi+ tpo+ sbd+ water+residential+ thertiary+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm09)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm09)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm09)#19572.95
----------------------------------------------------------------------------------------------
  #drop:tertiary
  dat.ppm10 <- ppm(Qz,trend = ~  clay+fpc+twi+ tpo+ sbd+ water+residential+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm10)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm10)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm10)#19572.81
-------------------------------------------------------------------------------------------------------------------
  #drop:fpc
  dat.ppm11 <- ppm(Qz,trend = ~  clay+twi+ tpo+ sbd+ water+residential+ hpop,
                   covariates = list(clay=clay.im, nitro=nitro.im, fpc= fpc.im,lot=lot_density.im, twi=twi.im, tpo = tpo.im,
                                     sbd = sbd.im,  water = water.im, residential=residential.im, secondary=secondry.im, thertiary=tertiary.im,
                                     hpop = hpop.im, temp = temp.im, elev = elev.im, suit1=suitable_1.im, suit2=suitable_2.im, suit3=suitable_3.im))
summary(dat.ppm11)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm11)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm11)#19573.48
----------------------------------------------------------------------------------------------------------------------
# The lurking variable plot gives you an idea of mis-specification of spatial trend:
windows(); diagnose.ppm(dat.ppm13, type = "raw")
savePlot(filename = "datppm13",type = c( "png"),device = dev.cur())

# Now add a Strauss spatial interaction term:
dat.ppm14 <- ppm(Qz, trend = ~ fpc +temp +habit02, interaction = Strauss(r = 6000), 
                 covariates = list(fpc=fpc.im, tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm14)$coefs.SE.CI
AIC(dat.ppm14) # 2459.733

windows(); plot(dat.ppm14, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm14, type = "raw")


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
dat.pred14 <- predict(dat.ppm13); pred <- dat.pred14

# Express predicted koala intensity in terms of the number of koalas per hectare:
pred$v <- pred$v * 1E04
summary(as.vector(pred$v))
windows(); hist(as.vector(pred$v))



# Predictions from the model:
breaks <- seq(from = 0, to = 4.3810 , length = 6)
col <- brewer.pal(n = 5, name = "Reds")

windows(); plot(x = xlim, y = ylim, main="Predictions from the model", type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "blue")

plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Predictions from the model14",type = c( "png"),device = dev.cur())

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
save.image("AU_Qld_koala_habitat_v08edited_distance.RData")

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
