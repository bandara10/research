library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):
#load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v08edited_distance_Ravi.RData")

mydatasighting <- read.csv("mydatasighting.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)

# Analyse data for 2011:
id <- mydata$yearnew >2010
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
head(mydata.n)
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
----------------------------------------------------------------------------
mydata.r <- raster(studyarea.shp)# creating a raster for cropping  study area.
res(mydata.r) <- 1000
r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset
# =================================================================================================================================
#1.total phosphorous:% Read in explanatory variable raster maps - 
fig.PTO_000_005 <- function() {
  PTO_000_005.r <- raster("raster_syn\\PTO_000_005.tif")
  windows(); plot(PTO_000_005.r)
    # Crop the raster to the shape fireferencele spatial extent:
  PTO_000_005.crop <- crop(x = PTO_000_005.r, y = r, snap = "near")
  windows(); plot(PTO_000_005.crop)
  writeRaster(PTO_000_005.r,"raster_crop\\PTO_000_005.TIF", overwrite=TRUE)
  # Plot to check:
  windows(); plot(PTO_000_005.r, main="PTO_000_005")
  points(acsel, cex = 0.4, col = 'red', pch = 15)
}
fig.PTO_000_005()
graphics.off()
----------------------------------------------------------------------------------------------------------
  #2.total phosphorous:% Read in explanatory variable raster maps - 
  fig.PTO_005_015 <- function() {
    PTO_005_015.r <- raster("raster_syn\\PTO_005_015.tif")
    windows(); plot(PTO_005_015.r, main="PTO_005_015")
    # Crop the raster to the shape fireferencele spatial extent:
    PTO_005_015.crop <- crop(x = PTO_005_015.r, y = r, snap = "near")
    windows(); plot(PTO_005_015.crop)
    writeRaster(PTO_005_015.r,"raster_crop\\PTO_005_015.TIF", overwrite=TRUE)
    # Plot to check:
   
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.PTO_005_015()
graphics.off()
----------------------------------------------------------------------------------------------------------
  #3.total phosphorous:% Read in explanatory variable raster maps - 
  fig.PTO_015_030 <- function() {
    PTO_015_030.r <- raster("raster_syn\\PTO_015_030.tif")
    windows(); plot(PTO_015_030.r, main="PTO_015_030")
    # Crop the raster to the shape fireferencele spatial extent:
    PTO_015_030.crop <- crop(x = PTO_015_030.r, y = r, snap = "near")
    windows(); plot(PTO_015_030.crop)
    writeRaster(PTO_015_030.r,"raster_crop\\PTO_015_030.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.PTO_015_030()
graphics.off()
----------------------------------------------------------------------------------------------------------
  #4.total phosphorous:% Read in explanatory variable raster maps - 
  fig.PTO_030_060 <- function() {
    PTO_030_060.r <- raster("raster_syn\\PTO_030_060.tif")
    windows(); plot(PTO_030_060.r, main="PTO_030_060")
    # Crop the raster to the shape fireferencele spatial extent:
    PTO_030_060.crop <- crop(x = PTO_030_060.r, y = r, snap = "near")
    windows(); plot(PTO_030_060.crop)
    writeRaster(PTO_030_060.r,"raster_crop\\PTO_030_060.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.PTO_030_060()
graphics.off()
----------------------------------------------------------------------------------------------------------
  #5.total phosphorous:% Read in explanatory variable raster maps - 
  fig.PTO_060_100 <- function() {
    PTO_060_100.r <- raster("raster_syn\\PTO_060_100.tif")
    windows(); plot(PTO_060_100.r, main="PTO_060_100")
    # Crop the raster to the shape fireferencele spatial extent:
    PTO_060_100.crop <- crop(x = PTO_060_100.r, y = r, snap = "near")
    windows(); plot(PTO_060_100.crop)
    writeRaster(PTO_060_100.r,"raster_crop\\PTO_060_100.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.PTO_060_100()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #6.total phosphorous:% Read in explanatory variable raster maps - 
  fig.PTO_100_200 <- function() {
    PTO_100_200.r <- raster("raster_syn\\PTO_100_200.tif")
    windows(); plot(PTO_100_200.r, main="PTO_100_200")
    # Crop the raster to the shape fireferencele spatial extent:
    PTO_100_200.crop <- crop(x = PTO_100_200.r, y = r, snap = "near")
    windows(); plot(PTO_100_200.crop)
    writeRaster(PTO_100_200.r,"raster_crop\\PTO_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.PTO_100_200()
graphics.off()
----------------------------------------------------------------------------------------------------------   
  #1.clay:% Read in explanatory variable raster maps - 
  fig.CLY_000_005 <- function() {
    CLY_000_005.r <- raster("raster_syn\\CLY_000_005.tif")
    windows(); plot(CLY_000_005.r, main="CLY_000_005")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_000_005.crop <- crop(x = CLY_000_005.r, y = r, snap = "near")
    windows(); plot(CLY_000_005.crop)
    writeRaster(CLY_000_005.r,"raster_crop\\CLY_000_005.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_000_005()
graphics.off()  
  
----------------------------------------------------------------------------------------------------------   
  #2.clay:% Read in explanatory variable raster maps - 
  fig.CLY_015_030 <- function() {
    CLY_015_030.r <- raster("raster_syn\\CLY_015_030.tif")
    windows(); plot(CLY_015_030.r, main="CLY_015_030")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_015_030.crop <- crop(x = CLY_015_030.r, y = r, snap = "near")
    windows(); plot(CLY_015_030.crop)
    writeRaster(CLY_015_030.r,"raster_crop\\CLY_015_030.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_015_030()
graphics.off()  
----------------------------------------------------------------------------------------------------------   
  #3.clay:% Read in explanatory variable raster maps - 
  fig.CLY_030_060 <- function() {
    CLY_030_060.r <- raster("raster_syn\\CLY_030_060.tif")
    windows(); plot(CLY_030_060.r, main="CLY_030_060")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_030_060.crop <- crop(x = CLY_030_060.r, y = r, snap = "near")
    windows(); plot(CLY_030_060.crop)
    writeRaster(CLY_030_060.r,"raster_crop\\CLY_030_060.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_030_060()
graphics.off()  
----------------------------------------------------------------------------------------------------------   
  #4.clay:% Read in explanatory variable raster maps - 
  fig.CLY_060_100 <- function() {
    CLY_060_100.r <- raster("raster_syn\\CLY_060_100.tif")
    windows(); plot(CLY_060_100.r, main="CLY_060_100")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_060_100.crop <- crop(x = CLY_060_100.r, y = r, snap = "near")
    windows(); plot(CLY_060_100.crop)
    writeRaster(CLY_060_100.r,"raster_crop\\CLY_060_100.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_060_100()
graphics.off()  
----------------------------------------------------------------------------------------------------------   
  #5.clay:% Read in explanatory variable raster maps - 
  fig.CLY_100_200 <- function() {
    CLY_100_200.r <- raster("raster_syn\\CLY_100_200.tif")
    windows(); plot(CLY_100_200.r, main="CLY_100_200")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_100_200.crop <- crop(x = CLY_100_200.r, y = r, snap = "near")
    windows(); plot(CLY_100_200.crop)
    writeRaster(CLY_100_200.r,"raster_crop\\CLY_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_100_200()
graphics.off()  
----------------------------------------------------------------------------------------------------------   
  #6.clay:% Read in explanatory variable raster maps - 
  fig.CLY_005_015 <- function() {
    CLY_005_015.r <- raster("raster_syn\\CLY_005_015.tif")
    windows(); plot(CLY_005_015.r, main="CLY_005_015")
    # Crop the raster to the shape fireferencele spatial extent:
    CLY_005_015.crop <- crop(x = CLY_005_015.r, y = r, snap = "near")
    windows(); plot(CLY_005_015.crop)
    writeRaster(CLY_005_015.r,"raster_crop\\CLY_005_015.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.CLY_005_015()
graphics.off()  

----------------------------------------------------------------------------------------------------------  
  #1.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_000_005 <- function() {
    BDW_000_005.r <- raster("raster_syn\\BDW_000_005.tif")
    windows(); plot(BDW_000_005.r, main="BDW_000_005")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_000_005.crop <- crop(x = BDW_000_005.r, y = r, snap = "near")
    windows(); plot(BDW_000_005.crop)
    writeRaster(BDW_000_005.r,"raster_crop\\BDW_000_005.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_000_005()
graphics.off()  
----------------------------------------------------------------------------------------------------------  
  #2.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_005_015 <- function() {
    BDW_005_015.r <- raster("raster_syn\\BDW_005_015.tif")
    windows(); plot(BDW_005_015.r, main="BDW_005_015")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_005_015.crop <- crop(x = BDW_005_015.r, y = r, snap = "near")
    windows(); plot(BDW_005_015.crop)
    writeRaster(BDW_005_015.r,"raster_crop\\BDW_005_015.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_005_015()
graphics.off() 
----------------------------------------------------------------------------------------------------------  
  #3.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_015_030 <- function() {
    BDW_015_030.r <- raster("raster_syn\\BDW_015_030.tif")
    windows(); plot(BDW_015_030.r, main="BDW_015_030")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_015_030.crop <- crop(x = BDW_015_030.r, y = r, snap = "near")
    windows(); plot(BDW_015_030.crop)
    writeRaster(BDW_015_030.r,"raster_crop\\BDW_015_030.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_015_030()
graphics.off() 
----------------------------------------------------------------------------------------------------------  
  #4.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_030_060 <- function() {
    BDW_030_060.r <- raster("raster_syn\\BDW_030_060.tif")
    windows(); plot(BDW_030_060.r, main="BDW_030_060.")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_030_060..crop <- crop(x = BDW_030_060.r, y = r, snap = "near")
    windows(); plot(BDW_030_060..crop)
    writeRaster(BDW_030_060.r,"raster_crop\\BDW_030_060.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_030_060()
graphics.off()  

----------------------------------------------------------------------------------------------------------  
  #5.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_060_100 <- function() {
    BDW_060_100.r <- raster("raster_syn\\BDW_060_100.tif")
    windows(); plot(BDW_060_100.r, main="BDW_060_100.")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_060_100.crop <- crop(x = BDW_060_100.r, y = r, snap = "near")
    windows(); plot(BDW_060_100.crop)
    writeRaster(BDW_060_100.r,"raster_crop\\BDW_060_100.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_060_100()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #6.BDW:% Read in explanatory variable raster maps - 
  fig.BDW_100_200 <- function() {
    BDW_100_200.r <- raster("raster_syn\\BDW_100_200.tif")
    windows(); plot(BDW_100_200.r, main="BDW_100_200.")
    # Crop the raster to the shape fireferencele spatial extent:
    BDW_100_200.crop <- crop(x = BDW_100_200.r, y = r, snap = "near")
    windows(); plot(BDW_100_200.crop)
    writeRaster(BDW_100_200.r,"raster_crop\\BDW_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.BDW_100_200()
graphics.off()

----------------------------------------------------------------------------------------------------------  
  #1.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_000_005 <- function() {
    NTO_000_005.r <- raster("raster_syn\\NTO_000_005.tif")
    windows(); plot(NTO_000_005.r, main="NTO_000_005.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_000_005.crop <- crop(x = NTO_000_005.r, y = r, snap = "near")
    windows(); plot(NTO_000_005.crop)
    writeRaster(NTO_000_005.r,"raster_crop\\NTO_000_005.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_000_005()
graphics.off()

----------------------------------------------------------------------------------------------------------  
  #2.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_015_030 <- function() {
    NTO_015_030.r <- raster("raster_syn\\NTO_015_030.tif")
    windows(); plot(NTO_015_030.r, main="NTO_015_030.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_015_030.crop <- crop(x = NTO_015_030.r, y = r, snap = "near")
    windows(); plot(NTO_015_030.crop)
    writeRaster(NTO_015_030.r,"raster_crop\\NTO_015_030.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_015_030()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #3.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_030_060 <- function() {
    NTO_030_060.r <- raster("raster_syn\\NTO_030_060.tif")
    windows(); plot(NTO_030_060.r, main="NTO_030_060.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_030_060.crop <- crop(x = NTO_030_060.r, y = r, snap = "near")
    windows(); plot(NTO_030_060.crop)
    writeRaster(NTO_030_060.r,"raster_crop\\NTO_030_060.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_030_060()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #4.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_060_100 <- function() {
    NTO_060_100.r <- raster("raster_syn\\NTO_060_100.tif")
    windows(); plot(NTO_060_100.r, main="NTO_060_100.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_060_100.crop <- crop(x = NTO_060_100.r, y = r, snap = "near")
    windows(); plot(NTO_060_100.crop)
    writeRaster(NTO_060_100.r,"raster_crop\\NTO_060_100.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_060_100()
graphics.off()

----------------------------------------------------------------------------------------------------------  
  #5.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_100_200 <- function() {
    NTO_100_200.r <- raster("raster_syn\\NTO_100_200.tif")
    windows(); plot(NTO_100_200.r, main="NTO_100_200.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_100_200.crop <- crop(x = NTO_100_200.r, y = r, snap = "near")
    windows(); plot(NTO_100_200.crop)
    writeRaster(NTO_100_200.r,"raster_crop\\NTO_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_100_200()
graphics.off()

----------------------------------------------------------------------------------------------------------  
  #1.NTO:% Read in explanatory variable raster maps - 
  fig.NTO_100_200 <- function() {
    NTO_100_200.r <- raster("raster_syn\\NTO_100_200.tif")
    windows(); plot(NTO_100_200.r, main="NTO_100_200.")
    # Crop the raster to the shape fireferencele spatial extent:
    NTO_100_200.crop <- crop(x = NTO_100_200.r, y = r, snap = "near")
    windows(); plot(NTO_100_200.crop)
    writeRaster(NTO_100_200.r,"raster_crop\\NTO_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.NTO_100_200()

graphics.off()
----------------------------------------------------------------------------------------------------------  
  #1.soil available water capacity:% Read in explanatory variable raster maps - 
  fig.SAWC_000_005 <- function() {
    SAWC_000_005.r <- raster("raster_syn\\SAWC_000_005.tif")
    windows(); plot(SAWC_000_005.r, main="SAWC_000_005.")
    # Crop the raster to the shape fireferencele spatial extent:
    SAWC_000_005.crop <- crop(x = SAWC_000_005.r, y = r, snap = "near")
    windows(); plot(SAWC_000_005.crop)
    writeRaster(SAWC_000_005.r,"raster_crop\\SAWC_000_005.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.SAWC_000_005()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #2.soil available water capacity:% Read in explanatory variable raster maps - 
  fig.SAWC_015_030 <- function() {
    SAWC_015_030.r <- raster("raster_syn\\SAWC_015_030.tif")
    windows(); plot(SAWC_015_030.r, main="SAWC_015_030.")
    # Crop the raster to the shape fireferencele spatial extent:
    SAWC_015_030.crop <- crop(x = SAWC_015_030.r, y = r, snap = "near")
    windows(); plot(SAWC_015_030.crop)
    writeRaster(SAWC_015_030.r,"raster_crop\\SAWC_015_030.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.SAWC_015_030()
graphics.off()
----------------------------------------------------------------------------------------------------

  #2.soil available water capacity:% Read in explanatory variable raster maps - 
  fig.SAWC_030_060 <- function() {
    SAWC_030_060.r <- raster("raster_syn\\SAWC_030_060.tif")
    windows(); plot(SAWC_030_060.r, main="SAWC_030_060.")
    # Crop the raster to the shape fireferencele spatial extent:
    SAWC_030_060.crop <- crop(x = SAWC_030_060.r, y = r, snap = "near")
    windows(); plot(SAWC_030_060.crop)
    writeRaster(SAWC_030_060.r,"raster_crop\\SAWC_030_060.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.SAWC_030_060()
graphics.off()

----------------------------------------------------------------------------------------------------------  
  #3.soil available water capacity:% Read in explanatory variable raster maps - 
  fig.SAWC_060_100 <- function() {
    SAWC_060_100.r <- raster("raster_syn\\SAWC_060_100.tif")
    windows(); plot(SAWC_060_100.r, main="SAWC_060_100.")
    # Crop the raster to the shape fireferencele spatial extent:
    SAWC_060_100.crop <- crop(x = SAWC_060_100.r, y = r, snap = "near")
    windows(); plot(SAWC_060_100.crop)
    writeRaster(SAWC_060_100.r,"raster_crop\\SAWC_060_100.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.SAWC_060_100()
graphics.off()
----------------------------------------------------------------------------------------------------------  
  #4.soil available water capacity:% Read in explanatory variable raster maps - 
  fig.SAWC_100_200 <- function() {
    SAWC_100_200.r <- raster("raster_syn\\SAWC_100_200.tif")
    windows(); plot(SAWC_100_200.r, main="SAWC_100_200.")
    # Crop the raster to the shape fireferencele spatial extent:
    SAWC_100_200.crop <- crop(x = SAWC_100_200.r, y = r, snap = "near")
    windows(); plot(SAWC_100_200.crop)
    writeRaster(SAWC_100_200.r,"raster_crop\\SAWC_100_200.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.SAWC_100_200()
graphics.off()
------------------------------------------------------------------------------------------------------
#habitat Distance to habitat followed by propotion of suitable habitat
  
  #1.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.habit1 <- function() {
    habit1.r <- raster("raster_syn\\habit1.tif")
    windows(); plot(habit1.r, main="habit1.")
    # Crop the raster to the shape fireferencele spatial extent:
    habit1.crop <- crop(x = habit1.r, y = r, snap = "near")
    windows(); plot(habit1.crop)
    writeRaster(habit1.r,"raster_crop\\habit1.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.habit1() 
graphics.off()
  
------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #2.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.habit2 <- function() {
    habit2.r <- raster("raster_syn\\habit2.tif")
    windows(); plot(habit2.r, main="habit2.")
    # Crop the raster to the shape fireferencele spatial extent:
    habit2.crop <- crop(x = habit2.r, y = r, snap = "near")
    windows(); plot(habit2.crop)
    writeRaster(habit2.r,"raster_crop\\habit2.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.habit2()
graphics.off()
 
------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #3.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.habit3 <- function() {
    habit3.r <- raster("raster_syn\\habit3.tif")
    windows(); plot(habit3.r, main="habit3.")
    # Crop the raster to the shape fireferencele spatial extent:
    habit3.crop <- crop(x = habit3.r, y = r, snap = "near")
    windows(); plot(habit3.crop)
    writeRaster(habit3.r,"raster_crop\\habit3.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.habit3()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  ------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #3.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.suitable_0 <- function() {
    suitable_0.r <- raster("raster_syn\\suitable_0.tif")
    windows(); plot(suitable_0.r, main="suitable_0.")
    # Crop the raster to the shape fireferencele spatial extent:
    suitable_0.crop <- crop(x = suitable_0.r, y = r, snap = "near")
    windows(); plot(suitable_0.crop)
    writeRaster(suitable_0.r,"raster_crop\\suitable_0.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.suitable_0()  
graphics.off()

------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #1.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.suitable_1 <- function() {
    suitable_1.r <- raster("raster_syn\\suitable_1.tif")
    windows(); plot(suitable_1.r, main="suitable_1.")
    # Crop the raster to the shape fireferencele spatial extent:
    suitable_1.crop <- crop(x = suitable_1.r, y = r, snap = "near")
    windows(); plot(suitable_1.crop)
    writeRaster(suitable_1.r,"raster_crop\\suitable_1.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.suitable_1()  
graphics.off()
 
------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #2.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.suitable_2 <- function() {
    suitable_2.r <- raster("raster_syn\\suitable_2.tif")
    windows(); plot(suitable_2.r, main="suitable_2.")
    # Crop the raster to the shape fireferencele spatial extent:
    suitable_2.crop <- crop(x = suitable_2.r, y = r, snap = "near")
    windows(); plot(suitable_2.crop)
    writeRaster(suitable_2.r,"raster_crop\\suitable_2.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.suitable_2()
graphics.off()
------------------------------------------------------------------------------------------------------
  #habitat Distance to habitat followed by propotion of suitable habitat
  
  #3.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.suitable_3<- function() {
    suitable_3.r <- raster("raster_syn\\suitable_3.tif")
    windows(); plot(suitable_3.r, main="suitable_3.")
    # Crop the raster to the shape fireferencele spatial extent:
    suitable_3.crop <- crop(x = suitable_3.r, y = r, snap = "near")
    windows(); plot(suitable_3.crop)
    writeRaster(suitable_3.r,"raster_crop\\suitable_3.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.suitable_3()
graphics.off()
-----------------------------------------------------------------------------------------------------------
  ------------------------------------------------------------------------------------------------------
  #habitat Distance to roads follwoed by road densities
  
  #1.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_bridleway<- function() {
    distance_bridleway.r <- raster("raster_syn\\distance_bridleway.tif")
    windows(); plot(distance_bridleway.r, main="distance_bridleway.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_bridleway.crop <- crop(x = distance_bridleway.r, y = r, snap = "near")
    windows(); plot(distance_bridleway.crop)
    writeRaster(distance_bridleway.r,"raster_crop\\distance_bridleway.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_bridleway() 
graphics.off()
------------------------------------------------------------------------------------------------------
  #habitat Distance to roads follwoed by road densities
  
  #2.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_cycleway <- function() {
    distance_cycleway.r <- raster("raster_syn\\distance_cycleway.tif")
    windows(); plot(distance_cycleway.r, main="distance_cycleway.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_cycleway.crop <- crop(x = distance_cycleway.r, y = r, snap = "near")
    windows(); plot(distance_cycleway.crop)
    writeRaster(distance_cycleway.r,"raster_crop\\distance_cycleway.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_cycleway() 
graphics.off()
------------------------------------------------------------------------------------------------------

  #3.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_footway <- function() {
    distance_footway.r <- raster("raster_syn\\distance_footway.tif")
    windows(); plot(distance_footway.r, main="distance_footway.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_footway.crop <- crop(x = distance_footway.r, y = r, snap = "near")
    windows(); plot(distance_footway.crop)
    writeRaster(distance_footway.r,"raster_crop\\distance_footway.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_footway() 
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #4.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_motorwayandlink <- function() {
    distance_motorwayandlink.r <- raster("raster_syn\\distance_motorwayandlink.tif")
    windows(); plot(distance_motorwayandlink.r, main="distance_motorwayandlink.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_motorwayandlink.crop <- crop(x = distance_motorwayandlink.r, y = r, snap = "near")
    windows(); plot(distance_motorwayandlink.crop)
    writeRaster(distance_motorwayandlink.r,"raster_crop\\distance_motorwayandlink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_motorwayandlink() 
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #5.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_path<- function() {
    distance_path.r <- raster("raster_syn\\distance_path.tif")
    windows(); plot(distance_path.r, main="distance_path.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_path.crop <- crop(x = distance_path.r, y = r, snap = "near")
    windows(); plot(distance_path.crop)
    writeRaster(distance_path.r,"raster_crop\\distance_path.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_path() 
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #6.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_pedestrian <- function() {
    distance_pedestrian.r <- raster("raster_syn\\distance_pedestrian.tif")
    windows(); plot(distance_pedestrian.r, main="distance_pedestrian.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_pedestrian.crop <- crop(x = distance_pedestrian.r, y = r, snap = "near")
    windows(); plot(distance_pedestrian.crop)
    writeRaster(distance_pedestrian.r,"raster_crop\\distance_pedestrian.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_pedestrian()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #7.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_primaryandlink <- function() {
    distance_primaryandlink.r <- raster("raster_syn\\distance_primaryandlink.tif")
    windows(); plot(distance_primaryandlink.r, main="distance_primaryandlink.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_primaryandlink.crop <- crop(x = distance_primaryandlink.r, y = r, snap = "near")
    windows(); plot(distance_primaryandlink.crop)
    writeRaster(distance_primaryandlink.r,"raster_crop\\distance_primaryandlink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_primaryandlink() 
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #8.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_residentil <- function() {
    distance_residentil.r <- raster("raster_syn\\distance_residentil.tif")
    windows(); plot(distance_residentil.r, main="distance_residentil.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_residentil.crop <- crop(x = distance_residentil.r, y = r, snap = "near")
    windows(); plot(distance_residentil.crop)
    writeRaster(distance_residentil.r,"raster_crop\\distance_residentil.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_residentil()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #9.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_secondaryandlink <- function() {
    distance_secondaryandlink.r <- raster("raster_syn\\distance_secondaryandlink.tif")
    windows(); plot(distance_secondaryandlink.r, main="distance_secondaryandlink.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_secondaryandlink.crop <- crop(x = distance_secondaryandlink.r, y = r, snap = "near")
    windows(); plot(distance_secondaryandlink.crop)
    writeRaster(distance_secondaryandlink.r,"raster_crop\\distance_secondaryandlink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_secondaryandlink()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #10.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_service <- function() {
    distance_service.r <- raster("raster_syn\\distance_service.tif")
    windows(); plot(distance_service.r, main="distance_service.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_service.crop <- crop(x = distance_service.r, y = r, snap = "near")
    windows(); plot(distance_service.crop)
    writeRaster(distance_service.r,"raster_crop\\distance_service.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_service()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #11.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_steps <- function() {
    distance_steps.r <- raster("raster_syn\\distance_steps.tif")
    windows(); plot(distance_steps.r, main="distance_steps.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_steps.crop <- crop(x = distance_steps.r, y = r, snap = "near")
    windows(); plot(distance_steps.crop)
    writeRaster(distance_steps.r,"raster_crop\\distance_steps.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_steps()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #12.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_tertiaryandlink <- function() {
    distance_tertiaryandlink.r <- raster("raster_syn\\distance_tertiaryandlink.tif")
    windows(); plot(distance_tertiaryandlink.r, main="distance_tertiaryandlink.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_tertiaryandlink.crop <- crop(x = distance_tertiaryandlink.r, y = r, snap = "near")
    windows(); plot(distance_tertiaryandlink.crop)
    writeRaster(distance_tertiaryandlink.r,"raster_crop\\distance_tertiaryandlink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_tertiaryandlink()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #13.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_track <- function() {
    distance_track.r <- raster("raster_syn\\distance_track.tif")
    windows(); plot(distance_track.r, main="distance_track.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_track.crop <- crop(x = distance_track.r, y = r, snap = "near")
    windows(); plot(distance_track.crop)
    writeRaster(distance_track.r,"raster_crop\\distance_track.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_track()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #14.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_trunkandlink <- function() {
    distance_trunkandlink.r <- raster("raster_syn\\distance_trunkandlink.tif")
    windows(); plot(distance_trunkandlink.r, main="distance_trunkandlink.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_trunkandlink.crop <- crop(x = distance_trunkandlink.r, y = r, snap = "near")
    windows(); plot(distance_trunkandlink.crop)
    writeRaster(distance_trunkandlink.r,"raster_crop\\distance_trunkandlink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_trunkandlink()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #14.distance to habitat catogories 1-3:% Read in explanatory variable raster maps - 
  fig.distance_unclassified <- function() {
    distance_unclassified.r <- raster("raster_syn\\distance_unclassified.tif")
    windows(); plot(distance_unclassified.r, main="distance_unclassified.")
    # Crop the raster to the shape fireferencele spatial extent:
    distance_unclassified.crop <- crop(x = distance_unclassified.r, y = r, snap = "near")
    windows(); plot(distance_unclassified.crop)
    writeRaster(distance_unclassified.r,"raster_crop\\distance_unclassified.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.distance_unclassified()
graphics.off()

------------------------------------------------------------------------------------------------------
  
  #1.densities road 
  fig.drains <- function() {
    drains.r <- raster("raster_syn\\drains.tif")
    windows(); plot(drains.r, main="drains.")
    # Crop the raster to the shape fireferencele spatial extent:
    drains.crop <- crop(x = drains.r, y = r, snap = "near")
    windows(); plot(drains.crop)
    writeRaster(drains.r,"raster_crop\\drains.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.drains()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #2.densities road 
  fig.footway <- function() {
    footway.r <- raster("raster_syn\\footway.tif")
    windows(); plot(footway.r, main="footway.")
    # Crop the raster to the shape fireferencele spatial extent:
    footway.crop <- crop(x = footway.r, y = r, snap = "near")
    windows(); plot(footway.crop)
    writeRaster(footway.r,"raster_crop\\footway.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.footway()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #3.densities road 
  fig.motorway <- function() {
    motorway.r <- raster("raster_syn\\motorway.tif")
    windows(); plot(motorway.r, main="motorway.")
    # Crop the raster to the shape fireferencele spatial extent:
    motorway.crop <- crop(x = motorway.r, y = r, snap = "near")
    windows(); plot(motorway.crop)
    writeRaster(motorway.r,"raster_crop\\motorway.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.motorway()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #4.densities road 
  fig.motorwaylink <- function() {
    motorwaylink.r <- raster("raster_syn\\motorwaylink.tif")
    windows(); plot(motorwaylink.r, main="motorwaylink.")
    # Crop the raster to the shape fireferencele spatial extent:
    motorwaylink.crop <- crop(x = motorwaylink.r, y = r, snap = "near")
    windows(); plot(motorwaylink.crop)
    writeRaster(motorwaylink.r,"raster_crop\\motorwaylink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.motorwaylink()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #5.densities road 
  fig.path <- function() {
    path.r <- raster("raster_syn\\path.tif")
    windows(); plot(path.r, main="path.")
    # Crop the raster to the shape fireferencele spatial extent:
    path.crop <- crop(x = path.r, y = r, snap = "near")
    windows(); plot(path.crop)
    writeRaster(path.r,"raster_crop\\path.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.path()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #5.densities road 
  fig.pedestrian <- function() {
    pedestrian.r <- raster("raster_syn\\pedestrian.tif")
    windows(); plot(pedestrian.r, main="pedestrian.")
    # Crop the raster to the shape fireferencele spatial extent:
    pedestrian.crop <- crop(x = pedestrian.r, y = r, snap = "near")
    windows(); plot(pedestrian.crop)
    writeRaster(pedestrian.r,"raster_crop\\pedestrian.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.pedestrian()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #6.densities road 
  fig.primary <- function() {
    primary.r <- raster("raster_syn\\primary.tif")
    windows(); plot(primary.r, main="primary.")
    # Crop the raster to the shape fireferencele spatial extent:
    primary.crop <- crop(x = primary.r, y = r, snap = "near")
    windows(); plot(primary.crop)
    writeRaster(primary.r,"raster_crop\\primary.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.primary()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #7.densities road 
  fig.residential <- function() {
    residential.r <- raster("raster_syn\\residential.tif")
    windows(); plot(residential.r, main="residential.")
    # Crop the raster to the shape fireferencele spatial extent:
    residential.crop <- crop(x = residential.r, y = r, snap = "near")
    windows(); plot(residential.crop)
    writeRaster(residential.r,"raster_crop\\residential.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.residential()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #8.densities road 
  fig.roads_motor <- function() {
    roads_motor.r <- raster("raster_syn\\roads_motor.tif")
    windows(); plot(roads_motor.r, main="roads_motor.")
    # Crop the raster to the shape fireferencele spatial extent:
    roads_motor.crop <- crop(x = roads_motor.r, y = r, snap = "near")
    windows(); plot(roads_motor.crop)
    writeRaster(roads_motor.r,"raster_crop\\roads_motor.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.roads_motor()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #9.densities road 
  fig.roads_other <- function() {
    roads_other.r <- raster("raster_syn\\roads_other.tif")
    windows(); plot(roads_other.r, main="roads_other.")
    # Crop the raster to the shape fireferencele spatial extent:
    roads_other.crop <- crop(x = roads_other.r, y = r, snap = "near")
    windows(); plot(roads_other.crop)
    writeRaster(roads_other.r,"raster_crop\\roads_other.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.roads_other()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #10.densities road 
  fig.secondry <- function() {
    secondry.r <- raster("raster_syn\\secondry.tif")
    windows(); plot(secondry.r, main="secondry.")
    # Crop the raster to the shape fireferencele spatial extent:
    secondry.crop <- crop(x = secondry.r, y = r, snap = "near")
    windows(); plot(secondry.crop)
    writeRaster(secondry.r,"raster_crop\\secondry.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.secondry()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #11.densities road 
  fig.secondry_link <- function() {
    secondry_link.r <- raster("raster_syn\\secondry_link.tif")
    windows(); plot(secondry_link.r, main="secondry_link.")
    # Crop the raster to the shape fireferencele spatial extent:
    secondry_link.crop <- crop(x = secondry_link.r, y = r, snap = "near")
    windows(); plot(secondry_link.crop)
    writeRaster(secondry_link.r,"raster_crop\\secondry_link.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.secondry_link()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #12.densities road 
  fig.tertiary <- function() {
    tertiary.r <- raster("raster_syn\\tertiary.tif")
    windows(); plot(tertiary.r, main="tertiary.")
    # Crop the raster to the shape fireferencele spatial extent:
    tertiary.crop <- crop(x = tertiary.r, y = r, snap = "near")
    windows(); plot(tertiary.crop)
    writeRaster(tertiary.r,"raster_crop\\tertiary.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.tertiary()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #12.densities road 
  fig.tertiary_liink <- function() {
    tertiary_liink.r <- raster("raster_syn\\tertiary_liink.tif")
    windows(); plot(tertiary_liink.r, main="tertiary_liink.")
    # Crop the raster to the shape fireferencele spatial extent:
    tertiary_liink.crop <- crop(x = tertiary_liink.r, y = r, snap = "near")
    windows(); plot(tertiary_liink.crop)
    writeRaster(tertiary_liink.r,"raster_crop\\tertiary_liink.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.tertiary_liink()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #13.densities road 
  fig.track <- function() {
    track.r <- raster("raster_syn\\track.tif")
    windows(); plot(track.r, main="track.")
    # Crop the raster to the shape fireferencele spatial extent:
    track.crop <- crop(x = track.r, y = r, snap = "near")
    windows(); plot(track.crop)
    writeRaster(track.r,"raster_crop\\track.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.track()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #14.densities road 
  fig.trunck <- function() {
    trunck.r <- raster("raster_syn\\trunck.tif")
    windows(); plot(trunck.r, main="trunck.")
    # Crop the raster to the shape fireferencele spatial extent:
    trunck.crop <- crop(x = trunck.r, y = r, snap = "near")
    windows(); plot(trunck.crop)
    writeRaster(trunck.r,"raster_crop\\trunck.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.trunck()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #15.densities road 
  fig.unclassified <- function() {
    unclassified.r <- raster("raster_syn\\unclassified.tif")
    windows(); plot(unclassified.r, main="unclassified.")
    # Crop the raster to the shape fireferencele spatial extent:
    unclassified.crop <- crop(x = unclassified.r, y = r, snap = "near")
    windows(); plot(unclassified.crop)
    writeRaster(unclassified.r,"raster_crop\\unclassified.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.unclassified()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #1.river densities
  fig.riverbank <- function() {
    riverbank.r <- raster("raster_syn\\riverbank.tif")
    windows(); plot(riverbank.r, main="riverbank.")
    # Crop the raster to the shape fireferencele spatial extent:
    riverbank.crop <- crop(x = riverbank.r, y = r, snap = "near")
    windows(); plot(riverbank.crop)
    writeRaster(riverbank.r,"raster_crop\\riverbank.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.riverbank()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #2.river densities
  fig.rivers <- function() {
    rivers.r <- raster("raster_syn\\rivers.tif")
    windows(); plot(rivers.r, main="rivers.")
    # Crop the raster to the shape fireferencele spatial extent:
    rivers.crop <- crop(x = rivers.r, y = r, snap = "near")
    windows(); plot(rivers.crop)
    writeRaster(rivers.r,"raster_crop\\rivers.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.rivers()
graphics.off()
------------------------------------------------------------------------------------------------------
  
  #3.river densities
  fig.stream <- function() {
    stream.r <- raster("raster_syn\\stream.tif")
    windows(); plot(stream.r, main="stream.")
    # Crop the raster to the shape fireferencele spatial extent:
    stream.crop <- crop(x = stream.r, y = r, snap = "near")
    windows(); plot(stream.crop)
    writeRaster(stream.r,"raster_crop\\stream.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.stream()
graphics.off()
------------------------------------------------------------------------------------------------------
#slopeSlope is always the ratio of change in vertical distance divided by the change in horizontal distance. 
  #In this ratio, the units of distance cancel out. Therefore, slope has no units
  #
  fig.Slope_percent_61 <- function() {
    Slope_percent_61.r <- raster("raster_syn\\Slope_percent_61.tif")
    windows(); plot(Slope_percent_61.r, main="Slope_percent_61")
    # Crop the raster to the shape fireferencele spatial extent:
    Slope_percent_61.crop <- crop(x = Slope_percent_61.r, y = r, snap = "near")
    windows(); plot(Slope_percent_61.crop)
    writeRaster(Slope_percent_61.r,"raster_crop\\Slope_percent_61.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.Slope_percent_61()
graphics.off()
------------------------------------------------------------------------------------------------------
#topographic wetness index
  #
  fig.TopoWI <- function() {
    TopoWI.r <- raster("raster_syn\\TopoWI.tif")
    windows(); plot(TopoWI.r, main="TopoWI")
    # Crop the raster to the shape fireferencele spatial extent:
    TopoWI.crop <- crop(x = TopoWI.r, y = r, snap = "near")
    windows(); plot(TopoWI.crop)
    writeRaster(TopoWI.r,"raster_crop\\TopoWI.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.TopoWI()
graphics.off()
------------------------------------------------------------------------------------------------------
  #topographic position index
  #
  fig.TopoPI <- function() {
    TopoPI.r <- raster("raster_syn\\TopoPI.tif")
    windows(); plot(TopoPI.r, main="TopoPI")
    # Crop the raster to the shape fireferencele spatial extent:
    TopoPI.crop <- crop(x = TopoPI.r, y = r, snap = "near")
    windows(); plot(TopoPI.crop)
    writeRaster(TopoPI.r,"raster_crop\\TopoPI.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.TopoPI()
graphics.off()
------------------------------------------------------------------------------------------------------
  #aspect
  #
  fig.aspect91 <- function() {
    aspect91.r <- raster("raster_syn\\aspect91.tif")
    windows(); plot(aspect91.r, main="aspect91")
    # Crop the raster to the shape fireferencele spatial extent:
    aspect91.crop <- crop(x = aspect91.r, y = r, snap = "near")
    windows(); plot(aspect91.crop)
    writeRaster(aspect91.r,"raster_crop\\aspect91.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.aspect91()
------------------------------------------------------------------------------------------------------
  #topographic position index
  #
  fig.dem <- function() {
    dem.r <- raster("raster_syn\\dem.tif")
    windows(); plot(dem.r, main="dem")
    # Crop the raster to the shape fireferencele spatial extent:
    dem.crop <- crop(x = dem.r, y = r, snap = "near")
    windows(); plot(dem.crop)
    writeRaster(dem.r,"raster_crop\\dem.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.dem()
graphics.off()
------------------------------------------------------------------------------------------------------
  #topographic position index
  #
  fig.hpop <- function() {
    hpop.r <- raster("raster_syn\\hpop.tif")
    windows(); plot(hpop.r, main="hpop")
    # Crop the raster to the shape fireferencele spatial extent:
    hpop.crop <- crop(x = hpop.r, y = r, snap = "near")
    windows(); plot(hpop.crop)
    writeRaster(hpop.r,"raster_crop\\hpop.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.hpop()
graphics.off()
------------------------------------------------------------------------------------------------------
  #lot density
  #
  fig.lot_density <- function() {
    lot_density.r <- raster("raster_syn\\lot_density.tif")
    windows(); plot(lot_density.r, main="lot_density")
    # Crop the raster to the shape fireferencele spatial extent:
    lot_density.crop <- crop(x = lot_density.r, y = r, snap = "near")
    windows(); plot(lot_density.crop)
    writeRaster(lot_density.r,"raster_crop\\lot_density.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.lot_density()
graphics.off()
------------------------------------------------------------------------------------------------------
  #temp
  #
  fig.temp_minimum1996_2005 <- function() {
    temp_minimum1996_2005.r <- raster("raster_syn\\temp_minimum1996_2005.tif")
    windows(); plot(temp_minimum1996_2005.r, main="temp_minimum1996_2005")
    # Crop the raster to the shape fireferencele spatial extent:
    temp_minimum1996_2005.crop <- crop(x = temp_minimum1996_2005.r, y = r, snap = "near")
    windows(); plot(temp_minimum1996_2005.crop)
    writeRaster(temp_minimum1996_2005.r,"raster_crop\\temp_minimum1996_2005.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.temp_minimum1996_2005()
graphics.off()
------------------------------------------------------------------------------------------------------
  #temp
  #
  fig.temp_maximum1990_2000 <- function() {
    temp_maximum1990_2000.r <- raster("raster_syn\\temp_maximum1990_2000.tif")
    windows(); plot(temp_maximum1990_2000.r, main="temp_maximum1990_2000")
    # Crop the raster to the shape fireferencele spatial extent:
    temp_maximum1990_2000.crop <- crop(x = temp_maximum1990_2000.r, y = r, snap = "near")
    windows(); plot(temp_maximum1990_2000.crop)
    writeRaster(temp_maximum1990_2000.r,"raster_crop\\temp_maximum1990_2000.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.temp_maximum1990_2000()
graphics.off()
------------------------------------------------------------------------------------------------------
  #temp
  #
  fig.tempmeanMGA65 <- function() {
    tempmeanMGA65.r <- raster("raster_syn\\tempmeanMGA65.tif")
    windows(); plot(tempmeanMGA65.r, main="tempmeanMGA65")
    # Crop the raster to the shape fireferencele spatial extent:
    tempmeanMGA65.crop <- crop(x = tempmeanMGA65.r, y = r, snap = "near")
    windows(); plot(tempmeanMGA65.crop)
    writeRaster(tempmeanMGA65.r,"raster_crop\\tempmeanMGA65.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.tempmeanMGA65()
graphics.off()
------------------------------------------------------------------------------------------------------
  #rainfall
  #rainfall
  #
  fig.rainfallmeanannual <- function() {
    rainfallmeanannual.r <- raster("raster_syn\\rainfallmeanannual.tif")
    windows(); plot(rainfallmeanannual.r, main="rainfallmeanannual")
    # Crop the raster to the shape fireferencele spatial extent:
    rainfallmeanannual.crop <- crop(x = rainfallmeanannual.r, y = r, snap = "near")
    windows(); plot(rainfallmeanannual.crop)
    writeRaster(rainfallmeanannual.r,"raster_crop\\rainfallmeanannual.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.rainfallmeanannual()
graphics.off()
------------------------------------------------------------------------------------------------------
  #rainfall
  
  #
  fig.foliagePC <- function() {
    foliagePC.r <- raster("raster_syn\\foliagePC.tif")
    windows(); plot(foliagePC.r, main="foliagePC")
    # Crop the raster to the shape fireferencele spatial extent:
    foliagePC.crop <- crop(x = foliagePC.r, y = r, snap = "near")
    windows(); plot(foliagePC.crop)
    writeRaster(foliagePC.r,"raster_crop\\foliagePC.TIF", overwrite=TRUE)
    # Plot to check:
    
    points(acsel, cex = 0.4, col = 'red', pch = 15)
  }
fig.foliagePC()
graphics.off()
----------------------------------
 
# Make a ppp object for spatstat: Kernel density
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
bandwidth <- 1500

# Select only those koalas within the owin:
# id <- inside.owin(x = mydata$x, y = mydata$y, w = dat.w)
# mydata <- mydata[id,]

windows(); plot(studyarea.shp, axes = TRUE)
points(x = acsel[,1], y = acsel[,2])

# Make a ppp object:
dat.ppp <- ppp(x = acsel[,1], y = acsel[,2], window = dat.w)

windows(); plot(dat.ppp, axes = TRUE)
---------------------------------------------------------------------
#gaussian kernel function. Its bandwidth defines the kernel’s window extent.
  dat.den <- density(dat.ppp,sigma=2000, diggle = TRUE,dimyx = dimyx)
windows();plot(dat.den, main=NULL)
contour(dat.den, add=TRUE)
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))
contour(dat.den, add=TRUE)
windows(); plot(dat.den)
---------------------------------------------------------------------
# Work out density for dat.ppp (using Diggle's edge correction):
dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)

# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre:
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))
contour(dat.den, add=TRUE)

windows(); plot(dat.den)

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 0.003, length = 6)
col <- brewer.pal(n = 5, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
points(acsel, cex = 0.4, col = 'red', pch = 15)
------------------------------------------------------------------------------------------------------------------------------------
  ## rhohat - total phosphorous:Computes a smoothing estimate of the intensity of a point process, as a function of a (continuous) spatial covariate.
  #Foliage projective cover
  # Convert tpo.r into a spatstat image object:
  fig.PTO_000_005.rho <- function() {
    PTO_000_005.im <- as.im(PTO_000_005.r)
    windows(); plot(PTO_000_005.im, axes = TRUE)
       # What is the nature of the association between koala sight locations and total phosphorous?
    PTO_000_005.rho <- rhohat(object = dat.ppp, covariate = PTO_000_005.im)
    windows(); plot(PTO_000_005.rho, xlab = "PTO_000_005 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\PTO_000_005",type = c( "png"),device = dev.cur())
      }
fig.PTO_000_005.rho() 
graphics.off()
# =================================================================================================================================
fig.PTO_005_015.rho <- function() {
  PTO_005_015.im <- as.im(PTO_005_015.r)
  windows(); plot(PTO_005_015.im, axes = TRUE)
  # What is the nature of the association between koala sight locations and total phosphorous?
  PTO_005_015.rho <- rhohat(object = dat.ppp, covariate = PTO_005_015.im)
  windows(); plot(PTO_005_015.rho, xlab = "PTO_005_015 (units)", main = "")
  savePlot(filename = "raster_crop\\rhohat\\PTO_005_015",type = c( "png"),device = dev.cur())
}
fig.PTO_005_015.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.PTO_015_030.rho <- function() {
    PTO_015_030.im <- as.im(PTO_015_030.r)
    windows(); plot(PTO_015_030.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    PTO_015_030.rho <- rhohat(object = dat.ppp, covariate = PTO_015_030.im)
    windows(); plot(PTO_015_030.rho, xlab = "PTO_015_030 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\PTO_015_030",type = c( "png"),device = dev.cur())
  }
fig.PTO_015_030.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.PTO_030_060.rho <- function() {
    PTO_030_060.im <- as.im(PTO_030_060.r)
    windows(); plot(PTO_030_060.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    PTO_030_060.rho <- rhohat(object = dat.ppp, covariate = PTO_030_060.im)
    windows(); plot(PTO_030_060.rho, xlab = "PTO_030_060 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\PTO_030_060",type = c( "png"),device = dev.cur())
  }
fig.PTO_030_060.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.PTO_060_100.rho <- function() {
    PTO_060_100.im <- as.im(PTO_060_100.r)
    windows(); plot(PTO_060_100.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    PTO_060_100.rho <- rhohat(object = dat.ppp, covariate = PTO_060_100.im)
    windows(); plot(PTO_060_100.rho, xlab = "PTO_060_100 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\PTO_060_100",type = c( "png"),device = dev.cur())
  }
fig.PTO_060_100.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.PTO_100_200.rho <- function() {
    PTO_100_200.im <- as.im(PTO_100_200.r)
    windows(); plot(PTO_100_200.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    PTO_100_200.rho <- rhohat(object = dat.ppp, covariate = PTO_100_200.im)
    windows(); plot(PTO_100_200.rho, xlab = "PTO_100_200 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\PTO_100_200",type = c( "png"),device = dev.cur())
  }
fig.PTO_100_200.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
fig.CLY_000_005.rho <- function() {
CLY_000_005.im <- as.im(CLY_000_005.r)
windows(); plot(CLY_000_005.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
CLY_000_005.rho <- rhohat(object = dat.ppp, covariate = CLY_000_005.im)
windows(); plot(CLY_000_005.rho, xlab = "CLY_000_005 (units)", main = "")
savePlot(filename = "raster_crop\\rhohat\\CLY_000_005",type = c( "png"),device = dev.cur())
}
fig.CLY_000_005.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.CLY_015_030.rho <- function() {
    CLY_015_030.im <- as.im(CLY_015_030.r)
    windows(); plot(CLY_015_030.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    CLY_015_030.rho <- rhohat(object = dat.ppp, covariate = CLY_015_030.im)
    windows(); plot(CLY_015_030.rho, xlab = "CLY_015_030 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\CLY_015_030",type = c( "png"),device = dev.cur())
  }
fig.CLY_015_030.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.CLY_030_060.rho <- function() {
    CLY_030_060.im <- as.im(CLY_030_060.r)
    windows(); plot(CLY_030_060.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    CLY_030_060.rho <- rhohat(object = dat.ppp, covariate = CLY_030_060.im)
    windows(); plot(CLY_030_060.rho, xlab = "CLY_030_060 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\CLY_030_060",type = c( "png"),device = dev.cur())
  }
fig.CLY_030_060.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.CLY_060_100.rho <- function() {
    CLY_060_100.im <- as.im(CLY_060_100.r)
    windows(); plot(CLY_060_100.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    CLY_060_100.rho <- rhohat(object = dat.ppp, covariate = CLY_060_100.im)
    windows(); plot(CLY_060_100.rho, xlab = "CLY_060_100 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\CLY_060_100",type = c( "png"),device = dev.cur())
  }
fig.CLY_060_100.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.CLY_100_200.rho <- function() {
    CLY_100_200.im <- as.im(CLY_100_200.r)
    windows(); plot(CLY_100_200.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    CLY_100_200.rho <- rhohat(object = dat.ppp, covariate = CLY_100_200.im)
    windows(); plot(CLY_100_200.rho, xlab = "CLY_100_200 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\CLY_100_200",type = c( "png"),device = dev.cur())
  }
fig.CLY_100_200.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.CLY_005_015.rho <- function() {
    CLY_005_015.im <- as.im(CLY_005_015.r)
    windows(); plot(CLY_005_015.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    CLY_005_015.rho <- rhohat(object = dat.ppp, covariate = CLY_005_015.im)
    windows(); plot(CLY_005_015.rho, xlab = "CLY_005_015 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\CLY_005_015",type = c( "png"),device = dev.cur())
  }
fig.CLY_005_015.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.BDW_000_005.rho <- function() {
    BDW_000_005.im <- as.im(BDW_000_005.r)
    windows(); plot(BDW_000_005.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    BDW_000_005.rho <- rhohat(object = dat.ppp, covariate = BDW_000_005.im)
    windows(); plot(BDW_000_005.rho, xlab = "BDW_000_005 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\BDW_000_005",type = c( "png"),device = dev.cur())
  }
fig.BDW_000_005.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.BDW_005_015.rho <- function() {
    BDW_005_015.im <- as.im(BDW_005_015.r)
    windows(); plot(BDW_005_015.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    BDW_005_015.rho <- rhohat(object = dat.ppp, covariate = BDW_005_015.im)
    windows(); plot(BDW_005_015.rho, xlab = "BDW_005_015 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\BDW_005_015",type = c( "png"),device = dev.cur())
  }
fig.BDW_005_015.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.BDW_015_030.rho <- function() {
    BDW_015_030.im <- as.im(BDW_015_030.r)
    windows(); plot(BDW_015_030.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    BDW_015_030.rho <- rhohat(object = dat.ppp, covariate = BDW_015_030.im)
    windows(); plot(BDW_015_030.rho, xlab = "BDW_015_030 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\BDW_015_030",type = c( "png"),device = dev.cur())
  }
fig.BDW_015_030.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.BDW_030_060.rho <- function() {
    BDW_030_060.im <- as.im(BDW_030_060.r)
    windows(); plot(BDW_030_060..im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    BDW_030_060.rho <- rhohat(object = dat.ppp, covariate = BDW_030_060.im)
    windows(); plot(BDW_030_060.rho, xlab = "BDW_030_060 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\BDW_030_060",type = c( "png"),device = dev.cur())
  }
fig.BDW_030_060..rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.BDW_060_100.rho <- function() {
    BDW_060_100.im <- as.im(BDW_060_100.r)
    windows(); plot(BDW_060_100.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    BDW_060_100.rho <- rhohat(object = dat.ppp, covariate = BDW_060_100.im)
    windows(); plot(BDW_060_100.rho, xlab = "BDW_060_100 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\BDW_060_100",type = c( "png"),device = dev.cur())
  }
fig.BDW_060_100.rho() 
graphics.off()

-----------------------------------------------------------------------------------------------------------------------------
  fig.NTO_000_005.rho <- function() {
    NTO_000_005.im <- as.im(NTO_000_005.r)
    windows(); plot(NTO_000_005.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    NTO_000_005.rho <- rhohat(object = dat.ppp, covariate = NTO_000_005.im)
    windows(); plot(NTO_000_005.rho, xlab = "NTO_000_005 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\NTO_000_005",type = c( "png"),device = dev.cur())
  }
fig.NTO_000_005.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.NTO_015_030.rho <- function() {
    NTO_015_030.im <- as.im(NTO_015_030.r)
    windows(); plot(NTO_015_030.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    NTO_015_030.rho <- rhohat(object = dat.ppp, covariate = NTO_015_030.im)
    windows(); plot(NTO_015_030.rho, xlab = "NTO_015_030 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\NTO_015_030",type = c( "png"),device = dev.cur())
  }
fig.NTO_015_030.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.NTO_030_060.rho <- function() {
    NTO_030_060.im <- as.im(NTO_030_060.r)
    windows(); plot(NTO_030_060.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    NTO_030_060.rho <- rhohat(object = dat.ppp, covariate = NTO_030_060.im)
    windows(); plot(NTO_030_060.rho, xlab = "NTO_030_060 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\NTO_030_060",type = c( "png"),device = dev.cur())
  }
fig.NTO_030_060.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.NTO_060_100.rho <- function() {
    NTO_060_100.im <- as.im(NTO_060_100.r)
    windows(); plot(NTO_060_100.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    NTO_060_100.rho <- rhohat(object = dat.ppp, covariate = NTO_060_100.im)
    windows(); plot(NTO_060_100.rho, xlab = "NTO_060_100 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\NTO_060_100",type = c( "png"),device = dev.cur())
  }
fig.NTO_060_100.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.NTO_100_200.rho <- function() {
    NTO_100_200.im <- as.im(NTO_100_200.r)
    windows(); plot(NTO_100_200.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    NTO_100_200.rho <- rhohat(object = dat.ppp, covariate = NTO_100_200.im)
    windows(); plot(NTO_100_200.rho, xlab = "NTO_100_200 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\NTO_100_200",type = c( "png"),device = dev.cur())
  }
fig.NTO_100_200.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.SAWC_000_005.rho <- function() {
    SAWC_000_005.im <- as.im(SAWC_000_005.r)
    windows(); plot(SAWC_000_005.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    SAWC_000_005.rho <- rhohat(object = dat.ppp, covariate = SAWC_000_005.im)
    windows(); plot(SAWC_000_005.rho, xlab = "SAWC_000_005 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\SAWC_000_005",type = c( "png"),device = dev.cur())
  }
fig.SAWC_000_005.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.SAWC_015_030.rho <- function() {
    SAWC_015_030.im <- as.im(SAWC_015_030.r)
    windows(); plot(SAWC_015_030.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    SAWC_015_030.rho <- rhohat(object = dat.ppp, covariate = SAWC_015_030.im)
    windows(); plot(SAWC_015_030.rho, xlab = "SAWC_015_030 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\SAWC_015_030",type = c( "png"),device = dev.cur())
  }
fig.SAWC_015_030.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.SAWC_030_060.rho <- function() {
    SAWC_030_060.im <- as.im(SAWC_030_060.r)
    windows(); plot(SAWC_030_060.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    SAWC_030_060.rho <- rhohat(object = dat.ppp, covariate = SAWC_030_060.im)
    windows(); plot(SAWC_030_060.rho, xlab = "SAWC_030_060 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\SAWC_030_060",type = c( "png"),device = dev.cur())
  }
fig.SAWC_030_060.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.SAWC_060_100.rho <- function() {
    SAWC_060_100.im <- as.im(SAWC_060_100.r)
    windows(); plot(SAWC_060_100.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    SAWC_060_100.rho <- rhohat(object = dat.ppp, covariate = SAWC_060_100.im)
    windows(); plot(SAWC_060_100.rho, xlab = "SAWC_060_100 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\SAWC_060_100",type = c( "png"),device = dev.cur())
  }
fig.SAWC_060_100.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.SAWC_100_200.rho <- function() {
    SAWC_100_200.im <- as.im(SAWC_100_200.r)
    windows(); plot(SAWC_100_200.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    SAWC_100_200.rho <- rhohat(object = dat.ppp, covariate = SAWC_100_200.im)
    windows(); plot(SAWC_100_200.rho, xlab = "SAWC_100_200 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\SAWC_100_200",type = c( "png"),device = dev.cur())
  }
fig.SAWC_100_200.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.habit1.rho <- function() {
    habit1.im <- as.im(habit1.r)
    windows(); plot(habit1.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    habit1.rho <- rhohat(object = dat.ppp, covariate = habit1.im)
    windows(); plot(habit1.rho, xlab = "habit1 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\habit1",type = c( "png"),device = dev.cur())
  }
fig.habit1.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.habit2.rho <- function() {
    habit2.im <- as.im(habit2.r)
    windows(); plot(habit2.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    habit2.rho <- rhohat(object = dat.ppp, covariate = habit2.im)
    windows(); plot(habit2.rho, xlab = "habit2 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\habit2",type = c( "png"),device = dev.cur())
  }
fig.habit2.rho() 
graphics.off()

-----------------------------------------------------------------------------------------------------------------------------
  fig.habit3.rho <- function() {
    habit3.im <- as.im(habit3.r)
    windows(); plot(habit3.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    habit3.rho <- rhohat(object = dat.ppp, covariate = habit3.im)
    windows(); plot(habit3.rho, xlab = "habit3 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\habit3",type = c( "png"),device = dev.cur())
  }
fig.habit3.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.suitable_0.rho <- function() {
    suitable_0.im <- as.im(suitable_0.r)
    windows(); plot(suitable_0.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    suitable_0.rho <- rhohat(object = dat.ppp, covariate = suitable_0.im)
    windows(); plot(suitable_0.rho, xlab = "suitable_0 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\suitable_0",type = c( "png"),device = dev.cur())
  }
fig.suitable_0.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.suitable_1.rho <- function() {
    suitable_1.im <- as.im(suitable_1.r)
    windows(); plot(suitable_1.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    suitable_1.rho <- rhohat(object = dat.ppp, covariate = suitable_1.im)
    windows(); plot(suitable_1.rho, xlab = "suitable_1 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\suitable_1",type = c( "png"),device = dev.cur())
  }
fig.suitable_1.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.suitable_2.rho <- function() {
    suitable_2.im <- as.im(suitable_2.r)
    windows(); plot(suitable_2.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    suitable_2.rho <- rhohat(object = dat.ppp, covariate = suitable_2.im)
    windows(); plot(suitable_2.rho, xlab = "suitable_2 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\suitable_2",type = c( "png"),device = dev.cur())
  }
fig.suitable_2.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.suitable_3.rho <- function() {
    suitable_3.im <- as.im(suitable_3.r)
    windows(); plot(suitable_3.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    suitable_3.rho <- rhohat(object = dat.ppp, covariate = suitable_3.im)
    windows(); plot(suitable_3.rho, xlab = "suitable_3 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\suitable_3",type = c( "png"),device = dev.cur())
  }
fig.suitable_3.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_bridleway.rho <- function() {
    distance_bridleway.im <- as.im(distance_bridleway.r)
    windows(); plot(distance_bridleway.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_bridleway.rho <- rhohat(object = dat.ppp, covariate = distance_bridleway.im)
    windows(); plot(distance_bridleway.rho, xlab = "distance_bridleway (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_bridleway",type = c( "png"),device = dev.cur())
  }
fig.distance_bridleway.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_cycleway.rho <- function() {
    distance_cycleway.im <- as.im(distance_cycleway.r)
    windows(); plot(distance_cycleway.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_cycleway.rho <- rhohat(object = dat.ppp, covariate = distance_cycleway.im)
    windows(); plot(distance_cycleway.rho, xlab = "distance_cycleway (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_cycleway",type = c( "png"),device = dev.cur())
  }
fig.distance_cycleway.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_footway.rho <- function() {
    distance_footway.im <- as.im(distance_footway.r)
    windows(); plot(distance_footway.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_footway.rho <- rhohat(object = dat.ppp, covariate = distance_footway.im)
    windows(); plot(distance_footway.rho, xlab = "distance_footway (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_footway",type = c( "png"),device = dev.cur())
  }
fig.distance_footway.rho() 
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_motorwayandlink.rho <- function() {
    distance_motorwayandlink.im <- as.im(distance_motorwayandlink.r)
    windows(); plot(distance_motorwayandlink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_motorwayandlink.rho <- rhohat(object = dat.ppp, covariate = distance_motorwayandlink.im)
    windows(); plot(distance_motorwayandlink.rho, xlab = "distance_motorwayandlink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_motorwayandlink",type = c( "png"),device = dev.cur())
  }
fig.distance_motorwayandlink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_path.rho <- function() {
    distance_path.im <- as.im(distance_path.r)
    windows(); plot(distance_path.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_path.rho <- rhohat(object = dat.ppp, covariate = distance_path.im)
    windows(); plot(distance_path.rho, xlab = "distance_path (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_path",type = c( "png"),device = dev.cur())
  }
fig.distance_path.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_pedestrian.rho <- function() {
    distance_pedestrian.im <- as.im(distance_pedestrian.r)
    windows(); plot(distance_pedestrian.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_pedestrian.rho <- rhohat(object = dat.ppp, covariate = distance_pedestrian.im)
    windows(); plot(distance_pedestrian.rho, xlab = "distance_pedestrian (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_pedestrian",type = c( "png"),device = dev.cur())
  }
fig.distance_pedestrian.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_primaryandlink.rho <- function() {
    distance_primaryandlink.im <- as.im(distance_primaryandlink.r)
    windows(); plot(distance_primaryandlink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_primaryandlink.rho <- rhohat(object = dat.ppp, covariate = distance_primaryandlink.im)
    windows(); plot(distance_primaryandlink.rho, xlab = "distance_primaryandlink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_primaryandlink",type = c( "png"),device = dev.cur())
  }
fig.distance_primaryandlink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_residentil.rho <- function() {
    distance_residentil.im <- as.im(distance_residentil.r)
    windows(); plot(distance_residentil.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_residentil.rho <- rhohat(object = dat.ppp, covariate = distance_residentil.im)
    windows(); plot(distance_residentil.rho, xlab = "distance_residentil (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_residentil",type = c( "png"),device = dev.cur())
  }
fig.distance_residentil.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_secondaryandlink.rho <- function() {
    distance_secondaryandlink.im <- as.im(distance_secondaryandlink.r)
    windows(); plot(distance_secondaryandlink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_secondaryandlink.rho <- rhohat(object = dat.ppp, covariate = distance_secondaryandlink.im)
    windows(); plot(distance_secondaryandlink.rho, xlab = "distance_secondaryandlink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_secondaryandlink",type = c( "png"),device = dev.cur())
  }
fig.distance_secondaryandlink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_service.rho <- function() {
    distance_service.im <- as.im(distance_service.r)
    windows(); plot(distance_service.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_service.rho <- rhohat(object = dat.ppp, covariate = distance_service.im)
    windows(); plot(distance_service.rho, xlab = "distance_service (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_service",type = c( "png"),device = dev.cur())
  }
fig.distance_service.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_steps.rho <- function() {
    distance_steps.im <- as.im(distance_steps.r)
    windows(); plot(distance_steps.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_steps.rho <- rhohat(object = dat.ppp, covariate = distance_steps.im)
    windows(); plot(distance_steps.rho, xlab = "distance_steps (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_steps",type = c( "png"),device = dev.cur())
  }
fig.distance_steps.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_tertiaryandlink.rho <- function() {
    distance_tertiaryandlink.im <- as.im(distance_tertiaryandlink.r)
    windows(); plot(distance_tertiaryandlink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_tertiaryandlink.rho <- rhohat(object = dat.ppp, covariate = distance_tertiaryandlink.im)
    windows(); plot(distance_tertiaryandlink.rho, xlab = "distance_tertiaryandlink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_tertiaryandlink",type = c( "png"),device = dev.cur())
  }
fig.distance_tertiaryandlink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_track.rho <- function() {
    distance_track.im <- as.im(distance_track.r)
    windows(); plot(distance_track.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_track.rho <- rhohat(object = dat.ppp, covariate = distance_track.im)
    windows(); plot(distance_track.rho, xlab = "distance_track (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_track",type = c( "png"),device = dev.cur())
  }
fig.distance_track.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_trunkandlink.rho <- function() {
    distance_trunkandlink.im <- as.im(distance_trunkandlink.r)
    windows(); plot(distance_trunkandlink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_trunkandlink.rho <- rhohat(object = dat.ppp, covariate = distance_trunkandlink.im)
    windows(); plot(distance_trunkandlink.rho, xlab = "distance_trunkandlink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_trunkandlink",type = c( "png"),device = dev.cur())
  }
fig.distance_trunkandlink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.distance_unclassified.rho <- function() {
    distance_unclassified.im <- as.im(distance_unclassified.r)
    windows(); plot(distance_unclassified.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    distance_unclassified.rho <- rhohat(object = dat.ppp, covariate = distance_unclassified.im)
    windows(); plot(distance_unclassified.rho, xlab = "distance_unclassified (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\distance_unclassified",type = c( "png"),device = dev.cur())
  }
fig.distance_unclassified.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.drains.rho <- function() {
    drains.im <- as.im(drains.r)
    windows(); plot(drains.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    drains.rho <- rhohat(object = dat.ppp, covariate = drains.im)
    windows(); plot(drains.rho, xlab = "drains (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\drains",type = c( "png"),device = dev.cur())
  }
fig.drains.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.footway.rho <- function() {
    footway.im <- as.im(footway.r)
    windows(); plot(footway.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    footway.rho <- rhohat(object = dat.ppp, covariate = footway.im)
    windows(); plot(footway.rho, xlab = "footway (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\footway",type = c( "png"),device = dev.cur())
  }
fig.footway.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.motorway.rho <- function() {
    motorway.im <- as.im(motorway.r)
    windows(); plot(motorway.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    motorway.rho <- rhohat(object = dat.ppp, covariate = motorway.im)
    windows(); plot(motorway.rho, xlab = "motorway (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\motorway",type = c( "png"),device = dev.cur())
  }
fig.motorway.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.path.rho <- function() {
    path.im <- as.im(path.r)
    windows(); plot(path.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    path.rho <- rhohat(object = dat.ppp, covariate = path.im)
    windows(); plot(path.rho, xlab = "path (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\path",type = c( "png"),device = dev.cur())
  }
fig.path.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.pedestrian.rho <- function() {
    pedestrian.im <- as.im(pedestrian.r)
    windows(); plot(pedestrian.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    pedestrian.rho <- rhohat(object = dat.ppp, covariate = pedestrian.im)
    windows(); plot(pedestrian.rho, xlab = "pedestrian (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\pedestrian",type = c( "png"),device = dev.cur())
  }
fig.pedestrian.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.primary.rho <- function() {
    primary.im <- as.im(primary.r)
    windows(); plot(primary.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    primary.rho <- rhohat(object = dat.ppp, covariate = primary.im)
    windows(); plot(primary.rho, xlab = "primary (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\primary",type = c( "png"),device = dev.cur())
  }
fig.primary.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.residential.rho <- function() {
    residential.im <- as.im(residential.r)
    windows(); plot(residential.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    residential.rho <- rhohat(object = dat.ppp, covariate = residential.im)
    windows(); plot(residential.rho, xlab = "residential (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\residential",type = c( "png"),device = dev.cur())
  }
fig.residential.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.roads_motor.rho <- function() {
    roads_motor.im <- as.im(roads_motor.r)
    windows(); plot(roads_motor.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    roads_motor.rho <- rhohat(object = dat.ppp, covariate = roads_motor.im)
    windows(); plot(roads_motor.rho, xlab = "roads_motor (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\roads_motor",type = c( "png"),device = dev.cur())
  }
fig.roads_motor.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.roads_other.rho <- function() {
    roads_other.im <- as.im(roads_other.r)
    windows(); plot(roads_other.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    roads_other.rho <- rhohat(object = dat.ppp, covariate = roads_other.im)
    windows(); plot(roads_other.rho, xlab = "roads_other (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\roads_other",type = c( "png"),device = dev.cur())
  }
fig.roads_other.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.secondry.rho <- function() {
    secondry.im <- as.im(secondry.r)
    windows(); plot(secondry.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    secondry.rho <- rhohat(object = dat.ppp, covariate = secondry.im)
    windows(); plot(secondry.rho, xlab = "secondry (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\secondry",type = c( "png"),device = dev.cur())
  }
fig.secondry.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.secondry_link.rho <- function() {
    secondry_link.im <- as.im(secondry_link.r)
    windows(); plot(secondry_link.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    secondry_link.rho <- rhohat(object = dat.ppp, covariate = secondry_link.im)
    windows(); plot(secondry_link.rho, xlab = "secondry_link (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\secondry_link",type = c( "png"),device = dev.cur())
  }
fig.secondry_link.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.tertiary.rho <- function() {
    tertiary.im <- as.im(tertiary.r)
    windows(); plot(tertiary.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    tertiary.rho <- rhohat(object = dat.ppp, covariate = tertiary.im)
    windows(); plot(tertiary.rho, xlab = "tertiary (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\tertiary",type = c( "png"),device = dev.cur())
  }
fig.tertiary.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.tertiary_liink.rho <- function() {
    tertiary_liink.im <- as.im(tertiary_liink.r)
    windows(); plot(tertiary_liink.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    tertiary_liink.rho <- rhohat(object = dat.ppp, covariate = tertiary_liink.im)
    windows(); plot(tertiary_liink.rho, xlab = "tertiary_liink (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\tertiary_liink",type = c( "png"),device = dev.cur())
  }
fig.tertiary_liink.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.track.rho <- function() {
    track.im <- as.im(track.r)
    windows(); plot(track.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    track.rho <- rhohat(object = dat.ppp, covariate = track.im)
    windows(); plot(track.rho, xlab = "track (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\track",type = c( "png"),device = dev.cur())
  }
fig.track.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.trunck.rho <- function() {
    trunck.im <- as.im(trunck.r)
    windows(); plot(trunck.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    trunck.rho <- rhohat(object = dat.ppp, covariate = trunck.im)
    windows(); plot(trunck.rho, xlab = "trunck (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\trunck",type = c( "png"),device = dev.cur())
  }
fig.trunck.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.unclassified.rho <- function() {
    unclassified.im <- as.im(unclassified.r)
    windows(); plot(unclassified.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    unclassified.rho <- rhohat(object = dat.ppp, covariate = unclassified.im)
    windows(); plot(unclassified.rho, xlab = "unclassified (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\unclassified",type = c( "png"),device = dev.cur())
  }
fig.unclassified.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.riverbank.rho <- function() {
    riverbank.im <- as.im(riverbank.r)
    windows(); plot(riverbank.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    riverbank.rho <- rhohat(object = dat.ppp, covariate = riverbank.im)
    windows(); plot(riverbank.rho, xlab = "riverbank (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\riverbank",type = c( "png"),device = dev.cur())
  }
fig.riverbank.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.rivers.rho <- function() {
    rivers.im <- as.im(rivers.r)
    windows(); plot(rivers.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    rivers.rho <- rhohat(object = dat.ppp, covariate = rivers.im)
    windows(); plot(rivers.rho, xlab = "rivers (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\rivers",type = c( "png"),device = dev.cur())
  }
fig.rivers.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.stream.rho <- function() {
    stream.im <- as.im(stream.r)
    windows(); plot(stream.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    stream.rho <- rhohat(object = dat.ppp, covariate = stream.im)
    windows(); plot(stream.rho, xlab = "stream (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\stream",type = c( "png"),device = dev.cur())
  }
fig.stream.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.Slope_percent_61.rho <- function() {
    Slope_percent_61.im <- as.im(Slope_percent_61.r)
    windows(); plot(Slope_percent_61.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    Slope_percent_61.rho <- rhohat(object = dat.ppp, covariate = Slope_percent_61.im)
    windows(); plot(Slope_percent_61.rho, xlab = "Slope_percent_61 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\Slope_percent_61",type = c( "png"),device = dev.cur())
  }
fig.Slope_percent_61.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.TopoWI.rho <- function() {
    TopoWI.im <- as.im(TopoWI.r)
    windows(); plot(TopoWI.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    TopoWI.rho <- rhohat(object = dat.ppp, covariate = TopoWI.im)
    windows(); plot(TopoWI.rho, xlab = "TopoWI (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\TopoWI",type = c( "png"),device = dev.cur())
  }
fig.TopoWI.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.TopoPI.rho <- function() {
    TopoPI.im <- as.im(TopoPI.r)
    windows(); plot(TopoPI.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    TopoPI.rho <- rhohat(object = dat.ppp, covariate = TopoPI.im)
    windows(); plot(TopoPI.rho, xlab = "TopoPI (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\TopoPI",type = c( "png"),device = dev.cur())
  }
fig.TopoPI.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.aspect91.rho <- function() {
    aspect91.im <- as.im(aspect91.r)
    windows(); plot(aspect91.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    aspect91.rho <- rhohat(object = dat.ppp, covariate = aspect91.im)
    windows(); plot(aspect91.rho, xlab = "aspect91 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\aspect91",type = c( "png"),device = dev.cur())
  }
fig.aspect91.rho()
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.dem.rho <- function() {
    dem.im <- as.im(dem.r)
    windows(); plot(dem.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    dem.rho <- rhohat(object = dat.ppp, covariate = dem.im)
    windows(); plot(dem.rho, xlab = "dem (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\dem",type = c( "png"),device = dev.cur())
  }
fig.dem.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.hpop.rho <- function() {
    hpop.im <- as.im(hpop.r)
    windows(); plot(hpop.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)
    windows(); plot(hpop.rho, xlab = "hpop (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\hpop",type = c( "png"),device = dev.cur())
  }
fig.hpop.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.lot_density.rho <- function() {
    lot_density.im <- as.im(lot_density.r)
    windows(); plot(lot_density.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    lot_density.rho <- rhohat(object = dat.ppp, covariate = lot_density.im)
    windows(); plot(lot_density.rho, xlab = "lot_density (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\lot_density",type = c( "png"),device = dev.cur())
  }
fig.lot_density.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.temp_minimum1996_2005.rho <- function() {
    temp_minimum1996_2005.im <- as.im(temp_minimum1996_2005.r)
    windows(); plot(temp_minimum1996_2005.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    temp_minimum1996_2005.rho <- rhohat(object = dat.ppp, covariate = temp_minimum1996_2005.im)
    windows(); plot(temp_minimum1996_2005.rho, xlab = "temp_minimum1996_2005 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\temp_minimum1996_2005",type = c( "png"),device = dev.cur())
  }
fig.temp_minimum1996_2005.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.temp_maximum1990_2000.rho <- function() {
    temp_maximum1990_2000.im <- as.im(temp_maximum1990_2000.r)
    windows(); plot(temp_maximum1990_2000.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    temp_maximum1990_2000.rho <- rhohat(object = dat.ppp, covariate = temp_maximum1990_2000.im)
    windows(); plot(temp_maximum1990_2000.rho, xlab = "temp_maximum1990_2000 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\temp_maximum1990_2000",type = c( "png"),device = dev.cur())
  }
fig.temp_maximum1990_2000.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.tempmeanMGA65.rho <- function() {
    tempmeanMGA65.im <- as.im(tempmeanMGA65.r)
    windows(); plot(tempmeanMGA65.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    tempmeanMGA65.rho <- rhohat(object = dat.ppp, covariate = tempmeanMGA65.im)
    windows(); plot(tempmeanMGA65.rho, xlab = "tempmeanMGA65 (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\tempmeanMGA65",type = c( "png"),device = dev.cur())
  }
fig.tempmeanMGA65.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.rainfallmeanannual.rho <- function() {
    rainfallmeanannual.im <- as.im(rainfallmeanannual.r)
    windows(); plot(rainfallmeanannual.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    rainfallmeanannual.rho <- rhohat(object = dat.ppp, covariate = rainfallmeanannual.im)
    windows(); plot(rainfallmeanannual.rho, xlab = "rainfallmeanannual (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\rainfallmeanannual",type = c( "png"),device = dev.cur())
  }
fig.rainfallmeanannual.rho
graphics.off()
-----------------------------------------------------------------------------------------------------------------------------
  fig.foliagePC.rho <- function() {
    foliagePC.im <- as.im(foliagePC.r)
    windows(); plot(foliagePC.im, axes = TRUE)
    # What is the nature of the association between koala sight locations and total phosphorous?
    foliagePC.rho <- rhohat(object = dat.ppp, covariate = foliagePC.im)
    windows(); plot(foliagePC.rho, xlab = "foliagePC (units)", main = "")
    savePlot(filename = "raster_crop\\rhohat\\foliagePC",type = c( "png"),device = dev.cur())
  }
fig.foliagePC.rho
graphics.off()
-----------------
---------------------------------------------------------------------------------------------------------------------------  
  save.image("AU_Qld_koala_habitat_v08edited_distance_Ravi15_03.RData")


# =================================================================================================================================
# Poisson point process model:

# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(aspect91= aspect91.im, CLY_005_015= CLY_005_015.im, dem=dem.im, distance_bridleway=distance_bridleway.im, distance_cycleway=distance_cycleway.im,
                                                         distance_footway=distance_footway.im, distance_motorwayandlink = distance_motorwayandlink.im, distance_pedestrian = distance_pedestrian.im,  distance_primaryandlink = distance_primaryandlink.im,
                                                         distance_residentil = distance_residentil.im, distance_secondaryandlink = distance_secondaryandlink.im, distance_steps = distance_steps.im, 
                                                         distance_tertiaryandlink = distance_tertiaryandlink.im, distance_track = distance_track.im,distance_trunkandlink = distance_trunkandlink.im, foliagePC = foliagePC.im, 
                                                         habit3 = habit3.im,habit1=habit1.im,habit2=habit2.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,PTO_000_005=PTO_000_005.im,
                                                         rainfallmeanannual = rainfallmeanannual.im, residential = residential.im, roads_motor=roads_motor.im, SAWC_030_060=SAWC_030_060.im, Slope_percent_61 = Slope_percent_61.im, suitable_0=suitable_0.im,
                                                         suitable_1=suitable_1.im, suitable_3=suitable_3.im,TopoWI=TopoWI.im, unclassified=unclassified.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409
windows();plot(dat.ppp)

# Backward stepwise variable selection:
dat.ppm01 <- ppm(Qz, trend =  ~ aspect91+CLY_005_015+dem+distance_bridleway+distance_cycleway+
                   distance_footway+distance_motorwayandlink +distance_pedestrian +distance_primaryandlink+
                   distance_residentil+distance_secondaryandlink+distance_steps +
                   distance_tertiaryandlink+distance_track+distance_trunkandlink+foliagePC+ 
                   habit3+habit1+habit2+hpop+lot_density+NTO_000_005+PTO_000_005+
                   rainfallmeanannual+residential+roads_motor+SAWC_030_060+Slope_percent_61 +suitable_0+
                   suitable_1+suitable_3+TopoWI+ unclassified,
                 covariates = list(aspect91= aspect91.im, CLY_005_015= CLY_005_015.im, dem=dem.im, distance_bridleway=distance_bridleway.im, distance_cycleway=distance_cycleway.im,
                      distance_footway=distance_footway.im, distance_motorwayandlink = distance_motorwayandlink.im, distance_pedestrian = distance_pedestrian.im,  distance_primaryandlink = distance_primaryandlink.im,
                      distance_residentil = distance_residentil.im, distance_secondaryandlink = distance_secondaryandlink.im, distance_steps = distance_steps.im, 
                      distance_tertiaryandlink = distance_tertiaryandlink.im, distance_track = distance_track.im,distance_trunkandlink = distance_trunkandlink.im, foliagePC = foliagePC.im, 
                      habit3 = habit3.im,habit1=habit1.im,habit2=habit2.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,PTO_000_005=PTO_000_005.im,
                      rainfallmeanannual = rainfallmeanannual.im, residential = residential.im, roads_motor=roads_motor.im, SAWC_030_060=SAWC_030_060.im, Slope_percent_61 = Slope_percent_61.im, suitable_0=suitable_0.im,
                      suitable_1=suitable_1.im, suitable_3=suitable_3.im,TopoWI=TopoWI.im, unclassified=unclassified.im))
-------------------------------------------------------------------------------------------------------------------------------------
  
# Values of the covariates 'roaden', 'temp' were NA or undefined at 0.5% (5 out of 1005) of the quadrature points. 
# See http://gis.stackexchange.com/questions/161893/creating-polygonal-windows-for-spatial-analysis-in-r for a work around.

# Extract the quadrature scheme from dat.ppm01:
Qz <- quad.ppm(dat.ppm01, drop = TRUE)
# This step will drop points within the quadrature scheme that had NA-values.

# Run the model again, but this time, use the corrected quadrature scheme:
dat.ppm02 <- ppm(Qz, trend =  ~ aspect91+CLY_005_015+dem+distance_bridleway+distance_cycleway+
                   distance_footway+distance_motorwayandlink +distance_pedestrian +distance_primaryandlink+
                   distance_residentil+distance_secondaryandlink+distance_steps +
                   distance_tertiaryandlink+distance_track+distance_trunkandlink+foliagePC+ 
                   habit3+habit1+habit2+hpop+lot_density+NTO_000_005+PTO_000_005+
                   rainfallmeanannual+residential+roads_motor+SAWC_030_060+Slope_percent_61 +suitable_0+
                   suitable_1+suitable_3+TopoWI+ unclassified,
                 covariates = list(aspect91= aspect91.im, CLY_005_015= CLY_005_015.im, dem=dem.im, distance_bridleway=distance_bridleway.im, distance_cycleway=distance_cycleway.im,
                                   distance_footway=distance_footway.im, distance_motorwayandlink = distance_motorwayandlink.im, distance_pedestrian = distance_pedestrian.im,  distance_primaryandlink = distance_primaryandlink.im,
                                   distance_residentil = distance_residentil.im, distance_secondaryandlink = distance_secondaryandlink.im, distance_steps = distance_steps.im, 
                                   distance_tertiaryandlink = distance_tertiaryandlink.im, distance_track = distance_track.im,distance_trunkandlink = distance_trunkandlink.im, foliagePC = foliagePC.im, 
                                   habit3 = habit3.im,habit1=habit1.im,habit2=habit2.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,PTO_000_005=PTO_000_005.im,
                                   rainfallmeanannual = rainfallmeanannual.im, residential = residential.im, roads_motor=roads_motor.im, SAWC_030_060=SAWC_030_060.im, Slope_percent_61 = Slope_percent_61.im, suitable_0=suitable_0.im,
                                   suitable_1=suitable_1.im, suitable_3=suitable_3.im,TopoWI=TopoWI.im, unclassified=unclassified.im))


summary(dat.ppm02)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm02)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm02) # 3413.045

step
# Drop tpo:roaden 
dat.ppm03 <- ppm(Qz, trend = ~ fpc +topo + slope + aspect + twi + tpo+hpop + temp + elev + habit01 + habit02 + habit03,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm03)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm03)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm03) # 3411.272


# Drop habitat03:tpo
dat.ppm04 <- ppm(Qz, trend = ~ fpc +topo + slope + aspect + twi + tpo+hpop + temp + elev + habit01 + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm04)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm04)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm04) # 3460.872


# Drop topo:
dat.ppm05 <- ppm(Qz, trend = ~ fpc + slope + aspect + twi + tpo+hpop + temp + elev + habit01 + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm05)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm05)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm05) # 3459.164


# Drop habit1:
dat.ppm06 <- ppm(Qz, trend = ~ fpc + slope + aspect + twi + tpo+hpop + temp + elev + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
tdat <- data.frame(summary(dat.ppm06)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm06) # 3457.573


# Drop topo:
dat.ppm07 <- ppm(Qz, trend = ~ fpc + slope + aspect + twi +hpop + temp + elev + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm07)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm07)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm07) # 3455.898


# Drop hpop:
dat.ppm08 <- ppm(Qz, trend = ~ fpc + slope + aspect + twi + temp + elev + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm08)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm08)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm08) # 3455.163


# Drop water
dat.ppm09 <- ppm(Qz, trend = ~ fpc + slope + twi + temp + elev + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm09)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm09)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm09) # 3453.953


# Drop twi
dat.ppm10 <- ppm(Qz, trend = ~ fpc + slope + twi + temp + elev + habit02,
                 covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                  water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                  habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm10)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm10)$coefs.SE.CI); tdat[order(tdat$Zval),]

AIC(dat.ppm10) # 3411.58
---------------------------------------------------------
  dat.ppm11 <- ppm(Qz, trend = ~ fpc + slope +temp + elev + habit02,
                   covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                    water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                    habit02 = habit02.im, habit03 = habit03.im))
tdat <- data.frame(summary(dat.ppm11)$coefs.SE.CI); tdat[order(tdat$Zval),]
summary(dat.ppm11)$coefs.SE.CI
AIC(dat.ppm11)
----------------------------------------------------------------------
  dat.ppm12 <- ppm(Qz, trend = ~ fpc +temp + elev + habit02,
                   covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                    water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                    habit02 = habit02.im, habit03 = habit03.im))
tdat <- data.frame(summary(dat.ppm12)$coefs.SE.CI); tdat[order(tdat$Zval),]
summary(dat.ppm12)$coefs.SE.CI
AIC(dat.ppm12)
---------------------------------------------------------------------------
  dat.ppm13 <- ppm(Qz, trend = ~ fpc +temp +habit02,
                   covariates =list(fpc= fpc.im, awc= awc.im, topo = topo.im, slope = slope.im, aspect = aspect.im, twi=twi.im, tpo = tpo.im, sbd = sbd.im,  
                                    water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im,
                                    habit02 = habit02.im, habit03 = habit03.im))
tdat <- data.frame(summary(dat.ppm13)$coefs.SE.CI); tdat[order(tdat$Zval),]
summary(dat.ppm13)$coefs.SE.CI
AIC(dat.ppm13)

windows(); plot(dat.ppm13, how = "image", se = TRUE)

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
dat.pred14 <- predict(dat.ppm02); pred <- dat.pred14

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
----------------------------------------------------------------------------
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
#Kernel density estmation
#Work out density for dat.ppp (using Diggle's edge correction):
dat.den <- density(dat.ppp, sigma = 1500, diggle = TRUE, dimyx = dimyx)

# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre
#:gaussian kernel function. Its bandwidth defines the kernel’s window extent.
dat.den <- density(dat.ppp,sigma=1500, diggle = TRUE,dimyx = dimyx)
windows();plot(dat.den, main=NULL)
contour(dat.den, add=TRUE)
head(dat.den)
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))
contour(dat.den, add=TRUE)
windows(); plot(dat.den)
# Plot the relationship
windows(); plot(effectfun(PPM1, "elev.im", se.fit=TRUE), main=NULL, cex.axis=0.6,cex.lab=0.6,legend=FALSE)

-------------------------------------------------------------------
  #Kernel Density Adjusted for Covariate
  # a kernel density map is generated using an elevation raster as a covariate. 
  #The outputs include a plot of  ρρ  vs. elevation and a raster map of  ρρ .
  
  library(spatstat)
# Compute rho using the ratio method
rho <- rhohat(dat.ppp, elev.im, method="ratio")
# Generate rho vs covariate plot
windows();plot(rho, main=NULL)
pred <- predict(rho)
pred <-pred / 1000
summary(as.vector(pred$v))
names(pred)
cl <-  interp.colours(c("yellow", "pink" ,"green"), 100) # Create color scheme
windows(); plot(pred, col=cl, las=2, main=NULL)

---------------------------------------------------------------------------------------------------
#Modeling intensity as a function of a covariate
  # Create the Poisson point process model
PPM1 <- ppm(dat.ppp ~ elev.im)
predict(PPM1, type = "cif", ngrid = 256)
coef(PPM1)
vcov(PPM1)
sqrt(diag(vcov(PPM1)))
round(vcov(PPM1, what = "corr"), 2)
windows();plot(PPM1, se=TRUE)#standard error of fitted intensity
M <- quadrat.test(PPM1, nx = 4, ny = 2)#Chi-squared test of fitted Poisson model fit using quadrat counts
windows(); plot(M, add = TRUE, cex = 1.5, col = "red")
#inspect the residual counts:observed number of points (top left), the predicted
#number of points according to the model (top right), and the Pearson residual (bottom)
windows();plot(dat.ppp, pch = ".")
plot(M, add = TRUE, cex = 1, col = "red")
#Residual measure
plot(predict(PPM1))
plot(dat.ppp, add = TRUE, pch = "+")
#smooth residuals
windows();diagnose.ppm(PPM1, which = "smooth")
###Lurking variable plot|If there is a spatial covariate Z(u) that plays an important role in the analysis,lurking variable plot of the residuals against Z.
lurking(PPM1, elev.im, type = "raw", main="plotshows approximate 5% significance bands for the cumulative residual")
lurking(PPM1, elev.im, type = "raw", cumulative = FALSE, main="display the derivative C′(z), which often indicates which values of z are
associated with a lack of fit.")
#four panel plot
diagnose.ppm(fit)
--------------------------------------------------------------------------------
#working on how to map predictions from model intensity AS a coariate.
  pre <- predict(PPM1)
head(pre)
pre$v <- pre$v * 1E04
summary(as.vector(pre$v))
windows(); hist(as.vector(pre$v))

# Predictions from the model:
breaks <- seq(from = 0, to = 3.3810 , length = 6)
col <- brewer.pal(n = 5, name = "Reds")

windows(); plot(x = xlim, y = ylim, main="Predictions from the model", type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pre$xcol, y = pre$yrow, z = t(pre$v), zlim = c(0, 0.015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "blue")

plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Predictions from the model14",type = c( "png"),device = dev.cur())
----------------------------------------------------------------------------------------------------