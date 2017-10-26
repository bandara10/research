setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
library(GISTools)
mydataSurveysighting <- read.csv("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\survey_sightings_latlongnew3.csv", header = TRUE)
mydata.s <- data.frame(mydataSurveysighting)
# Subset the data:
mydata.s <- mydata.s[c("LAT","LNG","X", "Y", "Year")]
names(mydata.s) <- tolower(names(mydata.s))
# Analyse data for 2011:
id <- mydata.s$year >2009
table(id)
mydata.s <- mydata.s[id,]
dim(mydata.s)
mydata.s <- mydata.s[c("lat","lng")]


# 354 rows

id <- mydata.s$lat > -30
mydata.s <- mydata.s[id,]
# STEP 2: Define a 3 km by 3 km grid and count up the number of koalas in each grid cell:

# Bring in MGA56 study area boundary map:
unzip("vector\\AU_Qld_study_area.zip")

studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 
windows(); plot(studyarea.shp)
# mydata in lat-lon:
coords <- SpatialPoints(mydata.s[, c("lng","lat")])
mydata.ll <- SpatialPoints(coords)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 
# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))

plot(mydata.mga)

# Add the MGA coordinates to mydata:
mydata.s$x <- coordinates(mydata.mga)[,1]
mydata.s$y <- coordinates(mydata.mga)[,2]
windows(); plot(studyarea.shp, axes = TRUE)


#Create a 3 km by 3 km grid:
mydata.r <- raster(studyarea.shp)
res(mydata.r) <- 1000
mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset


# Sample points from within each grid cell and create a SpatialPolygonsDataframe:
mydata.spolydf <- rasterToPolygons(mydata.r)
mydata.spolydf$layer <- 1:length(mydata.spolydf$layer)

names(mydata.spolydf)[1] <- "id"

windows(); plot(mydata.spolydf)
points(mydata.mga)
write.csv(mydata.mga, file="mydata.mga.csv", row.names=FALSE) 
# Count up the number of koala sightings in each polygon of the SpatialPolygonsDataframe:
rval <- poly.counts(pts = mydata.mga, polys = mydata.spolydf) # GIS tools
mydata.spolydf$n <- rval

writePolyShape(mydata.spolydf, "IE_study_area_boundary-ING")

#convert to a raster 


rr <- readOGR(getwd(), "IE_study_area_boundary-ING")
rrr <- raster(extent(rr))
res(rrr)=1000
koalan <- rasterize(rr, field="n", rrr)
proj4string(koalan) <- CRS("+init=epsg:28356") 
writeRaster(koalan,"koalan.TIF", overwrite=TRUE)



# --------------------------------------------------------------------------------------------------------------------
d_bridleway <- raster("raster_crop\\distance_bridleway.tif")
d_cycleway <- raster("raster_crop\\distance_cycleway.tif")
d_footway <- raster("raster_crop\\distance_footway.tif")
d_motorwayandlink <- raster("raster_crop\\distance_motorwayandlink.tif")
d_path <- raster("raster_crop\\distance_path.tif")
d_pedestrian <- raster("raster_crop\\distance_pedestrian.tif")
d_primaryandlink <- raster("raster_crop\\distance_primaryandlink.tif")
d_residentil <- raster("raster_crop\\distance_residentil.tif")
d_secondaryandlink <- raster("raster_crop\\distance_secondaryandlink.tif")
d_steps <- raster("raster_crop\\distance_steps.tif")
d_tertiaryandlink <- raster("raster_crop\\distance_tertiaryandlink.tif")
d_trunkandlink <- raster("raster_crop\\distance_trunkandlink.tif")
d_unclassified <- raster("raster_crop\\distance_unclassified.tif")
#---------- 
windows(); plot(d_bridleway)
myPredictors = c(
  "distance_bridleway","distance_cycleway", "distance_motorwayandlink", "distance_path", "distance_pedestrian", "distance_primaryandlink", "distance_residentil",
  "distance_secondaryandlink", "distance_steps", "distance_tertiaryandlink", "distance_trunkandlink", "distance_unclassified")
# add extensions
myPredictorsASC = paste(myPredictors,".TIF", sep = "")
myPredList = as.list(paste(myPredictors,".TIF", sep = ""))  
myStack = stack(myPredList)

surveydata <- extract(myStack, mydata.mga)

# STEP 4: Combine the outcome variable vector and the matrix of explanatory variables into a single data frame:
#b<- as.data.frame(stack.r,row.ames=null)
#write.csv(b, file="b.csv", row.names=FALSE) 
centroids <- getSpPPolygonsLabptSlots(rr)
colnames(centroids) <- c("x","y")
head(centroids)
extract.e <- extract(myStack, centroids)
extract.e = cbind(centroids,extract.e )
head(extract.e)
write.csv(extract.e, file="extract_e_survey.csv", row.names=FALSE) 
extract.e_survey <- read.csv("extract_e.csv", header = T)
# =================================================================================================================================
# STEP 5: Write the data out as a shape file:
coordinates(extract.e) <- ~x + y
proj4string(extract.e) <- CRS("+init=epsg:28356") 
class(extract.e)
shapefile(extract.e, "extract_e.shp")
plot(extract.e)
names(extract.e)

