library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)


setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\MARK_S")




near_roadtest.shp <- readShapeLines("near_roadstest.shp")
AB<-spsample(near_roadtest.shp, 1100, "regular")#lines to points


d = gDistance(p, habitatneartest.shp, byid=TRUE)
plot(d)

habitatneartest.shp <- readShapePoly("habitatneartest.shp")



AC<-spsample(habitatneartest.shp, 1100, "regular")# line to points 
plot(AC)
habitnear.r <- rasterize(x = habitatneartest.shp, y = mydata.r, field = "NEAR_DIST")
plot(near.r)
writeRaster(near.r,"near.r.TIF", overwrite=TRUE)
writeRaster(habitnear.r,"habitnear.r.r.TIF", overwrite=TRUE)
near.r <- rasterize(x = near_roadtest.shp, y = mydata.r, field = "NEAR_DIST")
dist2Line(habitatneartest.shp, near_roadtest.shp, distfun=distHaversine)
library(geosphere)
library (rgeos)
#calculate distance from roads to habitat

#GRID distance
roadsdist.r <- rasterize(x = near_roadtest.shp, y = mydata.r, field = "NEAR_DIST")
plot(roadsdist.r)
habitnear.r <- rasterize(x = habitatneartest.shp, y = mydata.r, field = "NEAR_DIST")
plot(habitnear.r)
#habitnear.r[]<-0
writeRaster(habitnear.r,"habitnear.r.TIF", overwrite=TRUE)
habitnear.r[!is.na(habitnear.r[])] <- 1# set not NA to 1
habitnear.r[is.na(habitnear.r[])] <- 0 # set NA to 0
plot(habitnear.r)
roadsdist.r <- rasterize(x = near_roadtest.shp, y = mydata.r, field = "NEAR_DIST")
roadsdist.r[!is.na(roadsdist.r[])] <- 2
roadsdist.r[is.na(roadsdist.r[])] <- 0
plot(roadsdist.r)
DD<-sum(roadsdist.r,habitnear.r)
plot(DD)
dx <- gridDistance(DD,origin=1) 
plot(dx)
writeRaster(DD,"distancefrom.r.TIF", overwrite=TRUE)
e <- gridDistance(DD,origin=2) 
writeRaster(e,"distancefromroads.r.TIF", overwrite=TRUE)
save.image(file="distance from roads to habitat")
######distance from uncassified roads to habtat suitable 3.
unclassified.r <- raster("roads_raster\\unclassified.tif")  #read unclassified roads raaster 
suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif") #read habitat suitable_3 high
plot(unclassified.r) # plot and see
plot(suitable_3.r)  # plot and see
unclassified.r [unclassified.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
plot(unclassified.r)
suitable_3.r [suitable_3.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
plot(suitable_3.r)
suitable3_unclassified<-sum(unclassified.r,suitable_3.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable3_unclassified) # 3`s are ovelappings. but can be reclassed to 2 or 1 again.
#writeRaster(suitable3_unclassified,"suitable3_unclassified.r.TIF", overwrite=TRUE)
suitable3_unclassified_dist<- gridDistance(suitable3_unclassified,origin=2, omit=NA) #get distance from 2 to 1.
plot(suitable3_unclassified_dist)
#writeRaster(suitable3_unclassified_dist,"suitable3_unclassified_dist.TIF", overwrite=TRUE)
#save.image(file="distance from roads to habitat")
#distance from unclassified roads to habitat2
unclassified.r <- raster("roads_raster\\unclassified.tif")  #read unclassified roads raaster 
suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif") #read habitat suitable_3 high
plot(suitable_2.r)
unclassified.r [unclassified.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_2.r [suitable_2.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable2_unclassified<-sum(unclassified.r,suitable_2.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable2_unclassified)
suitable2_unclassified_dist<- gridDistance(suitable2_unclassified,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable2_unclassified_dist)
##Distance from unclassified raods to habitat 1
#distance from unclassified roads to habitat2
unclassified.r <- raster("roads_raster\\unclassified.tif")  #read unclassified roads raaster 
suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif") #read habitat suitable_3 high
plot(suitable_1.r)
unclassified.r [unclassified.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_1.r [suitable_1.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable1_unclassified<-sum(unclassified.r,suitable_1.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable1_unclassified)
suitable1_unclassified_dist<- gridDistance(suitable1_unclassified,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable1_unclassified_dist)
####
##Distance from residential raods to habitat 1

residential.r <- raster("roads_raster\\residential.tif")  #read residential roads raaster 
plot(residential.r)
suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif") #read habitat suitable_3 high
plot(suitable_3.r)
residential.r [residential.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_3.r [suitable_3.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable3_residential<-sum(residential.r,suitable_3.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable3_residential)
suitable3_residential_dist<- gridDistance(suitable3_residential,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable3_residential_dist)
# distance from residential roads to habitat 2
residential.r <- raster("roads_raster\\residential.tif")  #read residential roads raaster 
suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif") #read habitat suitable_3 high
plot(suitable_2.r)
residential.r [residential.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_2.r [suitable_2.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable2_residential<-sum(residential.r,suitable_2.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable2_residential)
suitable2_residential_dist<- gridDistance(suitable2_residential,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable2_residential_dist)
#####
# distance from residential roads to habitat 3
residential.r <- raster("roads_raster\\residential.tif")  #read residential roads raaster 
suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif") #read habitat suitable_3 high
plot(suitable_1.r)
residential.r [residential.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_1.r [suitable_1.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable1_residential<-sum(residential.r,suitable_1.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable1_residential)
suitable1_residential_dist<- gridDistance(suitable1_residential,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable1_residential_dist)
########
# distance from secondry roads to habitat 1
secondry.r <- raster("roads_raster\\secondry.tif")  #read residential roads raaster 
suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif") #read habitat suitable_3 high
plot(suitable_1.r)
secondry.r [secondry.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_1.r [suitable_1.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable1_secondry<-sum(secondry.r,suitable_1.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable1_secondry)
suitable1_secondry_dist<- gridDistance(suitable1_secondry,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable1_secondry_dist)
# distance from secondry roads to habitat 2
secondry.r <- raster("roads_raster\\secondry.tif")  #read residential roads raaster
plot(secondry.r)
suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif") #read habitat suitable_3 high
plot(suitable_2.r)
secondry.r [secondry.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_2.r [suitable_2.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable2_secondry<-sum(secondry.r,suitable_2.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable2_secondry)
suitable2_secondry_dist<- gridDistance(suitable2_secondry,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable2_secondry_dist)

# distance from secondry roads to habitat 3
secondry.r <- raster("roads_raster\\secondry.tif")  #read residential roads raaster
plot(secondry.r)
suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif") #read habitat suitable_3 high
plot(suitable_3.r)
secondry.r [secondry.r >0] <- 1# set values greater than 0 to 1. Leave NA as it is.
suitable_3.r [suitable_3.r >0] <- 2# set values greater than 0 to 2 in this raster. leave NA as it is.
suitable3_secondry<-sum(secondry.r,suitable_3.r ) # add two raster to one. so 1 and 2 grid valuves
plot(suitable3_secondry)
suitable3_secondry_dist<- gridDistance(suitable3_secondry,origin=2, omit=NA) #get distance from 2 to1.
plot(suitable3_secondry_dist)



mystack<-stack(suitable3_unclassified_dist,suitable2_unclassified_dist,suitable1_unclassified_dist,
               suitable3_residential_dist,suitable2_residential_dist, suitable1_residential_dist,
               suitable3_secondry_dist, suitable2_secondry_dist, suitable1_secondry_dist)
save.image("distancetohabitat.RData")
plot(mystack)
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\distancetohabitat.RData")
