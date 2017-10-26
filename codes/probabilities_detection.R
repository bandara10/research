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
#mydata.p <- rasterToPolygons(mydata.r)
acsel<-as.data.frame(acsel,row.names = NULL)
colnames(acsel) <- c("x","y")# rename colomns to x, y to match negative data colomns
head(acsel)
----
windows(); plot(mydata.p, border = 'gray')
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(mydata.mga)
# Selected points in red:
points(acsel, cex = 0.2, col = 'red', pch = 15)# go to 103.

dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

windows(); plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(x = mydata$x, y = mydata$y,cex = 0.2, col = 'blue', pch = 15)
-----------
# kdate as a date:
# dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
# dat$kmonth <- format(dat$kdate, format = "%m")
# 2_Creates a random set of negatives
#############################################################edit this
windows();plot(mydata.p, border = 'gray')
# Sets the number of points to generate
n = 1000
#Generates random coordinates
writeRaster(mydata.r,"mydata.r.TIF", overwrite=TRUE)
myMask = raster("mydata.r.TIF") >= 0
x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
# Keep only points according to different conditions
myTDF = as.data.frame(cbind(x,y)) # negative data point
head(myTDF)
acsel2 <- gridSample(myTDF, myMask, n = 1)
points(x = acsel2$x, y = acsel2$y,cex = 0.2, col = 'blue', pch = 15)
#plot(acsel2)
# Merge the two datasets (positives and negatives
#############################################################
names(acsel) = c("x","y")
acsel$Haskoala = 1
acsel2$Haskoala = 0
acsel21 = rbind(acsel,acsel2)
# Remove potential duplicates falling in same pixel
#############################################################
source("Lib_DistEstimatesToItself.r")
acsel21$DistToItself = Lib_DistEstimatesToItself(acsel21$x,acsel21$y)
acsel21 = subset(acsel21, DistToItself > 100)
acsel21 = acsel21[,-4]
# displays the output
plot(mydata.p, border = 'gray')
points(acsel21$x,acsel21$y, pch = 15, col=c("blue","red")[acsel21$Haskoala+1])

# ---------------------------------------------------------------------------------------------------------------------------------

 "raster_crop\\roaden.TIF"
 
 "raster_crop\\hpop.TIF"
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
#-----
  mystack<-stack("raster_crop\\distance_bridleway.tif", "raster_crop\\distance_cycleway.tif", "raster_crop\\distance_footway.tif",
               "raster_crop\\distance_motorwayandlink.tif", "raster_crop\\distance_path.tif","raster_crop\\distance_pedestrian.tif",
               "raster_crop\\distance_primaryandlink.tif","raster_crop\\distance_residentil.tif","raster_crop\\distance_secondaryandlink.tif",
               "raster_crop\\distance_steps.tif", "raster_crop\\distance_tertiaryandlink.tif","raster_crop\\distance_trunkandlink.tif",
               "raster_crop\\distance_unclassified.tif")

acsel22 = cbind(acsel21,extract(mystack,acsel21[,1:2]))
str(acsel22)
acsel22 = na.omit(acsel22)

--------
  # Build the GLM object
  myStr = "Haskoala ~"
for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
print(myStr)
print("##########################################")
# Runs the GLM  
myGLM = glm(as.formula(myStr), data = acsel22, family = "binomial")
-------






myGLM = glm(Haskoala ~ distance_bridleway + distance_cycleway +distance_footway+distance_motorwayandlink+
            distance_path+ distance_pedestrian+distance_primaryandlink+
            distance_residentil+ distance_secondaryandlink+distance_steps+
            distance_tertiaryandlink+ distance_tertiaryandlink+ distance_trunkandlink+
            distance_unclassified, data = acsel22, family = "binomial")
summary(myGLM)
# Stores the residuals
acsel22$res = residuals(myGLM)
# Create a spatial correlogram of the residuals
myResCorr <- correlog(acsel22$x, acsel22$y, acsel22$res,na.rm=T, increment=1000, resamp=0, latlon = F)
#mySigVec = ifelse(myResCorr$p<0.05,1,0)  
#plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
#xlab="distance", ylab="Moran's I")
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")
abline(h=0)

# Crase approach to account for SA
#############################################################
myStack = stack(myPredList)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = raster("distance_bridleway.TIF") >= 0
plot(myMask)


# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

# Extract all predictors from raster stack
############################################################# 
myFD = cbind(acsel21,extract(myStack,acsel21[,1:2]))
str(myFD)
# removes lines with no data

myFD = na.omit(myFD)
# Build the GLM object
myStr = "Haskoala ~"
for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
print(myStr)
print("##########################################")
# Runs the GLM  
myGLM = glm(as.formula(myStr), data = myFD, family = "binomial")
#myGLM = glm(HasHPAI ~ brdn+ dsdn+dudn+hpdn+rddn+access+urban_r+dem+ncrop, data = myFD, family = "binomial")
# Show summary data
summary(myGLM)
# Stores the residuals
myFD$res = residuals(myGLM)
# Create a spatial correlogram of the residuals
myResCorr <- correlog(myFD$x, myFD$y, myFD$res,na.rm=T, increment=1000, resamp=0, latlon = F)
#mySigVec = ifelse(myResCorr$p<0.05,1,0)  
#plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
#xlab="distance", ylab="Moran's I")
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")
abline(h=0)

# Crase approach to account for SA
#############################################################
# Build the standard GLM object
myFD$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), myFD$x, myFD$y, 15000)
# Changes the formular string
myStrAR = paste(myStr,"+ AR1")
# Re-run the model  
myGLMR1 = glm(as.formula(myStrAR), data = myFD, family = "binomial")
summary(myGLMR1)
plot(myGLMR1)  



# Estimate correlogram of new residuals

myFD$R1 = residuals(myGLMR1)
Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=1000, resamp=0, latlon = F)              
#Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram
plot(Corr)              

# Estimates DAIC for each predictor
#############################################################
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "Haskoala ~"
  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
  myStrARN = paste(myStrN,"+ AR1")
  myGLMR2 = glm(as.formula(myStrARN), data = myFD, family = "binomial")
  mylrtest = lrtest(myGLMR1,myGLMR2)
  myChi2[i] = mylrtest[[4]]
  myPValue[i] = mylrtest[[6]]
}
myLRTable = data.frame(myChi2,myPValue)
myLRTable = cbind(myPredictors,myLRTable)
myLRTable

# ROC & predictions
#############################################################
# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

myPred = predict(myARStack, myGLMR1, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main="koala")
save.image("koalasightingprobability.RData")
writeRaster(myPred,"koalasightingprobability.TIF")
koalasightingprob.r<- raster("koalasightingprobability.TIF")

tpo.crop <- crop(x = koalasightingprob.r, y = r, snap = "near")
windows(); plot(tpo.crop)
writeRaster(tpo.crop,"koalasightingprobability2.TIF")
save.image("koalasightingprobability.RData")
library(swirl)
swirl()
ravi


