#1_Load Data
# Sets the general path
#############################################################
setwd("E:/My Documents/Dropbox/Data_Nepal")


# Loads data
#############################################################
#myD = read.dbf("case_series.dbf", as.is=TRUE)
#myD = read.dbf("case_seriesHPAI_Oct.dbf", as.is=TRUE)
#myD = read.dbf("alllab.dbf", as.is=TRUE)
#myD = read.dbf("Hninelabdata.dbf", as.is=TRUE)
#myD = read.dbf("Hfive.dbf", as.is=TRUE)
myD = read.dbf("H5N1_novem.dbf", as.is=TRUE)

#myD = subset(myD, TYPE == "OUT")

# Reproject the data
#############################################################
#OIE H5N1 dataset
#myND = myD[,c("Longitude","Latitude","Year","Start_Date")]
#Case_seriesHPAI_Oct
#myND = myD[,c("Latitude","Longitude","Year","Start_Date")]
#All laboratory data
#myND = myD[,c("Longitude","Latitude","SN","Date_Rec_1")]
#H9lab data
#myND = myD[,c("Longitude","Latitude","Year","Date_Rec_1")]
#H5labdata
#myND = myD[,c("Longitude","Latitude","Year","Date_Rec_1")]
#H5N1_novem
myND = myD[,c("Longitude","Latitude","Year","Start_Date")]
names(myND) = c("x_ll","y_ll","Year","Start_Date")

myInPtDF = SpatialPointsDataFrame(myND[c("x_ll","y_ll")], myND)
# Set the projection of the input spatial dataframe
myInProj = CRS("+proj=longlat +ellps=WGS84")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+init=epsg:32645") # UTM zone 45
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
# Add new projected coordinates to the table
myND$x_utm45 = myOutPtDF@coords[,1]
myND$y_utm45 = myOutPtDF@coords[,2]

# Loads State map  & map positives
#############################################################
myProvShape = readShapePoly("Shapefile/Nepal_DUTM")
# Plot Nepal maps
plot(myProvShape, axes = T)
# Add outbreak locations
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red")



###Loads the file list new: after removing climati cones, production system,distance to indian border,Elevation and forest:

myPredictors = c(
  # Poultry
  "ld_dn", "fw_dn", "lh_dn",
  
  # Anthropogenic
  
  # Topography
  
  #Distances and roads
  "rd_dn", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds",
  # Others
  "po_dn", "hh_dn","iw_pc", "pd_pc" )

# add extensions
myPredictorsASC = paste(myPredictors,".asc", sep = "")
myPredList = as.list(paste(myPredictors,".asc", sep = ""))                


myStack = stack(myPredList)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = raster("mask.asc") >= 0
plot(myMask)
plot(myProvShape, axes = T, add = T)

# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

# 2_Creates a random set of negatives
#############################################################
plot(myProvShape, axes = T)
# Sets the number of points to generate
n = 5000
#Generates random coordinates
x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
# Keep only points according to different conditions
myTDF = as.data.frame(cbind(x,y))
myTDF$IsLand = extract(myMask,myTDF[,1:2])
myTDF$DistToPos = Lib_DistEstimates(myTDF$x, myTDF$y, myND$x_utm45, myND$y_utm45)
myTDF = subset(myTDF, IsLand == 1 & myTDF$DistToPos > 1000)

# Merge the two datasets (positives and negatives
#############################################################
myFD1 = myND[,c("x_utm45","y_utm45")]
names(myFD1) = c("x","y")
myFD1$HasHPAI = 1
myFD2 = myTDF[,1:2]
names(myFD2) = c("x","y")
myFD2$HasHPAI = 0
myFD = rbind(myFD1,myFD2)

# Remove potential duplicates falling in same pixel
#############################################################
myFD$DistToItself = Lib_DistEstimatesToItself(myFD$x,myFD$y)
myFD = subset(myFD, DistToItself > 1000)
myFD = myFD[,-4]
# displays the output
plot(myProvShape, axes = T)
points(myFD$x, myFD$y, pch = 15, col=c("blue","red")[myFD$HasHPAI+1])

# Extract all predictors from raster stack
############################################################# 
myFD = cbind(myFD,extract(myStack,myFD[,1:2]))
str(myFD)
# removes lines with no data

myFD = na.omit(myFD)
