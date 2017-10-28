#1_Load Data
# Sets the general path
#############################################################
setwd("/Users/ravidissanayake/Documents/Dropbox/Lubies-FaoIndia/Data")


# Loads daya
#############################################################
myD = read.dbf("CaseSeriesBlockTotal/FinalCaseSeriesBlockWiseWBAT.dbf", as.is=TRUE)
#myD = subset(myD, TYPE == "OUT")

# Reproject the data
#############################################################
myND = myD[,c("LONGITUDE","LATITUDE","Year","Startdate")]
names(myND) = c("x_ll","y_ll","Year","Startdate")
myInPtDF = SpatialPointsDataFrame(myND[c("x_ll","y_ll")], myND)
# Set the projection of the input spatial dataframe
myInProj = CRS("+proj=longlat +ellps=WGS84")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+init=epsg:32646") # UTM zone 46
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
# Add new projected coordinates to the table
myND$x_utm46 = myOutPtDF@coords[,1]
myND$y_utm46 = myOutPtDF@coords[,2]

# Loads State map  & map positives
#############################################################
myProvShape = readShapePoly("Shapefile/WestBengalAssamTripura")
# Plot China maps
plot(myProvShape, axes = T)
# Add outbreak locations
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red")


# Loads predictors and masks
#############################################################
# loads the file list               

myPredictors = c(
                # Poultry
                "brdn","dsdn","dudn",
                # Anthropogenic
                "hpdn","access",
                # Topography
                "dem",
                # Others
                "ncrop")

# add extensions
myPredictorsASC = paste(myPredictors,".asc", sep = "")
myPredList = as.list(paste(myPredictors,".asc", sep = ""))                
# Creates a raster stack
myStack = stack(myPredList)
# Creates a mask for the analyses
myMask = raster("mask.asc")
plot(myMask)
plot(myStack)
#Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

###################################################################
#SAme session New commands
###################################################################
myD = read.dbf("CaseSeriesBlockTotal/FinalCaseSeriesBlockWiseWBAT.dbf")

#############################################################
myND = myD[,c("LONGITUDE","LATITUDE","Year","Startdate")]
names(myND) = c("x_ll","y_ll","Year","Startdate")
myInPtDF = SpatialPointsDataFrame(myND[c("x_ll","y_ll")], myND)
# Set the projection of the input spatial dataframe
myInProj = CRS("+proj=longlat +ellps=WGS84")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+init=epsg:32646") # UTM zone 46
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
# Add new projected coordinates to the table
myND$x_utm46 = myOutPtDF@coords[,1]
myND$y_utm46 = myOutPtDF@coords[,2]

# Loads State map  & map positives
#############################################################
myProvShape = readShapePoly("Shapefile/WestBengalAssamTripura")
# Plot maps
plot(myProvShape, axes = T)
# Add outbreak locations
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red")
myProvShape = readShapePoly("Shapefile/WestBengalAssamTripura")

# Loads predictors and masks
#############################################################
# loads the file list
myPredList = list("brdn.asc","dsdn.asc","dudn.asc","hpdn.asc","waterperpix.asc","access.asc","dem.asc","ncrop.asc")
# Creates a raster stack
myStack = stack(myPredList)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = raster("chdn.asc") >= 0
#plot(myMask)

plot(myStack[[7]])
plot(myStack[[6]])
plot(myStack[[5]])
plot(myStack[[4]])
plot(myStack[[3]])
plot(myStack[[2]])
plot(myStack[[1]])

points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red")
plot(myProvShape, axes = T, add = T)

######################################################################
########################################################################
# 2_Creates a random set of negatives
#############################################################
plot(myProvShape, axes = T)
# Sets the number of points to generate
n = 2000
#Generates random coordinates
x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
# Keep only points according to different conditions
myTDF = as.data.frame(cbind(x,y))
myTDF$IsLand = extract(myMask,myTDF[,1:2])
myTDF$DistToPos = Lib_DistEstimates(myTDF$x, myTDF$y, myND$x_utm46, myND$y_utm46)
myTDF = subset(myTDF, IsLand == 1 & myTDF$DistToPos > 3000)

# Merge the two datasets (positives and negatives
#############################################################
myFD1 = myND[,c("x_utm46","y_utm46")]
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
##########################################################################
######################################################################
#Boootstrap_Analysis
# Inputs to be checked
#############################################################
nBoot = 10
SensitivityTest = 0
#############################################################
#2_Plot Clusters
#################################################################
# Functions to draw circles and ellipses
#############################################################

myLib_DrawCircle = function(xc,yc,rc,nseg,mycol, mylty, mylwd){
  
  xx <- xc + rc*cos( seq(0,2*pi, length.out=nseg) )
  yy <- yc + rc*sin( seq(0,2*pi, length.out=nseg) )
  
  lines(xx,yy, col=mycol, lty = mylty, lwd = mylwd)
}

myLib_DrawEllipse = function(xc,yc,rl,rc,angle,nseg,mycol, mylty, mylwd){
  
  theta = seq(0, 2 * pi, length=(nseg))
  xx <- xc + rl * cos(theta) * cos(angle) - rc * sin(theta) * sin(angle)
  yy <- yc + rl * cos(theta) * sin(angle) + rc * sin(theta) * cos(angle)
  
  
  lines(xx,yy, col=mycol, lty = mylty, lwd = mylwd)
}



# Plot the outbreaks map
#############################################################
myProvShape = readShapePoly("Shapefile/WestBengalAssamTripura")
# Plot China maps
plot(myProvShape, axes = T, col = "white")
# Add outbreak locations
points(myND$x_utm46, myND$y_utm46, pch = 15, cex= 1, col = "red")

# Loads the clusters coordinates
myCC = read.table("Clusters.txt", header = T)
str(myCC)
# Get the UTM_46 coordinates
myInPtDF = SpatialPointsDataFrame(myCC[c("Longitude","Latitude")], myCC)
# Set the projection of the input spatial dataframe
myInProj = CRS("+proj=longlat +ellps=WGS84")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
# plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+init=epsg:32646") # UTM zone 46
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
# Add new projected coordinates to the table
myCC$x_utm46 = myOutPtDF@coords[,1]
myCC$y_utm46 = myOutPtDF@coords[,2]


##################################################################
# Subset to Spatial Cluster 30%
##################################################################
myMethod = "SP30"
myCC1 = subset(myCC, Method == myMethod)
myCC1
# Plot the background
plot(myProvShape, axes = T, col = "grey", main = myMethod)
# Add outbreak locations
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red")
for (i in 1:nrow(myCC1)){
  myLib_DrawCircle(myCC1$x_utm46[i],myCC1$y_utm46[i],myCC1$Radius[i]*1000,360,"blue", 1, 2)
}
#reads the ellipse data
myED = read.dbf("Spatial30.col.dbf", as.is = T)
myED = subset(myED, P_VALUE < 0.05)
for (i in 1:nrow(myED)){
  myLib_DrawEllipse(myED$X[i],myED$Y[i],myED$E_MAJOR[i],myED$E_MINOR[i],myED$E_ANGLE[i]*pi / 180,360,"green", 1, 2)
}
myED


myLRTable = as.data.frame(myPredictors)
myAUCVec = rep(0,nBoot)
# Create the mean raster
myMeanPred = myMask*0


for (i in 1:nBoot){
  
  
  # Creates a random set of negatives
  #############################################################
  #plot(myProvShape, axes = T)
  # Sets the number of points to generate
  n = 2000
  # Generates random coordinates
  x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
  y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
  # Keep only points according to different conditions
  myTDF = as.data.frame(cbind(x,y))
  myTDF$IsLand = extract(myMask,myTDF[,1:2])
  myTDF$DistToPos = Lib_DistEstimates(myTDF$x, myTDF$y, myND$x_utm46, myND$y_utm46)
  myTDF = subset(myTDF, IsLand == 1 & myTDF$DistToPos > 3000)
  
  # Merge the two datasets (positives and negatives
  #############################################################
  myFD1 = myND[,c("x_utm46","y_utm46")]
  names(myFD1) = c("x","y")
  myFD1$HasHPAI = 1
  myFD2 = myTDF[,1:2]
  names(myFD2) = c("x","y")
  myFD2$HasHPAI = 0
  if (SensitivityTest > 0){
    myVec = sample(floor(nrow(myFD2)*0.10), nrow(myFD2), replace = T)
    myFD2$HasHPAI[myVec] = 1
  }
  
  
  myFD = rbind(myFD1,myFD2)
  
  # Remove potential duplicates falling in same pixel
  #############################################################
  myFD$DistToItself = Lib_DistEstimatesToItself(myFD$x,myFD$y)
  myFD = subset(myFD, DistToItself > 1000)
  myFD = myFD[,-4]
  # displays the output
  #plot(myProvShape, axes = T)
  #points(myFD$x, myFD$y, pch = 15, col=c("blue","red")[myFD$HasHPAI+1])
  
  # Extract all predictors from raster stack
  ############################################################# 
  myFD = cbind(myFD,extract(myStack,myFD[,1:2]))
  # str(myFD)
  # removes lines with no data
  myFD = na.omit(myFD)
  
  # Runs a standard Logistic regression
  #############################################################
  # Build the GLM object
  myStr = "HasHPAI ~"
  for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
  # Runs the GLM  
  myGLM = glm(as.formula(myStr), data = myFD, family = "binomial")
  # Stores the residuals
  myFD$res = residuals(myGLM)
  # Crase approach to account for SA
  #############################################################
  # Build the standard GLM object
  myFD$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), myFD$x, myFD$y, 200000)
  # Changes the formular string
  myStrAR = paste(myStr,"+ AR1")
  # Re-run the model  
  myGLMR1 = glm(as.formula(myStrAR), data = myFD, family = "binomial")
  summary(myGLMR1)
  
  #points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
  # Plot the correlogram
  
  # Estimates DAIC for each predictor
  #############################################################
  myChi2 = rep(0,length(myPredictors))
  myPValue = rep(0,length(myPredictors))
  
  for (i in 1:length(myPredictors)){
    myShortPred = myPredictors[-i]
    myStrN = "HasHPAI ~"
    for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
    myStrARN = paste(myStrN,"+ AR1")
    myGLMR2 = glm(as.formula(myStrARN), data = myFD, family = "binomial")
    mylrtest = lrtest(myGLMR1,myGLMR2)
    myChi2[i] = mylrtest[[4]]
    myPValue[i] = mylrtest[[6]]
  }
  
  myLRTable = cbind(myLRTable,myChi2)
  
  
  # Predictions
  #############################################################
  myPred = predict(myARStack, myGLMR1, type = "response")
  myMeanPred = myMeanPred + myPred
  
  # ROC curve & GOF metrics
  myPred = prediction(predict(myGLM, type = "response"), myFD$HasHPAI)
  perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
  myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), myFD$HasHPAI, 0.5)
  myPredPerfs
  myAUCVec[i] = myPredPerfs[[1]]
}
myFTable = as.data.frame(rowMeans(myLRTable[,-1]))
myFTable = cbind(myPredictors,myFTable)
names(myFTable)=c("Predictors","Chisq")
myFTable$pvalue = pchisq(myFTable$Chisq,1,lower.tail = F)

myMeanPred = myMeanPred / nBoot

plot(myMeanPred)
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red", cex = 0.6)
myFTable

#writeRaster(myMeanPred, "c:/Temp/myPrediction.asc") 
#########################################################################
##########################################################################
#3_BRT_Analysis
###################
# Creates the BRT object
nVar = length(myPredictors)

myBRT <- gbm.step(myFD,gbm.x = 4:(nVar+3),gbm.y = 3,family = "bernoulli",tree.complexity = 5,learning.rate = 0.005,bag.fraction = 0.75)
summary(myBRT)
gbm.plot(myBRT)

myPDF = as.data.frame(values(myMask))

# Creates a dataframe with the predictors
for (i in 1:length(myPredictors)){
  myVec = values(raster(paste(myPredictors[i],".asc",sep="")))
  myPDF = cbind(myPDF,myVec)
}
names(myPDF) = c("Mask",myPredictors)
# Get the predicted values
myPredValues = predict.gbm(myBRT, n.trees = myBRT$gbm.call$best.trees,tree.complexity = 5, learning.rate = 0.005,newdata = myPDF, type = "response")
# Convert the predicted values to a raster
myPredR = setValues(myMask,myPredValues)*myMask
plot(myPredR)
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red", cex = 0.6)
writeRaster(myPredR,"BRTPrediction.asc")

############################################################
###########################################################
#Ellipse
############################################################
Elipse centered at (x0, y0) with axes _a_ and _b_:
  
  theta <- seq(0, 2 * pi, length=(big enough))
x <- x0 + a * cos(theta)
y <- y0 + b * sin(theta)
plot(x, y, type = "l")

:-)

Or, if you want to rotate it alpha degrees from the x-axis:
  
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
###################################################################
###################################################################
#3_Logistic_Regression
###################################################################
myPredictors = c(
  # Poultry
  "dudn","dsdn","brdn",
  # Anthropogenic
  "hpdn","access",
  # Topography
  "dem","ncrop")



# Runs a standard Logistic regression
#############################################################
# Build the GLM object
myStr = "HasHPAI ~"
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
myResCorr <- correlog(myFD$x, myFD$y, myFD$res,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(myResCorr$p<0.05,1,0)  
#plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
#xlab="distance", ylab="Moran's I")
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")
abline(h=0)

# Crase approach to account for SA
#############################################################
# Build the standard GLM object
myFD$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), myFD$x, myFD$y, 200000)
# Changes the formular string
myStrAR = paste(myStr,"+ AR1")
# Re-run the model  
myGLMR1 = glm(as.formula(myStrAR), data = myFD, family = "binomial")
summary(myGLMR1)
plot(myGLMR1)  

# Estimate correlogram of new residuals

myFD$R1 = residuals(myGLMR1)
Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram


# Estimates DAIC for each predictor
#############################################################
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "HasHPAI ~"
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
plot(myPred, xlab = "x-coord (UTM46)", ylab= "y-coord (UTM46)")
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red", cex = 0.6)

plot(myPred)
points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red", cex = 0.6)

# ROC curve & GOF metrics
myPred = prediction(predict(myGLM, type = "response"), myFD$HasHPAI)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), myFD$HasHPAI, 0.5)
myPredPerfs

#test = raster("C:\Dropbox\5_Supervision\Lubies-FaoIndia\Outputs\SatscanElliptical\






