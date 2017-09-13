#3_BRT_Analysis
###################
library(dismo)
library(gbm)
# Creates the BRT object
nVar = length(myfullstack)
 # select varibles 
newZTGLM5 <- ZTGLM.myFD5[c(1,2,53,8:23)]
myBRT <- gbm.step(newZTGLM5,gbm.x = 4:19,gbm.y = 3,family = "bernoulli",tree.complexity = 2,learning.rate = 0.01,bag.fraction = 0.75)
plot(myBRT)
par(mgp=c(3,1,0),mar=c(10,12,3,2)+0.1) # set mrgins for the plot
summary(myBRT, las=2, asp = 1)
dev.off()
par(mfrow=c(4,4))
gbm.plot(myBRT, n.plots = 16, write.title = FALSE)
myMask <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack\\mask\\mask.tif")
myPDF = as.data.frame(values(myMask))

predictions <- predict(myfullstack, myBRT, n.trees=myBRT$gbm.call$best.trees, type="response")
plot(predictions)

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
plot(myPredR,main="Predictive risk map for H5N1")
plot(myProvShape, axes = T, add = T)
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red", cex = 0.6)
writeRaster(myPredR,"BRTPredictionfao_18_06_14.asc")
gbm.plot.fits(myBRT)
find.int <- gbm.interactions(myBRT)
find.int$interactions
find.int$rank.list
#rainfall_paddy
gbm.perspec(myBRT,13, 11, y.range=c(3,15), z.range=c(0,1))
