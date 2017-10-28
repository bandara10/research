#3_BRT_Analysis
###################
# Creates the BRT object
nVar = length(myPredictors)

myBRT <- gbm.step(myFD,gbm.x = 4:(nVar+3),gbm.y = 3,family = "bernoulli",tree.complexity = 5,learning.rate = 0.004,bag.fraction = 0.75)
summary(myBRT)
gbm.plot(myBRT)
gbm.plot(myBRT,n.plots=21,write.title = FALSE)
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
plot(myPredR,main="Predictive risk map for H5N1")
plot(myProvShape, axes = T, add = T)
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red", cex = 0.6)
writeRaster(myPredR,"BRTPredictionfao_18_06_14.asc")
gbm.plot.fits(myBRT)
find.int <- gbm.interactions(myBRT)
find.int$interactions
find.int$rank.list
#rainfall_paddy
gbm.perspec(myBRT,17, 15, y.range=c(9,15), z.range=c(0,0.3))