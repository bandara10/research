
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



