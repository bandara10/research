# Inputs to be checked
#############################################################
nBoot = 10
SensitivityTest = 0
#############################################################


myLRTable = as.data.frame(myPredictors)
myAUCVec = rep(0,nBoot)
# Create the mean raster
myMeanPred = myMask*0
myBRTAUCVec = rep(0,nBoot)
myBRTOut = as.data.frame(1:length(myPredictors))
names(myBRTOut) = "Dummy"
row.names(myBRTOut) = myPredictors


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
#myStr = "HasHPAI ~"
#for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
## Runs the GLM  
#myGLM = glm(as.formula(myStr), data = myFD, family = "binomial")
## Stores the residuals
#myFD$res = residuals(myGLM)
## Crase approach to account for SA
#############################################################
# Build the standard GLM object
#myFD$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), myFD$x, myFD$y, 200000)
## Changes the formular string
#myStrAR = paste(myStr,"+ AR1")
## Re-run the model  
#myGLMR1 = glm(as.formula(myStrAR), data = myFD, family = "binomial")
#summary(myGLMR1)
#
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram

# Runs the BRT
#############################################################
myBRT <- gbm.step(myFD,gbm.x = 4:10,gbm.y = 3,family = "bernoulli",tree.complexity = 5,learning.rate = 0.005,bag.fraction = 0.75)
myBRTSum = summary(myBRT)
myVec = match(myPredictors, as.character(myBRTSum$var))
myBRTOut=cbind(myBRTOut,myBRTSum$rel.inf[myVec])

# Estimates DAIC for each predictor
##############################################################
#myChi2 = rep(0,length(myPredictors))
#myPValue = rep(0,length(myPredictors))
#
#for (i in 1:length(myPredictors)){
#  myShortPred = myPredictors[-i]
#  myStrN = "HasHPAI ~"
#  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
#  myStrARN = paste(myStrN,"+ AR1")
#  myGLMR2 = glm(as.formula(myStrARN), data = myFD, family = "binomial")
##  mylrtest = lrtest(myGLMR1,myGLMR2)
##  myChi2[i] = mylrtest[[4]]
##  myPValue[i] = mylrtest[[6]]
#  myChi2[i] = 0
#  myPValue[i] = 0
#
#}
#
#myLRTable = cbind(myLRTable,myChi2)


# Predictions
#############################################################
#myPred = predict(myARStack, myGLM, type = "response")
#myMeanPred = myMeanPred + myPred

# ROC curve & GOF metrics
#myPred = prediction(predict(myGLM, type = "response"), myFD$HasHPAI)
#perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
#myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), myFD$HasHPAI, 0.5)
myPredBRTPerfs = myLib.EstimateLogisticGOF(myBRT$fitted, myFD$HasHPAI, 0.5)

myAUCVec[i] = myPredPerfs[[1]]
myBRTAUCVec[i] = myPredBRTPerfs[[1]]
}


#myFTable = as.data.frame(rowMeans(myLRTable[,-1]))
#myFTable = cbind(myPredictors,myFTable)
#names(myFTable)=c("Predictors","Chisq")
#myFTable$pvalue = pchisq(myFTable$Chisq,1,lower.tail = F)

myBRTOutSum = rowMeans(myBRTOut[,-1])


#myMeanPred = myMeanPred / nBoot

#plot(myMeanPred)
#points(myND$x_utm46, myND$y_utm46, pch = 15, col = "red", cex = 0.6)
#myFTable



