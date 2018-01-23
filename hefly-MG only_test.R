Detection.data$res = residuals(Detection.model)

myResCorr <- correlog(Detection.data$x, Detection.data$y, Detection.data$res, na.rm=T, increment=10000, resamp=0, latlon = F)

plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,xlab="distance", ylab="Moran's I")

abline(h=0)



# Crase approach to account for SA
#############################################################
myStack = stack(habitat.rr)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = myStack$Annual_Precipitation
plot(myMask)
# Map predictions
AR1 = myMask * 0
plot(AR1)
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

# Crase approach to account for SA
#############################################################
# Build the standard GLM object
source("autoregresive2.r")
Detection.data$AR1 = myLib.AutoRegressiveMean(residuals(Detection.model), Detection.data$x, Detection.data$y, 15000)


Detection.model.2=glm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    + AR1
                    ,family= "binomial"
                    ,data=Detection.data)
summary(Detection.model.2)


plot(Detection.model.2)

Detection.data$R1 = residuals(Detection.model.2)
Corr <- correlog(Detection.data$x, Detection.data$y, Detection.data$R1,na.rm=T, increment=10000, resamp=0, latlon = F)              
#Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

#lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
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
  myGLMR2 = glm(as.formula(myStrARN), data = acsel22, family = "binomial")
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

# ROC curve & GOF metrics
myPred = prediction(predict(myGLM, type = "response"), acsel22$Haskoala)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), acsel22$Haskoala, 0.5)
myPredPerfs




