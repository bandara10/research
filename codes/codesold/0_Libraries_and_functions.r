# Loads required libraries
library(raster) #
library(rgdal)  
library(gstat)  #
library(ncf)    #
library(spatstat)
library(foreign)
library(maptools)
library(nlme)   
library(MASS)
library(ROCR)
library(vcd)
library(epicalc)
library(RColorBrewer) # 
library(classInt)






# Function to estimate the distance to the nearest point of the other vector
Lib_DistEstimates <- function (myX1Vec, myY1Vec, myX2Vec,myY2Vec){
MDist = myX1Vec * 0
for (i in 1:length(myX1Vec)) {
  myDistVec = (((myX1Vec[[i]] - myX2Vec)^2) + ((myY1Vec[[i]] - myY2Vec)^2))^0.5
  MDist[i] = min(myDistVec, na.rm = T)
}
return(MDist)
}

# Function to estimate the distance to the nearest point of the same vector
Lib_DistEstimatesToItself <- function (myXVec, myYVec){
myTemp_DF = as.data.frame(cbind(myXVec, myYVec))
myTemp_DF$MDist1 = 0
for (i in 1:nrow(myTemp_DF)) {
  myDistVec = (((myTemp_DF$myXVec[[i]] - myTemp_DF$myXVec)^2) + ((myTemp_DF$myYVec[[i]] - myTemp_DF$myYVec)^2))^0.5
  myDistVec = myDistVec[-i]
  myTemp_DF$MDist1[i] = min(myDistVec)
}
return(myTemp_DF$MDist1)
}

# Function to estimate the autoregressive term
myLib.AutoRegressiveQuantMod = function (myDepVec, myXVec, myYVec, myMaxDist) {
myARVec = myDepVec * 0
nRec =length(myARVec)
for (i in 1:nRec) {
  myDistVec = (((myXVec[i] - myXVec)^2)+((myYVec[i] - myYVec)^2))^0.5

  myDistVecS0 = subset(myDistVec, myDistVec != 0)
  myART = sum(1/myDistVecS0)

  mySGrid = cbind(myDistVec, myDepVec)
  mySSGrid = subset(mySGrid, myDistVec !=0 & myDistVec < myMaxDist)
  myAR = sum(mySSGrid[,2]/mySSGrid[,1])

  myARVec[i] = myAR/myART
myARVec = ifelse(is.na(myARVec),0,myARVec)
}
return(myARVec)
}

# Function to estimate the mean of values in a given neighbourhood
myLib.AutoRegressiveMean = function (myDepVec, myXVec, myYVec, myMaxDist) {
myARVec = myDepVec * 0
nRec =length(myARVec)
for (i in 1:nRec) {
  myDistVec = (((myXVec[i] - myXVec)^2)+((myYVec[i] - myYVec)^2))^0.5
  myTempDF = data.frame(cbind(myDistVec,myDepVec))

  myTempDF = subset(myTempDF, myDistVec != 0)
  myTempDF = subset(myTempDF, myDistVec < myMaxDist)

  myARVec[i] = mean(myTempDF$myDepVec, na.rm = T)
myARVec = ifelse(is.na(myARVec),0,myARVec)
}
return(myARVec)
}


myLib.EstimateLogisticGOF = function(myPred, myObs, myCutOff){
  # AUC
  #####################################################################
  myP = prediction(myPred,myObs)
  myPerf = performance(myP, measure = "auc")
  myAUC = myPerf@y.values[[1]]

  # Cohen's Kappa
  #####################################################################
  myPredClass = ifelse(myPred >= myCutOff,1,0)
  CorPos = sum(myPredClass * myObs)
  CorNeg = sum(ifelse(myPredClass == 1,0,1)*ifelse(myObs == 1,0,1))
  FalPos = sum(myPredClass * ifelse(myObs == 1,0,1))
  FalNeg = sum(ifelse(myPredClass == 1,0,1) * myObs)
  myConfM = matrix(c(CorPos,FalNeg,FalPos,CorNeg),2,2)
  myKappa = Kappa(myConfM)$Unweighted[[1]]

  # Proportion of explained deviance
  #####################################################################
  myDModel = -2*sum(ifelse(myObs == 0,log(1-myPred),log(myPred)))
  myProp = mean(myObs)
  myDNull = -2*sum(ifelse(myObs == 0,log(1-myProp),log(myProp)))
  myPropDev = 1-(myDModel/myDNull)

  myOutputs = c(myAUC,myKappa,myPropDev,myDModel)
  names(myOutputs) = c("AUC","Kappa","PDEV","DEV")
  return(myOutputs)
}
