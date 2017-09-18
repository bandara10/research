#3_BRT_Analysis
###################
library(dismo)
library(gbm)
library(spatial.tools)
library(VGAM)
library(mosaic)
library(spatstat)
library(faraway)
library(raster)
library(broom)
library(dismo)
library(randomForest)
library(forestFloor)
library(AUC)
library(rgl)
myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack = stack(myfullstack.a)
#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####

hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# plot koala presence locations
hefleydata.presence <-subset(hefleydata, presence==1)
coordinates(hefleydata.presence) <- ~x+y
plot(hefleydata.presence, add=TRUE)
# get only coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]
hefleydata.s = as.data.frame(hefleydata.s)
#select presence records

#### Step 3:  extract X=vector.boot from raster anad combine wth hefleydata, presence and group varibels.####
myFD = cbind(hefleydata.s,raster::extract(myfullstack,hefleydata.s))
# get presence and group vriables
myFD1 = cbind(myFD,hefleydata [6]) # get presence and group data. 
myFD1 = na.omit(myFD1) # remove all NA valuves

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 

#####Step 5: select only 1000 absences (monticarlo points as in hefleys method??)####
set.seed(12345)
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 
###
##### Step 6: now take a random sample of 80 and assign detected 1 non detected 0.####
train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = 80,replace=FALSE)
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0
##### Step 7: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected) 
Detection.data[53] <- lapply(Detection.data[53], as.numeric)
IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus no koalas data from myFD1.

ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.



#  use BRT to predict detectoin probabilities
#   https://cran.r-project.org/web/packages/gbm/gbm.pdf
# Creates the BRT object
# select varibles 
newZTGLM5 <- ZTGLM.myFD5[c(1,2,53,8:23, 37:45)]
myBRT <- gbm.step (newZTGLM5,gbm.x = 4:28,gbm.y = 3,family = "bernoulli", tree.complexity = 5, learning.rate = 0.01, bag.fraction = 0.5)
par(mgp=c(3,1,0),mar=c(10,12,3,2)+0.1) # set mrgins for the plot
summary(myBRT, las=1)
summary(myBRT)
gbm.plot(myBRT)
dev.off()
par(mfrow=c(4,4))
gbm.plot(myBRT)
#get predictions and plot full model.
predictions <- predict(myfullstack, myBRT, n.trees=myBRT$gbm.call$best.trees, type="response") #
plot(predictions)
#########Now fit the reduced model by selecting covariate most infulential
#         covariates selected from the dataset;newZTGLM5 : distance_trunkandlink         distance_trunkandlink 17.6012051
#distance_trunkandlink         distance_trunkandlink 17.6012051
#distance_steps                       distance_steps  9.9050443
#elev                                           elev  8.5460861
#distance_cycleway                 distance_cycleway  7.4069327
#distance_bridleway               distance_bridleway  7.3890307
set.seed(125)
myBRT.2 <- gbm.step(newZTGLM5, gbm.x = c(5,6,14,16,18), gbm.y = 3) #Build initial model #,family = "bernoulli",tree.complexity = 2,learning.rate = 0.01,bag.fraction = 0.75
gbm.plot(myBRT.2)
summary(myBRT.2, las=1)
predictions.2 <- predict(myfullstack, myBRT.2, n.trees=myBRT.2$gbm.call$best.trees, type="response")
#plot predictions
plot(predictions.2,main="Detection proabilities ")
plot(hefleydata.presence,  cex = 0.3,add=TRUE)
p.det = faraway::ilogit(predict(myBRT.2,n.trees=myBRT$gbm.call$best.trees, type="response",new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
hist(p.det, breaks=50)
ZTGLM.data$p.det=p.det
######
# This section for BRT for my interest
#
#writeRaster(predictions.2,"Detection proabilities_BRT.asc")
gbm.plot.fits(myBRT.2)
find.int <- gbm.interactions(myBRT.2)
find.int$interactions
DD<-find.int$rank.list
#rainfall_paddy
par(mgp=c(10,10,0),mar=c(0,3,3,2)+0.1)
gbm.perspec(myBRT,4, 3)
myBRT$cv.loss.matrix
######### Now we use this detection probabilities instead of detection probabilities derived from Hefley`s method.`

ZTGLM.data$p.det=p.det

######Step 5: - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected= glm(presence~twi + tpo + temp + aspect + elev+habit2pc+hpop+lot_density+sbd,
                   family="binomial",data=IPP.data)
summary(IPP.corrected)
# broom package: used tidy to get a table from model outputs. This doesnt work for VGLMs.
# get confidence intervals
confidenceintervals <- confint(IPP.corrected)
tidy(IPP.corrected,confidenceintervals)
####Step 6: Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####

#use only the significant covariates, tpo +hpop+lot_density+sbd
# VGAM: read about which family to use: https://www.r-project.org/doc/Rnews/Rnews_2008-2.pdf
ZTGLM.corrected = vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd
                       ,weights=1/p.det,family="pospoisson",data=ZTGLM.data) # zapoisson
summary(ZTGLM.corrected)
ZTGLM.corrected
qtplot(ZTGLM.corrected)
# step 7:  Map predictions
myPred = predict(myfullstack, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")
plot(hefleydata.presence, add=TRUE)
myPred2 = predict(myfullstack, IPP.corrected, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model-intensity of group or ??,")
plot(hefleydata.presence, add=TRUE)
myPred3.1 = predict(myfullstack, ZTGLM.corrected, type = "response")
plot(myPred3.1,  main="ZTGLM-Number of koalas in a grid - VGLM ")
plot(hefleydata.presence, add=TRUE)
writeRaster(myPred3, "ZTGLM.tif")
dev.off()






