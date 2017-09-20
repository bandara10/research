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
#library(dplyr)
library(usdm)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack <- stack(myfullstack.a)
#plot(myfullstack)
#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####
hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# plot koala presence locations
hefleydata.presence <-subset(hefleydata, presence==1)
coordinates(hefleydata.presence) <- ~x+y
plot(hefleydata.presence)

# get only coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]
hefleydata.s = as.data.frame(hefleydata.s)
#plot(hefleydata.s)
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
set.seed(1235)
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 
ZTGLM <- ZTGLM.myFD5

##### Step 6: now take a random sample of n% and assign detected 1 non detected 0.####

set.seed(1234)
train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = floor(0.50 * nrow(ZTGLM.myFD1)),replace=FALSE) # split the presence only dataset as detected and nondetected.
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0
##### Step 7: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected)  # This dataset length is same as presence dataset.
#Detection.data[53] <- lapply(Detection.data[53], as.numeric)
str(Detection.data)
#Detection.data <- Detection.data[c(-1,-2)]
#Detection.data <- subset(Detection.data, select=c(50, 1:49, 51))
IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus n=1000 deteected koalas data n=80 from ZTGLM.myFD4
#IPP.data <- IPP.data[c(-1,-2,-53)]
#IPP.data <- subset(IPP.data, select=c(50, 1:49))
ZTGLM.data=(detected)##ZTGLM.data# This is detected data randomly selected n= 80 .
#ZTGLM.data <- ZTGLM.data[c(-1,-2,-52)]
#ZTGLM.data <- subset(ZTGLM.data, select=c(50, 1:49))
##### Step 8:  analysis without detection correction factor#######
IPP.ignored=glm(presence~twi + tpo + temp + aspect + elev+habit2pc+hpop+lot_density+sbd,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
summary(IPP.ignored)
ZTGLM.ignored=vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,family="pospoisson", data=ZTGLM.data)
summary(ZTGLM.ignored)

#unclass(summary(Detection.model))
#Hefley method GLM
set.seed(123) # we create detection probabilities using two methods. glm, rf
#Detection model: steps as in Hefley`s code`

Detection.model=glm(presence~  distance_pedestrian + s1_residential_dist + distance_trunkandlink +
                      distance_tertiaryandlink+scale(group),family= "binomial", data=Detection.data)

summary(Detection.model)
#myPred = prediction(predict(Detection.model, type = "response"), Detection.data$presence)
#perf <- performance(myPred,measure = "tpr", x.measure = "fpr")
#plot(perf, colorize = T)
#summary(perf)
#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det = faraway::ilogit(predict(Detection.model, type="response",new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
#myBRT,n.trees=myBRT$gbm.call$best.trees ; add this as model in the above function.

hist(p.det, breaks=100)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det
## Stores the residuals and 
Detection.data$res = residuals(Detection.model) # library(ncf)
myResCorr <- correlog(Detection.data$x, Detection.data$y, Detection.data$res,na.rm=T, increment=10000, resamp=0, latlon = F)
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")

#############################################################
# Build the standard GLM object
Detection.data$Autoregres = myLib.AutoRegressiveMean(residuals(Detection.model), Detection.data$x, Detection.data$y, 1000)
# Re-run the model 
Detection.model.c=glm(presence~  distance_pedestrian + s1_residential_dist + distance_trunkandlink +
                      distance_tertiaryandlink+Autoregres+scale(group),family= "binomial", data=Detection.data)
summary(Detection.model.c)
# Estimate correlogram of new residuals
Detection.data$CraserRes = residuals(Detection.model.c)
CraserCorr <- correlog(Detection.data$x, Detection.data$y, Detection.data$CraserRes,na.rm=T, increment=2000, resamp=0, latlon = T)
# Plot the correlogram
plot(CraserCorr$mean.of.class[1:10], CraserCorr$correlation[1:10] ,type="b", pch=16, lwd=1.5, cex = 1.2,xlab="distance", ylab="Moran's I")
abline(h=0)
lines(CraserCorr$mean.of.class, CraserCorr$correlation, col = "red")
# Produces the spatial prediction

# Creates a dummy variable
Autoregres = myMask * 0
# Add it to the stack 
myfullstack = addLayer(myfullstack, Autoregres)
myfullstack@layers[length(myfullstack@layers)] = "Autoregres"
# Apply the prediction
myPred.new = predict(myfullstack, Detection.model.c, type = "response")
# plot the prediction
plot(myPred.new)
