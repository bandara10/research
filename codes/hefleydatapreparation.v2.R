library(spatial.tools)
library(VGAM)
library(mosaic)
library(spatstat)
library(faraway)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")


#below is how to create a raster stack from folder of rasters in two steps. set working directory to folder.
#setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack\\raster_stackfull")
#myfullstack <- list.files(pattern="\\.tif") # $ identify only tiff files. just ignore it for this folder with tifs only.

#stackfullnew <- stack(stackfull)
#go to mystakfull in line 103

#load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack\\hefley.RData")
hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)

# Subset the data:
hefleydata.s <- hefleydata[c("x","y")]
#str(hefleydata.s)
#coords <- SpatialPoints(hefleydata.s[, c("x","y")])
#hefleydata.s <- SpatialPoints(coords)
#p<-SpatialPointsDataFrame()
#plot(hefleydata.s)
#proj4string(hefleydata.s) <- CRS("+init=epsg:28356")

#setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
hpop <- raster("hpop.tif")
awc <- raster("awc.tif")
tpo <- raster("tpo.tif")
rain <- raster("rain.tif")
nitro <- raster("nitro.tif")
elev <- raster("elev.tif")
temp <- raster("temp.tif")
sbd <- raster("sbd.tif")
roaden <- raster("roadk.tif")
clay <- raster("clay.tif")
fpc <- raster("fpc.tif")
aspect <- raster("aspect.tif")
twi <- raster("twi.tif")
lot_density <- raster("lot_density.tif")
suitable_3 <- raster("suitable_3.tif")
suitable_2 <- raster("suitable_2.tif")
suitable_1 <- raster("suitable_1.tif")
habit3 <- raster("habit3.tif")
habit2 <- raster("habit2.tif")
habit1 <- raster("habit1.tif")
#myPredictors.new = c("lot_density", "suitable_3","suitable_2", "suitable_1","habit3", "habit2","habit1")
#myPredictors.new = as.list(paste(myPredictors.new,".TIF", sep = ""))
#myStack.habit = stack(myPredictors.new)
#myStack.habitc <- crop(myStack.habit,fpc) 
#myStack.habitc<-spatial_sync_raster(myStack.habitc,fpc, method = "ngb")
#writeRaster(myStack.habitc, filename=names(myStack.habitc), bylayer=TRUE,format="GTiff", overwrite=TRUE)


myPredictors.n = c("hpop","elev","fpc","rain","roadk", "temp", "clay","awc","tpo",  "nitro", "sbd","aspect", "twi","habit1", "habit2", "habit3","suitable_1","suitable_2","suitable_3","lot_density")
#myPredictors.na = paste(myPredictors.n,".TIF", sep = "")
myPredictors.na = as.list(paste(myPredictors.n,".TIF", sep = ""))  
myStack.hefely = stack(myPredictors.na)
myFD=as.data.frame(hefleydata.s)
myFD1 = cbind(myFD,extract(myStack.hefely,myFD)) # bind xy coordinates.
boxplot(myFD1[,3:22])
boxplot(myStack.hefely,las=2)
plot(myStack.hefely,asp=1)
hist(myStack.hefely, las=2)
boxplot(myStacksyn[,3:22],las=2)
plot(myStacksyn)
hist(myStacksyn, las=2)

dev.off()
#write.csv(myFD, file="myFD.csv", row.names=FALSE) 
# --------------------------------------------------------------------------------------------------------------------
#setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
d_bridleway <- raster("distance_bridleway.tif")
d_cycleway <- raster("distance_cycleway.tif")
d_footway <- raster("distance_footway.tif")
d_motorwayandlink <- raster("distance_motorwayandlink.tif")
d_path <- raster("distance_path.tif")
d_pedestrian <- raster("distance_pedestrian.tif")
d_primaryandlink <- raster("distance_primaryandlink.tif")
d_residentil <- raster("distance_residentil.tif")
d_secondaryandlink <- raster("distance_secondaryandlink.tif")
d_steps <- raster("distance_steps.tif")
d_tertiaryandlink <- raster("distance_tertiaryandlink.tif")
d_trunkandlink <- raster("distance_trunkandlink.tif")
d_unclassified <- raster("distance_unclassified.tif")

s1_residential_dist<- raster("s1_residential_dist.tif") 
s1_secondry_dist<- raster ("s1_secondry_dist.tif")     
s1_unclassified_dist<- raster("s1_unclassified_dist.tif")
s2_residential_dist<- raster("s2_residential_dist.tif")  
s2_secondry_dist<- raster("s2_secondry_dist.tif")    
s2_unclassified_dist<- raster("s2_unclassified_dist.tif") 
s3_residential_dist<- raster("s3_residential_dist.tif") 
s3_secondry_dist<- raster("s3_secondry_dist.tif")     
s3_unclassified_dist<- raster("s3_unclassified_dist.tif")



#---------- 
#windows(); plot(d_footway)
myPredictors = c(
  "distance_bridleway","distance_cycleway","distance_footway", "distance_motorwayandlink", "distance_path", "distance_pedestrian", "distance_primaryandlink", "distance_residentil",
  "distance_secondaryandlink", "distance_steps", "distance_tertiaryandlink", "distance_trunkandlink", "distance_unclassified",
  "s1_residential_dist","s1_secondry_dist","s1_unclassified_dist","s2_residential_dist","s2_secondry_dist","s2_unclassified_dist",
  "s3_residential_dist","s3_secondry_dist", "s3_unclassified_dist")
# add extensions
#myPredictorsASC = paste(myPredictors,".TIF", sep = "")
myPredList = as.list(paste(myPredictors,".TIF", sep = ""))  
myStack = stack(myPredList)
#myStack.c <- crop(myStack, extent(441829, 541829, 6901098, 7001098) ) # different extent so I cropped them to match previous dataset.
#myStacksyn<-spatial_sync_raster(myStack.c, myStack.hefely, method = "ngb") # synchronised for same sxtent
#writeRaster(myStacksyn, filename=names(myStacksyn), bylayer=TRUE,format="GTiff", overwrite=TRUE)
#plot(myStack)
#merge tow stacks tgether to one stack
myfullstack= stack(myStack.hefely,myStack)
plot(myfullstack[[1:16]], asp=1)
hist(myfullstack[[1:16]], las=2)
hist(myfullstack[[17:32]], las=2)
hist(myfullstack[[33:42]], las=2)
#myStackcc=unstack(myStacksyn)
#d_pedestrian <- subset(myStacksyn,6) # subset raster stack to get  new extent raster
#d_unclassified<-subset(myStacksyn,12) #subset raster stack to get  new extent raster
myFD=as.data.frame(hefleydata.s)
myFD3 = cbind(myFD,extract(myfullstack,myFD))
#write.csv(myFD3, file="myFD3.csv", row.names=FALSE) 

#myFD4 = cbind(myFD1,myFD3[,3:15]) #select covariates
myFD6 = cbind(myFD3,hefleydata [,6:7]) # get presence and group data. crop hefley data to study area.
#myFD6 = na.omit(myFD5)
ZTGLM.myFD6=myFD6[which(myFD6$presence==1),]
ZTGLM.myFD7=myFD6[which(myFD6$presence==0),]
set.seed(123)
ZTGLM.myFD8 <- sample(seq_len(nrow(ZTGLM.myFD7)), size = 1000)#use for ipp.data
ZTGLM.myFD9 <- ZTGLM.myFD7[ZTGLM.myFD8, ] #x.int data frame now
#IWLR data set similar to hefley.
ZTGLM.myFD10=rbind(ZTGLM.myFD6,ZTGLM.myFD9) # detected and montecarlo points hust a wild guess
#--------------------image saved till this poiint/data preparation
save.image("hefleynew.RData")
--------------------------
#select only 1000 observations.
#now take a random sample and assign detected 1 non detected 0.
#set.seed(123)
#sample=ZTGLM.myFD6[sample(nrow(ZTGLM.myFD6), 50), ]
#data radomly partitioned to detected and non detected
set.seed(123)
train_ind <- sample(seq_len(nrow(ZTGLM.myFD6)), size = 80)

detected <- ZTGLM.myFD6[train_ind, ]
notdetected <- ZTGLM.myFD6[-train_ind,] 
#not detected assigned valuve 0
notdetected$presence<- 0
#detection dataset is
Detection.data= rbind(detected,notdetected) # 156 ===1619x2

#IPP.data comes from detected data plus no koalas data from myFD6.
IPP.data=rbind(detected, ZTGLM.myFD9) # 1050 ====1461
#ZTGLM.data# get the appropriate 50======= #461x2
ZTGLM.data=(detected)
#####     analysis#######
IPP.ignored=glm(presence~s3_unclassified_dist,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
ZTGLM.ignored=vglm(group~awc+elev,family="pospoisson",data=ZTGLM.data)
###
Detection.model=glm(presence~s3_unclassified_dist ,family="binomial",data=Detection.data) #length 1619 #3 variables #reduce detection data to 200.

#Step 4 - Estimate the probability of detection for each presence-only location.

p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 variables 
hist(p.det,breaks=50)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD8)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det
##Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det
IPP.corrected=glm(presence~awc+elev,family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data) #length 1461. 3 variables.
summary(IPP.corrected)


##Step 6 - Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det
ZTGLM.corrected=vglm(group~awc+elev,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)


sink("modelsummary.txt")
ABC=list(IPP.ignored,ZTGLM.ignored, Detection.model, IPP.corrected, ZTGLM.corrected)
print(ABC)
colMeans(bootstrap.sample)[1]
colMeans(bootstrap.sample)[2]
colMeans(bootstrap.sample)[3]
colMeans(bootstrap.sample)[4]
colMeans(bootstrap.sample)[5]
colMeans(bootstrap.sample)[6]
sink()
sink()
sink()
dev.off()
plot(AR1)
AR1 = tpo * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(wac, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"
myPred = predict(myARStack, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")
myPred2 = predict(myARStack, IPP.corrected, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="IPP.model")
myPred2 = predict(myARStack, ZTGLM.corrected , type = "response")
plot(myPred, xlab = "x", ylab= "y",main="ZTGLM.corrected")
###############################################################################
###############################################################################
#	Calculating mean, standard deviation and 95%, equal-tailed confidence intervals
#	from the empirical distribution. See "Introduction to the Bootstrap" (Efron & Tibshirani 1994) 
#	for more details.
#
#     Estimates of the coefficient for the covariate "x" (0.9520575) of the inhomogeneous Poisson point 
#	process model and  estimate the intercept(0.9832485) and covariate "x" (0.531719) of the zero-truncated 
#	Poisson generalized linear model are close to the true values (1, 1, 0.5).   
#
#	Standard errors when detection bias is accounted for(0.1141398; 0.03696048; 0.0224879), and thus 95%  
#	confidence intervals, and standard errors (0.01185; 0.0099313; 0.0060487) are larger when the variability  
#	in the probability of detection is not accounted for. Note the estimated standard errors(0.01185; 0.0099313; 0.0060487)
#	are obtained from steps 5 & 6 above.	
##############################################################################
set.seed(1234)

tpnbs=function()	{
  bss=resample(ZTGLM.data$presence:dim(ZTGLM.data)[ZTGLM.data$presence])
  IPP.data.bss=rbind(IPP.data[bss,1:47],IPP.data[which(IPP.data$presence==0),1:47])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(presence~s3_unclassified_dist+scale(group),family="binomial",data=Detection.data.bss)
  
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(ZTGLM.myFD8))) 
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  
  IPP.model=glm(presence~awc+elev,family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data.bss) 
  ZTGLM.corrected=vglm(group~awc+elev,weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
  c(coef(IPP.model),coef(ZTGLM.corrected))
}
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps) *tpnbs() # Ravi changed the moasic::do and tpnbs()
str(bootstrap.sample)
bootstrap.sample=data.matrix(bootstrap.sample) #Ravi created a matrix from the dataframe
save.image("Hefleyfull.RData")
######## each step gives mean for intercept or covariate valuve 
colMeans(bootstrap.sample)[1] 
sd(bootstrap.sample)[1]
qdata(c=(.025),bootstrap.sample[,1]) # ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[2]
sd(bootstrap.sample)[2]
qdata(c=(.025),bootstrap.sample[,2])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[3]
sd(bootstrap.sample)[3]
qdata(c=(.025),bootstrap.sample[,3])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[4]
sd(bootstrap.sample)[4]
qdata(c=(.025),bootstrap.sample[,4])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[5]
sd(bootstrap.sample)[5]
qdata(c=(.025),bootstrap.sample[,5])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[6]
sd(bootstrap.sample)[6]
qdata(c=(.025),bootstrap.sample[,6])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[7]
colMeans(bootstrap.sample)[8]
colMeans(bootstrap.sample)[9]
colMeans(bootstrap.sample)[10]
colMeans(bootstrap.sample)[11]
colMeans(bootstrap.sample)[12]
colMeans(bootstrap.sample)[13]
colMeans(bootstrap.sample)[14]
colMeans(bootstrap.sample)[15]
colMeans(bootstrap.sample)[16]
colMeans(bootstrap.sample)[17]
colMeans(bootstrap.sample)[18]
