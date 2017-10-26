setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
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

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
hpop <- raster("hpop.tif")
awc <- raster("awc.tif")
tpo <- raster("tpo.tif")
rain <- raster("rain.tif")
nitro <- raster("nitro.tif")
elev <- raster("elev.tif")
temp <- raster("temp.tif")
sbd <- raster("sbd.tif")
roaden <- raster("roaden.tif")
clay <- raster("clay.tif")
fpc <- raster("fpc.tif")
aspect <- raster("aspect.tif")
twi <- raster("twi.tif")
#plot(fpc )
myPredictors.n = c("hpop","awc","tpo", "rain", "nitro", "elev", "temp", "sbd","roaden","clay","fpc","aspect", "twi")
#myPredictors.na = paste(myPredictors.n,".TIF", sep = "")
myPredictors.na = as.list(paste(myPredictors.n,".TIF", sep = ""))  
myStack.hefely = stack(myPredictors.na)

myFD=as.data.frame(hefleydata.s)
myFD1 = cbind(myFD,extract(myStack.hefely,myFD)) # bind xy coordinates.

#write.csv(myFD, file="myFD.csv", row.names=FALSE) 
# --------------------------------------------------------------------------------------------------------------------
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
d_bridleway <- raster("raster_crop\\distance_bridleway.tif")
d_cycleway <- raster("raster_crop\\distance_cycleway.tif")
d_footway <- raster("raster_crop\\distance_footway.tif")
d_motorwayandlink <- raster("raster_crop\\distance_motorwayandlink.tif")
d_path <- raster("raster_crop\\distance_path.tif")
d_pedestrian <- raster("raster_crop\\distance_pedestrian.tif")
d_primaryandlink <- raster("raster_crop\\distance_primaryandlink.tif")
d_residentil <- raster("raster_crop\\distance_residentil.tif")
d_secondaryandlink <- raster("raster_crop\\distance_secondaryandlink.tif")
d_steps <- raster("raster_crop\\distance_steps.tif")
d_tertiaryandlink <- raster("raster_crop\\distance_tertiaryandlink.tif")
d_trunkandlink <- raster("raster_crop\\distance_trunkandlink.tif")
d_unclassified <- raster("raster_crop\\distance_unclassified.tif")
#---------- 
#windows(); plot(d_bridleway)
myPredictors = c(
  "distance_bridleway","distance_cycleway", "distance_motorwayandlink", "distance_path", "distance_pedestrian", "distance_primaryandlink", "distance_residentil",
  "distance_secondaryandlink", "distance_steps", "distance_tertiaryandlink", "distance_trunkandlink", "distance_unclassified")
# add extensions
myPredictorsASC = paste(myPredictors,".TIF", sep = "")
myPredList = as.list(paste(myPredictors,".TIF", sep = ""))  
myStack = stack(myPredList)
myStack.c <- crop(myStack, fpc )
myFD=as.data.frame(hefleydata.s)
myFD3 = cbind(myFD,extract(myStack.c,myFD))
#write.csv(myFD3, file="myFD3.csv", row.names=FALSE) 

myFD4 = cbind(myFD1,myFD3[,3:14]) #select covariates
myFD6 = cbind(myFD4,hefleydata [,6:7]) # select xy ?covaiates of two stacks joined
myFD6 = na.omit(myFD5)
ZTGLM.myFD6=myFD6[which(myFD6$presence==1),]
ZTGLM.myFD7=myFD6[which(myFD6$presence==0),]
set.seed(123)
ZTGLM.myFD8 <- sample(seq_len(nrow(ZTGLM.myFD7)), size = 1000)#use for ipp.data
ZTGLM.myFD9 <- ZTGLM.myFD7[ZTGLM.myFD8, ] #x.int data frame now
#IWLR data set similar to hefley.
ZTGLM.myFD10=rbind(ZTGLM.myFD6,ZTGLM.myFD9) # detected and montecarlo points hust a wild guess
#select only 1000 observations.
#now take a random sample and assign detected 1 non detected 0.
#set.seed(123)
#sample=ZTGLM.myFD6[sample(nrow(ZTGLM.myFD6), 50), ]
#data radomly partitioned to detected and non detected
set.seed(123)
train_ind <- sample(seq_len(nrow(ZTGLM.myFD6)), size = 50)

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
IPP.ignored=glm(presence~fpc,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
ZTGLM.ignored=vglm(group~fpc,family="pospoisson",data=ZTGLM.data)
###
Detection.model=glm(presence~fpc,family="binomial",data=Detection.data) #length 1619 #3 variables #reduce detection data to 200.

#Step 4 - Estimate the probability of detection for each presence-only location.

p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 variables 
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD8)))
IPP.data$p.det=p.det   # IPP.data number of obserarions=1461

ZTGLM.data$p.det=p.det
##Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det
IPP.corrected=glm(presence~fpc,family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data) #length 1461. 3 variables.
summary(IPP.corrected)


##Step 6 - Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det
ZTGLM.corrected=vglm(group~fpc,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)
###############################################################################
set.seed(1234)
tpnbs=function()	{
  bss=resample(ZTGLM.data$presence:dim(ZTGLM.data)[ZTGLM.data$presence])
  IPP.data.bss=rbind(IPP.data[bss,1:30],IPP.data[which(IPP.data$presence==0),1:30])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(presence~fpc+scale(group),family="binomial",data=Detection.data.bss)
  
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=p.det.bss
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  
  IPP.model=glm(presence~fpc,family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data.bss) 
  ZTGLM.corrected=vglm(group~fpc,weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
  c(coef(IPP.model)[2],coef(ZTGLM.corrected))
}
n.bootstraps=500
bootstrap.sample = mosaic::do(n.bootstraps) *tpnbs() # Ravi changed the moasic::do and tpnbs()
str(bootstrap.sample)
bootstrap.sample=data.matrix(bootstrap.sample) #Ravi created a matrix from the dataframe
save.image("Hefley.RData")
######## each step gives mean for intercept or covariate valuve 
colMeans(bootstrap.sample)[1] 
sd(bootstrap.sample)[1]
qdata(c=(.025),bootstrap.sample[,1]) # ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[2]
sd(bootstrap.sample)[2]
qdata(c=(.025),bootstrap.sample[,2])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[3]
sd(bootstrap.sample)[3]
qdata(c=(.025),bootstrap.sample2[,3])#ravi:removed after .975 after .025, add c=
#######
