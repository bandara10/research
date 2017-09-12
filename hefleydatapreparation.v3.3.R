##### load libraries ####
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
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

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

##### Step 6: now take a random sample of 80 and assign detected 1 non detected 0.####
train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = 80,replace=FALSE)
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0
##### Step 7: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected) 
Detection.data[53] <- lapply(Detection.data[53], as.numeric)
str(Detection.data)
Detection.data <- Detection.data[c(-1,-2)]
Detection.data <- subset(Detection.data, select=c(50, 1:49, 51))
IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus no koalas data from myFD1.
IPP.data <- IPP.data[c(-1,-2,-53)]
IPP.data <- subset(IPP.data, select=c(50, 1:49))
ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.
ZTGLM.data <- ZTGLM.data[c(-1,-2,-52)]
ZTGLM.data <- subset(ZTGLM.data, select=c(50, 1:49))
##### Step 8:  analysis without detection correction factor#######
IPP.ignored=glm(presence~twi + tpo + temp + aspect + elev+habit2pc+hpop+lot_density+sbd,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
summary(IPP.ignored)
ZTGLM.ignored=vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,family="pospoisson", data=ZTGLM.data)
summary(ZTGLM.ignored)

#############  .
set.seed(123) # we create detection probabilities using two methods. glm, rf
#Detection model: steps as in Hefley`s code`
Detection.model=glm(presence~  distance_pedestrian + s1_residential_dist + distance_trunkandlink +
                      distance_tertiaryandlink+ scale(group),family= "binomial", data=Detection.data)
#### Random  forest detection model
model <- presence~  distance_pedestrian + s1_residential_dist + distance_trunkandlink +
  distance_tertiaryandlink
rf1 <- randomForest(model, data=Detection.data)
plot(rf1)
str(rf1)
pr <- predict(myfullstack, rf1)
#Random floor method
 # Now use random floor for modelling

y= Detection.data$presence
X= Detection.data
X = X[,!names(X)=="presence"]
# regression doesnt work. lets fix this by telling rf to do classification instead of regression as we have fewer valuves.
# this to be done for the detection model as well.
Detection.data$presence <- as.character(Detection.data$presence)
Detection.data$presence <- as.factor(Detection.data$presence)

rf.default = randomForest::randomForest(X,y,ntree=50)
rf.robust  = randomForest:: randomForest(X,y,sampsize=25,ntree=5000, mtry=4, keep.inbag = y,keep.forest = T) 

#Roc curves doesnt work for this data set. because message; are you sure want to do regresison. 

plot(AUC::roc(rf.default$votes[,2],y) ,main="ROC: default black, robust is red")
plot(AUC::roc(rf.robust$votes[,2],y),col=2,add = y)
plot.new()
print(AUC:: auc(AUC::roc(rf.default$votes[,2],y)))
print(AUC:: auc(AUC::roc(rf.robust$votes[,2],y)))

#computing feature contributions and visualization  
  
ff = forestFloor(rf.robust,X,binary_reg = T,calc_np=T)
Col = fcol(ff,cols=1,outlier.lim = 2.5)
plot(ff,col=Col,plot_GOF = T)  
#A clear interaction effect can be observed  : troubeshoot"# reinstall rgdl and run rgl.close() to close all and then rerun above show3d from randomforestfloor.
#Fitted surface(grey) estimate explained variance. 
forestFloor::show3d(ff,c(1,5),5,col=Col,plot_GOF = T)


#####

 ### to evalvuate the model create presence location and absence location datasets from hefleydata.s.
hefleydata.presence <-subset(hefleydata, presence==1)
# model evaluvation using all presence records and absence records 
absence.selected <- ZTGLM.myFD2[ZTGLM.myFD3, ]
(rf1.evauvate <-evaluate(ZTGLM.myFD1, absence.selected, rf1))   
plot(pr, main= "Detection probabilities- random forest")
tr <- threshold(rf1.evauvate, "spec_sens") #  set a threshold valuve fordetected and non detected. 
plot(pr > tr, main= "presence/absenec- random")

### random forest detection probabilities for next steps
p.det.rf = faraway::ilogit(predict(rf1,new=ZTGLM.data))# 
hist(p.det.rf, breaks= 100)
#unclass(summary(Detection.model))



#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det = faraway::ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
hist(p.det, breaks=100)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
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
ZTGLM.corrected = vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd
                       ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)
ZTGLM.corrected
# step 7:  Map predictions
myPred = predict(myfullstack, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")
plot(hefleydata.presence, add=TRUE)
myPred2 = predict(myfullstack, IPP.corrected, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model-intensity of group or ??,")
plot(hefleydata.presence, add=TRUE)
myPred3 = predict(myfullstack, ZTGLM.corrected, type = "response")
plot(myPred3,  main="ZTGLM-Number of koalas in a grid - VGLM ")
plot(hefleydata.presence, add=TRUE)
#writeRaster(myPred3, "ZTGLM.tif")
dev.off()
####

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

####### step 8: Two phase algoritham#######################################################################
set.seed(1234)
tpnbs=function()	{
  bss=resample(1:dim(ZTGLM.data)[1])
  IPP.data.bss=rbind(IPP.data[bss,1:53],IPP.data[which(IPP.data$presence==0),1:53])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                        distance_tertiaryandlink+scale(group),family="binomial",data=Detection.data.bss)
  
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(ZTGLM.myFD3))) 
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  IPP.model= glm(presence~ twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd, 
                 family = "binomial", weights = (1/p.det)*10000^(1-presence), data = IPP.data.bss)
  
  ZTGLM.corrected=vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,
                       weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
  c(coef(IPP.model),coef(ZTGLM.corrected))
}

#### step 9: boostrapping ####
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps)*tpnbs() # Ravi changed the moasic::do and tpnbs()
bootstrap.sample <- data.frame(bootstrap.sample)
#bootstrap.sample=as.data.frame (bootstrap.sample) #Ravi created a matrix from the dataframe
######## each step gives mean of coefficients for intercept or covariates

#### step 10:  Get means and sd ####
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
colMeans(bootstrap.sample)[19]
colMeans(bootstrap.sample)[20]
colMeans(bootstrap.sample)[21]
colMeans(bootstrap.sample)[22]
colMeans(bootstrap.sample)[23]
colMeans(bootstrap.sample)[24]
colMeans(bootstrap.sample)[25]
colMeans(bootstrap.sample)[26]

 
#change resolution. by 4. get points in each cell. treat them as groups and size.
#<- disaggregate(meuse.raster, fact=4)

