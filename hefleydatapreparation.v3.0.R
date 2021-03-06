##### load libraries ####
library(spatial.tools)
library(VGAM)
library(mosaic)
library(spatstat)
library(faraway)
library(raster)
library(dplyr)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

myfullstack.a <- list.files(pattern="\\.tif") 
myfullstack = stack(myfullstack.a)
#check the coliniarity of variables 
#vif(myfullstack)# check for coliniarity
#vifstep(myfullstack, th=10) # select variables which have VIF less than 10.
#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####
hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# get only coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]
#### Step 3:  extract variables from raster anad combine wth hefleydata, presence and group varibels.####
myFD = cbind(hefleydata.s,raster::extract(myfullstack,hefleydata.s))
# get presence and group vriables
myFD1 = cbind(myFD,hefleydata [,6:7]) # get presence and group data. 
myFD1 = na.omit(myFD1) # remove all NA valuves
#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 
#####Step 5:slect only 1000 absences (monticarlo points as in hefleys method??)####
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 
##### Step 6: now take a random sample of 80 and assign detected 1 non detected 0.####
train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = 80,replace=FALSE)
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence<- 0
##### Step 7: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected) 
#glimpse(Detection.data, n=10)
IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus no koalas data from myFD1.
ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.
##### Step 8:  analysis without detection correction factor#######
IPP.ignored=glm(presence~twi+tpo+temp +aspect+awc+clay+elev+fpc+habit2pc+hpop+
                  lot_density+nitro+roadk+sbd+suitable_1+suitable_3,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
summary(IPP.ignored)
ZTGLM.ignored=vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,family="pospoisson", data=ZTGLM.data)
summary(ZTGLM.ignored)

#############  detection model and select the best model.
#Detection model: steps as in Hefley`s code`
set.seed(123)
##use step function here then use significant variable in the next model. I asume this is correct way to do it.
#So the correct model is the second model in this step. or elase run the line 62.
Detection.model.1=step(glm(presence~Dis_habitat_suitable_1+Dis_habitat_suitable_2+Dis_habitat_suitable_3+
                             distance_bridleway+distance_motorwayandlink+distance_path+distance_pedestrian+
                             distance_primaryandlink+distance_residentil+distance_secondaryandlink+distance_tertiaryandlink+
                             distance_trunkandlink+ distance_unclassified+s1_residential_dist+s1_unclassified_dist+s2_residential_dist+
                             s2_unclassified_dist+s3_residential_dist +scale(group), family= "binomial",data=Detection.data))  #length 1619 #3 variables #reduce detection data to 200.
summary(Detection.model.1)

# Significant covariates are 1.distance_pedestrian, 2.s1_residential_dist, 3.
# distance_trunkandlink, 4.distance_tertiaryandlink
Detection.model=glm(presence~ distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                      distance_tertiaryandlink+scale(group), family= "binomial",data=Detection.data)

summary(Detection.model)
#confint(Detection.model); cov2cor(vcov(Detection.model))
#####Step 4 - Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 variables 
hist(p.det, breaks=70)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det
######Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected.1= step(glm(presence~twi+tpo+temp +aspect+awc+clay+elev+fpc+habit2pc+hpop+
                            lot_density+nitro+roadk+sbd+suitable_1+suitable_3,
                          family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data)) 
summary(IPP.corrected.1)
# use significant covariates in the model. They are twi,tpo,temp,aspect,elev,habit2pc,hpop,lot_density,sbd
IPP.corrected= glm(presence~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,
                   family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data)
summary(IPP.corrected)
####Step 6 - Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####
ZTGLM.corrected.1 = vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd
                         ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected.1)
#use only the significant covariates, tpo +hpop+lot_density+sbd
ZTGLM.corrected = vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd
                       ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)

summary(ZTGLM.corrected)
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
####### step:7 Two phase algoritham#######################################################################



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
#### step:8 boostrapping ####
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps) *tpnbs() # Ravi changed the moasic::do and tpnbs()

bootstrap.sample=as.data.frame (bootstrap.sample) #Ravi created a matrix from the dataframe

means= as.matrix(apply(bootstrap.sample,2, mean))
sds= as.matrix(apply(bootstrap.sample,2, sd))
means.sds= cbind(means, sds )
####
#for (i in 1:ncol(bootstrap.sample)) {

hist(bootstrap.sample[,i], breaks=50)
}

###

######## each step gives mean of coefficients for intercept or covariates

#### step"9 Get means and sd ####
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

#### step:10 Map predictions# the above algoritham give colom means, How to pam
#### that prediction becuse ut is not directly from a model####
##Below  map is not from above algo but from the models run before that.
myPred = predict(myfullstack, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")

myPred2 = predict(myfullstack, IPP.corrected, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model-intensity of group or ??")
myPred3 = predict(myfullstack, ZTGLM.corrected , type = "response")
plot(myPred3,main="ZTGLM-Number of koalas in a grid model") 

