#test <- readShapeSpatial("1994_2014fullstudyarea.shp")
# read shape file fullgrid1997-2014
#create a dataf frame
# create a colomn and store sum
#create a colomn divide sum by sum to get presence
#test <- as(test, "data.frame")
#test$group <- rowSums(test[18:22])
#test <- mutate(test, presence=group/group)
#hefleydata <- test[,26:29]
#hefleydata[is.na(hefleydata)] <- 0
#write.csv(hefleydata, file="hefley_fishnet_rastermatch2010_2014.csv")
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
library(dplyr)
library(usdm)
library(ROCR)
library(maptools)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\renamed_vars")
# variables renamed and selected based on vif.
#Datapreparation:raster; data not trasformed. Koala data for from 2010-2015.
# spatialEco pachage can trasformdata to norm", "rstd", "std", "stretch", "nl", "slog", "sr".
# use #vif for varibale selection.  step function to select the best model. Then fit models ipp ignored and corrected models.
#If you want to formally test the hypothesis of Complete Spatial Randomness you can do this using the chi-squared 
# chi-squared test based on quadrat counts (quadrat.test) from spatstat package.

myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack <- stack(myfullstack.a)
#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####
hefleydata <- read.csv("hefley_fishnet_rastermatch2010_2014.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# plot koala presence locations
hefleydata.presence <-subset(hefleydata, presence==1)
coordinates(hefleydata.presence) <- ~x+y
plot(hefleydata.presence)

# get coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]
hefleydata.s = as.data.frame(hefleydata.s)
#plot(hefleydata.s)
#select presence records
#### Step 3:  extract X=vector.boot from raster anad combine wth hefleydata, presence and group varibels.####
myFD = cbind(hefleydata.s,raster::extract(myfullstack,hefleydata.s))
# get presence and group vriables
myFD1 = cbind(myFD,hefleydata [,3:4]) # get presence and group data. 
myFD1 = na.omit(myFD1) # remove all NA valuves

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 
#####Step 5: select only 1000 absences (monticarlo points as in hefleys method??)####
set.seed(12356)
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 
ZTGLM <- ZTGLM.myFD5
##### Step 6: now take a random sample of n% and assign detected 1 non detected 0.####

set.seed(12345)
train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = floor(0.80 * nrow(ZTGLM.myFD1)),replace=FALSE) # split the presence only dataset as detected and nondetected.
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0
##### Step 7: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected)  # This dataset length is same as presence dataset.
IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus n=1000 deteected koalas data n=80 from ZTGLM.myFD4
ZTGLM.data=(detected)##ZTGLM.data# This is detected data randomly selected n= 50% .

# Selection of explanatory varibaels based on VIF#
vifstep (myfullstack, th=10) # select variables which have Varience nflation Factor less than 10.

##### Step 8:  analysis without detection correction factor####### mode selected from below method:
IPP.ignored = glm(presence ~ distance_secondaryandlink,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
summary(IPP.ignored)
ZTGLM.ignored=vglm(group ~ rain_mean_annual+habit3percent,family="pospoisson", data=ZTGLM.data)
summary(ZTGLM.ignored)
#unclass(summary(Detection.model))

##### perform forward selection by specifying a starting model and the range of models which we want to examine in the search.#### 
null=glm(presence~ 1, family= "binomial", data=Detection.data) # null model
full=glm(presence ~ city_dis+dis_habit1+dis_habit3+
           Dis_habitat_suitable_2+Dis_habitat_suitable_3+distance_motorwayandlink+distance_primaryandlink
         +distance_residentil+distance_secondaryandlink+
           distance_tertiaryandlink+visitorarea_dis+ scale(group),
         family= "binomial", data=Detection.data) # full set of explanatory varibales.
# perform forward selection using the command step
step(null, scope=list(lower=null, upper=full), direction="forward")
### model slection in this way did not give valid results in predictions and maps. hwo to use selected variables.
###
#Variable selection VGLM model
null=glm(presence~ 1, family= "binomial", data=Detection.data)
full=glm (presence~DP_QLD_FPC20141 +elevation+habit1percent +habit2percent+habit3percent+hpop_2015+lot_density+NTO_000_005+PTO_000_005+
  rain_mean_annual+SAWC_000_005+Slope_percent_61+temp_maximum1990_2000+temp_minimum1996_2005+TopoWI+
    scale(group),family= "binomial", data=Detection.data)
step(null, scope=list(lower=null, upper=full), direction="forward")
####   #####
logistic
#Hefley method GLM
set.seed(1238) # we create detection probabilities using two methods. glm, rf
#Detection model: steps as in Hefley`s code`
Detection.model=glm(presence ~ city_dis+dis_tertiaryandlink+s3_residential_dist+visitorarea_dis, family= "binomial", data=Detection.data)

summary(Detection.model)
# check the prediction map right here. #use only the selected variabels as
myPred1 = predict(myfullstack, Detection.model, type = "response")

plot(myPred1, xlab = "x", ylab= "y",main="detection model")
plot(hefleydata.presence,pch=1, add=TRUE)
###### #######
myPred = prediction(predict(Detection.model, type = "response"), Detection.data$presence)
perf <- performance(myPred,measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)


##### Stores the residuals  plot corellagram ####
Detection.data$res = residuals(Detection.model) # library(ncf)
myResCorr <- correlog(Detection.data$x, Detection.data$y, Detection.data$res,na.rm=T, increment=1000, resamp=0, latlon = F)
plot(myResCorr$mean.of.class[1:100], myResCorr$correlation[1:100] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")

#rerun the modelwith residuals: looks awfull. inst it. use crase approach for auto corelation. have to refer his paper for this. 
Detection.model.1=glm(presence~  distance_pedestrian + s1_residential_dist + distance_trunkandlink +
                        distance_tertiaryandlink + res + scale(group),family= "binomial", data=Detection.data)

summary(Detection.model.1) # not a good model.
##### #####
#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det = faraway::ilogit(predict(Detection.model, type="response", new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
#myBRT,n.trees=myBRT$gbm.call$best.trees ; add this as model in the above function.

hist(p.det, breaks=100)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
ZTGLM.data$p.det=p.det

######Step 5: - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####

IPP.corrected= glm(presence~twi + tpo + temp + aspect + dem+habit2pc+hpop+lot_density+sbd,
                   family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data)

summary(IPP.corrected)
# broom package: used tidy to get a table from model outputs. This doesnt work for VGLMs.
# get confidence intervals and compare ignored and corrected models.
confidenceintervals <- confint(IPP.ignored)
tidy(IPP.ignored,confidenceintervals)
summary(IPP.ignored)

#visualise response vs covariates eg. lot density.
histogram(IPP.data$lot_density ,fitted(IPP.corrected),xlab='lot',
          ylab='Occurrence', add = TRUE, col="red")
####Step 6: Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####

#use only the significant covariates, tpo +hpop+lot_density+sbd
# VGAM: read about which family to use: https://www.r-project.org/doc/Rnews/Rnews_2008-2.pdf
ZTGLM.corrected = vglm(group~twi + tpo + temp + aspect + elev+habit2pc+hpop+lot_density+sbd
                       ,weights=1/p.det,family="pospoisson",data=ZTGLM.data) # zapoisson

summary(ZTGLM.corrected)
# step 7:  Map predictions
myPred1 = predict(myfullstack, Detection.model, type = "response")
plot(myPred1, xlab = "x", ylab= "y",main="detection model")
plot(hefleydata.presence, pch=1,add=TRUE)
myPred2 = predict(myfullstack, IPP.corrected, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model-intensity of group")
plot(hefleydata.presence,pch=1, add=TRUE)
myPred3.1 = predict(myfullstack, ZTGLM.corrected, type = "response")
plot(myPred3.1,  main="ZTGLM-Number of koalas in a grid - VGLM ")
plot(hefleydata.presence,pch=1, add=TRUE)
#writeRaster(myPred3, "ZTGLM.tif")
dev.off()
####
unlist(ZTGLM.corrected@predictors)> .02
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
  IPP.data.bss=rbind(IPP.data[bss,1:54],IPP.data[which(IPP.data$presence==0),1:54])
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


######### discussion points ########
# koala detection map: we corrected for detection erros. People report koalas from places where they went and found koalas.
#This doesnt mean that it is a random sample from where koalas are present and people detected some and reported some. 
# Detected kaoals are from a biased sample as people visit only limited areas. So our detection bias correction is 
# should be interpreted as correcting detection bias in places where detection was was carried out. basically, we corrected
# for detection erros in detected areas but not in undetected areas our model resilved this. IF we use more spatially spread 
# dataset with fw years of data then this model may behave well.

#####
# In the field what could happe
#1. koala is found in suitabel habitat only. This iclude biotic and abitic enviroment is suitable.
#   However, but not in all suitable habitats due to unsuitable abiotic enviroment. orcompetitio or threat.
#2. habitat suitable - koala have acccess to suitabel habitat- koala is present-  people go to that habitat- people see koala- some poeple report koala. It is a detection. but can invlove an error. as detection is not always prefect.
#                                                                                                          - people do not see koala due to detection error/lack of experience. or due to influnce of covariates like forest type, time of the day. 
#                                                        - koala present but moved to an another location when Pople go to that location -  So do not see koalas because it is not there. if people observe this location repetedly thye might find koalas once they return back.                                                                             
#                                                        - koala present - people do not go to that habitat or no access conserved areas- no koala sighting records.
#3. habitat is suitable- but koalas do not have access to that suitable habitat due to phisical barriers. So Koalas are not present. People go to suitable habitat  areas and do not see koals as koals are not present.
#4. habitat is not suitable- koalas are not present. 

#Key features of koala habitat or threshhold levels: eg. Temperature, rainfall, type of forest. etc.

# Detection model/ GLM model: Probability of detecting koalas if present and people have equall access to all areas. 
# IPP model/ GLM model= Intensity of koala groups in each grid cell. 
# and they have a size. These are number of groups and they have a size. Better approach would be to increase number of 
#sightings adding more sightings from other years. Then instead of 1000 k grids, create smaller grids say 1/4 th of 1000k, count the number of kolas in each small grid and consider them as groups.
#Then ZTGLM model which will be modelled at 1000k grid will give us the total count of koalas @ 1000k grid. 
# ZTGLM model: or VGLM model :   Total number of koalas in each grid. Groups have sizes and this is the total of each group size.
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1                   aspect91 1.019162
# 2                BDW_000_005 3.454197
# 3                   city_dis 2.914351
# 4                CLY_000_005 4.508156
# 5                den_primary 1.171345
# 6            den_residential 4.237767
# 7          den_secondry_link 1.183622
# 8         den_tertiary_liink 1.217534
# 9                 dis_habit1 2.192523
# 10                dis_habit3 5.726589
# 11    Dis_habitat_suitable_1 2.285179
# 12    Dis_habitat_suitable_2 3.119045
# 13    Dis_habitat_suitable_3 5.706824
# 14  distance_motorwayandlink 5.998391
# 15   distance_primaryandlink 2.415658
# 16       distance_residentil 2.569289
# 17 distance_secondaryandlink 2.855088
# 18  distance_tertiaryandlink 2.500948
# 19           DP_QLD_FPC20141 1.675547
# 20                 elevation 6.074694
# 21             habit1percent 1.325542
# 22             habit2percent 1.340830
# 23             habit3percent 1.182947
# 24                 hpop_2015 4.353555
# 25               lot_density 6.383154
# 26               NTO_000_005 3.538900
# 27               PTO_000_005 3.225335
# 28          rain_mean_annual 7.864017
# 29              SAWC_000_005 4.291834
# 30          Slope_percent_61 2.163289
# 31     temp_maximum1990_2000 2.722183
# 32     temp_minimum1996_2005 8.656966
# 33                    TopoWI 1.553572
# 34           visitorarea_dis 4.193013