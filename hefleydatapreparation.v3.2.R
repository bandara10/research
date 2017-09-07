##### load libraries ####
library(spatial.tools)
library(VGAM)
library(mosaic)
library(spatstat)
library(faraway)
library(raster)
#library(dplyr)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

myfullstack.a <- list.files(pattern="\\.tif") 
myfullstack = stack(myfullstack.a)

#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####

hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# get only coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]

#### Step 3:  extract X=vector.boot from raster anad combine wth hefleydata, presence and group varibels.####
myFD3 = cbind(hefleydata.s,raster::extract(myfullstack,hefleydata.s))
# get presence and group vriables
myFD6 = cbind(myFD3,hefleydata [,6:7]) # get presence and group data. 
myFD6 = na.omit(myFD6) # remove all NA valuves

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD6=myFD6[which(myFD6$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD7=myFD6[which(myFD6$presence==0),] # select all absence data 

#####Step 6:slect only 1000 absences (monticarlo points as in hefleys method??)####
ZTGLM.myFD8 <- sample(seq_len(nrow(ZTGLM.myFD7)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD9 <- ZTGLM.myFD7[ZTGLM.myFD8, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD10=rbind(ZTGLM.myFD6,ZTGLM.myFD9) 

##### Step 7: now take a random sample of 80 and assign detected 1 non detected 0.####
train <- sample(seq_len(nrow(ZTGLM.myFD6)), size = 80,replace=FALSE)
detected <- ZTGLM.myFD6[train, ]
notdetected <- ZTGLM.myFD6[-train,] 
#not detected assigned valuve 0
notdetected$presence<- 0

##### Step 8: Create the final data sets for the analysis####
Detection.data= rbind(detected,notdetected) 
#glimpse(Detection.data, n=10)
IPP.data=rbind(detected, ZTGLM.myFD9) #IPP.data comes from detected data plus no koalas data from myFD6.
ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.

##### Step 9:  analysis without detection correction factor#######
IPP.ignored=glm(presence~s3_unclassified_dist,family="binomial",weights=10000^(1-presence),data=IPP.data) # IPP.data2 added.
summary(IPP.ignored)
ZTGLM.ignored=vglm(group~awc+elev,family="pospoisson", data=ZTGLM.data)
summary(ZTGLM.ignored)

#############  detection model and select the best model.
set.seed(123)
##use step function here then use significant variable in the next model. I asume this is correct way to do it.
#So the correct model is the second model in this step. or elase run the line 62.
#Detection.model.1=step(glm(presence~Dis_habitat_suitable_1+Dis_habitat_suitable_2+Dis_habitat_suitable_3+
#                          distance_bridleway+distance_motorwayandlink+distance_path+distance_pedestrian+
#                          distance_primaryandlink+distance_residentil+distance_secondaryandlink+distance_tertiaryandlink+
#                          distance_trunkandlink+ distance_unclassified+s1_residential_dist+s1_unclassified_dist+s2_residential_dist+
#                          s2_unclassified_dist+s3_residential_dist +scale(group), family= "binomial",data=Detection.data))  #length 1619 #3 X=vector.boot #reduce detection data to 200.
#summary(Detection.model.1)

# Significant covariates are 1.distance_pedestrian, 2.s1_residential_dist, 3.
# distance_trunkandlink, 4.distance_tertiaryandlink
Detection.model=glm(presence~ distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                      distance_tertiaryandlink+scale(group), family= "binomial",data=Detection.data)

summary(Detection.model)
#confint(Detection.model); cov2cor(vcov(Detection.model))

#####Step 4 - Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
hist(p.det, breaks=70)
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD8)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det

######Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
#IPP.corrected.1= step(glm(presence~twi+tpo+temp +aspect+awc+clay+elev+fpc+habit2pc+hpop+
#                    lot_density+nitro+roadk+sbd+suitable_1+suitable_3,
#                  family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data)) 
#summary(IPP.corrected.1)
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


####
set.seed(1234)
tpnbs=function()	{
  bss=resample(1:dim(ZTGLM.data)[1])
  IPP.data.bss=rbind(IPP.data[bss,1:53],IPP.data[which(IPP.data$presence==0),1:53])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                        distance_tertiaryandlink+scale(group),family="binomial",data=Detection.data.bss)
  
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(ZTGLM.myFD8))) 
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

myPred2 = predict(myfullstack, IPP.ignored, type = "response")
plot(myPred2, col = rainbow(100), xlab = "x", ylab= "y",main=" IPP model-intensity of group or ??,")
myPred3 = predict(myfullstack, ZTGLM.corrected, type = "response")
plot(myPred3,main="ZTGLM-Number of koalas in a grid model") 




-----------------------
  #means= as.matrix(apply(bootstrap.sample,2, mean))
  #sds= as.matrix(apply(bootstrap.sample,2, sd))
  #means.sds= cbind(means, sds )
  ####
  #
  x=bootstrap.sample 
for (i in 1:ncol()) {
  
  hist(x[,i], breaks=50)
}

bsr <- function(dat=IPP.data, ) {
  
}
#####http://www.davidbroadstock.org/bootstrap
#https://stats.stackexchange.com/questions/63652/how-does-bootstrapping-in-r-actually-work
#step : 1 Data based resampling
bs <- function(formula,data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- glm(formula, data=d) 
  return(coef(fit)) 
} 
bootstrap_regression<-boot(data=Detection.data, statistic=bs, 
                           R=999, formula=presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                             distance_tertiaryandlink+scale(group))
#step: 2 Residual based resampling
fit <- fitted(lm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                   distance_tertiaryandlink+scale(group) , data=Detection.data))
e <- residuals(lm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                    distance_tertiaryandlink+scale(group)  , data=Detection.data))
X <- model.matrix(lm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                       distance_tertiaryandlink+scale(group)  , data=Detection.data))

bs <- function(data, indices) {
  y.star <- fit + e[indices]
  mod <- lm(y.star ~ X -1)
  coefficients(mod)
} 

bootstrap_regression_residual<-boot(data=Detection.data, statistic=bs, 
                                    R=999)
bootstrap_regression_residual = predict(myfullstack, bootstrap_regression_residual, type = "response")
plot(bootstrap_regression_residual)
#Step 3: Plot the empirical distribution functions

layout(rbind(cbind(1,2),cbind(3,4)))

hist(bootstrap_regression$t[,1],main="Intercept - data based resampling",
     xlab="(Dashed line: mean of bootstrap replications)", ylab="")

abline(v=mean(bootstrap_regression$t[,1]), lwd=3,lty=2)

hist(bootstrap_regression$t[,2],main="Slope - data based resampling",
     xlab="(Dashed line: mean of bootstrap replications)", ylab="")

abline(v=mean(bootstrap_regression$t[,2]), lwd=3,lty=2)

hist(bootstrap_regression_residual$t[,1],main="Intercept - residual based resampling",
     xlab="(Dashed line: mean of bootstrap replications)", ylab="")

abline(v=mean(bootstrap_regression_residual$t[,1]), lwd=3,lty=2)

hist(bootstrap_regression_residual$t[,2],main="Slope - residual based resampling",
     xlab="(Dashed line: mean of bootstrap replications)", ylab="")

abline(v=mean(bootstrap_regression_residual$t[,2]), lwd=3,lty=2)


#Step 4: Extract the bootstrap confidence intervals

print(boot.ci(bootstrap_regression, type="bca", index=1))
print(boot.ci(bootstrap_regression, type="bca", index=2))

print(boot.ci(bootstrap_regression_residual, type="bca", index=1))
print(boot.ci(bootstrap_regression_residual, type="bca", index=2))




#####
observe.lm <- glm(presence~awc+hpop, data=IPP.data)
vif(observ.lm)
confidenceEllipse(observ.lm)
N=500 # may not work for some catogorical varibales
bootstrapRandom<- function(dat=IPP.data, mod.formula= formula(presence~awc+hpop)) {
  dat.boot<-  dat[sample(x = NROW(dat), size = NROW(dat), replace = T),]
  boot.lm<- lm(mod.formula, dat=dat.boot)
  coef(boot.lm)
}
vector.boot<- t(replicate(N, bootstrapRandom()))
#standarad erros
apply(vector.boot, MARGIN = 2, sd)
#precentile CIs (trasposed so it is oriented like output of confint)
t(apply(vector.boot, MARGIN = 2, quantile, probs= c(0.025, 0.975)))
#symptomatic normal CIs
confint(observ.lm)
# covarience among parameters
cov2cor(vcov(observ.lm))
par(mfrow=c(2,2))
multiplehisto <- function(X=vector.boot) {
  for (i in 1:NCOL(X)){
    hist(X[,i], freq = F, 
         main = colnames(X)[i],
         xlab=colnames(X)[i])
  }
  
}


multiplehisto()
#look at corelation among parameter estimates frm bootstrap
pairs(vector.boot)
# The residula /fixed efffect bootstrap method
# extract residuals from observe.lm model.
resid.model.1<- resid(observe.lm)
plot(density(resid.model.1, bw=0.5))
#now check 
par(mfrow= c(1,2))
plot(resid.model.1 ~ IPP.data$awc)
plot(resid.model.1 ~ IPP.data$hpop)

# write a function to perfomr bootstrapping
bootstrapFromResidulas<- function(mod.object= observe.lm, dat=IPP.data) {
  resids= mod.object$resid # extract residulas from model
  fittedValues=mod.object$fitted #extract fitted valuves.
  matr <- model.matrix(mod.object)
  # generating new values for eeach y[i], by adding bootstraped residulas to fitted values
  Y <- fittedValues+sample(resids, length(resids), replace = T)
  #using model matrix for the predictors
  model.boot<- glm(Y~ 0 + matr, dat=IPP.data) # refit the model with new Y values
  coef(model.boot) # extract the coefficients.
}
residual.boot.N <-t(replicate(N, bootstrapFromResidulas()))
par(mfrow=c(2,2))
multiplehisto(X=residual.boot.N)
pairs(residual.boot.N)
t(apply(residual.boot.N, MARGIN = 2, quantile, probs=c(0.025, 0.975)))
confint(observ.lm) sontompe.
# for comparison , generate confidence intervals using Monte carlo simulations(parametric bootstrap)
SimulationUnderModel <- function(model=observe.lm){
  #extract design matrix
  matr <- model.matrix(model)
  rse <- summary(model)$sigma
  df=model$df
  
  # incorporate uncertanity in RSE
  rse.sim <- rse* sqrt ( df/rchisq(1, df = df))
  # simulate data (response) conditional on the simulated RSE.
  y.sim <- rnorm(n = nrow(matr),
                 mean = matr%*% coef(model), sd, rse.sim)
  #0+ design matrix (since the intercept is already in the design matrix)
  lm.sim <- glm(y.sim~0+ matr) # fit the model with simulated response
  coef(lm.sim)}
sim.coef <- t(replicate(N, SimulationUnderModel()))
t(apply(sim.coef, MARGIN = 2, quantile, probs=c(0.25, 0.975)))
