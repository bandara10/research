library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(lmtest);library(spatial.tools);library(VGAM);library(mosaic);library(faraway);library(gstat)  #
library(ncf);library(foreign);library(nlme)   ;library(MASS);library(ROCR);library(vcd)
library(RColorBrewer);library(classInt);library(ppmlasso);library(usdm) ; library(ncf); library(epicalc)

##############

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")

##########Step 1: load koala data frm tow sources. BoalaBASE and Wildnet#####
## step 1 and two are common to all methods.

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata=mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
all.locations= subset(mydata,yearnew >1998:2015, select=X:Y)

# Select records based on distance

source("Lib_DistEstimatesToItself.r")## set a minimum distance betweek koalas
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 200)
selected.locations = select.locations[,1:2]

# Get wildnet data

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata.RData") # wildnet data.

selected.locations=rbind(selected.locations, wildnetdata)

######Step 2: load rasters and crop to my extent of interest######

myfullstack.a <- list.files(pattern="\\.tif$",full.names=TRUE)
myfullstack.a = scale(stack(myfullstack.a))

#select my extent / myextent=drawExtent()

myextent=c(387900, 553100, 6862400, 7113600) 
myfullstack=crop(myfullstack.a, myextent, snap="near") ###crop to my area of interest
habitat.r<- subset(myfullstack, c(1,2,4,6,13,16,23, 24, 25,26,27,28,29,39,40,41)) # select my variables

#select koalas from grids

selected.locations=gridSample(selected.locations, habitat.r, n = 100)
loc=SpatialPoints((selected.locations))
plot(loc, add=TRUE)

#plot the stack and data for visualization.
points=SpatialPoints(select.locations)
fun <- function() {
  plot(LGA.shp, add = TRUE, col = "red", pch = 1)
}

#plot rasters only##

plot(habitat.r,  nc = 4, nr =3)

# plot raster and overlay koala locations

plot(habitat.r,addfun = fun,  nc =4, nr =3)
#--------------------------------------------------------------------------------------
# ######## All climate rasters .

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = scale(stack(myfullstack.b))
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.b, myextent, snap="near")

habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,33,39)) # select my variables

#------------------------------------------------------------------------------------------

######Step 3. Create quadrature points and species xy data. Remove NA. Keep xy colomns for preditions and rasterize. ######

big.grp <- as.data.frame(habitat.rr, xy=TRUE) # to be used later in Hefley method.

#

bigquad <- as.data.frame(habitat.rr, xy=TRUE,na.rm=TRUE)
xydata <- bigquad[c(1,2)] ### use to get predictions and rasterise.

# ####optional step to create a new varibale. fpc buffer area. 2 km buffer. Then make xy in to km.######
# ##fpc buffer variable.
# bigquadxy=bigquad[c(1,2)]
# fpcnew=myfullstack$fpcnew
# buff=extract(fpcnew,bigquadxy, fun=max,buffer=2,df=TRUE)
# #rename this vvariable as buffer
# colnames(buff)[2] <- "fpcnew_buff"
# buff=buff[c(2)]
# bigquad.t=cbind(buff,bigquad)
# bigquad=bigquad.t[c(2:19,1)]

bigquad =cbind((bigquad[,1:2]/(1000)),bigquad[c(-1,-2)]) # xy convert to km and integers. decimals take days to run..
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'  
xydatan <- as.data.frame(lapply(bigquad[c(1,2)], as.integer)) 
dd=bigquad[c(-1,-2)]
bigquad <- cbind(xydatan,dd )   # Quadrature points XY as integers.

##### Step 4:Create quadrature points with values set to minimum of distance variables 
#             use for bias corretion model
#             model based control of observer bias based on same level of distance

bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 
bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)

####### Step 5.Koala XY data.


#select locations over study area as  dataset is for the full area..

testlayer <- habitat.rr[[1]]   # get a raster

selected.locations.t=cbind(selected.locations,(extract(testlayer,selected.locations))) # locations are NA for outside study area.

# only data from study area. Remove NA rows and keep xy.

selected.locations<- as.data.frame(na.omit(selected.locations.t))

selected.locations=selected.locations[c(1,2)]


###
sp.xy = data.frame(selected.locations)

loc=SpatialPoints(sp.xy)

sp.xy=(sp.xy[,1:2]/(1000))
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

##Check kerneal smoothing maps with spatstat. my.owin=as.owin(my.ppp)
loca=SpatialPoints(sp.xy)
my.ppp=as.ppp(loca)
my.density=density(my.ppp, edge=TRUE) # sigma=1, 
plot(my.density)
plot(loca, add=TRUE)
# plot  LGA 10 map.
LGA.shp <- readShapePoly("LGA10new.shp")
plot(LGA.shp, add=TRUE)


#####step 5: Models: Enviromental only. Test using different env variables.
#                   Enviromental and distance variables
#                   Bias corrected.

#Model 1:  Enviromental varibales only
#ppm.1 = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)
ppm.1 = ~ poly(clay
               ,sbd
               ,tpo
               ,AnnualMeanTemperature
               ,awc
               ,habit1decimal
               ,habit2decimal
               ,habit3decimal
               ,fpcnew
               ,AnnualPrecipitation
               ,degree = 2, raw = TRUE)


ppm.1 = ~ poly(Annual_Mean_Temperature
               ,Annual_Precipitation
               ,Max_Temperature_of_Warmest_Month
               ,Min_Temperature_of_Coldest_Month
               ,fpcnew
               ,degree = 2, raw = TRUE)



scales = c( 0.5, 1, 2, 4, 8,16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.1)#find the spatial resolution best for analysis

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, 
                    criterion = "nlgcv", n.fits = 100)#Fitting a regularisation path of point process models.#a LASSO penalty that optimises non-linear GCV

ppmFit.1$beta  

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "blockCV", n.blocks = 50,block.size = 50, sp.scale = 1, n.fits = 100)
ppmFit.1$beta

#use interaction terms; family = "area.inter", r = 1    |     interaction = Strauss(r = 6)


ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "blockCV", n.blocks = 5,
                    block.size = 10, sp.scale = 1, n.fits = 50, family = "area.inter", r = .1)
ppmFit.1$beta

#Check residuals]

diagnose.ppmlasso(ppmFit.1)
resid.plot = diagnose(ppmFit.1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.1, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.1, which = "y", type = "Pearson", compute.sd = TRUE)

#Check level of interaction

kenv = envelope(ppmFit.1,fun = Kinhom) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

#Predict and plot

pred.fit.e = predict.ppmlasso(ppmFit.1, newdata=bigquad) #get prdictions.
predi.fit.e <- cbind(xydata, pred.fit.e)     #predictions anad xy locations.
pred.model.1<- rasterFromXYZ(as.data.frame(predi.fit.e )[, c("x", "y", "pred.fit.e")])
plot(pred.model.1, main=" koala density-warton method/ env only")
plot(loc, add=TRUE)
(total_koala=cellStats(pred.model.1,sum))        #check total koal in the area


#Model :2 Enviromental and distance covariate which are of original scale (bigquad).

ppm.2 = ~ poly(clay
               ,sbd
               ,tpo
               ,AnnualMeanTemperature
               ,awc
               ,habit1decimal
               ,habit2decimal
               ,habit3decimal
               ,fpcnew
               ,AnnualPrecipitation,degree = 2 
               ,raw = TRUE)+poly(distance_primaryandlink
                                ,distance_motorwayandlink
                                ,degree = 2, raw = TRUE)
### raster data set 2

ppm.2 = ~ poly(Annual_Mean_Temperature
               ,Annual_Precipitation
               ,Max_Temperature_of_Warmest_Month
               ,fpcnew
               ,degree = 2, raw = TRUE)+poly(distance_primaryandlink
                                             ,distance_motorwayandlink
                                             ,degree = 2, raw = TRUE)


ppmFit.2 = ppmlasso(ppm.2, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 2, n.fits = 100)
diagnose.ppmlasso(ppmFit.2)
#             Check residuals
resid.plot = diagnose(ppmFit.2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.2, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.2, which = "y", type = "Pearson", compute.sd = TRUE)

#check level of interaction

kenv = envelope(ppmFit.2, fun = Kest) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")
# Make predictions

pred.fit.dis = predict.ppmlasso(ppmFit.2, newdata=bigquad)
predi.fit.dis <- cbind(xydata, pred.fit.dis) 
pred.model.2<- rasterFromXYZ(as.data.frame(predi.fit.dis)[, c("x", "y", "pred.fit.dis")])
plot(pred.model.2, main=" koala density-warton method/ env_dis")
(total_koala=cellStats(pred.model.2,sum))


# Make predictions to new data. Bias correction step:now predict using quadrature (biguad.2) with minimum distances 

pred.fit.crt = predict.ppmlasso(ppmFit.2, newdata=bigquad.2)
predi.fit.crt <- cbind(xydata, pred.fit.crt) # xydatan was chnaged to xydata.
pred.model.crt  <- rasterFromXYZ(as.data.frame(predi.fit.crt )[, c("x", "y", "pred.fit.crt")])
plot(pred.model.crt , main=" koala density-warton method/ bias corrected")
(total_koala=cellStats(pred.model.crt, sum))

#####step 6.#plot all with same legend######

par(mfrow=c(2,2),oma=c(1,1,1,1))
plot(pred.model.1, zlim = c(0, 7),main="env only")
plot(pred.model.2, zlim = c(0, 7),main="env_dis")
plot(pred.model.crt, zlim = c(0, 7), main="bias corrected")
LGA.shp <- readShapePoly("LGA10new.shp")
plot(LGA.shp, add=TRUE)



##### Detail way of fitting the model:

#A regularisation path of Poisson point process models

quad.1k = sample.quad(bigquad, 1)
ppm.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, criterion = "nlgcv")
diagnose(ppm.fit)

#4.3 Block cross-validation
#block cross-validation as a method for choosing the LASSO penalty
#area interaction model with radius 2k and lasso penalty by5-fold cross validation.

final.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, 
                     criterion = "blockCV", n.blocks = 5, block.size = 10)
#Predict and plot

pred.final.fit = predict.ppmlasso(final.fit, newdata=bigquad.2)

predictions.final.fit <- cbind(xydata, pred.final.fit) # xydatan was chnaged to xydata.
pred.final<- rasterFromXYZ(as.data.frame(predictions.final.fit )[, c("X", "Y", "pred.final.fit")])
plot(pred.final, main=" Koala density-warton method bias corrected blockCV")

#diagnostic plots:residuals and lurking variable plots.

diagnose.ppmlasso(final.fit)
diagnose.ppmlasso(final.fit, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(final.fit, which = "y", type = "Pearson", compute.sd = TRUE)
resid.plot = diagnose(final.fit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas biascorrected interaction model")

########## For large maps.

opar <- par() #make a copy of current settings

#par(opar)          # restore original settings

mypar <- par(mar=c(1.5,1.5,1.5,1.5), oma=c(0.5,0.5,0.5,0.5))


######## =====================================############
#             First run the above code the get objects data to run belwo code.
#
#             IWLR |DWP
#       use same koala and raster dataset prepared above for comparison.
#
########======================================#############

# selected.locations = koala data over the study area.

# habitat.rr = my raster subset for the analysis.

acsel210 = selected.locations
acsel210$Haskoala <- 1
myfullstack.sub=habitat.rr

selected.loc=SpatialPoints(selected.locations)

#subset raster as df

myfullstack.subdf= as.data.frame(myfullstack.sub) 


#############

# get all presence absences. Put all presences and mark them as 1 and all others zero. Then get 0/1.
# This is the response data in iwlr and dwpr. Generally need presence locations only for other methods but
#this is approximating ppm with logistic regression. So need 0/1.

r <- raster(myfullstack.sub, layer=4)   # get a raster layer

dfull <- as.data.frame(r, xy = TRUE)

dpart = cbind(acsel210,extract(r,acsel210[,1:2]))

dpart <- subset(acsel210, Haskoala==1)

#check minimum valuve of raster r
r

values(r )[values(r) > -1.414443] = 0
rspec <- r

plot(r)
cells <- cellFromXY(rspec, as.matrix(dpart[, c("X", "Y")]))

rspec[cells] <- dpart$Haskoala

plot(rspec)
text(dpart$X, dpart$Y, lab = dpart$Haskoala)

Press.g<- as.data.frame(rspec, xy=TRUE)   # Has snome NA valuves.

#### Press.g need this for Hefley method for group variable.

#rename variable awc

colnames(Press.g)[3] <- "koala"
ptest=cbind(Press.g,myfullstack.subdf)   ### now remove na
Press= na.omit(ptest)  # same as sp.at=Press ;  Pres=koala

# save this for future use in comapringlikelihood.

Press.my=Press

# koala presence locations

Pres <- Press[,3]   # this create a vector, Press[c(3)] create a sublist.

### now get design matrix

X.des=Press[,4:10]
X.des=as.matrix(X.des)  # quadrature points.

####      iwlr      #######==============================
         


up.wt = (1.e6)^(1 - Pres)  # up.wt = (10^6)^(1 - Pres)
iwlr = glm(Pres ~ X.des, family = binomial(), weights = up.wt)

## check coefficents
iwlr
# 

dd1 <- as.data.frame(X.des) # get coordinates of the design matrix for predictions. similr to warton method.
pred.iwlr = predict(iwlr, newdata=dd1)

#r <- raster(myfullstack.sub, layer=2) 

dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only

pred.iwlr <- cbind(xydatan, pred.iwlr)
pred.iwlreg <- rasterFromXYZ(as.data.frame(pred.iwlr)[, c("x", "y", "pred.iwlr")])
plot(pred.iwlreg,asp=1)
plot(selected.loc, add=TRUE)
###       DWPR     ######===================================

p.wt = rep(1.e-6, length(Pres))
p.wt[Pres == 0] = 8310/sum(Pres == 0)
dwpr = glm(Pres/p.wt ~ X.des, family = poisson(), weights = p.wt)

# check coefficents 
dwpr
#

dd <- as.data.frame(X.des)
pred.dwpr = predict(dwpr, newdata=dd)

dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg,asp=1)

# plot lcoations over the map
plot(selected.loc, add=TRUE)

############## STOP here ###############

# 5.3 Assessing the variability in likelihood for different numbers quadrature points.
# to generate different sie quadrature points.

quad=Press.my[c(-3)]
quad=as.data.frame((quad))
#quad <- as.data.frame(habitrasters, xy=TRUE) 
quad <- na.omit(quad)   ### question
colnames(quad)[1] <- 'X'; colnames(quad)[2] <- 'Y'

#load("Quad100m.RData") #xy and variables.

# another way is to use Pres.my to generate points to generate Quadraatures and use them seperatly for iwlr|dwpr.

n.quad = c(1000, 2000, 5000, 8310) # number of quadrature poiints.
quad.inc = sample(1:dim(quad)[1], 1000)
assign(paste("quad.", n.quad[1], sep = ""), quad[quad.inc[1:n.quad[1]],])
for (i in 2:length(n.quad)){
  quad.inc = c(quad.inc, sample(setdiff(1:dim(quad)[1], quad.inc),
                                (n.quad[i] - n.quad[i - 1])))
  assign(paste("quad.", n.quad[i], sep = ""), quad[quad.inc[1:n.quad[i]],])
}

#3# use them to iwlr and dwpr. 



#compare the likelihood of PPMs fitted using downweighted Poisson regression:
#create species data
#Press has x y koala 0 1  and covariates
# get koala 1 and select variables excluding koala

sp.dat.1=subset(Press.my, koala==1)
sp.dat=sp.dat.1[c(-3)]
colnames(sp.dat)[1] <- 'X'; colnames(sp.dat)[2] <- 'Y'
#to compare likelihoods from iwlr , dwpr and ppmlasso. recall : same as sp.at=Press ;  Pres=koala

# remember to chnage variable names.

sp.dat$Pres = 1
loglik = rep(NA, length(n.quad))

for (i in 1:length(n.quad)){
  quad = get(paste("quad.", n.quad[i], sep = ""))
  quad$Pres = 1
  all.dat = na.omit(data.frame(rbind(sp.dat, quad)))
  X.des = as.matrix(cbind(poly(all.dat$AnnualMeanTemperature, all.dat$twi, all.dat$tpo,
                               all.dat$distance_trunkandlink, all.dat$habit3decimal,  
                               all.dat$nitro,all.dat$roadk, degree = 2, raw = TRUE)))   
  p.wt = rep(1.e-8, dim(all.dat)[1])
  p.wt[all.dat$Pres == 0] = 10000/n.quad[i]
  z = all.dat$Pres/p.wt
  dwpr = glm(z ~ X.des, family = poisson(), weights = p.wt)
  mu = dwpr$fitted
  loglik[i] = sum(p.wt*(z*log(mu) - mu))
}
plot(n.quad, loglik, log = "x", type = "o")

#coudn`t get the Renners plot with all lines.


#######   Hefley METHOD                              #########
#                                                     #
#            Use same dataset                         #
#              Press.g
#               habitat.rr from climatic data folder. ###


#data required are : Koala location as 0 | over the entire study area.
#                  :   raster data for above koala 0 | locations.
# we saved data set "Press.my" has this variables.
#  renamed object 

# Press.g   # x, y, koala
# big.grp # x,y variables
#rename koala colomn to presence:


#### This method is to correct group size so model. same covariate possible.
# create a group variable and make it 1all locations. Otehrwise aggregate koalas in a grid and consider them as one unit or group..

# r is an empty raster creaated above.

group.r=rasterize(selected.locations, r) # 

plot(group.r)

group.r[is.na(group.r[])] <- 0 

group.rr=group.r+r

plot(group.rr)

group=getValues(group.rr)

# remove xy from big.grp  data

myFDc=cbind(Press.g,group,big.grp)
myFDc=myFDc[c(-1,-2)]

myFD1=na.omit(myFDc)
colnames(myFD1)[1] <- "presence"

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 

#####Step 5: select only 1000 absences (monticarlo points as in hefleys method??)####
set.seed(123)                                   # size data comes from absence locations.
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 

##### Step 6: now take a random sample of 80 and assign detected 1 non detected 0.####

train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = 300,replace=FALSE)
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0

##### Step 7: Create the final data sets for the analysis####

Detection.data= rbind(detected,notdetected) 
str(Detection.data)

IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus no koalas data from myFD1.

ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.


IPP.ignored=glm(presence~Annual_Mean_Temperature
                + Annual_Precipitation
                + Max_Temperature_of_Warmest_Month
                + Min_Temperature_of_Coldest_Month
                + fpcnew
                ,family="binomial",weights=100000^(1-presence),data=IPP.data) # IPP.data2 added.



summary(IPP.ignored)
detecPress.gk=subset(Detection.data, select=c(-2))
Press.gk=as.data.frame(Press.gk)

# # ROC curve & GOF metrics
# myPred = prediction(predict(Detection.model, type = "response"), Press.gk$koala)
# perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
# plot(perf, colorize = T)
# myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), Press.gk$koala, 0.5)
# myPredPerfs


#############  
set.seed(125)
#Detection model: steps as in Hefley`s code`
Detection.model=glm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    ,family= "binomial"
                    , data=Detection.data)

unclass(summary(Detection.model))

###################### Improve the model autocorrelation check and residual approach.#####
#  got to step 4 if want to ignore this section
######## Detection modelimproved by adding 
# 1. Correlogrm prepared, Moran`s` I
# 2. spatial autocorelation corrected using residulas.
# 3. ROC curve 
# 4. LRtest 
###               ###               ###           ### 
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
library(epicalc)
myPredictors=c("distance_motorwayandlink","distance_primaryandlink","fpcnew")
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "presence ~"
  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
  myStrARN = paste(myStrN,"+ AR1")
  Detection.model.2 = glm(as.formula(myStrARN), data = Detection.data, family = "binomial")
  mylrtest = lrtest(Detection.model,Detection.model.2)
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

myPred = predict(myARStack, Detection.model.2, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main=" detection koala")

# ROC curve & GOF metrics
myPred = prediction(predict(Detection.model.2, type = "response"), Detection.data$presence)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(Detection.model.2, type = "response"), Detection.data$presence, 0.5)
myPredPerfs



#####################

#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(Detection.model.2,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 

# set very small valuves to 0.1
p.det <- ifelse(p.det<0.1,0.1,p.det)

hist(p.det, breaks=70)


# p.det 1 or very low valuves create convergencce issues.
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det

######Step 5: - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected=glm(presence~Annual_Mean_Temperature
                + Annual_Precipitation
                + Max_Temperature_of_Warmest_Month
                + Min_Temperature_of_Coldest_Month
                + fpcnew
                ,weights=(1/p.det)*10000^(1-presence),family="binomial",  data=IPP.data)


summary(IPP.corrected)
####Step 6: Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####

#use only the significant covariates, tpo +hpop+lot_density+sbd
ZTGLM.corrected=vglm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)

# step 7:  Detion model Map predictions

myPred = predict(habitat.rr, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")

#        IPP ignored model Map predictions

myPred2 = predict(habitat.rr, IPP.ignored, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model |detection ignored- number of koalas")

#        IPP corrected model Map predictions

myPred2 = predict(habitat.rr, IPP.corrected, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model |detection considered- number of koalas")

myPred3 = predict(habitat.rr, ZTGLM.corrected, type = "response")
plot(myPred3,  main="ZTGLM-group size koalas in a grid model")


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

#### chnage below variable list based on s=identified variables in each model.

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
save.image(file="hefleydatapreparation.v3.3.RData")
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





######   Drorazio method  ##############

          #use same data 
####  ####  ####    ###  ####  ####  ####    ###




