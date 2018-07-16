library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(lmtest);library(spatial.tools);library(VGAM);library(mosaic);library(faraway);library(gstat)  #
library(ncf);library(foreign);library(nlme)   ;library(MASS);library(ROCR);library(vcd)
library(RColorBrewer);library(classInt);library(ppmlasso);library(usdm) ; library(ncf); library(epicalc)
library(optiRum)
##############

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")

##### import rasters uncorrelated set.
setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/myvars")
myfullstack.v <- list.files(pattern="\\.tif$",full.names=TRUE)
myfullstack.v = stack(myfullstack.v)
names(myfullstack.v)

#### set extent 
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.v, myextent, snap="out") # myextent <- SpatialPoints(selected.locations)
plot(habitat.rr)

# myextent <- readOGR("Chull.shp") ### code to create Chull is in research.
# habitat.rr <- mask(myfullstack.v, myextent)
# plot(habitat.rr)

### Emanuele data.
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Emanuele_data//koala1997_2013month.RData")
myalldata <- koala1997_2013
myalldata <- myalldata[c(1:3)]
names(myalldata )[1]<-"X"
names(myalldata )[2]<-"Y"
selected.locations <- myalldata
#selected.locations= subset(myalldata,yearnew >2010 & yearnew < 2013)
selected.locations = selected.locations[,1:2]
### remove duplicate
selected.locations <- selected.locations[duplicated(selected.locations$X & selected.locations$Y), ]
selected.locations=SpatialPoints(selected.locations)
selected.locations= na.omit(selected.locations) 

############ Create araster:
mydata.r <-  subset(habitat.rr, c(1))
res(mydata.r) <- 500
#extent <- extend(selected.locations, extent(selected.locations) + 3000)
# Sample points from within each grid cell:
selected.locations <- gridSample(selected.locations, mydata.r, n = 1)
plot(selected.locations)
selected.locations <- as.data.frame(selected.locations)



# # Select records based on distance
# source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas.
# all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
# select.locations = subset(all.locations, all.locations$disttoitself >=1) #3 350
# selected.locations = select.locations[,1:2] # selected locations are used again in line 94
# loc=SpatialPoints(selected.locations)
# plot(loc)

######Step 2: load rasters and crop to my extent of interest######

# ######## All climate rasters . Here  have scaled rasters,

# myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
# myfullstack.b = scale(stack(myfullstack.b))
# myfullstack.b = stack(myfullstack.b) # 63 variables # check VIF.
# 
# #plot(myfullstack.b)
# #myextent<-extent(min(selected.locations$X)-2000,max(selected.locations$X)+2000
# #,min(selected.locations$Y)-20000, max(selected.locations$Y)+20000)
# #myextent<-extent(min(selected.locations$X)-3000,max(selected.locations$X)+3000,min(selected.locations$Y)-3000,max(selected.locations$Y)+3000)
# home_range <- raster("home_range.tif") # for emanuel models
# myextent=c(387900, 553100, 6862400, 7113600) 
# habitat.rr=crop(myfullstack.b, myextent, snap="near")
# 
# habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,33, 38)) # select my variables
# names(habitat.rr)
# plot(habitat.rr)



# habitat.rr = na.omit(habitat.rr)
# (vifstep(habitat.rr, th=3))
# 
# ### correlation pearson 
# habitat.rr<- subset(myfullstack.v, c(34,33,30,29,27,23,20,11,10,9,8,7,6,4,3))
# 
# writeRaster(habitat.rr, filename= names(habitat.rr), bylayer= TRUE, format="GTiff")
######Step 3. Create quadrature points and species xy data. Remove NA. Keep xy colomns for preditions and rasterize. ###### starting with Warton method.

big.grp <- as.data.frame(habitat.rr, xy=TRUE) # to be used later in Hefley method.

#

bigquad <- as.data.frame(habitat.rr, xy=TRUE,na.rm=TRUE)
xydata <- bigquad[c(1,2)] ### use to get predictions and rasterise.



bigquad =cbind((bigquad[,1:2]/(1000)),bigquad[c(-1,-2)]) # xy convert to km and integers. decimals take days to run..
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'  
xydatan <- as.data.frame(lapply(bigquad[c(1,2)], as.integer)) 
dd=bigquad[c(-1,-2)]
bigquad <- cbind(xydatan,dd )   # Quadrature points XY as integers.

##### Step 4:Create quadrature points with values set to minimum of distance variables ####
#             use for bias corretion model
#             model based control of observer bias based on same level of distance

bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) # roads_dist
bigquad.2$distance_tertiaryandlink = min(bigquad.2$distance_tertiaryandlink) #distance_unclassified

####### Step 5.Koala XY data.####


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



#####step 5: Models: Warton method: Enviromental only. Test using different env variables#######
#                   Enviromental and distance variables
#                   Bias corrected.

#Model 1:  Enviromental varibales only#####
#ppm.1 = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)


ppm.1 = ~ poly(twi
               ,tpo
               ,Precipitation_Seasonality
               ,Precipitation_of_Driest_Quarter
               ,nitro
               ,hpop
               ,habit3decimal
               ,habit2decimal
               ,fpcnew
               ,elev
               ,awc
               ,aspect91
               ,roads_other
               ,degree = 2, raw = TRUE)



scales = c( 0.5, 1, 2, 4, 8, 16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.1)#find the spatial resolution best for analysis



ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "nlgcv", n.blocks = 5,
                    block.size = 10, sp.scale = 1, n.fits = 100) # family = "area.inter", r = 2 # 5, 10blockCV
ppmFit.1$beta

#Check residuals]

diagnose.ppmlasso(ppmFit.1)
resid.plot = diagnose(ppmFit.1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.1, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.1, which = "y", type = "Pearson", compute.sd = TRUE)

#Check level of interaction
#Calculates simulation envelopes for goodness-of-fit

kenv = envelope.ppmlasso(ppmFit.1,fun = Kinhom) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope", shadecol= "yellow")
plot(kenv, sqrt(./pi) ~ r,  shadecol= "yellow")#shade=NULL,
#Predict and plot

pred.fit.e = predict.ppmlasso(ppmFit.1, newdata=bigquad) #get prdictions.
predi.fit.e <- cbind(xydata, pred.fit.e)     #predictions anad xy locations.
pred.model.1<- rasterFromXYZ(as.data.frame(predi.fit.e )[, c("x", "y", "pred.fit.e")])
plot(pred.model.1, main=" koala density-warton method/ env only")
plot(loc, add=TRUE)
(total_koala=cellStats(pred.model.1,sum))        #check total koal in the area


#Model :2 Enviromental and distance covariate which are of original scale (bigquad).#####

ppm.2 = ~ poly(twi
               ,tpo
               ,Precipitation_Seasonality
               ,Precipitation_of_Driest_Quarter
               ,nitro
               ,hpop
               ,habit3decimal
               ,habit2decimal
               ,fpcnew
               ,elev
               ,awc
               ,aspect91
               ,roads_other
               ,degree = 2, raw = TRUE)+ poly(distance_primaryandlink,
                                              distance_tertiaryandlink,degree = 2, raw = TRUE)

# interactions allowed between variables but not between climatic and distance variables.
# linear, quadratic and first order interactions.

ppmFit.2 = ppmlasso(ppm.2, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "blockCV", n.blocks = 5,
                    block.size = 10, sp.scale = 1, n.fits = 100)#, family = "area.inter", r = 2
diagnose.ppmlasso(ppmFit.2)
#             Check residuals
resid.plot = diagnose(ppmFit.2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.2, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.2, which = "y", type = "Pearson", compute.sd = TRUE)

#check level of interaction
kenv = envelope.ppmlasso(ppmFit.2, fun = Kest, nsim = 20) # #check k envelop and see what level of interaction occurs.
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

plot(loc,add=TRUE)

##### Select the correct folder for writing.
jpeg("warton3.jpeg",width=12,height=8.4,units="in",res=300) # 8:4 for three mpas. one map 8:10
# png(width=3,height=3,units="in",res=600)

par(mfrow=c(1,2),oma=c(1,4,1,2))
#plot(pred.model.1, zlim = c(0, 5),main="env only", legend.width = 1.5)
plot(pred.model.2, zlim = c(0, 5),main="env_dis", legend.width = 1.5)
plot(pred.model.crt, zlim = c(0, 5), main="bias corrected", legend.width = 1.5)

graphics.off()

#############Model 2. Data preparation for DWPR ####
# selected.locations = koala data over the study area.

# habitat.rr = my raster subset for the analysis.

acsel210 = selected.locations   # same as sp.xy in the above model
acsel210$Haskoala <- 1
myfullstack.sub=habitat.rr

selected.loc=SpatialPoints(selected.locations)

#subset raster as df

myfullstack.subdf= as.data.frame(myfullstack.sub) 



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

Press.g<- as.data.frame(rspec, xy=TRUE)   

#### Press.g need this for Hefley method for group variable.

#rename variable awc

colnames(Press.g)[3] <- "koala"
ptest=cbind(Press.g,myfullstack.subdf)   ### now remove na
Press= na.omit(ptest)  # same as sp.at=Press ;  Pres=koala

Press.my = Press# save this for future use in comapringlikelihood.

# koala presence locations
# set.seed(123)
# smp_size <- floor(0.75 * nrow(Press))
# train_ind <- sample(seq_len(nrow(Press)), size = smp_size)
# 
# train <- Press[train_ind, ]
# test <- Press[-train_ind, ]
#Press=train ## after data partitioning
Pres <- Press[,3]   # this create a vector, Press[c(3)] create a sublist.
X.des=Press[,4:17]

# Split the data into training and test set
## 75% of the sample size


## set the seed to make your partition reproducible

# X.des.2=as.matrix(X.des.2)  # quadrature points with common level of distance to all locations.
# X.des=as.matrix(X.des)  # quadrature points.

## ENV and distance 
X.des.ed=X.des 
X.des.ed <- as.matrix(X.des.ed)   #variables together

## get distance cvariates only

X.des.d <- X.des[c(-1,-2,-5,-6,-7)]
X.des.d <- as.matrix(X.des.d)

##Env variables only
X.des.e <- X.des[c(-3,-4)]
X.des.e <- as.matrix(X.des.e)

#### for bias correction 
X.des.ed <- as.data.frame(X.des.ed)
X.des.2 <- X.des.ed
X.des.2$distance_primaryandlink = min(X.des.2$distance_primaryandlink) 
X.des.2$distance_tertiaryandlink = min(X.des.2$distance_tertiaryandlink)#distance_motorwayandlink
X.des.bc <- as.matrix(X.des.2)

###       DWPR     ######===================================
# QUesiton : how to evaluvate this model? AUC and TSS.
# 

p.wt = rep(1.e-6, length(Pres))
#X.des <- X.des[,-6]
p.wt[Pres == 0] = 44415/sum(Pres == 0)




##### 1. Env only model X.des is a matrix

dwpr= glm(Pres/p.wt ~ poly(X.des.e, degree = 2, raw = TRUE ), family = poisson(), weights = p.wt)
dd <- as.data.frame(X.des.e) # Press.my
pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)

#cv.glm(as.data.frame(cbind(Pres/p.wt, X.des.e)), dwpr, K=50)

## predict to full dataset before partitioning.
dd <- as.data.frame(X.des.e) # Press.my
pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)

pred.dwpr=logit.prob(pred.dwpr)
dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)] # coordinates

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg.e <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg.e,asp=1, main = "dwpr_bias not corrected")


###### 2. env distance model  ## seperate Env+ distance. not together as one dataframe.

dwpr = glm(Pres/p.wt ~ poly(X.des.e,degree = 2, raw = TRUE ) + poly(X.des.d,degree = 2, raw = TRUE ), family = poisson(), weights = p.wt)

pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)
pred.dwpr=logit.prob(pred.dwpr)

xydatan <- Press[c(1,2)] # coordinates

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg.e.d <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg.e.d,asp=1, main = "dwpr_bias not corrected")

#### 3. bias correction ; not working properly
dd.2 <- as.data.frame(X.des.bc)
pred.dwpr = predict(dwpr, newdata=dd.2, response=TRUE)
pred.dwpr=logit.prob(pred.dwpr)
xydatan <- Press[c(1,2)] # coordinates

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg.e.d <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg.e.d,asp=1, main = "dwpr_bias not corrected")


# 
# ###### env distance model  together as one dataframe.
# dwpr = glm(Pres/p.wt ~ poly(X.des.ed,degree = 2, raw = TRUE ) , family = poisson(), weights = p.wt)
# dd <- as.data.frame(X.des.ed)
# # dd=poly(X.des.ed,degree = 2, raw = TRUE )
# # dd <- as.data.frame(dd)
# pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)
# pred.dwpr=logit.prob(pred.dwpr)
# dfull <- as.data.frame(r, xy = TRUE)
# 
# xydatan <- Press[c(1,2)] # coordinates
# 
# pred.dwpr <- cbind(xydatan, pred.dwpr)
# pred.dwpreg <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
# plot(pred.dwpreg,asp=1, main = "dwpr_bias not corrected")


###========
###### bias corrected model with the matrix. but it should be the same previous model to new data that didnt work.
X.des.bc <- as.matrix(X.des.bc)
dwpr = glm(Pres/p.wt ~ poly(X.des.bc, degree = 2, raw = TRUE ) , family = poisson(), weights = p.wt)
dd <- as.data.frame(X.des.bc)
# dd=poly(X.des.ed,degree = 2, raw = TRUE )
# dd <- as.data.frame(dd)
pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)
pred.dwpr=logit.prob(pred.dwpr)
dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)] # coordinates

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg1 <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg1,asp=1, main = "dwpr_bias not corrected")

##### Select the correct folder for writing.
jpeg("dwpr.jpeg",width=12,height=8.4,units="in",res=300) # 8:4 for three mpas. one map 8:10
# png(width=3,height=3,units="in",res=600)

par(mfrow=c(1,2),oma=c(1,4,1,2))
#plot(pred.model.1, zlim = c(0, 5),main="env only", legend.width = 1.5)
plot(pred.dwpreg.e.d, main="env_dis", legend.width = 1.5)
plot(pred.dwpreg1, main="bias corrected", legend.width = 1.5)

graphics.off()
#############



prob=predict(dwpr,type=c("response"))
Press$prob=prob
library(pROC)
g <- roc(koala ~ prob, data = Press)
g
plot(g, main= "dwpr")

library(boot)
(cv.err <- cv.glm(koala.x, dwprr)$delta)

##ok. lets create a complete dataframe.
koala <- Pres/p.wt
koala.x <- cbind(koala, X.des)
koala.x <- as.data.frame(koala.x)
dwprr <- glm(koala ~ poly(X.des, degree= 2, raw=TRUE),family = poisson(), weights = p.wt)


cv.glm(koala.x, dwpr, K=50)

# $K
# [1] 50
# 
# $delta
# [1] 49029683700 49029683700
###=======
# bias correction map, for bias correction new data set`s distance variables are set to a common level as in ppmlasso method.`
# X.de.2 has this common level.
# #  is the data set for bias correction by setting distance to minimum valuve of distance.

## for bias correction env + min distance
# set same distance level to distance covariate as in lasso method.
X.des.2 <- X.des
X.des.2$distance_primaryandlink = min(X.des.2$distance_primaryandlink) 
X.des.2$distance_tertiaryandlink = min(X.des.2$distance_tertiaryandlink)#distance_motorwayandlink
X.des.bc <- as.matrix(X.des.2)


X.des.bc <- poly(X.des.bc,degree = 2, raw = TRUE )
X.des.bc <- as.data.frame(X.des.bc)

pred.dwpr.2 = predict(dwpr, newdata=X.des.bc,response=TRUE) # bias corrected check distance variables in  dd2 dataset.
pred.dwpr.2=logit.prob(pred.dwpr.2)

pred.dwpr.2 <- cbind(xydatan, pred.dwpr.2)
pred.dwpreg.2 <- rasterFromXYZ(as.data.frame(pred.dwpr.2)[, c("x", "y", "pred.dwpr.2")])
plot(pred.dwpreg.2, main = " dwpr bias corrected")
plot(selected.loc, add=TRUE)

#######                #########
###      Maxent model  ####
###                    ####
###                    ####
##### Maxent model===============
## we use same datasets and criteria. import data and select in the same manner we did for previous analysis.
## Maxnet is an R package that for fitting Maxent species distribution models using the R package glmnet.
##### Sstep 1. 
setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/myvars")
myfullstack.v <- list.files(pattern="\\.tif$",full.names=TRUE)
myfullstack.v = stack(myfullstack.v)
names(myfullstack.v)
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.v, myextent, snap="near")
# writeRaster(habitat.rr, filename = names(habitat.rr), bylayer = TRUE, format= "ascii", overwrite=TRUE)
habitat.rr

### Emanuele data.
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Emanuele_data//koala1997_2013month.RData")
myalldata <- koala1997_2013
myalldata <- myalldata[c(1:3)]
names(myalldata )[1]<-"X"
names(myalldata )[2]<-"Y"
selected.locations <- myalldata
#selected.locations= subset(myalldata,yearnew >2010 & yearnew < 2013)
selected.locations = selected.locations[,1:2]
### remove duplicate
selected.locations <- selected.locations[duplicated(selected.locations$X & selected.locations$Y), ]
selected.locations=SpatialPoints(selected.locations)
selected.locations= na.omit(selected.locations) 

############ Create araster:
mydata.r <-  subset(habitat.rr, c(1))
res(mydata.r) <- 500
#extent <- extend(selected.locations, extent(selected.locations) + 3000)
# Sample points from within each grid cell:
selected.locations <- gridSample(selected.locations, mydata.r, n = 1)
plot(selected.locations)
selected.locations <- as.data.frame(selected.locations)
selected.loc <- selected.locations
# # loc = SpatialPoints(selected.loc)
# # plot(loc, add = TRUE)
# # remove duplicates
# selected.loc.dups=duplicated(selected.loc[, c("x", "y")])
# selected.loc <-selected.loc[!selected.loc.dups, ]
# 
# 
# # Step 3. Get raster data
# 
# myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
# myfullstack.b = stack(myfullstack.b)
# myextent<-extent(min(selected.loc$x)-3000,max(selected.loc$x)+3000,min(selected.loc$y)-3000,max(selected.loc$y)+3000)
# 
# habitat.rr=crop(myfullstack.b, myextent, snap="near")
# 
# habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,34, 41)) # select my variables
# plot(habitat.rr)
# 
# habitat.rr=scale(habitat.rr)
# plot(habitat.rr)

#### step 4.#####   selected.loc <-  selected.locations
fold <- kfold(selected.loc, k=10) # add an index that makes five random groups of observations
selected.loctest <- selected.loc[fold == 1, ] # hold out one fifth as test data
selected.loctrain <- selected.loc[fold != 1, ] # the other four fifths are training data

habitat.rr <- scale(habitat.rr)


TrainEnv <- extract(habitat.rr,selected.loctrain)

# Background points are for model validation.No partitioning. 
set.seed(0)
backgr <- randomPoints(habitat.rr, 20000)
absvals <- extract(habitat.rr, backgr)

#we make a presence/absence collumn where 1=presence for our occurance points, and 0=absence for our background points
presabs <- c(rep(1, nrow(TrainEnv)), rep(0, nrow(absvals)))

sdmdata <- data.frame(cbind(presabs, rbind(TrainEnv, absvals)))

sdmdata = na.omit(sdmdata)
#we make a subset of that dataset, without the presence and absence values
data <- sdmdata[,-1]
data = na.omit(data)

#we run the maxnet function to fit the SDM
#maxnet fits using glmnet
loc = sdmdata$presabs    # presence absence data set. data has env other variables only.
library(maxent); library(maxnet)
koala.model<-maxnet(loc, data,maxnet.formula(loc, data, classes="lh")) #, regmult = 2

# Default choices for feature classes and regularization.
#  The second call overrides the default feature classes, using only linear and
# quadratic features for the continuous predictor variables.
####
#fits Maxent models using the same feature classes
#(linear, quadratic, hinge, etc.) and regularization
#mod <- maxnet(loc, data, maxnet.formula(loc, data, classes="default"))
summary(koala.model)
##three types of response plots
plot(koala.model, type = "exponential")
plot(koala.model, type = "cloglog")
plot(koala.model, type = "logistic",common.scale = T)

plot(koala.model, "fpcnew")
plot(koala.model, "Annual_Mean_Temperature")

gg=predict(habitat.rr,koala.model, response=TRUE)
# ps <- unscale(gg,habitat.rr)
# pp=unscale(gg, center = NULL, scale = NULL)
# plot(ps)

gg=predict(habitat.rr, koala.model,clamp=T, response=TRUE)
ggg<- logit.prob(gg)
plot(ggg)

#### what if I set distance raster to minimum valuve and stack them, then use to predict.

names(habitat.rr)
plot(habitat.rr, 3)
distance_primaryandlink <- subset(habitat.rr, 3)
distance_tertiaryandlink <- subset(habitat.rr, 4)
distance_primaryandlink <- (distance_primaryandlink*0)+-1.332977
distance_tertiaryandlink <- (distance_tertiaryandlink*0)+-1.141844

habit <- subset(habitat.rr, c(1,2,5,6,7,8,9,10,11,12,13,14,15,16,17))
habitat.r.r <- stack(habit,distance_tertiaryandlink,distance_primaryandlink )

gga=predict(habitat.r.r, koala.model,clamp=T, response=TRUE)
ggga<- logit.prob(gga)
plot(ggga)


# library(boot)
# (cv.err <- cv.glm(koala.x, dwprr)$delta)
## plot
jpeg("maxent.jpeg",width=12,height=8.4,units="in",res=300) # 8:4 for three mpas. one map 8:10
# png(width=3,height=3,units="in",res=600)

par(mfrow=c(1,2),oma=c(1,4,1,2))
plot(ggg,zlim = c(0,0.0025), main="Maxent model bias not corrected", legend.width=2)
plot(ggga,zlim = c(0,0.0025),main="Maxent model bias corrected", legend.width=2)


graphics.off()



## Maxent from software  to get AUC. This step is not reqired if AUC is getable xy can be used provided that env varibales in the same format.
#http://biodiversityinformatics.amnh.org/open_source/maxent/Maxent_tutorial2017.pdf
selected.location <- as.data.frame(selected.locations)
selected.location$species <- "koala"
selected.location <- selected.location[c(3,1,2)]
names(selected.location) <- tolower(names(selected.location))
write.csv(selected.location, "select.maxent.csv",row.names=FALSE)
#### step 2. create rasters in ascii format.#####

setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn/warton_data/warton_data_allclimate/myvars")
myfullstack.v <- list.files(pattern="\\.tif$",full.names=TRUE)
myfullstack.v = stack(myfullstack.v)
names(myfullstack.v)
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.v, myextent, snap="near")
writeRaster(habitat.rr, filename = names(habitat.rr), bylayer = TRUE, format= "ascii", overwrite=TRUE)
habitat.rr

### run maxent from programe. Output asc file is read and mapped.

koala=raster("koala.asc")
plot(koala, main = "bias not corrected")

##### step 3. bias correction needs distance rasters set to minimum distance. #####
distance_primaryandlink= raster("distance_primaryandlink.asc")
d=as.data.frame(distance_primaryandlink)
distance_primaryandlink1=(distance_primaryandlink*0)+-1.412223    # is this distance correct.
plot(distance_primaryandlink1)

writeRaster(distance_primaryandlink1, "distance_primaryandlink1.asc")

distance_motorwayandlink= raster("distance_motorwayandlink.asc")
d=as.data.frame(distance_motorwayandlink)
distance_motorwayandlink1=(distance_motorwayandlink*0)+-1.435722   # is this distance correct.
plot(distance_motorwayandlink1)
writeRaster(distance_motorwayandlink1, "distance_motorwayandlink1.asc")

koala.unbias=raster("koala_1.asc")
plot(koala.unbias, main = " bias corrected")

par(mfrow=c(1,2))

##### step 4. Get outputs from maxent to preare roc curve.
# install.packages("ROCR", dependencies=TRUE)
# install.packages("vcd", dependencies=TRUE)
library(ROCR)
library(vcd)
library(boot)

setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/maxent/results")
presence <- read.csv("koala_samplePredictions.csv")
background <- read.csv("koala_backgroundPredictions.csv")
pp <- presence$Logistic.prediction  # get the column of predictions

testpp <- pp[presence$Test.or.train=="test"] # select only test points
trainpp <- pp[presence$Test.or.train=="train"] # select only test points
bb <- background$Logistic    # porbability of back ground points.

###now ROC from maxent tu
combined <- c(testpp, bb)

label <- c(rep(1,length(testpp)),rep(0,length(bb))) # labels: 1=present, 0=random

pred <- prediction(combined, label)

perf <- performance(pred, "tpr", "fpr")

plot(perf, colorize=TRUE)

# combine into a single vector
# labeled predictions
# True / false positives, for ROC curve
# Show the ROC curve
performance(pred, "auc")@y.values[[1]] # Calculate the AUC
#0.9084246

## step 5. make a bootstrap estimate of the standard deviation of the AUC.
AUC <- function(p,ind) {
  pres <- p[ind]
  combined <- c(pres, bb)
  label <- c(rep(1,length(pres)),rep(0,length(bb)))
  predic <- prediction(combined, label)
  return(performance(predic, "auc")@y.values[[1]])
}
b1 <- boot(testpp, AUC, 100) # do 100 bootstrap AUC calculations

# # Bootstrap Statistics :
#   original       bias    std. error
# t1* 0.9084246 0.0001803703 0.004921075
b1 # gives estimates of standard error and bias

#Bootstrap results can also be used to determine confidence intervals for the AUC:
boot.ci(b1)

#### Cohen Kappa
# calculates Kappa for the threshold given by the minimum presence prediction
confusion <- function(thresh) {
  return(cbind(c(length(testpp[testpp>=thresh]), length(testpp[testpp<thresh])),
               c(length(bb[bb>=thresh]), length(bb[bb<thresh]))))
}
mykappa <- function(thresh) {
  return(Kappa(confusion(thresh)))
}
(mykappa(min(trainpp))) # kappa value
# value       ASE     z Pr(>|z|)
# Unweighted 0.00377 0.0002549 14.79 1.73e-49
# Weighted   0.00377 0.0002549 14.79 1.73e-49

# ## If we want to use the threshold that minimizes the sum of sensitivity and
# specificity on the test data, we can do the following, using the true positive rate and false positive rate
# values from the "performance" object used above to plot the ROC curve:
fpr = perf@x.values[[1]]
tpr = perf@y.values[[1]]
sum = tpr + (1-fpr)
index = which.max(sum)
cutoff = perf@alpha.values[[1]][[index]]
(mykappa(cutoff)) # kappa value 
# 
# value   ASE     z   Pr(>|z|)
# Unweighted 0.2514 0.011 22.86 1.286e-115
# Weighted   0.2514 0.011 22.86 1.286e-115



##  determine binomial probabilities for these two threshold values
mybinomial <- function(thresh) {
  conf <- confusion(thresh)
  trials <- length(testpp)
  return(binom.test(conf[[1]][[1]], trials, conf[[1,2]] / length(bb), "greater"))
}

(mybinomial(min(trainpp)))

# number of successes = 451, number of trials = 451, p-value = 3.943e-09
# alternative hypothesis: true probability of success is greater than 0.958
# 95 percent confidence interval:
#   0.9933796 1.0000000

(mybinomial(cutoff))
# number of successes = 391, number of trials = 451, p-value < 2.2e-16
# alternative hypothesis: true probability of success is greater than 0.1725
# 95 percent confidence interval:
#   0.8377475 1.0000000
# sample estimates:
#   probability of success 
# 0.8669623 

# ##### additional steps: gird sample is not required.####
# #selected.locations=gridSample(selected.locations, myfullstack.b, n = 1)
# selected.locations$species <- "koala"
# colnames(selected.locations)[2] <- 'latitude'
# colnames(selected.locations)[1] <- 'longitude'
# selected.locations <- selected.locations[c(3,1,2)]
# write.csv(selected.locations, "selected.locations.csv",row.names=FALSE )
# koala=raster("koala_9.asc")
# plot(koala)



#######   Hefley METHOD                              #########
#                                                     #
#            Use same dataset  ru IWLR |DWPR firt to get data. #
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

group.r=rasterize(selected.locations, r) #  rcomes from line 356.

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
Press.gk=as.data.frame(detecPress.gk)

# # ROC curve & GOF metrics
# myPred = prediction(predict(Detection.model, type = "response"), Press.gk$koala)
# perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
# plot(perf, colorize = T)
# myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), Press.gk$koala, 0.5)
# myPredPerfs


###
set.seed(125)
#Detection model: steps as in Hefley`s code`
Detection.model=glm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    ,family= "binomial"
                    , data=Detection.data)

unclass(summary(Detection.model))

#Plot models.

myPred0 = predict(habitat.rr, Detection.model, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred0, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main=" detection koala")






#### Sstep 8 Improve the model autocorrelation check and residual approach.#####
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
###
myStack = stack(habitat.rr)
# Plot the predictors
#plot(myStack)
# Creates a mask layer
myMask = myStack$Annual_Precipitation >=-2

plot(myMask)
# Map predictions
AR1 = myMask * 0
plot(AR1)
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

# Crase approach to account for SA
###
# Build the standard GLM object
source("autoregresive2.r")
Detection.data$AR1 = myLib.AutoRegressiveMean(residuals(Detection.model), Detection.data$x, Detection.data$y, 18000)
AR1

Detection.model.2=glm(presence~ distance_primaryandlink
                      +distance_motorwayandlink
                      + fpcnew
                      +AR1 ,family= "binomial"
                      ,data=Detection.data)
summary(Detection.model.2)


#plot(Detection.model.2)

Detection.data$R1 = residuals(Detection.model.2)
Corr <- correlog(Detection.data$x, Detection.data$y, Detection.data$R1,na.rm=T, increment=18000, resamp=0, latlon = F)              

plot(Corr) 
# Estimates DAIC for each predictor
####
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
####
# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"


myPred2 = predict(myARStack, Detection.model.2, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred2, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main=" detection koala")

# ROC curve & GOF metrics
myPred = prediction(predict(Detection.model.2, type = "response"), Detection.data$presence)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(Detection.model.2, type = "response"), Detection.data$presence, 0.5)
myPredPerfs

###

###
# three  additional varibales, res, AR, R1 created for Detection model.2.  hence, they area not available in IPPdata, ZTGLM.data...
# t=extract(myPred2, ZTGLM.data$x &ZTGLM.data$y)
# To=ilogit(t)

#####Step 9: Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 

# set very small valuves to 0.1
p.det <- ifelse(p.det<0.1,0.1,p.det)

hist(p.det, breaks=70)


# p.det 1 or very low valuves create convergencce issues.
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det

######Step 10: - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected=glm(presence~Annual_Mean_Temperature
                  + Annual_Precipitation
                  + Max_Temperature_of_Warmest_Month
                  + Min_Temperature_of_Coldest_Month
                  + fpcnew
                  ,weights=(1/p.det)*10000^(1-presence),family="binomial",  data=IPP.data)


summary(IPP.corrected)
####Step 11: Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####

#use only the significant covariates, tpo +hpop+lot_density+sbd
ZTGLM.corrected=vglm(presence~ distance_primaryandlink
                     +distance_motorwayandlink
                     + fpcnew
                     ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)

# step 7:  Detecion model Map predictions # 

myPred = predict(habitat.rr, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")

#        IPP ignored model Map predictions

myPred2 = predict(habitat.rr, IPP.ignored, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model |detection ignored- number of koalas")

#        IPP corrected model Map predictions

myPred2.1 = predict(habitat.rr, IPP.corrected, type = "response")
plot(myPred2.1, xlab = "x", ylab= "y",main=" IPP model |detection considered- number of koalas")


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

####### step 12: Two phase algoritham#######################################################################

#### chnage below variable list based on s=identified variables in each model.




set.seed(1234)
tpnbs=function()	{
  bss=resample(1:dim(ZTGLM.data)[1])
  IPP.data.bss=rbind(IPP.data[bss,1:12],IPP.data[which(IPP.data$presence==0),1:12])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(presence ~ distance_primaryandlink+distance_motorwayandlink+fpcnew + scale(group),family="binomial",data=Detection.data.bss)
  p.det.bss=ilogit(predict(Detection.model, new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(ZTGLM.myFD3))) 
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  
  IPP.model= glm(presence~ Annual_Mean_Temperature+Annual_Precipitation+ Max_Temperature_of_Warmest_Month
                 +Min_Temperature_of_Coldest_Month+fpcnew, family = "binomial", weights = (1/p.det)*10000^(1-presence), data = IPP.data.bss)
  
  ZTGLM.corrected=vglm(group~distance_primaryandlink+distance_motorwayandlink+fpcnew,weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
  
  c(coef(IPP.model),coef(ZTGLM.corrected))
}

#### step 9: boostrapping 
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps)*tpnbs() # Ravi changed the moasic::do and tpnbs()
bootstrap.sample <- data.frame(bootstrap.sample)
#bootstrap.sample=as.data.frame (bootstrap.sample) #Ravi created a matrix from the dataframe
######## each step gives mean of coefficients for intercept or covariates

#### step 10:  Get means and sd 
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

#### Step 1. read function preapre data####
source("functions2.r")

## recall rasters being used. "habitat.rr"
## we read rasters directly from the folder.

x.files=c("aspect91.tif"
          ,"awc.tif"
          ,"elev.tif"
          ,"fpcnew.tif"
          ,"habit1decimal.tif"
          ,"habit2decimal.tif"
          ,"habit3decimal.tif"
          ,"hpop.tif"
          ,"nitro.tif"
          ,"Precipitation_of_Driest_Quarter.tif"
          , "Precipitation_Seasonality.tif"
          ,"sbd.tif"
          ,"tpo.tif"
          ,"twi.tif"
          ,"roads_other.tif")

d <- stack(x.files)

s.occupancy <- scale(d,scale=TRUE,center = TRUE)
plot(s.occupancy)


w.files=c("distance_primaryandlink.tif", "distance_tertiaryandlink.tif") 

dd <- stack(w.files)

s.detection <- scale(dd, scale=TRUE,center = TRUE)
plot(s.detection)

### Emanuele data.
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Emanuele_data//koala1997_2013month.RData")
myalldata <- koala1997_2013
myalldata <- myalldata[c(1:3)]
names(myalldata )[1]<-"X"
names(myalldata )[2]<-"Y"
selected.locations <- myalldata
#selected.locations= subset(myalldata,yearnew >2010 & yearnew < 2013)
selected.locations = selected.locations[,1:2]
### remove duplicate
selected.locations <- selected.locations[duplicated(selected.locations$X & selected.locations$Y), ]
selected.locations=SpatialPoints(selected.locations)
selected.locations= na.omit(selected.locations) 

############ Create araster:
mydata.r <- raster("elev.tif")
mydata.r <-  subset(habitat.rr, c(1))
res(mydata.r) <- 500
#extent <- extend(selected.locations, extent(selected.locations) + 3000)
# Sample points from within each grid cell:
selected.locations <- gridSample(selected.locations, mydata.r, n = 1)
plot(selected.locations)
selected.locations <- as.data.frame(selected.locations)

## set extent 

#### set extent 
myextent=c(387900, 553100, 6862400, 7113600) 
#crop base on extent of koala locations.
s.occupancy=crop(s.occupancy,myextent)
s.detection= crop(s.detection,myextent)

plot(s.occupancy)
plot(s.detection)
###

pb=selected.locations

pb.loc=SpatialPoints(pb) # over 6629 beyond study area. 
# get locations over rasters only.
pb.occupancy=extract(s.occupancy,pb.loc) 
pb.detection=extract(s.detection,pb.loc) 
# only data from study area.
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,] #3663#  distance covariates extract from rasters. only one here
pb.occupancy=pb.occupancy[is.complete.pb,]#3663# env covariates extract from  raster


#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection)) #s.detection is raster stack.

# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell


s.area=area.back*nrow(X.back) #study area



# adding column of ones - po locations # Is this pb.detection?????
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # pb.occupancy is all presence locations. 
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)


# Step 2. Analising the data ========================================================

#Analyzing Presence-Only data
(pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)) 
#estimates are obtained but not se.
# If any of the eigenvalues are zero or negative, there is obviously a problem, 
# perhaps stemming from under identification or boundary conditions
# One thing you might try is to refit the model from its solution, to see 
# if the gradients get smaller, and the NaN's clear up.
#
### Step 3. Trasfom data to et predictions.####
coef <- pb.fit$coefs
coeff <- coef[c(1,2)]

#occupancy rasters

f <- s.occupancy

# get predictions for each raster

f2 <- subset(f,c(1))*coeff$Value[[2]]
f3 <- subset(f,c(2))*coeff$Value[[3]]
f4 <- subset(f,c(3))*coeff$Value[[4]]
f5 <- subset(f,c(4))*coeff$Value[[5]]
f6 <- subset(f,c(5))*coeff$Value[[6]]



fn <- exp(f2+f3+f4+f5+f6+coeff$Value[[1]])

plot(fn)

plot(fn, main= "intensity 2000 data")


#detection rasters

ff <- s.detection

# get predictions for each raster

ff13 <- subset(ff,c(1))*coeff$Value[[7]] 
ff14 <- subset(ff,c(2))*coeff$Value[[8]]
fn2 <-  (ff13+ff14)+coeff$Value[[7]]   

fnn <- logit.prob(fn2)

plot(fnn, main =" probability of detection")

#final intensity

fn3 <- fn*fnn

plot(fn3, main = "bias corrected intensity")


# plot intensity baed on ecological varibales and detection based on observer bias vribales.

par(mfrow=c(1,3))

plot(fn3); plot(fnn) ; plot(fn)


###### Standard erros

coef <- pb.fit$coefs
se<- coef[c(1,3)]

#occupancy rasters

f <- s.occupancy

# get predictions for each raster

fs2 <- subset(f,c(1))*se$`Standard error`[[2]]
fs3 <- subset(f,c(2))*se$`Standard error`[[3]]
fs4 <- subset(f,c(3))*se$`Standard error`[[4]]
fs5 <- subset(f,c(4))*se$`Standard error`[[5]]
fs6 <- subset(f,c(5))*se$`Standard error`[[6]]



fns <- exp(fs2+fs3+fs4+fs5+fs6+coeff$Value[[1]])

plot(fns)

plot(fns, main= "Standard error")


#detection rasters

ff <- s.detection

# get predictions for each raster

ffs13 <- subset(ff,c(1))**se$`Standard error`[[7]] 
ffs14 <- subset(ff,c(2))**se$`Standard error`[[8]]
fns2 <-  (ffs13+ffs14)+coeff$Value[[7]]   

fnns <- logit.prob(fns2)

plot(fnns, main =" probability of detection")

#final intensity

fns3 <- fns*fnns

plot(fns3, main = "bias corrected intensity")


# plot intensity baed on ecological varibales and detection based on observer bias vribales.

par(mfrow=c(1,3))

plot(fns3); plot(fnns) ; plot(fns)




save.image(file = "models_ravi.Rdata")
