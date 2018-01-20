library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(lmtest)
library(spatial.tools)
library(VGAM)
library(mosaic)
library(faraway)
library(gstat)  #
library(ncf)    #
library(foreign)
library(nlme)   
library(MASS)
library(ROCR)
library(vcd)
library(RColorBrewer) # 
library(classInt)
library(ppmlasso)
library(usdm)

##############
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]

# Prepare data @km.
all.locations= subset(mydata,yearnew >1998:2015, select=X:Y)
# all.locations=all.locations/1000
# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 100)
selected.locations = select.locations[,1:2]

###combine wildnet data 
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata.RData")
wildnet <- wildnetdata[c("X","Y")] 
mydata=rbind(selected.locations, wildnet)


#select koalas in the study region
(b=SpatialPoints(all.locations))#chextent then keep a small buffer area.
selected.locations2<- subset(selected.locations, X > 401117 & X < 554668)
selected.locations <- subset(selected.locations2, Y > 6883463 & Y < 7084959) # xy only within the study area.
#points(selected.locations$x, selected.locations$y, pch=15, col="red")


#prepare rasters removing na in rivers for better maps
# myfullstack.a <- list.files(pattern="\\.tif$")
# myfullstack = scale(stack(myfullstack.a))
# #rasters have na valuves along brisbane river. change this to 0.
# myfullstack[is.na(myfullstack[])] <- 0
# large_studyarea.shp <- readShapePoly("large_studyarea.shp")
# myfullstack=mask(myfullstack,large_studyarea.shp)
----------------

## now write all these rasters to a new folder with layer name intact.
# ####writeRaster((myfullstack+ mask), names((myfullstack+ mask)), bylayer=TRUE, driver='GTIFF')
----------------------
#myfullstack=crop(myfullstack,lga10.shp)
# extent(myfullstack) <- extent(c(xmin(myfullstack), xmax(myfullstack), ymin(myfullstack), ymax(myfullstack))/1000)
habitat.r<- subset(myfullstack, c(1,2,4,6,13,16,23, 24, 25,26,27,28,29,39,40,41))
habitat.r=crop(habitat.r,b)


#plot(habitat.r,1)
#Quadrature points. stt[is.na(stt)] <- 0
bigquad <- as.data.frame(habitat.r, xy=TRUE,na.rm=T)
bigquadxy=bigquad[c(1,2)]
fpcnew=myfullstack$fpcnew
buff=extract(fpcnew,bigquadxy, fun=max,buffer=2,df=TRUE)
#rename this vvariable as buffer
colnames(buff)[2] <- "fpcnew_buff"
buff=buff[c(2)]
bigquad.t=cbind(buff,bigquad)
bigquad=bigquad.t[c(2:19,1)]
##### in kilometers.
bigquad =cbind((bigquad[,1:2]/(1000)),bigquad[c(-1,-2)])
#loc=SpatialPoints(bigquad)
#koala dataxy
sp.xy = data.frame(selected.locations)
sp.xy=(sp.xy[,1:2]/(1000))

sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

##### to predict using model based control of observer bias based on same level of distance
bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 
bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)
#stt <- na.omit(stt)
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'
# stt[is.na(stt)] <- 0
xydatan <- bigquad[c(1,2)]
#  xy requires as integers otherwise takes years to run the model..
xydata <- as.data.frame(lapply(xydatan, as.integer)) # stt[] <- lapply(stt, as.integer)#this line edited on 04/01# make only xy integer in line with dadta shared with Mark S. 
dd=bigquad[c(-1,-2)]
bigquad <- cbind(xydata,dd )

ddd=bigquad.2[c(-1,-2)]
bigquad.2 <- cbind(xydata,ddd )

xydatan <- bigquad[c(1,2)]
  

#ppm.1 = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)
ppm.1 = ~ poly(AnnualMeanTemperature,habit1decimal,fpcnew_buff,AnnualPrecipitation,degree = 1, raw = TRUE)
scales = c( 0.5, 1, 2, 4, 8,16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.1)

#Model1:
ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
diagnose.ppmlasso(ppmFit.1)
resid.plot = diagnose(ppmFit.1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")

##########
opar <- par() #make a copy of current settings
#par(opar)          # restore original settings
mypar <- par(mar=c(1.5,1.5,1.5,1.5), oma=c(0.5,0.5,0.5,0.5))

#Predict and plot
pred.fit.e = predict.ppmlasso(ppmFit.1, newdata=bigquad)

predictions.fit.e <- cbind(xydata, pred.fit.e) # xydatan was chnaged to xydata.
pred.final0.e<- rasterFromXYZ(as.data.frame(predictions.fit.e )[, c("X", "Y", "pred.fit.e")])
pred.final0.e
plot(pred.final0.e, main=" koala density-warton method/ env only")

##### check total koal in the area
(total_koala=cellStats(pred.final0.e,sum))
lga10.shp <- readShapePoly("LGA10new.shp")
plot(lga10.shp, add=TRUE)
pred.crop=crop(pred.final0.e,lga10.shp)
e=mask(pred.crop,lga10.shp)
plot(e)
plot(lga10.shp, add=TRUE)
points(selected.locations$x, selected.locations$y, pch=15, col="red")

#Model :2a

#### Env and distance both
ppm.form = ~ poly(awc,clay,elev,fpcnew, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi, habit1decimal,degree = 2, raw = TRUE)+poly(distance_primaryandlink,distance_motorwayandlink,degree = 2, raw = TRUE)


#4.2 Fitting a regularisation path of point process models
#a LASSO penalty that optimises non-linear GCV
ppmFit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
diagnose.ppmlasso(ppmFit)

#bias not corrected
pred.fit.correctn = predict.ppmlasso(ppmFit, newdata=bigquad)
predictions.fit.correctn <- cbind(xydata, pred.fit.correctn) # xydatan was chnaged to xydata.
pred.final0.correctn<- rasterFromXYZ(as.data.frame(predictions.fit.correctn )[, c("X", "Y", "pred.fit.correctn")])
plot(pred.final0.correctn, main=" koala density-warton method/ env_dis")
(total_koala=cellStats(pred.final0.correctn,sum))

#overlay LGA
plot(lga10.shp, add=TRUE)
pred.crop=crop(pred.final0.correctn,lga10.shp)
e=mask(pred.crop,lga10.shp)
plot(e)
plot(lga10.shp, add=TRUE)

#Model:2b Bias corrected.
# now correct for bias.
pred.fit.correct = predict.ppmlasso(ppmFit, newdata=bigquad.2)
predictions.fit.correct <- cbind(xydata, pred.fit.correct) # xydatan was chnaged to xydata.
pred.final0.correct<- rasterFromXYZ(as.data.frame(predictions.fit.correct )[, c("X", "Y", "pred.fit.correct")])
plot(pred.final0.correct, main=" koala density-warton method/ bias corrected")
(total_koala=cellStats(pred.final0.correct,sum))
#plot all with same legend
par(mfrow=c(2,2),oma=c(1,1,1,1))
plot(pred.final0.e, zlim = c(0, 5),main="env only")

plot(pred.final0.correctn, zlim = c(0, 5),main="env_dis")

plot(pred.final0.correct,zlim = c(0, 5), main="bias corrected")
plot(lga10.shp, add=TRUE)
#overlay LGA
plot(lga10.shp, add=TRUE)
pred.crop=crop(pred.final0.correct,lga10.shp)
e=mask(pred.crop,lga10.shp)
plot(e)
plot(lga10.shp, add=TRUE)

### residuals:
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit)
#K-envelop
kenv = envelope(ppmFit, fun = Kinhom) # simulated envelop for summary function 
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")


#######
##lurking varibale plots and residulas seperatly.
#diagnostic plots:residuals and lurking variable plots.
diagnose.ppmlasso(ppmFit)
diagnose.ppmlasso(ppmFit, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit, which = "y", type = "Pearson", compute.sd = TRUE)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas biascorrected interaction model")

# area interaction model with details.
#family = "area.inter", r = 2
### another way to fit area interaction model. can add: criterion = "blockCV", n.blocks = 5, block.size = 10 
final.fit = ppmlasso(ppmFit,sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 5,
                      block.size = 10,sp.scale = 1, n.fits = 100, family = "area.inter", r = 1)
diagnose(final.fit, which = "smooth", type = "Pearson")
pred.interaction = predict(final.fit, newdata=bigquad.2)
pred.inter.action <- cbind(xydata, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter)
diagnose.ppmlasso(final.fit)
(total_koala=cellStats(pred.ct.inter,sum))

opar






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
###### runs this code without error.

### another way to fit area interaction model. can add: criterion = "blockCV", n.blocks = 5, block.size = 10 
final.fita = ppmlasso(final.fit,sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 5,
                      block.size = 10,sp.scale = 1, n.fits = 100, family = "area.inter", r = .5)
diagnose(final.fita, which = "smooth", type = "Pearson")
pred.interaction = predict(final.fita, newdata=bigquad.2)
pred.inter.action <- cbind(xydata, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter )
diagnose.ppmlasso(final.fita)









