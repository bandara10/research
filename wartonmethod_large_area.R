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
all.locations=all.locations/1000

# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > .5)
selected.locations = select.locations[,1:2]

#select koalas in the study region
(b=SpatialPoints(selected.locations))#chextent then keep a small buffer area.
selected.locations2<- subset(selected.locations, X > 438.9503 & X < 553.9503)
selected.locations <- subset(selected.locations2, Y > 6873.579 & Y < 7108.579) # xy only within the study area.
points(selected.locations$x, selected.locations$y, pch=15, col="red")

#prepare rasters,
myfullstack.a <- list.files(pattern="\\.tif$")
myfullstack = scale(stack(myfullstack.a)) 
extent(myfullstack) <- extent(c(xmin(myfullstack), xmax(myfullstack), ymin(myfullstack), ymax(myfullstack))/1000)
habitat.r<- subset(myfullstack, c(1,2,4,6,13,16,23, 24, 25,26,27,28,29,39,40,41))

#plot(habitat.r,1)
#Quadrature points.
bigquad <- as.data.frame(habitat.r, xy=TRUE,na.rm=T)

#koala dataxy
sp.xy = data.frame(selected.locations)
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

##### to predict using model based control of observer bias based on same level of distance
bigquad[, c(1,2)] <- sapply(bigquad[, c(1,2)], as.integer)
bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 
bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)
# stt[is.na(stt)] <- 0
xydatan <- bigquad[c(1,2)]
# stt requires xy as integers.
xydata <- as.data.frame(lapply(xydatan, as.integer))  


  #gt species selected data as a dataframe



#ppm.form.e = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)
ppm.form.e = ~ poly(AnnualMeanTemperature,habit1decimal,fpcnew,AnnualPrecipitation,degree = 1, raw = TRUE)
scales = c( 0.5, 1, 2, 4, 8,16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form.e)

ppmFit.e = ppmlasso(ppm.form.e, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
diagnose.ppmlasso(ppmFit.e)
resid.plot = diagnose(ppmFit.e, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")

##########
opar <- par() #make a copy of current settings
#par(opar)          # restore original settings
mypar <- par(mar=c(1.5,1.5,1.5,1.5), oma=c(0.5,0.5,0.5,0.5))
###########

#Predict and plot
pred.fit.e = predict.ppmlasso(ppmFit.e, newdata=bigquad)

predictions.fit.e <- cbind(xydata, pred.fit.e) # xydatan was chnaged to xydata.
pred.final0.e<- rasterFromXYZ(as.data.frame(predictions.fit.e )[, c("X", "Y", "pred.fit.e")])
pred.final0.e
extent(pred.final0.e) <- extent(c(xmin(pred.final0.e), xmax(pred.final0.e), ymin(pred.final0.e), ymax(pred.final0.e))*1000)
plot(pred.final0.e, main=" koala density-warton method/ env only")

lga10.shp <- readShapePoly("LGA10new.shp")
plot(lga10.shp, add=TRUE)
pred.crop=crop(pred.final0.e,lga10.shp)
e=mask(pred.crop,lga10.shp)
plot(e)
plot(lga10.shp, add=TRUE)
points(selected.locations$x, selected.locations$y, pch=15, col="red")
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
plot(pred.final0.correctn, main=" koala density-warton method/ bias corrected")

extent(pred.final0.correctn) <- extent(c(xmin(pred.final0.correctn), xmax(pred.final0.correctn), ymin(pred.final0.correctn), ymax(pred.final0.correctn))*1000)
plot(pred.final0.correctn, main=" koala density-warton method/ env_dis")

#overlay LGA
plot(lga10.shp, add=TRUE)
pred.crop=crop(pred.final0.correctn,lga10.shp)
e=mask(pred.crop,lga10.shp)
plot(e)
plot(lga10.shp, add=TRUE)

# now correct for bias.
pred.fit.correct = predict.ppmlasso(ppmFit, newdata=bigquad.2)
predictions.fit.correct <- cbind(xydata, pred.fit.correct) # xydatan was chnaged to xydata.
pred.final0.correct<- rasterFromXYZ(as.data.frame(predictions.fit.correct )[, c("X", "Y", "pred.fit.correct")])
plot(pred.final0.correct, main=" koala density-warton method/ bias corrected")
#overlay LGA
extent(pred.final0.correct) <- extent(c(xmin(pred.final0.correct), xmax(pred.final0.correct), ymin(pred.final0.correct), ymax(pred.final0.correct))*1000)
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
kenv = envelope(ppmFit, fun = Kinhom, nsim=39) # simulated envelop for summary function 
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")


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









