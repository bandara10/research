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
##########Step 1: load koala data frm tow sources. BoalaBASE and Wildnet#####
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
all.locations= subset(mydata,yearnew >1998:2015, select=X:Y)

source("Lib_DistEstimatesToItself.r")## set a minimum distance betweek koalas
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 300)
selected.locations = select.locations[,1:2]

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata.RData") # wildnet data.
wildnet <- wildnetdata[c("X","Y")] 
mydata=rbind(selected.locations, wildnet)

######Step 2: load rasters and crop to my extent of interest######
myfullstack.a <- list.files(pattern="\\.tif$")
myfullstack = scale(stack(myfullstack.a))

myextent=c(387900, 553100, 6862400, 7113600) 
myfullstack=crop(myfullstack, myextent, snap="near") ###crop to my area of interest
# ########
# myfullstack.a <- list.files("path\\warton_data_allclimate",pattern="\\.tif$")
# myfullstack = scale(stack(myfullstack.a))

habitat.r<- subset(myfullstack, c(1,2,4,6,13,16,23, 24, 25,26,27,28,29,39,40,41)) # select my variables

######Step 3. Create quadrature points and species xy data. Remove NA. Keep xy colomns for preditions and rasterize. ######
bigquad <- as.data.frame(habitat.r, xy=TRUE,na.rm=T)

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
sp.xy = data.frame(selected.locations)
sp.xy=(sp.xy[,1:2]/(1000))
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

#####step 5: Models: Enviromental only. Test using different env variables.
#                   Enviromental and distance variables
#                   Bias corrected.

#Model 1:  Enviromental varibales only
#ppm.1 = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)
ppm.1 = ~ poly(clay,sbd,tpo,AnnualMeanTemperature,awc,habit1decimal,habit2decimal,habit3decimal,hpop,lot_density ,fpcnew,AnnualPrecipitation,degree = 1, raw = TRUE)

scales = c( 0.5, 1, 2, 4, 8,16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.1)#find the spatial resolution best for analysis

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, criterion = "nlgcv",n.fits = 100)#Fitting a regularisation path of point process models.#a LASSO penalty that optimises non-linear GCV
#           
ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 50,block.size = 50, sp.scale = 1, n.fits = 100)

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 5,
                    block.size = 10, sp.scale = 1, n.fits = 50, family = "area.inter", r = 1)

#             Check residuals
diagnose.ppmlasso(ppmFit.1)
resid.plot = diagnose(ppmFit.1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.1, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.1, which = "y", type = "Pearson", compute.sd = TRUE)

#             Check level of interaction
kenv = envelope(ppmFit.1,fun = Kinhom) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

#Predict and plot
pred.fit.e = predict.ppmlasso(ppmFit.1, newdata=bigquad) #get prdictions.
predi.fit.e <- cbind(xydata, pred.fit.e)     #predictions anad xy locations.
pred.model.1<- rasterFromXYZ(as.data.frame(predi.fit.e )[, c("x", "y", "pred.fit.e")])
plot(pred.model.1, main=" koala density-warton method/ env only")

(total_koala=cellStats(pred.model.1,sum))        #check total koal in the area


#Model :2 Enviromental and distance covariate which are of original scale (bigquad).

ppm.2 = ~ poly(AnnualMeanTemperature,habit1decimal,hpop,lot_density ,fpcnew,AnnualPrecipitation,degree = 2, 
                  raw = TRUE)+poly(distance_primaryandlink,distance_motorwayandlink,degree = 2, raw = TRUE)

ppmFit.2 = ppmlasso(ppm.2, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 50)
diagnose.ppmlasso(ppmFit.2)
#             Check residuals
resid.plot = diagnose(ppmFit.2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.2, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.2, which = "y", type = "Pearson", compute.sd = TRUE)

#             check level of interaction


kenv = envelope(ppmFit.2, fun = Kest) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

pred.fit.dis = predict.ppmlasso(ppmFit.2, newdata=bigquad)
predi.fit.dis <- cbind(xydata, pred.fit.dis) 
pred.model.2<- rasterFromXYZ(as.data.frame(predi.fit.dis)[, c("x", "y", "pred.fit.dis")])
plot(pred.model.2, main=" koala density-warton method/ env_dis")
(total_koala=cellStats(pred.model.2,sum))


# Bias correction step:now predict using quadrature (biguad.2) with minimum distances 
pred.fit.crt = predict.ppmlasso(ppmFit.2, newdata=bigquad.2)
predi.fit.crt <- cbind(xydata, pred.fit.crt) # xydatan was chnaged to xydata.
pred.model.crt  <- rasterFromXYZ(as.data.frame(predi.fit.crt )[, c("x", "y", "pred.fit.crt")])
plot(pred.model.crt , main=" koala density-warton method/ bias corrected")
(total_koala=cellStats(pred.model.crt, sum))

#####step 6.#plot all with same legend######
par(mfrow=c(2,2),oma=c(1,1,1,1))
plot(pred.model.1, zlim = c(0, 5),main="env only")
plot(pred.model.2, zlim = c(0, 5),main="env_dis")
plot(pred.model.crt, zlim = c(0, 5), main="bias corrected")
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



