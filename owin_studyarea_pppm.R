#### Load libraries ####
library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatstat.utils)

###Objectives ####
# 1. Identify which variables to select for the model by Rhohat
# 2. Estimate density Kernel smoothing with Diggle's edge correction
# 3. PPM without interaction term
# 4. PPM with interaction term Geyer(r = 5000, sat = 20)|nteraction = Strauss(r = 6000)
# 5. Density from Adaptive kernel smoothing using sparr.dividing surfaces to make relative risk density map:
# 5.1 Calculate a pilot bandwidth (using leave-one-out least-squares cross-validation) and a global bandwidth
# 5.2 (using maximal smoothing) and then feed these values to the bivariate.density function:

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")
##########Step 1: load koala data from BoalaBASE, Wildnet and Gold coast#####
## step 1 and two are common to all methods.

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata=mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]

# Get wildnet data

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata_old.RData") # wildnet data.
mydata2 <- wildnetdata[c("X","Y","yearnew")]
# 
# mydata3=rbind(mydata, mydata2)

# Get gold coast data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//koala_gold_coast.RData")
mydata3 <- koala_gold_coast[c("X","Y","yearnew")]
myalldata=rbind(mydata, mydata2,mydata3)
table(myalldata$yearnew)

#myalldata[!duplicated(myalldata$X &myalldata$Y), ]
all.locations= subset(myalldata,yearnew >2005,select=X:Y) #, select=X:Y REmove xy for ppp density
all.locations <- unique(all.locations)

pp=SpatialPoints(all.locations)
plot(pp)

#########=======================
## for marks in ppp density only
selected.locations= all.locations
selected.locations=unique(selected.locations)
##### Stratified sampling #####
# #We split a data.frame into color groups. From each such a group, we sample 200 rows.
# df2 <- lapply(split(myalldata, myalldata$yearnew),
#               function(subdf) subdf[sample(1:nrow(subdf), 200),])
# 
# d=do.call('rbind', df2) #merged into 1 data.frame

####### Select records based on distance###### distance points are only selected.
###better way is to select one from grid.

# source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas
# all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
# select.locations = subset(all.locations, all.locations$disttoitself > 1000)#2000
# selected.locations = select.locations[,1:2]
# # # selected.locations$yearnew=as.factor(selected.locations$yearnew)
# # selected.locations=SpatialPoints(selected.locations)
# # proj4string(selected.locations) <- CRS("+init=epsg:28356") 

######### Bring in MGA56 square study area boundary map:######

aus.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\australia.shp")#AU_Qld_detail_studyarea_outline-MGA56
proj4string(aus.shp) <- CRS("+init=epsg:4326") 
aus.shp <- spTransform(aus.shp, crs(studyarea.shp))
plot(aus.shp)

studyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 

# detail study area

qld.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\qld.shp")#AU_Qld_detail_studyarea_outline-MGA56
proj4string(qld.shp) <- CRS("+init=epsg:4326") 
qld.shp <- spTransform(qld.shp, crs(studyarea.shp))

#### selected area for model
dstudyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_detail_studyarea_outline-MGA56.shp")#
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 


plot(qld.shp, axes = TRUE)
plot(dstudyarea.shp, axes = TRUE, col="blue")
plot(studyarea.shp, add = TRUE,col="yellow")
points(x = selected.locations$X, y = selected.locations$Y)

#points(x=all.locations$X, y= all.locations$Y)
#### select one from grid.

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = stack(myfullstack.b)
names(myfullstack.b)
#plot(myfullstack.b,1)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.b, myextent, snap="near")
habitat.rr <- scale(habitat.rr)
library(usdm)
vifstep(habitat.rr, th=10)
# mydata.r <- raster(studyarea.shp)
# res(mydata.r) <- 1000
# mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset
selected.locations=gridSample(all.locations, habitat.rr, n = 2)# disaggregate a raster factor 4 and get one point.
selected.locations <- as.data.frame(selected.locations)
plot(selected.locations)
ss=SpatialPoints(selected.locations)
###
# Make a ppp object for spatstat:=================================================================================================================================
# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:

source("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\owin2sp_source.r")

# studyarea.shp= readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp") # warton data folder.
# #dstudyarea.shp <- dstudyarea.shp[, -(1:3)] 

dat.w <- as(as(studyarea.shp, "SpatialPolygons"), "owin") # warton data climate
dat.spw <- owin2SP(dat.w)

# Set projection of the owin as GDA94 / SA Lambert:
proj4string(dat.spw) <- CRS("+init=epsg:28356") 
dat.spw <- raster::intersect(x = dstudyarea.shp, y = dat.spw)

###### Convert the sp object into a window:#####
dat.w <- as.owin(dat.spw)
plot(dat.w)

# Select only those koalas within the owin:
id <- inside.owin(x = selected.locations$X, y = selected.locations$Y, w = dat.w)
selected.locations <- selected.locations[id,]

# acsel=selected.locations
plot(studyarea.shp, axes = TRUE)
points(x = selected.locations[,1], y = selected.locations[,2])


# Make a ppp object:
dat.ppp <- ppp(x = selected.locations[,1], y = selected.locations[,2], window = dat.w)

unitname(dat.ppp) <- c("meter","meter")
plot(dat.ppp, axes = TRUE)

### add marks and plot density# select year..#####
marks(dat.ppp) <- selected.locations$yearnew
summary(dat.ppp)


#split based on year whihc is a factor
split_dat.ppp <- split(dat.ppp)
plot(split(dat.ppp))
#as.matrix(lapply(split_dat.ppp,npoints),ncol=1)
plot(density(split(dat.ppp)), ribbon = TRUE)

# b=density(split(dat.ppp), ribbon = TRUE)
# 
# plot(rhohat(dat.ppp), b)

#estimation of the density
dens_all <- density(split_dat.ppp)
plot(dens_all*1E06)

contour(density(split_dat.ppp), axes = F)


dat.ppp <- unmark(dat.ppp)
densit <- density.ppp(dat.ppp)
plot(densit*1E06)

####use the density map to look at how the intensity of events (number of points per unit area) relates to the density of sighting image

plot(rhohat(unmark(dat.ppp),densit ))
#Quadrat test
#p-value < 2.2e-16| null hypothesis of the point pattern being generated by complete spatial random process is rejected, 
# we have some evidence that the point pattern is inhomogenous or nonstationary.
quadrat.test(dat.ppp)

#If the point pattern follow Complete Spatial Randomness (CSR) then there is a
#known relationship between this count number (K) and the distance considered (r).

en=envelope(dat.ppp,fun=Kest,funargs=list(correction="border"),global=TRUE)
plot(en)
####points are more dispersed than expected
#by using Kest we assume that the point pattern was generated by one homogenous 
#intensity function characterized by the average intensity of our point pattern
###The way Kinhom works is by deriving an intensity estimate from the data (similar to density.ppp) 
##and by weighting each point based on their estimated intensity
ee_inhom1 <- envelope(dat.ppp,fun=Kinhom,global = TRUE)
plot(ee_inhom)
### Bring in rasters #####

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = stack(myfullstack.b)
names(myfullstack.b)
#plot(myfullstack.b,1)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.b, myextent, snap="near")
habitat.rr <- scale(habitat.rr)
#habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,33,38)) # select my variables
#plot(habitat.rr,2)


# Convert rasters into a spatstat image object=================================================================================================================================
##Rhohat fpcnew ####
fpc.im <- as.im(habitat.rr$fpcnew )
amt.im <- as.im(habitat.rr$AnnualMeanTemperature)
apt.im <- as.im(habitat.rr$AnnualPrecipitation)
mtwm.im <- as.im(habitat.rr$Max_Temperature_of_Warmest_Month)
mtcm.im <- as.im(habitat.rr$Min_Temperature_of_Coldest_Month)
pdm.im <- as.im(habitat.rr$Precipitation_of_Driest_Month)
pwm.im <- as.im(habitat.rr$Precipitation_of_Wettest_Month)
awc.im <- as.im(habitat.rr$awc)
elev.im <- as.im(habitat.rr$elev)
clay <- as.im(habitat.rr$clay)
nitro.im  <- as.im(habitat.rr$nitro)
citydis.im <- as.im(habitat.rr$city_dis)
sbd.im  <- as.im(habitat.rr$sbd)
roughness.im <- as.im(habitat.rr$roughness)
tpo.im  <- as.im(habitat.rr$tpo)
twi.im  <- as.im(habitat.rr$twi)
hpop.im <- as.im(habitat.rr$hpop)
habit3.im  <- as.im(habitat.rr$habit3decimal)
habit2.im  <- as.im(habitat.rr$habit2decimal)
habit1.im  <- as.im(habitat.rr$habit1decimal)
habitdis2.im <- as.im(habitat.rr$Dis_habitat_suitable_2)
roadsM.im  <- as.im(habitat.rr$roads_motor)
roadsother.im  <- as.im(habitat.rr$roads_other)
dmwl.im <- as.im(habitat.rr$distance_motorwayandlink)
dprl.im <- as.im(habitat.rr$distance_primaryandlink)
#Poisson point process model =================================================================================================================================
predList=list(fpc = fpc.im
     ,mtcm = mtcm.im
     ,mtwm = mtwm.im 
     ,pdm=pdm.im
     ,pwm=pwm.im
     ,apt = apt.im
     ,awc =awc.im
     ,elev= elev.im
     , roughness= roughness.im
     ,nitro= nitro.im
     ,sbd= sbd.im
     ,tpo= tpo.im
     ,twi = twi.im
     ,hpop= hpop.im
     ,citydis= citydis.im
     ,habit3 = habit3.im
     ,habit2 = habit2.im
     ,habit1 = habit1.im
     ,habitdis2 = habitdis2.im
     ,roadsM = roadsM.im
     ,roadsother= roadsother.im
     ,dprl=dprl.im
     ,dmwl= dmwl.im)
# Null model:
(dat.ppm00 <- ppm(dat.ppp, trend = ~ 1)) #,covariates = predList
    diagnose.ppm(dat.ppm00)                                                   
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409

######dat.ppm01 <- ppm(dat.ppp ~ polynom(elev,mtwm,2),covariates = predList)


###### Extract the quadrature scheme from dat.ppm01:#####
Qz <- quad.ppm(dat.ppm00, drop = TRUE)
plot(Qz)


# This step will drop points within the quadrature scheme that had NA-values.

####### Saturated model:#### 
dat.ppm01 <- step(ppm(Qz, trend = ~ roughness+roadsother+fpc+mtcm + mtwm+pwm+pdm
                      + awc +citydis+ elev+ nitro+ sbd + tpo+ twi+ roadsother+ hpop+habitdis2+habit2, covariates = predList)) #+ roadsM + roadsother
# m1 <- ppm(Qz ~ polynom(elev,mtwm ,2),covariates = predList)

diagnose.ppm(dat.ppm01)

eem1 <- envelope(dat.ppp,Kest,funargs = list(lambda=dat.ppm01),global=TRUE)
plot(eem1)

##model summary
summary(dat.ppm01)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm01)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm01) # 3466.415

###plot
pred <- predict(dat.ppm01, type="trend") # ngrids smooth the map.
pred=pred*1E04
pred.grid <- as(pred,"SpatialGridDataFrame")
pred.r=raster(pred.grid)
plot(pred.r)
writeRaster(pred.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred01.tif",overwrite=TRUE)
######==========



plot(pred* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer

pred1 <- pred* 1E06
max(pred1)
coll <- colourmap(terrain.colors(100,1), range=c(0,3.306665))
plot(coll,col.ticks="black")
plot(pred1, col=coll)


#draw a cool perspective map
# 
# pred1 <- pred* 1E06
# persp(pred1,box=FALSE,visible=TRUE)
# persp(pred1, colin=pred1, box=FALSE,visible=TRUE,theta=1,phi=100,expand=2)

diagnose.ppm(dat.ppm01,which = "smooth")### ,which = "smooth"
# use the fitted intensity in the Kinhom function to see if the observed point pattern 
# is more or less clustered than expected from the model fit:

eem <- envelope(dat.ppp,Kinhom,funargs = list(lambda=dat.ppm01),global=TRUE)
# plot(eem)
#clustered poisson point process models
dat.ppm02 <- kppm (Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ roadsother+ hpop+habitdis2 
                             ,covariates = predList) #+ roadsM + roadsother
coef(dat.ppm02)

test=predict.kppm(dat.ppm02,ngrid=200)
test.t=test* 1E06
plot(test.t)



######## write raster to plot in arcgis
test.grid <- as(test.t,"SpatialGridDataFrame")
test.r=raster(test.grid)
writeRaster(test.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/test.tif",overwrite=TRUE)
######==========

plot(dat.ppm02,what="statistic",pause=FALSE)
pred2 <- predict(dat.ppm02, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred* 1E06)


###### Now add a Geyer spatial interaction term:#####
dat.ppm03 <- ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ roadsother+ hpop+habitdis2
                 ,covariates = predList) #interaction = Geyer(r = 2000, sat = 20),interaction = Strauss(r = 1000) ,
diagnose.ppm(dat.ppm03)
summary(dat.ppm03)
coef(dat.ppm03)
summary(dat.ppm03)$coefs.SE.CI

AIC(dat.ppm03) 
plot(dat.ppm03,ngrid=200)
diagnose.ppm(dat.ppm03)


######## QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions##
##re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
dat.qq03 <- qqplot.ppm(fit = dat.ppm07, nsim = 10)
plot(dat.qq03, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# Predict=================================================================================================================================
# Model with spatial interaction term:
pred3 <- predict(dat.ppm03, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred3* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer
persp(pred3)
pred3 <- pred3* 1E06
max(pred3)
coll <- colourmap(terrain.colors(100,1), range=c(0,4.908895))
plot(coll,col.ticks="black")
plot(pred1, col=coll)

#### Express predicted koala intensity in terms of the number of koalas per square km:####
pred3$v <- pred3$v * 1E06
summary(as.vector(pred3$v))
hist(as.vector(pred3$v))
max(pred3)

##### step 2. Model with distance covariates

dat.ppm04 <- step(ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ 
                        habit3 + habit2 + habit1 + roadsM + roadsother+dprl.im+dmwl.im, covariates = predList))
diagnose.ppm(dat.ppm04)
# see if the observed point pattern is more or less clustered than expected from the model fit:


####interaction = Strauss(r = 2000) 
dat.ppm04 <- step(ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ 
                        habit3 + habit2 + habit1 + roadsM + roadsother+dprl.im+dmwl.im,interaction = Strauss(r = 1000), covariates = predList))
diagnose.ppm(dat.ppm04)

summary(dat.ppm04)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm04)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm04) # 3466.415
# Predict=================================================================================================================================
# Model with spatial interaction term:
pred4 <- predict(dat.ppm04, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred4* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer
persp(pred4)
pred4 <- pred4* 1E06
max(pred4)
coll <- colourmap(topo.colors(200),reverse = TRUE, range=c(0,0.01884608))

coll <- colourmap(c("yellow", "blue", "green", "red"), breaks=c(0,0.01,0.02,0.03,0.04))

plot(coll,col.ticks="black")
plot(pred1, col=coll)



###### Map Predictions from the model:#### Another color scheme
breaks <- seq(from = 0, to = 0.1792389, length = 10)
col <- brewer.pal(n = 9, name = "Greens")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred4$xcol, y = pred4$yrow, z = t(pred4$v), zlim = c(0, 0.31), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "black")
plot(studyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, 
      lab = breaks, cols = col, shift = 0, cex = 0.75)




##### Actual kernel smoothed data:#####
###### Work out density for dat.ppp (using Diggle's edge correction):#####
# set the dimensions for the grid and the bandwidth:
dimyx = c(200,200)
bandwidth <- 1500

dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)
plot(dat.den*1E06)
# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre:
dat.den$v <- dat.den$v *1E06
summary(as.vector(dat.den$v))
max(dat.den)
# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]

x.points <- seq(from = 387950.3, to = 552950.3, by = 1000); x.lab <- x.points / 1000
y.points <- seq(from = 6862579, to = 7113579, by = 1000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 9.5, length = 10)
col <- brewer.pal(n = 9, name = "Blues")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


breaks <- seq(from = 0, to = 9.447769, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), zlim = c(0, 0.3007344), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)




 par(pin = c(1 * 5, ratio * 5), omi = c(0.5,0,0,0))
plot(xylims[1,], xylims[2,], type = "n", xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.15), col = col, breaks = breaks, add = TRUE)
points(cas.ppp, col = "red", pch = 16)
plot(bnd.shp, col = "transparent", border = "black", add = TRUE)
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = -25000, yb = 8700000, xr = 0, yt = 9000000, lab = breaks, cols = col, shift = 0, cex = 0.75)
legend(x = "topleft", legend = "FMD-outbreak wards", pch = 16, col = "red", cex = 0.75, bty = "n")
epi.saveplot("ZM_fmd_ppm05_predictions")


#Adaptive kernel smoothing using sparr =================================================================================================================================
# Adaptive smoothing and dividing surfaces to make relative risk density map:

library(spatialkernel); library(sparr); library(rgl); library(spatialkernel)

# Calculate a pilot bandwidth (using leave-one-out least-squares cross-validation) and a global bandwidth (using maximal smoothing) and then feed these values to the bivariate.density function:

dat.pilot <- LSCV.density(dat.ppp)
dat.global <- OS(dat.ppp, nstar = NULL)
dat.sden <- bivariate.density(data = dat.ppp, pilotH = dat.pilot, globalH = dat.global, adaptive = TRUE, edgeCorrect = TRUE)

dat.sden$Zm <- dat.sden$Zm / 0.000001 / 0.000001
hist(as.vector(dat.sden$Zm)) 

summary(as.vector(dat.sden$Zm))

xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 100, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

 plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.sden$xcol, y = dat.sden$yrow, z = t(dat.sden$v), col = col, breaks = breaks, add = TRUE)
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


# =================================================================================================================================
# Save the workspace because these analyses take forever to run:
save.image("AU_Qld_koala_habitat_v01.RData")


