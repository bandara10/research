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
all.locations= subset(myalldata,yearnew ==2010, select=X:Y)
pp=SpatialPoints(all.locations)
plot(pp)
##### Stratified sampling #####
#We split a data.frame into color groups. From each such a group, we sample 200 rows.
df2 <- lapply(split(myalldata, myalldata$yearnew),
              function(subdf) subdf[sample(1:nrow(subdf), 200),])

d=do.call('rbind', df2) #merged into 1 data.frame


####### Select records based on distance######

source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 50)
selected.locations = select.locations[,1:2]
mydata.ll=SpatialPoints(selected.locations)
proj4string(mydata.ll) <- CRS("+init=epsg:28356") 

######### Bring in MGA56 square study area boundary map:######

studyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 

# detail study area

dstudyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(x = mydata$X, y = mydata$Y)


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

plot(dat.ppp, axes = TRUE)


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

amt.im <- as.im(habitat.rr$Annual_Mean_Temperature)

apt.im <- as.im(habitat.rr$Annual_Precipitation)

mtwm.im <- as.im(habitat.rr$Max_Temperature_of_Warmest_Month)

mtcm.im <- as.im(habitat.rr$Min_Temperature_of_Coldest_Month)

awc.im <- as.im(habitat.rr$awc)

elev.im <- as.im(habitat.rr$elev)

nitro.im  <- as.im(habitat.rr$nitro)

sbd.im  <- as.im(habitat.rr$sbd)

tpo.im  <- as.im(habitat.rr$tpo)

twi.im  <- as.im(habitat.rr$twi)

habit3.im  <- as.im(habitat.rr$habit3decimal)
habit2.im  <- as.im(habitat.rr$habit2decimal)
habit1.im  <- as.im(habitat.rr$habit1decimal)

roadsM.im  <- as.im(habitat.rr$roads_motor)

roadsother.im  <- as.im(habitat.rr$roads_other)

dmwl.im <- as.im(habitat.rr$distance_motorwayandlink)

dprl.im <- as.im(habitat.rr$distance_primaryandlink)
#Poisson point process model =================================================================================================================================
predList=list(fpc = fpc.im
     ,mtcm = mtcm.im
     ,mtwm = mtwm.im 
     ,apt = apt.im
     ,awc =awc.im
     ,elev= elev.im
     ,nitro= nitro.im
     ,sbd= sbd.im
     ,tpo= tpo.im
     ,twi = twi.im
     ,habit3 = habit3.im
     ,habit2 = habit2.im
     ,habit1 = habit1.im
     ,roadsM = roadsM.im
     ,roadsother= roadsother.im
     ,dprl=dprl.im
     ,dmwl= dmwl.im)
# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1,covariates = predList) 
                                                       
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409

###### Extract the quadrature scheme from dat.ppm01:#####
Qz <- quad.ppm(dat.ppm00, drop = TRUE)
plot(Qz)
# This step will drop points within the quadrature scheme that had NA-values.

####### Saturated model:#### ??(x, y)=exp(??0 +??1x+??2y)
dat.ppm01 <- step(ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ 
                        habit3 + habit2 + habit1 + roadsM + roadsother, covariates = predList))
summary(dat.ppm01)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm01)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm01) # 3466.415
pred <- predict(dat.ppm01, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer
persp(pred)
pred1 <- pred* 1E06
max(pred1)
coll <- colourmap(terrain.colors(100,1), range=c(0,4.908895))
plot(coll,col.ticks="black")
plot(pred1, col=coll)


###### Now add a Geyer spatial interaction term:#####
dat.ppm07 <- ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+habit3 + habit2 + habit1 + roadsM + roadsother
                 , interaction = Geyer(r = 2000, sat = 20),covariates = predList)
summary(dat.ppm07)
coef(dat.ppm07)
summary(dat.ppm07)$coefs.SE.CI

AIC(dat.ppm07) 
plot(dat.ppm07,ngrid=200)
diagnose.ppm(dat.ppm07)


######## QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions##
##re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
dat.qq07 <- qqplot.ppm(fit = dat.ppm07, nsim = 10)
plot(dat.qq07, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# Predict=================================================================================================================================
# Model with spatial interaction term:
pred <- predict(dat.ppm07, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer
persp(pred)
pred1 <- pred* 1E06
max(pred1)
coll <- colourmap(terrain.colors(100,1), range=c(0,4.908895))
plot(coll,col.ticks="black")
plot(pred1, col=coll)

#### Express predicted koala intensity in terms of the number of koalas per square km:####
pred$v <- pred$v * 1E06
summary(as.vector(pred$v))
hist(as.vector(pred$v))
max(pred)

##### step 2. Model with distance covariates

dat.ppm01 <- step(ppm(Qz, trend = ~ fpc +  mtcm + mtwm + apt + awc + elev+ nitro+ sbd + tpo+ twi+ 
                        habit3 + habit2 + habit1 + roadsM + roadsother+dprl.im+dmwl.im, covariates = predList))
summary(dat.ppm01)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm01)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm01) # 3466.415
# Predict=================================================================================================================================
# Model with spatial interaction term:
pred <- predict(dat.ppm01, type="trend",ngrid=200) # ngrids smooth the map.
plot(pred* 1E06) #Express predicted koala intensity in terms of the number of koalas per square kilometer
persp(pred)
pred1 <- pred* 1E06
max(pred1)
coll <- colourmap(terrain.colors(200,1), range=c(0,4.072643))
plot(coll,col.ticks="black")
plot(pred1, col=coll)



###### Map Predictions from the model:#### Another color scheme
breaks <- seq(from = 0, to = 0.31, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.31), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "black")
plot(studyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)




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


