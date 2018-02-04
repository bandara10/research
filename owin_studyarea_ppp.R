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
all.locations= subset(myalldata,yearnew >1998, select=X:Y)
pp=SpatialPoints(all.locations)

##### Stratified sampling #####
#We split a data.frame into color groups. From each such a group, we sample 200 rows.
df2 <- lapply(split(myalldata, myalldata$yearnew),
              function(subdf) subdf[sample(1:nrow(subdf), 200),])

d=do.call('rbind', df2) #merged into 1 data.frame


####### Select records based on distance######

source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 500)
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

###### Work out density for dat.ppp (using Diggle's edge correction):#####
# set the dimensions for the grid and the bandwidth:
dimyx = c(200,200)
bandwidth <- 1500

dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)
plot(dat.den)
# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre:
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]

x.points <- seq(from = 387950.3, to = 552950.3, by = 1000); x.lab <- x.points / 1000
y.points <- seq(from = 6862579, to = 7113579, by = 1000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = .007, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)

### Bring in rasters #####

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = scale(stack(myfullstack.b))
plot(myfullstack.b,2)
# myextent=c(387900, 553100, 6862400, 7113600) 
myextent=studyarea.shp
habitat.rr=crop(myfullstack.b, myextent, snap="near")

habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,33,39)) # select my variables
plot(habitat.rr,2)


# Convert rasters into a spatstat image object=================================================================================================================================
##Rhohat fpcnew ####
fpc.im <- as.im(habitat.rr$fpcnew )
plot(fpc.im )

# What is the nature of the association between koala sight locations and total phosphorous?
fpc.rho <- rhohat(object = dat.ppp, covariate = fpc.im)

plot(fpc.rho, xlab = "Foliage projective cover (units)", main = "")

#Annual_Mean_Temperature==================================================================================================================================
amt.im <- as.im(habitat.rr$Annual_Mean_Temperature)
plot(amt.im )

# What is the nature of the association between koala sight locations and total phosphorous?
amt.rho <- rhohat(object = dat.ppp, covariate = amt.im)

plot(amt.rho, xlab = "Annual mean temperature (units)", main = "")


#Annual_Precipitation==================================================================================================================================
apt.im <- as.im(habitat.rr$Annual_Precipitation)
plot(apt.im )

# What is the nature of the association between koala sight locations and total phosphorous?
apt.rho <- rhohat(object = dat.ppp, covariate = apt.im)

plot(apt.rho, xlab = "Annual precipitation (units)", main = "")

#Max_Temperature_of_Warmest_Month==================================================================================================================================
mtwm.im <- as.im(habitat.rr$Max_Temperature_of_Warmest_Month)
plot(mtwm.im )

# What is the nature of the association between koala sight locations and total phosphorous?
mtwm.rho <- rhohat(object = dat.ppp, covariate = mtwm.im)

plot(mtwm.rho, xlab = "Max_Temperature_of_Warmest_Month (units)", main = "")

#Min_Temperature_of_Coldest_Month==================================================================================================================================
mtcm.im <- as.im(habitat.rr$Min_Temperature_of_Coldest_Month)
plot(mtcm.im )

# What is the nature of the association between koala sight locations and total phosphorous?
mtcm.rho <- rhohat(object = dat.ppp, covariate = mtcm.im)

plot(mtcm.rho, xlab = "Max_Temperature_of_Coldest_Month (units)", main = "")

#distance_motorwayandlink==================================================================================================================================
dmwl.im <- as.im(habitat.rr$distance_motorwayandlink)
plot(dmwl.im )

# What is the nature of the association between koala sight locations and total phosphorous?
dmwl.rho <- rhohat(object = dat.ppp, covariate = dmwl.im)

plot(dmwl.rho, xlab = "distance_motorwayandlink (units)", main = "")

#distance_primaryandlink==================================================================================================================================
dprl.im <- as.im(habitat.rr$distance_primaryandlink)
plot(dprl.im )

# What is the nature of the association between koala sight locations and total phosphorous?
dprl.rho <- rhohat(object = dat.ppp, covariate = dprl.im)

plot(dprl.rho, xlab = "distance_primaryandlink (units)", main = "")


#Poisson point process model =================================================================================================================================
predList=list(fpc = fpc.im
     ,dprl = dprl.im 
     ,dmwl = dmwl.im
     ,mtcm = mtcm.im
     ,mtwm = mtwm.im 
     ,apt = apt.im)
# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1,covariates = predList) 
                                                          
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409

###### Extract the quadrature scheme from dat.ppm01:#####
Qz <- quad.ppm(dat.ppm00, drop = TRUE)
plot(Qz)
# This step will drop points within the quadrature scheme that had NA-values.

####### Saturated model:####
dat.ppm01 <- step(ppm(Qz, trend = ~ fpc*apt + dprl + dmwl + mtcm + mtwm + apt , covariates = predList))
summary(dat.ppm01)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm01)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm01) # 3466.415


###### Now add a Geyer spatial interaction term:#####
dat.ppm07 <- ppm(Qz, trend = ~ fpc+dprl+ dmwl + mtcm + mtwm + apt, interaction = Geyer(r = 5000, sat = 20),covariates = predList)
summary(dat.ppm07)$coefs.SE.CI

AIC(dat.ppm07) 
plot(dat.ppm07,ngrid=200)
diagnose.ppm(dat.ppm07)


######## QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions##
##re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
dat.qq07 <- qqplot.ppm(fit = dat.ppm07, nsim = 10)
plot(dat.qq05, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# Predict=================================================================================================================================
# Model with spatial interaction term:
trend1 <- predict(dat.ppm07, type="trend",ngrid=200) # ngrids smooth the map.
plot(trend1)
persp(trend1)

#### Express predicted koala intensity in terms of the number of koalas per hectare:####
pred$v <- pred$v * 1E04
summary(as.vector(pred$v))
hist(as.vector(pred$v))


###### Map Predictions from the model:####
breaks <- seq(from = 0, to = .000003, length = 5)
col <- brewer.pal(n = 4, name = "Reds")

plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "black")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)

##### Actual kernel smoothed data:#####
breaks <- seq(from = 0, to = 0.0015, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)




windows(); par(pin = c(1 * 5, ratio * 5), omi = c(0.5,0,0,0))
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

pool.pilot <- LSCV.density(dat.ppp)
pool.global <- OS(dat.ppp, nstar = NULL)
dat.pool <- bivariate.density(data = dat.ppp, pilotH = pool.pilot, globalH = pool.global, WIN = dat.w)

dat.pool$Zm <- dat.pool$Zm / 0.000001 / 0.000001
hist(as.vector(dat.pool$Zm)) 

summary(as.vector(dat.pool$Zm))

# Save the workspace because these analyses take forever to run:
save.image("AU_Qld_koala_habitat_v01.RData")

xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 100, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.sden$xcol, y = dat.sden$yrow, z = t(dat.sden$v), col = col, breaks = breaks, add = TRUE)
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


# =================================================================================================================================



