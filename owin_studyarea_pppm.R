#### Load libraries ####
library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatstat.utils);library(usdm)

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
all.locations= subset(myalldata,yearnew >1997,select=X:Y) #, select=X:Y REmove xy for ppp density
all.locations <- unique(all.locations)

pp=SpatialPoints(all.locations)
plot(pp, add=TRUE)

#### step 2select one from grid.========================================================================================================================

myfullstack.b <- list.files(path="spatstatmodel",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = stack(myfullstack.b)
names(myfullstack.b)

# myextent=c(387900, 553100, 6862400, 7113600) 
studyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data\\warton_data_allclimate\\ss_area.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356")

myextent=c(387900, 553100, 6862400, 7113600)
#myextent=studyarea.shp
habitat.rr=crop(myfullstack.b, myextent, snap="near") # use pp
#habitat.rr <- scale(habitat.rr)
set.seed(123)
selected.locations=gridSample(all.locations, habitat.rr, n = 1)# disaggregate a raster factor 4 and get one point.
selected.locations <- as.data.frame(selected.locations)
plot(selected.locations)
ss=SpatialPoints(selected.locations)
write.csv(selected.locations, "selected.locations.csv")




####
vifstep(habitat.rr, th=10)
# VIFs of the remained variables -------- 
#   Variables      VIF
# 1                            awc 3.926543
# 2                           clay 2.362773
# 3                           elev 2.825733
# 4                         fpcnew 1.503247
# 5                  habit1decimal 1.171915
# 6                  habit2decimal 1.162066
# 7                  habit3decimal 1.075773
# 8                           hpop 4.073797
# 9                    lot_density 4.170613
# 10                         nitro 2.978881
# 11 Precipitation_of_Driest_Month 2.657214
# 12     Precipitation_Seasonality 3.712386
# 13                           sbd 3.201861
# 14       Temperature_Seasonality 5.820057
# 15                           tpo 2.627350
# 16                           twi 1.512787
# mydata.r <- raster(studyarea.shp)
# res(mydata.r) <- 1000
# mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset


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

######### Bring in MGA56 square study area boundary map:######

aus.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\australia.shp")#AU_Qld_detail_studyarea_outline-MGA56
proj4string(aus.shp) <- CRS("+init=epsg:4326") 
aus.shp <- spTransform(aus.shp, crs(studyarea.shp))
#plot(aus.shp)

 # detail study area

qld.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\qld.shp")#AU_Qld_detail_studyarea_outline-MGA56
proj4string(qld.shp) <- CRS("+init=epsg:4326") 
qld.shp <- spTransform(qld.shp, crs(studyarea.shp))

#### selected area for model
dstudyarea.shp <- readShapePoly("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_detail_studyarea_outline-MGA56.shp")#
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

studyarea.shp <- crop(dstudyarea.shp,myextent)


plot(qld.shp, axes = TRUE)
plot(dstudyarea.shp, axes = TRUE, col="blue")
plot(studyarea.shp, add = TRUE,col="yellow")
points(x = selected.locations$X, y = selected.locations$Y)

#points(x=all.locations$X, y= all.locations$Y)

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

# dat.ppp=rescale(dat.ppp, 1000, unitname="km")
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

# Convert rasters into a spatstat image object=================================================================================================================================
##Rhohat fpcnew ####
fpc.im <- as.im(habitat.rr$fpcnew )
fpc.rho <- rhohat(object = dat.ppp, covariate = fpc.im)
plot(fpc.rho)
plot(fpc.rho ,xlim = c(0, 100),ylim = c(0, 6e+10), xlab = "lot_density(units)", main = "") 


amt.im <- as.im(habitat.rr$AnnualMeanTemperature)
apt.im <- as.im(habitat.rr$AnnualPrecipitation)
mtwm.im <- as.im(habitat.rr$Max_Temperature_of_Warmest_Month)
mtcm.im <- as.im(habitat.rr$Min_Temperature_of_Coldest_Month)
pdm.im <- as.im(habitat.rr$Precipitation_of_Driest_Month)
pwm.im <- as.im(habitat.rr$Precipitation_of_Wettest_Month)
awc.im <- as.im(habitat.rr$awc)
elev.im <- as.im(habitat.rr$elev)
clay.im <- as.im(habitat.rr$clay)
nitro.im  <- as.im(habitat.rr$nitro)
citydis.im <- as.im(habitat.rr$city_dis)
sbd.im  <- as.im(habitat.rr$sbd)
#roughness.im <- as.im(habitat.rr$roughness)
tpo.im  <- as.im(habitat.rr$tpo)
twi.im  <- as.im(habitat.rr$twi)
hpop.im <- as.im(habitat.rr$hpop)
lot.im <- as.im(habitat.rr$lot_density)
habit3.im  <- as.im(habitat.rr$habit3decimal)
habit2.im  <- as.im(habitat.rr$habit2decimal)
habit1.im  <- as.im(habitat.rr$habit1decimal)
habitdis2.im <- as.im(habitat.rr$Dis_habitat_suitable_2)
roadsM.im  <- as.im(habitat.rr$roads_motor)
roadsother.im  <- as.im(habitat.rr$roads_other)
dmwl.im <- as.im(habitat.rr$distance_motorwayandlink)
dprl.im <- as.im(habitat.rr$distance_primaryandlink)
#Poisson point process model =================================================================================================================================
predList=list(
              amt=amt.im 
              ,apt=apt.im 
              ,mtwm=mtwm.im 
              ,mtcm=mtcm.im 
              ,pdm=pdm.im 
              ,pwm=pwm.im 
              ,awc=awc.im 
              ,elev=elev.im 
              ,clay=clay.im 
              ,nitro=nitro.im  
              ,sbd=sbd.im 
              ,tpo=tpo.im  
              ,twi=twi.im  
              ,hpop=hpop.im 
              ,lot= lot.im
              ,habit3=habit3.im  
              ,habit2=habit2.im  
              ,habit1=habit1.im  
            
              )
# Null model:fpc=fpc.im 
(dat.ppm00 <- ppm(dat.ppp, trend = ~ 1)) #,covariates = predList
    diagnose.ppm(dat.ppm00)                                                   
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409

###### Extract the quadrature scheme from dat.ppm01:#####
Qz <- quad.ppm(dat.ppm00, drop = TRUE)
plot(Qz)


# This step will drop points within the quadrature scheme that had NA-values.

####### Saturated model:#### 
dat.ppm01 <- step(ppm(dat.ppp, trend = ~ fpc
                      + awc + elev+ nitro+ tpo+ hpop*habit2+habit2+habit3,covariates = predList)) #+ roadsM + roadsother
# Estimate        S.E.       CI95.lo     CI95.hi Ztest        Zval
# (Intercept) -17.89636094 0.061275733 -1.801646e+01 -17.7762627   *** -292.062780
# elev         -1.69757834 0.072997905 -1.840652e+00  -1.5545051   ***  -23.255165
# dprl         -0.81895297 0.053805661 -9.244101e-01  -0.7134958   ***  -15.220573
# citydis      -0.27893751 0.032338545 -3.423199e-01  -0.2155551   ***   -8.625543
# tpo          -0.10470027 0.028988172 -1.615160e-01  -0.0478845   ***   -3.611827
# fpc           0.05501796 0.027654066  8.169829e-04   0.1092189     *    1.989507
# awc           0.09811207 0.048502497  3.048919e-03   0.1931752     *    2.022825
# habit1        0.06683203 0.030034659  7.965183e-03   0.1256989     *    2.225164
# habit3        0.09571611 0.022103409  5.239422e-02   0.1390380   ***    4.330378
# habit2        0.16583757 0.022129926  1.224637e-01   0.2092114   ***    7.493815
# nitro         0.31975717 0.036986334  2.472653e-01   0.3922490   ***    8.645279
# hpop          0.12275546 0.008701996  1.056999e-01   0.1398111   ***   14.106587

dat.ppm011 <- ppm(Qz, trend = ~ awc+clay+fpc+habit1+habit2+habit2+hpop+lot+nitro+pdm+tpo+twi,covariates = predList) #+ roadsM
diagnose.ppm(dat.ppm011)

#drop nitro

dat.ppm012 <- ppm(Qz, trend = ~ awc+clay+fpc+habit1+habit2+habit2+hpop+lot+pdm+tpo+twi,covariates = predList) #+ roadsM
diagnose.ppm(dat.ppm012)
AIC(dat.ppm012)
dat.ppm012=update(dat.ppm012, .~. ,interaction =  Strauss(2000))
diagnose.ppm(dat.ppm012)

pred012 <- predict(dat.ppm012, type="trend") # ngrids smooth the map.
pred012=pred012*1E06
pred012.grid <- as(pred012,"SpatialGridDataFrame")

pred012.r=raster(pred012.grid)

pred012.r=pred012.r*2.526440

plot(pred012.r)
plot(pp, add=TRUE,cex=.5, pch=1)

writeRaster(pred012.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred0122.r.tif",overwrite=TRUE)





dat.ppm011 <- ppm(Qz, trend = ~ 
fpc+apt+mtwm+mtcm+awc+elev+lot+clay+nitro+citydis+sbd+roughness+tpo+habit3+habit1+roadsM,interaction =  Strauss(2000),covariates = predList) #+ roadsM 

#dat.ppm011<- ppm(Qz ~ polynom(mtcm,mtwm, 2)+ polynom(pwm,pdm,2)+polynom(habit1,lot,2)+polynom(elev,clay,2),covariates = predList)


diagnose.ppm(dat.ppm01)
###drop
dat.ppm011 <- ppm(Qz, trend = ~ 
fpc+apt+mtwm+mtcm+elev+clay+nitro+dprl+tpo+habit3+lot*habit3+habit1,covariates = predList) #+ roadsM 

diagnose.ppm(dat.ppm011)


eem1 <- envelope(dat.ppp,Kest,funargs = list(lambda=dat.ppm01),global=TRUE)
plot(eem1)

##model summary
summary(dat.ppm011)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm011)$coefs.SE.CI)
tdat1[order(tdat1$Zval),] 
AIC(dat.ppm01) # 3466.415

###plot
pred <- predict(dat.ppm011, type="trend") # ngrids smooth the map.
pred=pred*1E06
pred.grid <- as(pred,"SpatialGridDataFrame")
pred.r=raster(pred.grid)
plot(pred.r)
writeRaster(pred.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred014.tif",overwrite=TRUE)


### standard erros
se <- predict.ppm(dat.ppm01, se=TRUE, interval = c("confidence"))
sse=se$se
sse=sse*1E06
sse.grid <- as(sse,"SpatialGridDataFrame")
sse.r=raster(sse.grid)
plot(sse.r)
writeRaster(sse.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred014sse.tif",overwrite=TRUE)


plot(sse)
######CI Upper==========
ciu <- predict.ppm(dat.ppm01,interval = c("confidence"), level = 0.95)
ciiu=ci$`97.5%`
ciiu=ciiu* 1E06
ciiu.grid <- as(ciiu,"SpatialGridDataFrame")
ciiu.r=raster(ciiu.grid)
plot(ciiu.r)
writeRaster(ciiu.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred01ciiu.tif",overwrite=TRUE)

###### CI Lower
ciL <- predict.ppm(dat.ppm01,interval = c("confidence"), level = 0.95)
ciiL=ciL$`2.5%`
ciiL=ciiL* 1E06
ciiL.grid <- as(ciiL,"SpatialGridDataFrame")
ciiL.r=raster(ciiL.grid)
plot(ciiL.r)
writeRaster(ciiL.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/pred01ciiL.tif",overwrite=TRUE)


#########


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
dat.ppm02 <- step(kppm (Qz, trend = ~ fpc+citydis+ awc + elev+ nitro+ sbd + tpo+ twi+ hpop+habit1+habit2+habit3+dprl 
                   ,clusters = "Thomas",covariates = predList)) #+ roadsM + roadsother

##model summary
summary(dat.ppm02)$coefs.SE.CI
tdat2<- data.frame(summary(dat.ppm02)$coefs.SE.CI)
tdat2[order(tdat2$Zval),] 
AIC(dat.ppm02) # 3466.415

predk=predict.kppm(dat.ppm02,ngrid=200)
predk=predk* 1E06
plot(predk)

######## write raster to plot in arcgis
predk.grid <- as(predk,"SpatialGridDataFrame")
predk.grid.r=raster(predk.grid)
plot(predk.grid.r)
writeRaster(predk.grid.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/predkppm.tif",overwrite=TRUE)


### standard erros
se2 <- predict(dat.ppm02, se=TRUE)
sse2=se2$se
sse2=sse2*1E06
sse2.grid <- as(sse2,"SpatialGridDataFrame")
sse2.r=raster(sse2.grid)
plot(sse2.r)
writeRaster(sse2.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/predkkpmse.tif",overwrite=TRUE)


plot(sse)

######================

# plot(dat.ppm02,what="statistic",pause=FALSE)
# pred2 <- predict(dat.ppm02, type="trend",ngrid=200) # ngrids smooth the map.
# plot(pred* 1E06)


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
dat.qq03 <- qqplot.ppm(fit =dat.ppm01, nsim = 10)
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

##### Actual kernel smoothed data:#####
###### Work out density for dat.ppp (using Diggle's edge correction):#####
# set the dimensions for the grid and the bandwidth:
dimyx = c(200,200)
bandwidth <- 2000

dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)
plot(dat.den*1E06)

dat.den=dat.den*1E06
dat.den.grid <- as(dat.den,"SpatialGridDataFrame")
dat.den.r=raster(dat.den.grid)
plot(dat.den.r)
writeRaster(dat.den.r, "C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/preddat_denr.tif",overwrite=TRUE)




#Adaptive kernel smoothing using sparr =================================================================================================================================
# Adaptive smoothing and dividing surfaces to make relative risk density map:

library(spatialkernel); library(sparr); library(rgl); library(spatialkernel)

# Calculate a pilot bandwidth (using leave-one-out least-squares cross-validation) and a global bandwidth (using maximal smoothing) and then feed these values to the bivariate.density function:

dat.pilot <- LSCV.density(dat.ppp)
dat.global <- OS(dat.ppp, nstar = NULL)
dats.den <- bivariate.density(data = dat.ppp, pilotH = dat.pilot, globalH = dat.global, adaptive = TRUE, 
                              edgeCorrect = TRUE)
p=as.im(dats.den)
dat.sden=p*1E06

dat.sden.grid <- as(dat.sden,"SpatialGridDataFrame")
dat.sden.r=raster(dat.sden.grid)
plot(dat.sden.r)

plot(p)







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


