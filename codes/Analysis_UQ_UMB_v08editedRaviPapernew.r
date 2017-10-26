library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):
#load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v07.RData")
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v08edited_distance_Ravi.RData")

mydatasighting <- read.csv("sightingdata\\mydatasighting_cleaned.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)

# Analyse data for 2011:
id <- mydata$yearnew <2000
table(id)
mydata <- mydata[id,]
dim(mydata)
names(mydata)


# 354 rows

# Create a new variable called "caller". If home phone number is missing, use the work number. If work number is missing we use the caller name:
# mydata$caller <- ifelse(mydata$homephone == " ", as.character(mydata$workphone), as.character(mydata$homephone))
# If either home or work phone number is missing, use the callername:
# mydata$caller <- ifelse(mydata$caller == " ", as.character(mydata$callername), as.character(mydata$caller))
# sort(table(mydata$caller))

# Drop those records without a valid caller identifier (n = 16):
# id <- mydata$caller != " "
# mydata <- mydata[id,]

# Create a unique 1 to n caller identifier:
# tcaller <- unique(mydata$caller)
# tid <- 1:length(tcaller)
# lookup <- data.frame(id = tid, caller = tcaller)

# Add the unique caller identifier to the mydata data frame:
# mydata$uniq <- lookup$id[match(mydata$caller, lookup$caller)]
# sort(table(mydata$uniq))

# Drop the multiple reporters:
# id <- mydata$uniq != 182 & mydata$uniq != 19 & mydata$uniq != 64 & mydata$uniq != 65 & mydata$uniq != 9 & mydata$uniq != 62 & mydata$uniq != 23 & mydata$uniq != 18 & mydata$uniq != 88 
# mydata <- mydata[id,]


# Bring in MGA56 study area boundary map:
unzip("vector\\AU_Qld_study_area.zip")
studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 
plot(studyarea.shp)
id <- mydata$lat > -30 # remove record with wrong coordinate
mydata <- mydata[id,]
head(mydata)
# mydata in lat-lon:
coords <- SpatialPoints(mydata[, c("lng", "lat")])
mydata.ll <- SpatialPointsDataFrame(coords, mydata)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))
plot(mydata.mga)
# Add the MGA coordinates to mydata:
mydata$x <- coordinates(mydata.mga)[,1]
mydata$y <- coordinates(mydata.mga)[,2]

windows(); plot(mydata.mga, cex = 0.3, col = 'blue', pch = 15)

# Create a 3 km by 3 km grid:
mydata.r <- raster(studyarea.shp)
res(mydata.r) <- 1000
mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#extenet should be similar to the selected area not the full dataset
#writeRaster(mydata.r,"testing.asc")
 #windows(); plot(mydata.r)
# create the disance matrix
#####Go to 109 for grid based sampling. START Distaance based selection using the function
# Function to estimate the distance to the nearest point of the same vector
source("Lib_DistEstimatesToItself.r")
# start distance based selection method
mydata.n <- mydata.mga[c(3,4)]
mydata.n = as.data.frame(mydata.n)
myInPtDF = SpatialPointsDataFrame(mydata.n[c("x","y")], mydata.n)
windows();plot(myInPtDF, axes = T)
myInPtDF$disttoitself = Lib_DistEstimatesToItself(myInPtDF$x, myInPtDF$y)
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 100)
windows();plot(myInPtDF2, axes = T)
acsel<-myInPtDF2[,1:2]
#distance based selection. the next steps are to veryfy only.
acsel = myInPtDF[,1:2]
acsel$all = 1
myFD2 = myInPtDF2[,1:2]
myFD2$all = 0
myFD = rbind(acsel,myFD2)
windows();plot(myFD, axes = T,col=c("red","blue")[myFD$all+1])

acsel<-as.data.frame(myFD$all==0)#myFD3
windows(); plot(acsel)
#########END of distance based selection#####
# Sample points from within each grid cell:
# kdate as a date:
# dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
# dat$kmonth <- format(dat$kdate, format = "%m")

# =================================================================================================================================
#1.total phosphorous:% Read in explanatory variable raster maps - 

tpo.r <- raster("raster\\AU_Qld_soil_total_phosphorous_000_005cm-wgs84.tif")
projection(tpo.r) <- CRS("+init=epsg:4326")
# Reproject to MGA-56:
tpo.r <- projectRaster(from = tpo.r, crs = CRS("+init=epsg:28356"))
windows(); plot(tpo.r, main="tpo")
# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
#chnaged reolution 34*34 to 100 *102
r <- raster(ncol = 100, nrow = 100)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
tpo.crop <- crop(x = tpo.r, y = r, snap = "near")
windows(); plot(tpo.crop)
# Resample:okay. lets keep this resolution for all rasters. why do we need to change this match r which is our grid 
tpo.r <- resample(x =tpo.crop , y = r)
writeRaster(tpo.r,"raster_crop\\tpo.TIF", overwrite=TRUE)
# Plot to check:
windows(); plot(tpo.r, main="tpo")
points(acsel, cex = 0.4, col = 'red', pch = 15)


---------------------------------------------------------------------------------
--------------------------------------------------------------
# Make a ppp object for spatstat:

# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:
source("owin2sp_source.r")

dat.w <- as(as(studyarea.shp, "SpatialPolygons"), "owin")
dat.spw <- owin2SP(dat.w)

# Set projection of the owin as GDA94 / SA Lambert:
proj4string(dat.spw)
proj4string(dat.spw) <- CRS("+init=epsg:28356") 
dat.spw <- raster::intersect(x = dstudyarea.shp, y = dat.spw)

# Convert the sp object into a window:
dat.w <- as.owin(dat.spw)

# set the dimensions for the grid and the bandwidth:
dimyx = c(200,200)
bandwidth <- 2000

# Select only those koalas within the owin:
# id <- inside.owin(x = mydata$x, y = mydata$y, w = dat.w)
# mydata <- mydata[id,]

windows(); plot(studyarea.shp, axes = TRUE)
points(x = acsel[,1], y = acsel[,2])
windows();plot(acsel)
# Make a ppp object:
acsel1<-acsel$x
acsel2<-acsel$y
dat.ppp <- ppp(x = acsel1, y = acsel2, window = dat.w)

windows(); plot(dat.ppp, axes = TRUE)

# Work out density for dat.ppp (using Diggle's edge correction):
dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)

# Express density as the number of koalas per hectare (100 m * 100 m) cf number of koalas per square metre:
dat.den$v <- dat.den$v / 1e-04
summary(as.vector(dat.den$v))

windows(); plot(dat.den)

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]
x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 0.2, length =8)
col <- brewer.pal(n = 7, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
points(acsel, cex = 0.4, col = 'red', pch = 15)

## rhohat - total phosphorous:Computes a smoothing estimate of the intensity of a point process, as a function of a (continuous) spatial covariate.
-----------------------------------------------------------------------------------------------------------------------------------
  #read all cropped rasters and make a stack for easy reading
  
  stack.n<- list.files(pattern="\\.tif$")
stack <- stack(stack.n)
summary(stack)
stack.c<- crop(x = stack, y = mydata.r, snap = "near")
head(stack)
writeRaster(stack.c, names(stack.n), bylayer=TRUE, driver='GTIFF')

-------------------------------------------------------------------------------------------------------------------


# =================================================================================================================================
# Poisson point process model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(aspect = aspect.im, residential = residential.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409
windows();plot(dat.ppp)
# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409
plot(dat.ppp)

# Backward stepwise variable selection:
dat.ppm01 <- ppm(dat.ppp, trend = ~ tpo + sbd + water + nitro + roaden + hpop + temp + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
# Values of the covariates 'roaden', 'temp' were NA or undefined at 0.5% (5 out of 1005) of the quadrature points. 
# See http://gis.stackexchange.com/questions/161893/creating-polygonal-windows-for-spatial-analysis-in-r for a work around.

# Extract the quadrature scheme from dat.ppm01:
Qz <- quad.ppm(dat.ppm01, drop = TRUE)
# This step will drop points within the quadrature scheme that had NA-values.

# Run the model again, but this time, use the corrected quadrature scheme:
dat.ppm02 <- ppm(Qz, trend = ~ tpo + sbd + water + nitro + roaden + hpop + temp + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm02)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm02)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm02) # 3413.045

# Drop tpo:temp +
dat.ppm03 <- ppm(Qz, trend = ~ sbd +tpo+ water + nitro + roaden +  hpop + elev + habit01 + habit02 + habit03, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm03)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm03)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm03) # 3411.272


# Drop habitat03:tpo
dat.ppm04 <- ppm(Qz, trend = ~ sbd + roaden+water + nitro + habit03+ hpop + elev + habit01 + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm04)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm04)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm04) # 3460.872


# Drop habit1:
dat.ppm05 <- ppm(Qz, trend = ~ sbd + roaden+water + nitro + habit03+ hpop + elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm05)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm05)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm05) # 3459.164


# Drop roadden:
dat.ppm06 <- ppm(Qz, trend = ~ sbd + water + nitro + habit03+ hpop + elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm06)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm06)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm06) # 3457.573


# Drop hpop:
dat.ppm07 <- ppm(Qz, trend = ~ sbd + water + nitro + habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm07)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm07)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm07) # 3455.898


# Drop nitro:
dat.ppm08 <- ppm(Qz, trend = ~ sbd + water +  habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm08)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm08)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm08) # 3455.163


# Drop water
dat.ppm09 <- ppm(Qz, trend = ~ sbd + habit03+ elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm09)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm09)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm09) # 3453.953


# Drop habit03
dat.ppm10 <- ppm(Qz, trend = ~ sbd +  elev  + habit02, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm10)$coefs.SE.CI
AIC(dat.ppm10) # 3411.58

windows(); plot(dat.ppm10, how = "image", se = TRUE)

# The lurking variable plot gives you an idea of mis-specification of spatial trend:
windows(); diagnose.ppm(dat.ppm10, type = "raw")
savePlot(filename = "datppm10",type = c( "png"),device = dev.cur())

# Now add a Strauss spatial interaction term:
dat.ppm11 <- ppm(Qz, trend = ~ sbd +  elev  + habit02, interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm11)$coefs.SE.CI
AIC(dat.ppm11) # 2459.733

windows(); plot(dat.ppm11, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm11, type = "raw")


# Drop +  elev:
dat.ppm12 <- ppm(Qz, trend = ~ sbd   + habit02, interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm12)$coefs.SE.CI
AIC(dat.ppm12) # 2459.799

windows(); plot(dat.ppm12, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm12, type = "raw")


# Drop + habit02
dat.ppm13 <- ppm(Qz, trend = ~ sbd   , interaction = Strauss(r = 6000), 
   covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im, habit01 = habit01.im, habit02 = habit02.im, habit03 = habit03.im))
summary(dat.ppm13)$coefs.SE.CI
AIC(dat.ppm13) # 2459.806

windows(); plot(dat.ppm13, how = "image", se = TRUE)
windows(); diagnose.ppm(dat.ppm13, type = "raw")


# QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
# dat.qq11 <- qqplot.ppm(fit = dat.ppm11, nsim = 10)

# windows(); plot(dat.qq11, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# =================================================================================================================================
# Model with spatial interaction term:
dat.pred13 <- predict(dat.ppm13); pred <- dat.pred13

# Express predicted koala intensity in terms of the number of koalas per hectare:
pred$v <- pred$v * 1E04
summary(as.vector(pred$v))
windows(); hist(as.vector(pred$v))



# Predictions from the model:
breaks <- seq(from = 0, to = 0.00013, length = 6)
col <- brewer.pal(n = 5, name = "Reds")

windows(); plot(x = xlim, y = ylim, main="Predictions from the model", type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "blue")

plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Predictions from the model",type = c( "png"),device = dev.cur())

# Actual kernel smoothed data:
breaks <- seq(from = 0, to = 0.0015, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n",main="Actual kernel smoothed data", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Actual kernel smoothed data",type = c( "png"),device = dev.cur())


# =================================================================================================================================


# Save the workspace because these analyses take forever to run:
save.image("AU_Qld_koala_habitat_v08.RData")

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
