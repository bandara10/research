library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v01.RData")

mydatasighting <- read.csv("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\mydatasighting.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)

# Analyse data for 2011:
id <- mydata$yearnew == 2011
table(id)
mydata <- mydata[id,]
dim(mydata)
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
unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\vector\\AU_Qld_study_area.zip")

studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 

id <- mydata$lat > -30
mydata <- mydata[id,]

# mydata in lat-lon:
coords <- SpatialPoints(mydata[, c("lng", "lat")])
mydata.ll <- SpatialPointsDataFrame(coords, mydata)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))

# Add the MGA coordinates to mydata:
mydata$x <- coordinates(mydata.mga)[,1]
mydata$y <- coordinates(mydata.mga)[,2]

windows(); plot(mydata.mga, cex = 0.3, col = 'blue', pch = 15)

# Create a 3 km by 3 km grid:
mydata.r <- raster(mydata.mga)
res(mydata.r) <- 3000
mydata.r <- extend(mydata.r, extent(mydata.r) + 3000)

# windows(); plot(mydata.r)

# Sample points from within each grid cell:
acsel <- gridSample(mydata.mga, mydata.r, n = 1)
mydata.p <- rasterToPolygons(mydata.r)

windows(); plot(mydata.p, border = 'gray')
points(mydata.mga)
# Selected points in red:
points(acsel, cex = 0.2, col = 'red', pch = 15)

dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

windows(); plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(x = mydata$x, y = mydata$y)

# kdate as a date:
# dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
# dat$kmonth <- format(dat$kdate, format = "%m")


# =================================================================================================================================
# Read in explanatory variable raster maps - total phosphorous:

tpo.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_soil_total_phosphorous_000_005cm-wgs84.tif")
projection(tpo.r) <- CRS("+init=epsg:4326")
windows(); plot(tpo.r)

# Reproject to MGA-56:
tpo.r <- projectRaster(from = tpo.r, crs = CRS("+init=epsg:28356"))

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
tpo.crop <- crop(x = tpo.r, y = r, snap = "near")
windows(); plot(tpo.crop)

# Resample:
tpo.r <- resample(x = tpo.crop, y = r)
writeRaster(tpo.r, filename = "tpo.r",overwrite=TRUE)
# Plot to check:
windows(); plot(tpo.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Read in explanatory variable raster maps - Soil bulk density:

sbd.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_soil_bulk_density_000_005cm-wgs84.tif")
projection(sbd.r) <- CRS("+init=epsg:4326")
windows(); plot(sbd.r)

# Reproject to MGA-56:
sbd.r <- projectRaster(from = sbd.r, crs = CRS("+init=epsg:28356"))

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
sbd.crop <- crop(x = sbd.r, y = r, snap = "near")
windows(); plot(sbd.crop)

# Resample:
sbd.r <- resample(x = sbd.crop, y = r)

# Plot to check:
windows(); plot(sbd.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Soil clay:

clay.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_soil_clay_000_005cm-wgs84.tif")
projection(clay.r) <- CRS("+init=epsg:4326")
windows(); plot(clay.r)

# Reproject to MGA-56:
clay.r <- projectRaster(from = clay.r, crs = CRS("+init=epsg:28356"))

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
clay.crop <- crop(x = clay.r, y = r, snap = "near")
windows(); plot(clay.crop)

# Resample:
clay.r <- resample(x = clay.crop, y = r)
writeRaster(clay.r, filename = "clay.r",overwrite=TRUE)
# Plot to check:
windows(); plot(clay.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Soil water capacity:

water.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_soil_available_water_capacity_000_005cm-wgs84.tif")
projection(water.r) <- CRS("+init=epsg:4326")
windows(); plot(water.r)

# Reproject to MGA-56:
water.r <- projectRaster(from = water.r, crs = CRS("+init=epsg:28356"))

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
water.crop <- crop(x = water.r, y = r, snap = "near")
windows(); plot(water.crop)

# Resample:
water.r <- resample(x = water.crop, y = r)
writeRaster(water.r, filename = "water.r",overwrite=TRUE)
# Plot to check:
windows(); plot(water.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Soil nitrogen:

nitro.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_soil_nitrogen_000_005cm-wgs84.tif")
projection(nitro.r) <- CRS("+init=epsg:4326")
windows(); plot(nitro.r)

# Reproject to MGA-56:
nitro.r <- projectRaster(from = nitro.r, crs = CRS("+init=epsg:28356"))

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
nitro.crop <- crop(x = nitro.r, y = r, snap = "near")
windows(); plot(nitro.crop)

# Resample:
nitro.r <- resample(x = nitro.crop, y = r)
writeRaster(nitro.r, filename = "nitro.r",overwrite=TRUE)
# Plot to check:
windows(); plot(nitro.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Road network (expressed as the number of metres of road per square km grid):

unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\vector\\AU_Qld_road_length.zip")

roaden.shp <- readShapePoly("AU_Qld_road_length-MGA65.shp")
proj4string(roaden.shp) <- CRS("+init=epsg:28356") 

# Vector to raster conversion:
roaden.r <- rasterize(x = roaden.shp, y = mydata.r, field = "Sum_lengh")
writeRaster(roaden.r, filename = "roaden.r",overwrite=TRUE)
windows(); plot(roaden.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Human population density

hpop.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_hpop-MGA65.asc")
projection(hpop.r) <- CRS("+init=epsg:28356")
windows(); plot(hpop.r)

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
hpop.crop <- crop(x = hpop.r, y = r, snap = "near")
windows(); plot(hpop.crop)

# Resample:
hpop.r <- resample(x = hpop.crop, y = r)
writeRaster(hpop.r, filename = "hpop.r",overwrite=TRUE)
# Plot to check:
windows(); plot(hpop.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Rainfall

rain.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_rainfall-MGA65.asc")
projection(rain.r) <- CRS("+init=epsg:28356")

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
rain.crop <- crop(x = rain.r, y = r, snap = "near")
windows(); plot(rain.crop)

# Resample:
rain.r <- resample(x = rain.crop, y = r)
writeRaster(rain.r, filename = "rain.r",overwrite=TRUE)
# Plot to check:
windows(); plot(rain.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Temperature

temp.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_temp-MGA65.asc")
projection(temp.r) <- CRS("+init=epsg:28356")

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
temp.crop <- crop(x = temp.r, y = r, snap = "near")
windows(); plot(temp.crop)

# Resample:
temp.r <- resample(x = temp.crop, y = r)
writeRaster(temp.r, filename = "temp.r",overwrite=TRUE)
# Plot to check:
windows(); plot(temp.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Elevation

elev.r <- raster("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster\\AU_Qld_dem-MGA65.asc")
projection(elev.r) <- CRS("+init=epsg:28356")

windows(); plot(elev.r, axes = TRUE)

# Create an empty raster of appropriate dimensions. Remind ourselves of the dimensions of mydata.r which we created earlier:
mydata.r
r <- raster(ncol = 52, nrow = 56)
extent(r) <- extent(mydata.r)

# Crop the raster to the shape file spatial extent:
elev.crop <- crop(x = elev.r, y = r, snap = "near")
windows(); plot(elev.crop)

# Resample:
elev.r <- resample(x = elev.crop, y = r)
writeRaster(elev.r, filename = "elev.r",overwrite=TRUE)
# Plot to check:
windows(); plot(elev.r)
plot(mydata.p, add = TRUE)


# ---------------------------------------------------------------------------------------------------------------------------------
# Habitat
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v01.RData")
unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\vector\\AU_Qld_koala_habitat_suitability.zip")
habit.shp <- readShapePoly("AU_Qld_koala_habitat_suitability_MGA56.shp")
proj4string(habit.shp) <- CRS("+init=epsg:28356") 

# Vector to raster conversion:
habit.r <- rasterize(x = habit.shp, y = mydata.r, field = "suitabilit")

windows(); plot(habit.r)
plot(mydata.p, add = TRUE)


# =================================================================================================================================
# Make a ppp object for spatstat:

# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:
source("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\owin2sp_source.r")

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
bandwidth <- 1500

# Select only those koalas within the owin:
# id <- inside.owin(x = mydata$x, y = mydata$y, w = dat.w)
# mydata <- mydata[id,]

windows(); plot(studyarea.shp, axes = TRUE)
points(x = acsel[,1], y = acsel[,2])

# Make a ppp object:
dat.ppp <- ppp(x = acsel[,1], y = acsel[,2], window = dat.w)

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

breaks <- seq(from = 0, to = 0.002, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


# =================================================================================================================================
# rhohat - total phosphorous:

# Convert tpo.r into a spatstat image object:
tpo.im <- as.im(tpo.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
tpo.rho <- rhohat(object = dat.ppp, covariate = tpo.im)

windows(); plot(tpo.rho, xlab = "Total phosphorous (units)", main = "")


------------------------------------------------------------------------------------------------------------------------------
# rhohat - total nitrogen:
  
# Convert tpo.r into a spatstat image object:
nitro.im <- as.im(nitro.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total Nitrogen?
nitro.rho <- rhohat(object = dat.ppp, covariate = nitro.im)

windows(); plot(nitro.rho, xlab = "Total soil nitrogen (units)", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat - soil bulk density:

# Convert tpo.r into a spatstat image object:
sbd.im <- as.im(sbd.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations a soil bulk density?
sbd.rho <- rhohat(object = dat.ppp, covariate = sbd.im)

windows(); plot(sbd.rho, xlab = "Soil bulk density (units)", main = "")


-----------------------------------------------------------------------------------------------------------------------------------
# rhohat - clay:

# Convert tpo.r into a spatstat image object:
clay.im <- as.im(clay.r)
# windows(); plot(tpo.im, axes = TRUE)

# What is the nature of the association between koala sight locations and clay?
clay.rho <- rhohat(object = dat.ppp, covariate = clay.im)

windows(); plot(clay.rho, xlab = "clay (units)", main = "")


------------------------------------------------------------------------------------------------------------------------------
# rhohat - elevation:

# Convert .ph.r into a spatstat image object:
elev.im <- as.im(elev.r)

windows(); plot(elev.im, axes = TRUE)

# What is the nature of the association between koala sight locations and elevation?
elev.rho <- rhohat(object = dat.ppp, covariate = elev.im)

windows(); plot(elev.rho, xlim = c(0, 600), xlab = "Elevation (metres)", main = "")

# Most koalas sighted at elevations of less than 100 metres. Plot raster map to show areas less than 100 m:
elev100.r <- elev.r < 100

windows(); plot(elev100.r, axes = TRUE)
points(x = mydata$x, y = mydata$y)


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - human population density:

# Convert .ph.r into a spatstat image object:
hpop.im <- as.im(hpop.r)

windows(); plot(hpop.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)

windows(); plot(hpop.rho, xlim = c(0, 3000), ylim = c(0, 0.1e-06), xlab = "Human population density", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - annual rainfall (mm):

# Convert .ph.r into a spatstat image object:
rainfall.im <- as.im(rainfall.r)

windows(); plot(rainfall.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
rainfall.rho <- rhohat(object = dat.ppp, covariate = rainfall.im)

windows(); plot(rainfall.rho, xlab = "Annual rainfall (mm)", main = "")

# Strong association with rainfall between 1200 and 1400 mm per year:
crainfall.r <- rainfall.r >1200 & rainfall.r < 1400
windows(); plot(crainfall.r)

crainfall.im <- as.im(crainfall.r)


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - roads:

roaden.im <- as.im(roaden.r)

windows(); plot(roaden.im, axes = TRUE)

roaden.rho <- rhohat(object = dat.ppp, covariate = roaden.im)

windows(); plot(roaden.rho, xlab = "Road density (metres per square km)", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - habitat:

unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\vector\\AU_Qld_koala_habitat_suitability.zip")
habit.shp <- readShapePoly("AU_Qld_koala_habitat_suitability_MGA56.shp")
proj4string(habit.shp) <- CRS("+init=epsg:28356") 

names(habit.shp)
table(habit.shp$suitabilit)

# Take habitat suitability area 3 only:
id <- habit.shp$suitabilit == 3
habit.shp <- habit.shp[id,]

windows(); plot(habit.shp)

# Work out distance to the nearest habitat 3 area:
habit.w <- as(as(habit.shp, "SpatialPolygons"), "owin")
windows(); plot(habit.w)

habit.im <- distmap(habit.w)
windows(); plot(habit.im, axes = TRUE)

# Turn habit.im into a raster and re-sample to match the dimensions of the other rasters:

# Raster to im (again!):

# What is the nature of the association between koala sight locations and habitat 3 areas?
habit.rho <- rhohat(object = dat.ppp, covariate = habit.im)
windows(); plot(habit.rho, xlab = "Distance to the nearest habitat area 3", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - water:

water.im <- as.im(water.r)
windows(); plot(water.im, axes = TRUE)

water.rho <- rhohat(object = dat.ppp, covariate = water.im)
windows(); plot(water.rho, xlab = "Available water capacity (units)", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - temp:

temp.im <- as.im(temp.r)
windows(); plot(temp.im, axes = TRUE)

temp.rho <- rhohat(object = dat.ppp, covariate = temp.im)
windows(); plot(temp.rho, xlab = "Average annual temperature (degrees Celcius)", main = "")


# =================================================================================================================================
# Poisson point process model:

# Null model:
dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409


# Saturated model:
dat.ppm01 <- ppm(dat.ppp, trend = ~ tpo + sbd + water + nitro + roaden + hpop + temp + elev, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm01)$coefs.SE.CI
tdat1 <- data.frame(summary(dat.ppm01)$coefs.SE.CI)
tdat1[order(tdat$Zval),] 
AIC(dat.ppm01) # 3466.415


# Drop water:
dat.ppm03 <- ppm(dat.ppp, trend = ~ tpo + sbd + nitro + roaden + hpop + temp + elev, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm03)$coefs.SE.CI
tdat2 <- data.frame(summary(dat.ppm03)$coefs.SE.CI)
tdat2[order(tdat2$Zval),] 
AIC(dat.ppm03) # 3457.997


# Drop nitor:

dat.ppm04 <- ppm(dat.ppp, trend = ~ tpo + sbd + roaden + hpop + temp + elev, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm04)$coefs.SE.CI
tdat4 <- data.frame(summary(dat.ppm04)$coefs.SE.CI)
tdat4[order(tdat4$Zval),] 
AIC(dat.ppm04) # 3456.439


# Drop tpo:

dat.ppm05 <- ppm(dat.ppp, trend = ~ sbd + roaden + hpop + temp + elev, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm05)$coefs.SE.CI
tdat5 <- data.frame(summary(dat.ppm05)$coefs.SE.CI)
tdat5[order(tdat5$Zval),] 
AIC(dat.ppm05) # 3455.217


# Drop elev:

dat.ppm06 <- ppm(dat.ppp, trend = ~ sbd + roaden + hpop + temp, covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm06)$coefs.SE.CI
tdat6 <- data.frame(summary(dat.ppm06)$coefs.SE.CI)
tdat6[order(tdat6$Zval),] 
AIC(dat.ppm06) # 3455.217

# The lurking variable plot gives you an idea of mis-specification of spatial trend:
windows(); diagnose.ppm(dat.ppm06)


# Now add a Geyer spatial interaction term:
dat.ppm07 <- ppm(dat.ppp, trend = ~ sbd + roaden + hpop + temp, interaction = Geyer(r = 5000, sat = 10),covariates = list(tpo = tpo.im, sbd = sbd.im, water = water.im, nitro = nitro.im, roaden = roaden.im, hpop = hpop.im, temp = temp.im, elev = elev.im))
summary(dat.ppm07)$coefs.SE.CI

AIC(dat.ppm07) # 2040.951

windows(); diagnose.ppm(dat.ppm05)


# QQ plot of residuals from fitted model. Useful for checking distributional assumptions, particularly assumptions re spatial dependence. Careful - this takes ages to run. To reduce time, set nsim to 10.
dat.qq05 <- qqplot.ppm(fit = dat.ppm05, nsim = 10)
windows(); plot(dat.qq05, xlab = "Raw residuals: mean quantile of simulations", ylab = "Data quantile")
# epi.saveplot("ZM_fmd_ppm04_qqplot")
# Observed values don't completely lie within the boundaries of the simulation envelopes.


# =================================================================================================================================

# Model with spatial interaction term:
dat.pred05 <- predict(dat.ppm05); pred <- dat.pred05

# Express predicted koala intensity in terms of the number of koalas per hectare:
pred$v <- pred$v * 1E04
summary(as.vector(pred$v))
windows(); hist(as.vector(pred$v))


# Predictions from the model:
breaks <- seq(from = 0, to = 0.0015, length = 6)
col <- brewer.pal(n = 5, name = "Reds")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = pred$xcol, y = pred$yrow, z = t(pred$v), zlim = c(0, 0.0015), col = col, breaks = breaks, add = TRUE)
points(x = acsel[,1], y = acsel[,2], pch = 16, cex = 0.75, col = "black")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)

# Actual kernel smoothed data:
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


# =================================================================================================================================
# Adaptive kernel smoothing using sparr

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
