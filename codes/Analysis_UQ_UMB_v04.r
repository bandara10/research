library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer)

setwd("D:\\TEMP")

load(file = "D:\\Contracts\\Australia\\UQ\\Koalas\\data\\AU_Qld_koala_habitat_v01.RData")

mydatasighting <- read.csv("D:\\Contracts\\Australia\\UQ\\Koalas\\data\\mydatasighting.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

table(mydata$fateid)

# Bring in MGA56 study area boundary map:
unzip("D:\\Contracts\\Australia\\UQ\\Koalas\\maps\\AU_Qld_study_area.zip")

studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 

dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

windows(); plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(x = mydata$x, y = mydata$y)

# kdate as a date:
# dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
# dat$kmonth <- format(dat$kdate, format = "%m")


# =================================================================================================================================

# Create an observation window which is an intersection of the square study area boundaries and the detailed study area:
source("D:\\Contracts\\Australia\\UoA\\Rib_fractures\\code\\owin2sp_source.r")

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
id <- inside.owin(x = mydata$x, y = mydata$y, w = dat.w)
mydata <- mydata[id,]

windows(); plot(studyarea.shp, axes = TRUE)
points(x = mydata$x, y = mydata$y)

# Make a ppp object:
dat.ppp <- ppp(x = mydata$x, y = mydata$y, window = dat.w)

windows(); 
plot(dat.ppp, axes = TRUE)

# Work out density for dat.ppp (using Diggle's edge correction):
dat.den <- density(dat.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)

# Express density as the number of koalas per square kilometre (cf number of koalas per square metre):
dat.den$v <- dat.den$v / 0.000001

windows(); plot(dat.den)

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]

x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 100, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)


# =================================================================================================================================
# rhohat analyses - elevation:

elev.r <- raster("D:\\Contracts\\Australia\\UQ\\Koalas\\maps\\elevation.asc")
projection(elev.r) <- CRS("+init=epsg:28356")

windows(); plot(elev.r, axes = TRUE)

# Crop the elevation raster to the study area extent:
elev.r <- crop(elev.r, extent(studyarea.shp), snap = "out")
windows(); plot(elev.r)

windows(); plot(elev.r, col = topo.colors(n = 15), axes = TRUE)

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

hpop.r <- raster("D:\\Contracts\\Australia\\UQ\\Koalas\\maps\\hpop.asc")
projection(hpop.r) <- CRS("+init=epsg:28356")

windows(); plot(hpop.r, axes = TRUE)

# Crop the elevation raster to the study area extent:
hpop.r <- crop(hpop.r, extent(studyarea.shp), snap = "out")
windows(); plot(hpop.r)

# Convert .ph.r into a spatstat image object:
hpop.im <- as.im(hpop.r)

windows(); plot(hpop.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
hpop.rho <- rhohat(object = dat.ppp, covariate = hpop.im)

windows(); plot(hpop.rho, xlim = c(0, 6000), xlab = "Human population density", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - annual rainfall (mm):

rainfall.r <- raster("D:\\Contracts\\Australia\\UQ\\Koalas\\maps\\rainfall.asc")
projection(rainfall.r) <- CRS("+init=epsg:28356")

windows(); plot(rainfall.r, axes = TRUE)

# Crop the elevation raster to the study area extent:
rainfall.r <- crop(rainfall.r, extent(studyarea.shp), snap = "out")
windows(); plot(rainfall.r)

# Convert .ph.r into a spatstat image object:
rainfall.im <- as.im(rainfall.r)

windows(); plot(rainfall.im, axes = TRUE)

# What is the nature of the association between koala sight locations and human population density?
rainfall.rho <- rhohat(object = dat.ppp, covariate = rainfall.im)

windows(); plot(rainfall.rho, xlab = "Annual rainfall (mm)", main = "")



# ---------------------------------------------------------------------------------------------------------------------------------
# rhohat analyses - roads:

unzip("D:\\Contracts\\Australia\\UQ\\Koalas\\vector\\AU_Qld_study_area_mask_roads.zip")
roads.shp <- readShapeLines("AU_Qld_study_area_mask_lines-MGA56.shp")
proj4string(roads.shp) <- CRS("+init=epsg:28356") 

names(roads.shp)
table(roads.shp$highway)

# Take primary, secondary and tertiary roads:
id <- roads.shp$highway == "primary" | roads.shp$highway == "secondary" | roads.shp$highway == "tertiary"
roads.shp <- roads.shp[id,]

# Crop the roads to match the extent of the study area:
roads.shp <- roads.shp[studyarea.shp, ]

# Work out distance to the nearest road:
roads.psp <- as.psp(roads.shp, W = dat.w, check = FALSE)
roads.im <- distmap(roads.psp)

windows(); plot(roads.im, axes = TRUE)

# What is the nature of the association between koala sight locations and distance to the nearest road?
roads.rho <- rhohat(object = dat.ppp, covariate = roads.im)

windows(); plot(roads.rho, xlab = "Distance to the nearest road (m)", main = "")

# Save everything we've done so far:
save.image(file = "AU_Qld_koala_habitat_v01.RData")


# =================================================================================================================================
# Adaptive kernel smoothing using sparr

# Adaptive smoothing and dividing surfaces to make relative risk density map:
library(spatialkernel); library(sparr); library(rgl); library(spatialkernel)


n1 <- sum(PBC$data$ID)
n2 <- nrow(PBC$data) - n1

pool.pilot <- CV.sm(PBC)
pool.global <- OS(PBC$data[,1:2])
pbc.pool <- bivariate.density(data = PBC$data[,1:2], pilotH = pool.pilot, globalH = pool.global, WIN = PBC$owin)
f.pilot <- CV.sm(PBC$data[PBC$data$ID == 1, 1:2])
g.pilot <- CV.sm(PBC$data[PBC$data$ID == 0, 1:2])
f.global <- g.global <- OS(PBC$data[,1:2], nstar = sqrt(n1 * n2))
pbc.case <- bivariate.density(data = PBC$data, ID = 1, pilotH = f.pilot, globalH = f.global, WIN = PBC$owin, gamma = pbc.pool$gamma)



# Calculate a pilot bandwidth (using leave-one-out least-squares cross-validation) and a global bandwidth (using maximal smoothing) and then feed these values to the bivariate.density function:

dat.pilot <- LSCV.density(dat.ppp)
dat.global <- OS(dat.ppp, nstar = NULL)
dat.sden <- bivariate.density(data = dat.ppp, pilotH = dat.pilot, globalH = dat.global, adaptive = TRUE, edgeCorrect = TRUE)

dat.sden$Zm <- dat.sden$Zm / 0.000001 / 0.000001
hist(as.vector(dat.sden$Zm)) 

summary(as.vector(dat.sden$Zm))

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
