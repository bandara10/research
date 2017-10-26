library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer)

# setwd("D:\\TEMP")

mydatasighting <- read.csv("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\mydatasighting.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)

# Subset the data:
mydata <- mydata[c("LAT","LNG","X","Y", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 1999, select = c("x","y"))

# Bring in MGA56 study area boundary map:
unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_study_area.zip")

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
bandwidth <- 2500

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

windows();
plot(dat.den)

# Set x and ylims for plot:
xlim <- bbox(studyarea.shp)[1,]; ylim <- bbox(studyarea.shp)[2,]

x.points <- seq(from = 440000, to = 540000, by = 20000); x.lab <- x.points / 1000
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n",)
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = terrain.colors(n = 10), add = TRUE)
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)


# =================================================================================================================================
# rhohat analyses

# Soil pH raster:
unzip("F:\\Maps\\Oceania\\Australia\\raster\\AU_CSIRO_soil_pH.zip")
ph.r <- raster("AU_CSIRO_ph_0_30-wgs84.tif")
projection(ph.r) <- CRS("+init=epsg:4326")

windows(); plot(ph.r, axes = TRUE)

# Crop the raster to South Australia (WGS84) extent:
ph.r <- crop(ph.r, extent(sawgs84.shp), snap = "out")
windows(); plot(ph.r)

# Re-project the WGS84 raster map to EPSG 3107 (takes a while to run):
ph.r <- projectRaster(ph.r, crs =  CRS("+init=epsg:3107"))

windows(); plot(ph.r, col = topo.colors(n = 15), axes = TRUE)

# Resample the pH raster:
dim(ph.r)

tcas.denr <- raster(tcas.den)
.ph.r <- resample(x = ph.r, y = tcas.denr, method = "ngb")

windows(); plot(sa.shp, axes = TRUE)
plot(.ph.r, add = TRUE)

# Convert .ph.r into a spatstat image object:
ph.im <- as.im(.ph.r)

windows(); plot(ph.im, axes = TRUE)
plot(tdat.w, add = TRUE)

# What is the nature of the association between fracture-positive farms ans soil pH? 
# Note use of baseline to adjust for the spatial distribution of pop herds. See page 83 -- 85 Baddeley and page 181 of Baddeley et al (2015):

# Soil pH and case farms (only):
ph.rho <- rhohat(object = tcas.ppp, covariate = ph.im)

windows(); plot(ph.rho, xlim = c(4,10), xlab = "Soil pH", main = "")

# Soil pH and cases farms, adjusting for the spatial distribution of the sheep farm population at risk:
ph.rho <- rhohat(object = tcas.ppp, covariate = ph.im, baseline = tpop.den, method = "transform")

windows(); plot(ph.rho, xlim = c(4,10), xlab = "Soil pH", main = "")


# ---------------------------------------------------------------------------------------------------------------------------------
# Soil phosphorous:

# https://data.csiro.au/dap/landingpage?execution=e1s6

# Soil pH raster:
unzip("F:\\Maps\\Oceania\\Australia\\raster\\AU_CSIRO_soil_pH.zip")
p.r <- raster("PTO_000_005_05_N_P_AU_NAT_C_20140801.tif")
projection(p.r) <- CRS("+init=epsg:4326")

windows(); plot(p.r, axes = TRUE)

# Crop the raster to South Australia (WGS84) extent:
p.r <- crop(p.r, extent(sawgs84.shp), snap = "out")
windows(); plot(p.r)

# Re-project the WGS84 raster map to EPSG 3107 (takes a while to run):
p.r <- projectRaster(p.r, crs =  CRS("+init=epsg:3107"))

windows(); plot(p.r, col = topo.colors(n = 15), axes = TRUE)

# Resample the pH raster:
dim(p.r)

tcas.denr <- raster(tcas.den)
.p.r <- resample(x = p.r, y = tcas.denr, method = "ngb")

windows(); plot(sa.shp, axes = TRUE)
plot(.p.r, add = TRUE)

# Convert .ph.r into a spatstat image object:
p.im <- as.im(.p.r)

windows(); plot(p.im, axes = TRUE)
plot(tdat.w, add = TRUE)

# What is the nature of the association between fracture-positive farms ans soil pH? 
# Note use of baseline to adjust for the spatial distribution of pop herds. See page 83 -- 85 Baddeley and page 181 of Baddeley et al (2015):

# Soil phosphorous and case farms (only):
p.rho <- rhohat(object = tcas.ppp, covariate = p.im)

windows(); plot(p.rho, xlab = "Soil phosporous", main = "")

# Soil pH and cases farms, adjusting for the spatial distribution of the sheep farm population at risk:
p.rho <- rhohat(object = tcas.ppp, covariate = p.im, baseline = tpop.den, method = "transform")

windows(); plot(p.rho, xlab = "Soil phosporous", main = "")


# =================================================================================================================================
# Adaptive kernel smoothing using sparr

# Load the R workspace saved earlier:
load("D:\\Contracts\\Australia\\UoA\\Rib_fractures\\data\\EAS_farm_loc_mapped_day1_end.RData")

x.points <- seq(from = 400000, to = 1600000, by = 200000); x.lab <- x.points / 1000
y.points <- seq(from = 1400000, to = 2600000, by = 200000); y.lab <- y.points / 1000

windows(); par(mfrow = c(1,1))
plot(sa.shp, col = "transparent", border = "black", xaxt = "n", yaxt = "n", xlab = "Easting (m)", ylab = "Northing (m)", axes = TRUE); box()
points(tpop.ppp, pch = 16, col = "grey", cex = 0.5)
points(tcas.ppp, pch = 16, col="black", cex = 0.5)
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
legend(x = "bottomleft", legend = c("Positive", "Negative"), pch = 16, col = c("black", "grey"), bty = "n")


# Adaptive smoothing and dividing surfaces to make relative risk density map:
library(spatialkernel); library(sparr); library(rgl); library(spatialkernel)

ncas <- tcas.ppp$n
npop <- tpop.ppp$n

# Calculate a pilot bandwidth (using leave-one-out least-squares cross-validation) and a global bandwidth (using maximal smoothing) and then feed these values to the bivariate.density function:

pop.pilot <- LSCV.density(tpop.ppp)
pop.global <- OS(tpop.ppp, nstar = NULL)
pop.den <- bivariate.density(data = tpop.ppp, pilotH = pop.pilot, globalH = pop.global)

cas.pilot <- LSCV.density(tcas.ppp)
cas.global <- OS(tcas.ppp, nstar = NULL)
cas.den <- bivariate.density(tcas.ppp, pilotH = cas.pilot, globalH = cas.global, gamma = pop.den$gamma)

summary(as.vector(cas.den$Zm))
cas.den$Zm

cas.den$Zm <- cas.den$Zm / 0.000001 / 0.000001
hist(as.vector(cas.den$Zm)) 

pop.den$Zm <- pop.den$Zm / 0.000001 / 0.000001
hist(as.vector(pop.den$Zm)) 

summary(as.vector(pop.den$Zm))

# Take the ratio of case.den to pop.den to return a prevalence surface. Use the tolerance function to identify areas of significantly elevated prevalence:
risk.den <- risk(f = cas.den, g = pop.den, log = TRUE, plotit = FALSE)
summary(risk.den)

risk.tol <- tolerance(rs = risk.den, pooled = pop.den, test = "upper")
# Image plot showing the spatial distribution of prevalence.
# Contour lines delineate areas of significantly elevated log prevalence:

# Save the workspace because these analyses take forever to run:
save.image("EAS_farm_loc_mapped_day1_end.RData")


summary(as.numeric(risk.den$rsM))
dim(risk.den$rsM)

breaks <- seq(from = -2, to = 2, by = 0.5)
col <- rev(brewer.pal(n = 8, name = "RdBu"))

windows(); par(pty = "s")
plot(x = xlim, y = ylim, type = "n", axes = TRUE, xaxt = "n", yaxt = "n", xlab = "Easting (km)", ylab = "Northing (km)")
axis(1, at = xlabs * 1000, labels = xlabs)
axis(2, at = ylabs * 1000, labels = ylabs)
image(x = risk.den$f$X, y = risk.den$f$Y, z = risk.den$rsM, zlim = c(-2, 2), col = col, breaks = breaks, add = TRUE)
plot(sa.shp, col = "transparent", border = "black", add = TRUE); box()
# points(tpop.ppp, pch = 16, col = "grey", cex = 0.5)
# points(tcas.ppp, pch = 16, col = "black", cex = 0.5)
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
contour(x = risk.tol$X, y = risk.tol$Y, z = risk.tol$P, levels = 0.05, drawlabels = FALSE, add = TRUE, lwd = 1, lty = 2, col = "black")
plot(dat.proj.w3, add = TRUE)
metre(xl = xlim[1], yb = ylim[1], xr = xlim[1] + 50000, yt = ylim[1] + 500000, lab = breaks, cols = col, shift = 0, cex = 0.75)

# spatial model
# http://www.soilquality.org.au/factsheets/soil-ph-south-austral


# =================================================================================================================================
