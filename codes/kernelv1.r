library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)

# setwd("D:\\TEMP")

# from Geoscience Australia GeoData Topo layers (in GDA94 = )
# make list of projections and find GDA94 and SA Lambert

# Load EPSG codes:
#EPSG <- make_EPSG()
#EPSG[grep("GDA94", EPSG$note), 1:2]

# Bring in WGS84 map of Australian states (we use this later when we are dealing with the pH map):
unzip("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_study_area.zip")

studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
plot(studyarea.shp)
mydatasighting<-read.csv("mydatasighting.csv", header = TRUE)
mydatasightingN <- mydata[c("LAT","LNG","X","Y", "Sighting","SexId", "FateId", "LGA" , "monthnew","yearnew","Suburb","PostCode")]
myND <- subset(mydatasightingN, yearnew <=1999, select=c("X","Y"))
names(mydatasightingN)


# mainlands.shp <- readShapePoly("mainlands.shp")
mainlands.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")

proj4string(mainlands.shp)
# allocated the CRS as GDA94
proj4string(mainlands.shp) <- CRS("+init=epsg:4283") 
proj4string(mainlands.shp)

# Project to GDA94 / SA Lambert (EPSG 3107):
mainlands.shp <- spTransform(mainlands.shp,  CRS("+init=epsg:3107"))

windows(); plot(mainlands.shp, axes = TRUE)

# Create a SA boundary map:
names(mainlands.shp)

mainlands.shp$STATE
id <- mainlands.shp$STATE == "SA"
sa.shp <- mainlands.shp[id,]

windows(); plot(sa.shp, col = "red")

# Read in fracture data:
dat <- read.csv("SA_rib_fractures_2007-2015.csv", header = TRUE, sep = ",")

# kdate as a date:
dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
dat$kmonth <- format(dat$kdate, format = "%m")


# =================================================================================================================================
# The data we've just read in is for 2007 to 2015 (inclusive). Select a year of interest and then subset the data:

# tdat <- subset(dat, kyear == 2015)
tdat <- dat

id <- complete.cases(tdat)
tdat <- tdat[id,]

# Reproject point data from lat long to projected SA Lambert 
coords <- SpatialPoints(tdat[, c("lon", "lat")])
tdat.ll <- SpatialPointsDataFrame(coords, tdat)
proj4string(tdat.ll) <- CRS("+init=epsg:4326")

# Re-project to EPSG 3107 and add coordinates to tdat.proj:
tdat.gda <- spTransform(tdat.ll,  CRS("+init=epsg:3107"))
coord <- coordinates(tdat.gda)
tdat$x <- coord[,1]; tdat$y <- coord[,2]

# Create an observation window and a marked ppp object. Plot the data:
tdat.w <- convexhull.xy(x = c(tdat$x, 1577766.2, 710874.4), y = c(tdat$y, 1232747, 1260488))

# Dilate the owin by 20 km:
tdat.w <- dilation(w = tdat.w, r = 20000)

windows(); plot(sa.shp)
points(x = tdat$x, y = tdat$y, pch = 16, cex = 0.1, col = "red")
plot(tdat.w, add = TRUE)

# Alternative. Read a polygon shape file into R and then using that shape file to create an owin object. For example:
# dat.shp <- readShapePoly("myshapefile.shp")
# dat.w <- as(as(dat.shp, "SpatialPolygons"), "owin")

# Intersect from raster package
source("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\owin2sp_source.r")
tdat.spw <- owin2SP(tdat.w)

# Set projection of the owin as GDA94 / SA Lambert:
proj4string(tdat.spw)
proj4string(tdat.spw) <- CRS("+init=epsg:3107") 
tdat.spw <- raster::intersect(x = sa.shp, y = tdat.spw)

# Convert the sp object into a window:
tdat.w <- as.owin(tdat.spw)

windows(); plot(sa.shp)
plot(tdat.w, border = "blue", lwd = 2, add = TRUE)

# writePolyShape(tdat.spw, "AU_SA_convex_hull")
# write.table(tdat, 'sheep_dat.csv', row.names = FALSE, col.names = TRUE, sep = ",")

# tdat.spw <- readShapePoly("AU_SA_convex_hull_edit")
# tdat.w <- as.owin(tdat.spw)


# =================================================================================================================================
# Spatial kernel density surfaces
library(spatialkernel); library(splancs); library(RColorBrewer)

# dimyx = c(200, 200)
dimyx = c(200,200)
bandwidth <- 10000 # 10 kilometres

windows(); plot(sa.shp, axes = TRUE)
points(x = tdat$x[tdat$nfrac > 0], y = tdat$y[tdat$nfrac > 0], pch = 16, cex = 1)

# Case farms:
tcas.ppp <- ppp(x = tdat$x[tdat$nfrac > 0], y = tdat$y[tdat$nfrac > 0], window = tdat.w)

# Pop farms:
tpop.ppp <- ppp(x = tdat$x, y = tdat$y, window = tdat.w)

windows(); par(mfrow = c(1,2))
plot(tcas.ppp, axes = TRUE); plot(tpop.ppp, axes = TRUE)

# Work out density for case.ppp and pop.ppp (using Diggle's edge correction). Increase the bandwidth for the population at risk by a factor of 1.5 (because rhohat doesn't work):
tcas.den <- density(tcas.ppp, sigma = bandwidth, diggle = TRUE, dimyx = dimyx)
tpop.den <- density(tpop.ppp, sigma = bandwidth * 1.5, diggle = TRUE, dimyx = dimyx)

windows(); par(mfrow = c(1,2))
plot(tcas.den); plot(tpop.den)

# Risk density surface - create a matrix of fracture risk by dividing cas.den by pop.den:
tz <- (t(tcas.den$v) / t(tpop.den$v))
summary(as.vector(tz))
hist(as.vector(tz))

# tz <- ifelse(is.infinite(tz), 0, tz)
# tz <- ifelse(tz < 0, 0, tz)
# tz <- ifelse(tz > 0.10, NA, tz)

# Perspective plot - good for reality checking:
persp(tz, zlim = c(0,1))
image(tz)

# Set x and ylims for plot:
xlim <- bbox(sa.shp)[1,]; ylim <- bbox(sa.shp)[2,]

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)")
image(x = tcas.den$xcol, y = tcas.den$yrow, z = tz, zlim = c(0,1), col = terrain.colors(n = 10), add = TRUE)
plot(sa.shp, col = "transparent", border = "black", add = TRUE); box()


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
