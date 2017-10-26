library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):

mydatasighting <- read.csv("sightingdata\\mydatasighting_cleaned.csv", header = TRUE)
mydata <- data.frame(mydatasighting)
names(mydata)
# Subset the data: for three year catogories.
mydata <- mydata[c("LAT","LNG","X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))
mydata99 <- subset(mydata, yearnew < 2000) #select 1997-1999/ then yearnew>=2000|yearnew<2003
names(mydata99) <- tolower(names(newdata99))
head(mydata99)
#2000 to 2002 data
mydata02 <- subset(mydata, yearnew > 1999 & yearnew < 2003) 
# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)
head(mydata02)
mydata05 <- subset(mydata, yearnew > 2002 & yearnew < 2006) 
head(mydata05)
mydata08 <- subset(mydata, yearnew > 2005 & yearnew < 2009) 
head(mydata08)
mydata11 <- subset(mydata, yearnew > 2008 & yearnew < 2012) 
head(mydata11)
mydata14 <- subset(mydata, yearnew > 2011) 
head(mydata14)

--------------------------------------------------------------------------------------------------------------------------
# Analyse data for 1997-1999:
  # Bring in MGA56 study area boundary map:
  unzip("vector\\AU_Qld_study_area.zip")
studyareaa.shp <- readShapePoly("vector\\newstudyarea56.shp")
proj4string(studyareaa.shp) <- CRS("+init=epsg:28356") 
windows();plot(studyareaa.shp)

id <- mydata99$lat > -30 # remove record with wrong coordinate
mydata99 <- mydata99[id,]
head(mydata99)
# mydata in lat-lon:
coords <- SpatialPoints(mydata99[, c("lng", "lat")])

mydata.99 <- SpatialPointsDataFrame(coords, mydata99)
proj4string(mydata.99) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.99, CRS("+init=epsg:28356"))
windows();plot(mydata.mga)
# Add the MGA coordinates to mydata:
mydata99$x <- coordinates(mydata.mga)[,1]
mydata99$y <- coordinates(mydata.mga)[,2]

windows(); plot(mydata.mga, cex = 0.3, col = 'blue', pch = 15)
--------------------------------------------------------------------------
  mydata.r <- raster(studyareaa.shp)
windows();plot(studyarea.shp)
res(mydata.r) <- 1000
mydata.r <- extend(mydata.r, extent(mydata.r) + 1000)#
----------------------------------------------------------------------------------------------------------------------
source("Lib_DistEstimatesToItself.r")
# start distance based selection method
mydata.n <- mydata.mga[c(3,4)]
mydata.n = as.data.frame(mydata.n)
myInPtDF = SpatialPointsDataFrame(mydata.n[c("x","y")], mydata.n)
windows();plot(myInPtDF, axes = T)
myInPtDF$disttoitself = Lib_DistEstimatesToItself(myInPtDF$x, myInPtDF$y)
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 200)
windows();plot(myInPtDF2, axes = T)
acsel= myInPtDF2[,1:2]

-------------------------------------------------------------------------------------
  # Sample points from within each grid cell:
  acsel <- gridSample(acsel, mydata.r, n = 1)
mydata.p <- rasterToPolygons(mydata.r)
windows(); plot(mydata.p, border = 'gray')
# Selected points in red:
points(acsel, cex = 0.5, col = 'red', pch = 15)

-----------------------------------------------------------
# SP <- as(dstudyarea.shp, "SpatialPolygons")
#W <- as(SP, "owin")
#plot(W)

mydata99 <- na.omit(mydata99)

mydata99 <-unique(mydata99)
any(duplicated(mydata99))
mydata99$yearnew <- as.factor(mydata99$yearnew)
east<- mydata99$x
north<- mydata99$y
yearnew<-mydata99$yearnew
head(mydata99)
x <- ppp(east,north, c(445607,531146), c(6911470, 7003109),marks=yearnew)
unitname(x) <- c("metre", "metres")
#ABC<-ppp(east, north, window = W)
#windows();plot(ABC)
plot(x)
C<-plot(split(x))
data1997<- split(x)$"1997"#Full dataset can be splitted like this for the analysis
plot(data1997)
D <- density(x,sigma=1000)
E <- density(x,sigma=1000)
plot(D, col=grey(seq(1, 0, length = 512)))
contour(density(x, 10), axes = FALSE)
#adaptive density
A <- by(x, FUN = adaptive.density)
plot(A)
######
DDD<-as.matrix(D)
persp(D,theta = 30)
###
D[locator(2)]
logZ <- D(log(C))
---------------------------------------------------
dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp.shp")

windows(); plot(dstudyarea.shp, axes = TRUE)
studyareaa.shp <- readShapePoly("vector\\newstudyarea56.shp")
projection(studyareaa.shp) <- CRS("+init=epsg:28356")
# Reproject to MGA-56:

plot(studyareaa.shp, add = TRUE)
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(x = mydata$x, y = mydata$y,cex = 0.2, col = 'blue', pch = 15)

---------------------------------------------------------------------------------------------------------------------------
  # Create an observation window which is an intersection of the square study area boundaries and the detailed study area:
  source("owin2sp_source.r")

dat.w <- as(as(studyareaa.shp, "SpatialPolygons"), "owin")
#observation window is a map. try this windw.
#SP <- as(dstudyarea.shp, "SpatialPolygons")
#W <- as(SP, "owin")
#dat.spw <- owin2SP(dat.w)

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

windows(); plot(studyareaa.shp, axes = TRUE)
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
xlim <- bbox(studyareaa.shp)[1,]; ylim <- bbox(studyareaa.shp)[2,]
x.points <- seq(from = 300000, to = 540000, by = 20000); x.lab <- x.points / 100
y.points <- seq(from = 6900000, to = 7200000, by = 20000); y.lab <- y.points / 1000

breaks <- seq(from = 0, to = 0.004, length = 6)
col <- brewer.pal(n = 5, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
points(acsel, cex = 0.4, col = 'red', pch = 15)
---------------------------------------------------------------------------------------------------------------------------------------------
  # Actual kernel smoothed data:
  breaks <- seq(from = 0, to = 0.004, length = 5)
col <- brewer.pal(n = 4, name = "Blues")

windows(); plot(x = xlim, y = ylim, type = "n",main="Actual kernel smoothed data", axes = TRUE, xlab = "Easting (km)", ylab = "Northing (km)", xaxt = "n", yaxt = "n")
image(x = dat.den$xcol, y = dat.den$yrow, z = t(dat.den$v), zlim = c(0, 0.004), col = col, breaks = breaks, add = TRUE)
# points(x = mydata$x, y = mydata$y, pch = 1, cex = 0.25, col = "grey")
plot(dstudyarea.shp, col = "transparent", border = "black", add = TRUE); box()
axis(side = 1, at = x.points, labels = x.lab, tick = TRUE)
axis(side = 2, at = y.points, labels = y.lab, tick = TRUE)
metre(xl = xlim[1] + 5000, yb = ylim[1], xr = xlim[1] + 7500, yt = ylim[1] + 40000, lab = breaks, cols = col, shift = 0, cex = 0.75)
savePlot(filename = "Actual kernel smoothed data",type = c( "png"),device = dev.cur())
-----------------------------------------------------
  #To plot average distance as a function of neighbor order for the first 100 closest neighbors:
  #The bottom axis shows the neighbor order number and the left axis shows the average distance.
  ANN <- apply(nndist(dat.ppp, k=1:150),2,FUN=mean)
plot(ANN ~ eval(1:150), type="b", main=NULL )
#k test
p<-dat.ppp

K <- Kest(p)
plot(K, main=NULL)
K
# 
E<-envelope(,Kest,nsim=3)#didnt converge
plot(E)
#L function
L <- Lest(p, main=NULL)
plot(L)
plot(L, . -r ~ r, main="L function with L expected")
########################
plot(elev.r)
line <- drawLine()
alt <- extract(elev.r, line, along = TRUE)
plot(alt[[1]], type = 'l', ylab = "Altitude (m)")
library(geosphere)
## Calculate great circle distance between the two ends of the line

