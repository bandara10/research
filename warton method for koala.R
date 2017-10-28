library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files( path="wartondata", pattern="\\.tif$", full.names = TRUE) 

myenv.stack <- stack(myenv)
plot(myenv.stack)
rasterxy <- as.data.frame(myenv.stack, xy=TRUE)
coordinates(rasterxy) <- ~x+y
quad.1 = sample.quad(env.grid = rasterxy, sp.scale = 1, file = "Quad")
#extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)

#lapply outputs a list of rasters
# m is data matrix from getvalues stack.
#r is an empry raster
ll <- getValues(myenv.stack)
lll <- as.data.frame(ll)
lll <- lll[c(2,3)]

#create an empty raster
start_y <- 441829 #100
start_x <- 6901098 #200
griddf <- expand.grid(X = seq(from = start_y, by = 1, l = 10000),
                      Y = seq(from = start_x, by = 1, l = 10000))

A <- rnorm(10000, 0, 0)

envs<- cbind(griddf,A)

xy.rr <- rasterFromXYZ(as.data.frame(envs)[, c("X", "Y", "A")])
#proj4string(xy.rr) <- CRS("+init=epsg:28356")
#xy.rr[is.na(xy.rr[])] <- 0 # this raster is fine but includs the sea.

plot(xy.rr)
####
l <- lapply(1:ncol(lll), function(i) {
  setValues(xy.rr,lll[ , i])
} )

# stack the list
s <- stack(l)
plot(s[[1]])

ss <- as.data.frame(s, xy=TRUE)
#ss[is.na(ss)] <- 0
colnames(ss)[1] <- "X"
colnames(ss)[2] <- "Y"
sss <- na.omit(ss)
xydatan <- sss[c(1,2)] # to create a raster
coordinates(kolaxy) <- X+Y
plot(kolaxy, add= TRUE)


#koala data
kolaxy <- read.csv("wartondata\\kolaxy.csv", header = TRUE)
kolaxy <- kolaxy[sample(nrow(kolaxy), 80), ]
xydatan <- na.omit(xydatan)
#kolaxy <- kolaxy[,1:ncol(kolaxy)]/1000 # selects every row and 2nd to last columns
#coordinates(kolaxy)= ~ X+Y
#proj4string(kolaxy) <- CRS("+init=epsg:28356")
#plot(kolaxy)
## Pre-standardise observer bias variables
#stand.A.1=scale.default(ss$A.1, center = TRUE, scale = TRUE) #standarise
#stand.A.2=scale.default(ss$A.2, center = TRUE, scale = TRUE) #standarise
#ss$A.1 = stand.A.1
#ss$A.2 = stand.A.2

#########+ +poly (A.1, degree = 2)
ppmForm = ~  poly(A.1, A.2, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = kolaxy, env.grid = sss, sp.scale = 1)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = sss
#newEnv$A.1 = min(stand.A.1)
pred.biasCorrect = predict(ppmFit, newdata=newEnv)
#newEnv$A.1 = min(stand.A.1)
pred.biasCorrect.1 = predict(ppmFit, newdata=sss)
predictions <- cbind(xydatan, pred.biasCorrect.1)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = kolaxy,
        env.grid = sss, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")