library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files( path="wartondata", pattern="\\.tif$", full.names = TRUE) 
myenv.stack <- stack(myenv)
### crop to extent to get square area
plot(myenv.stack[[1]])
newex <- drawExtent()
stack.c <- crop(myenv.stack, newex)
plot(stack.c)
myenv.stack2 <- stack.c

d <- (myenv.stack2[[1]])
e <- (myenv.stack2[[2]])
f <- (myenv.stack2[[3]])

dd <- rasterToPoints(d, spatial = TRUE)
ee <- rasterToPoints(e, spatial = TRUE)
ff <- rasterToPoints(f, spatial = TRUE)

dd <- as.data.frame(dd)
ee<- as.data.frame(ee)
ff<- as.data.frame(ff)
rasterXY <- cbind(ee$x, ee$y, dd$distance_tertiaryandlink, ee$elev, ff$temp)
rasterXYa <- as.data.frame(rasterXY)
colnames(rasterXYa)[1] <- "X"
colnames(rasterXYa)[2] <- "Y"
colnames(rasterXYa)[3] <- "DIS_TER"
colnames(rasterXYa)[4] <- "ELEV"
colnames(rasterXYa)[5] <- "TEMP
 erXYa <-  my rasterXYa$X/1000
 asterXYa$y/1000
 asterXY <- rasterxy [c(6,7,3,4,5)]




# quadrature points not working.
quad.1 = sample.quad(env.grid =rasterXY , sp.scale = 1, file = "Quad")


#koala data
kolaxy <- read.csv("wartondata\\kolaxy.csv", header = TRUE)
kolaxy <- subset(kolaxy, y> 6901098 & y < 7000000, select=x:y) # limit y selection

# minimum distace between kola 1000m
kolaxy$disttoitself = Lib_DistEstimatesToItself(kolaxy$x, kolaxy$y)
kolaxy2 = subset(kolaxy, kolaxy$disttoitself > 400)
kxy <- as.data.frame(kolaxy2[c(1,2)])
kxy$X <- kxy$x/1000 ; kxy$Y <- kxy$y/1000
kXY <- kxy [c(3,4)]

write.csv(kXY, "koalaxy.csv", row.names=FALSE)
plot(kxy, add=TRUE)

#read data from folder directly
koalaxy <- read.csv("wartondata\\koalaxy.csv", header = TRUE)

env <- read.csv("wartondata\\env.csv", header = TRUE)

#########
ppmForm = ~  poly(elev,temp, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = koalaxy, env.grid = env, sp.scale = s1)
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

#####
