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
rasterxy <- na.omit(rasterxy)
write.csv(rasterxy, "wartonmatch.csv")
# convert xy to km
rasterxy$X <- rasterxy$x/1000 ; rasterxy$Y <- rasterxy$y/1000
rasterXY <- rasterxy [c(5,6,3,4)]
colnames(rasterXY)[3] <- "ELEV"
colnames(rasterXY)[4] <- "TEMP"
write.csv(rasterXY, "wartonmatch.csv")
#koala data
kolaxy <- read.csv("wartondata\\kolaxy.csv", header = TRUE)
kolaxy <- subset(kolaxy, y> 6901098 & y < 7000000, select=x:y) # limit y selection

# minimum distace between kola 1000m
kolaxy$disttoitself = Lib_DistEstimatesToItself(kolaxy$x, kolaxy$y)
kolaxy2 = subset(kolaxy, kolaxy$disttoitself > 400)
kxy <- as.data.frame(kolaxy2[c(1,2)])
kxy$X <- kxy$x/1000 ; kxy$Y <- kxy$y/1000
kXY <- kxy [c(3,4)]

write.csv(kXY, "koalaxy.csv")
plot(kxy, add=TRUE)

#########
ppmForm = ~  poly(elev,temp, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = kXY, env.grid = rasterxy, sp.scale = 1)
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