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
# plot(myenv.stack[[1]])
# newex <- drawExtent()
# stack.c <- crop(myenv.stack, newex)
# plot(stack.c)
# myenv.stack2 <- stack.c

d <- (myenv.stack[[1]])
k <- (myenv.stack[[2]])
f <- (myenv.stack[[3]])

e <- as.data.frame(d, xy=TRUE)
ee <- e[c(1,2)]
ee$X  <- ee$x/1000 
ee$Y <- ee$y/1000
eee <- ee [c(3,4)]
#bind this with raster valuves

# sss <- as.data.frame(myenv.stack, centroids= TRUE,  xy = TRUE)
# a <- sss[c(1,2,4)]
# ar <- rasterFromXYZ(as.data.frame(a)[, c("x", "y", "elev")])

dd <- as.data.frame(d)
kk<- as.data.frame(k)
ff<- as.data.frame(f)
rasterXY <- cbind(eee$X, eee$Y, dd$distance_tertiaryandlink, kk$elev, ff$temp)
rasterXYa <- as.data.frame(rasterXY)
colnames(rasterXYa)[1] <- "X"
colnames(rasterXYa)[2] <- "Y"
colnames(rasterXYa)[3] <- "ELEV"
colnames(rasterXYa)[4] <- "TEMP"

#rasterXYaX <- rasterXYa$X/1000 
#rasterXYa$Y <- rasterXYa$Y/1000
arr <- rasterFromXYZ(as.data.frame(rasterXYa)[, c("X", "Y", "TEMP")], res=1)
arrr <- rasterFromXYZ(as.data.frame(rasterXYa)[, c("X", "Y", "ELEV")], res=1)
st <- stack(arr, arrr)
res(st) <- .5
plot(st[[1]])

stt <- as.data.frame(st, xy=TRUE, na.rm=T)
colnames(stt)[1] <- "X"
colnames(stt)[2] <- "Y"
stt[is.na(stt)] <- 0
xydatan <- stt[c(1,2)]
# quadrature points not working.
quad.1 = sample.quad(env.grid =stt , sp.scale = 1, file = "Quad")


#koala data
kolaxy <- read.csv("wartondata\\kolaxy.csv", header = TRUE)
#kxy <- subset(kolaxy, y> 6951098 & y < 6991098, select=x:y) # limit y selection

kolaxy <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom

# # minimum distace between kola 1000m
# kolaxy$disttoitself = Lib_DistEstimatesToItself(kolaxy$x, kolaxy$y)
# kolaxy2 = subset(kolaxy, kolaxy$disttoitself > 100)
# kxy <- as.data.frame(kolaxy2[c(1,2)])
kxy$X <- kxy$x/1000 ; kxy$Y <- kxy$y/1000
kXY <- kxy [c(3,4)]

#kXY1<- kXY[sample(nrow(kXY), 100), ]
# select points fall within the raster. rasterise a variable(arr). extract valuves with points. remove na records.
abc <- extract(arr, kXY)
abc1 <- as.data.frame(abc)
abc2 <- cbind(abc1, kXY)
abc2 <- na.omit(abc2)
koalaxyT <- abc2[c(2,3)]
# #read data from folder directly
#koalaxy <- read.csv("wartondata\\koalaxy.csv", header = TRUE)
# 
# env <- read.csv("wartondata\\env.csv", header = TRUE)

#########
ppmForm = ~  poly(TEMP,ELEV, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = kolaxy, env.grid = stt, sp.scale = 1)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = stt
#newEnv$A.1 = min(stand.A.1)
pred.biasCorrect = predict(ppmFit, newdata=newEnv)
#newEnv$A.1 = min(stand.A.1)
#pred.biasCorrect.1 = predict(ppmFit, newdata=sss)
predictions <- cbind(xydatan, pred.biasCorrect)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = koalaxyT,
        env.grid = stt, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")

#####
