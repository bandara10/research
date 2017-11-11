library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
list<- list.files( path="wartondata\\wd", pattern="\\.tif$", full.names = TRUE) 
stack.r <- stack(list)
# setting resolution to .5 to this stack results in 200000 rows and colomns
valuves.rr <- as.data.frame(stack.r, xy= TRUE, na.rm=T)
valuves.rr$X<- valuves.rr$x/1000
valuves.rr$Y<- valuves.rr$y/1000

valuves.rrT <- valuves.rr[c(9,10 ,6,7,8,3,4,5)]

a<- valuves.rrT[c(1,2,3)]
b <- valuves.rrT[c(1,2,4)]
c <- valuves.rrT[c(1,2,5)]
d<- valuves.rrT[c(1,2,6)]
e <- valuves.rrT[c(1,2,7)]
f <- valuves.rrT[c(1,2,8)]
#rasterise each varibale with resolution 1 and then set resolution to .5
habit3pc <- rasterFromXYZ(as.data.frame(a)[, c("X", "Y", "habit3pc")], res=1)
sbd <- rasterFromXYZ(as.data.frame(b)[, c("X", "Y", "sbd")], res=1)
temp <- rasterFromXYZ(as.data.frame(c)[, c("X", "Y", "temp")], res=1)
clay <- rasterFromXYZ(as.data.frame(d)[, c("X", "Y", "clay")], res=1)
distance_tertiaryandlink <- rasterFromXYZ(as.data.frame(e)[, c("X", "Y", "distance_tertiaryandlink")], res=1)
elev <- rasterFromXYZ(as.data.frame(f)[, c("X", "Y", "elev")], res=1)
 # create a stack and set resolution to .5
raster.s <- stack(habit3pc,sbd , temp, clay , distance_tertiaryandlink, elev )
res(raster.s) <- .5
plot(raster.s[[3]])
#quad.1 = sample.quad(env.grid = env, sp.scale = 1, file = "Quad")

#load("TestPPM.RData")
#load("SpEnvData.RData")
# create a dataframe from one raster
env <- as.data.frame(raster.s, xy=TRUE, na.rm=T)
colnames(env )[1] <- "X"
colnames(env )[2] <- "Y"
#env [is.na(env )] <- 0

# # extract valuves from raster stack to this dataframe
 env1 <- env[c(1,2)]
# #for one raster
 rt <- raster.s[[2]]
 res(rt) <- .5
 
 presvals <- extract(raster.s, env1)
test <- as.data.frame(presvals)
test1 <- cbind(test, env1)
# #xydatan <- env[c(1,2)]
# # raster from xy data
sbd <- rasterFromXYZ(as.data.frame(test1)[, c("X", "Y", "temp")])
# # Get koala data XY
kolaxy <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom
kolaxy2 <- subset(kolaxy, X > 442 & X < 540)
kolaxy3 <- subset(kolaxy2, Y > 6902 & Y < 7000)




# xy divide by 1000
#kolaxy4<- kolaxy3[sample(nrow(kolaxy3), 250), ]

stand.distance_tertiaryandlink = scale.default(env $distance_tertiaryandlink, center = TRUE, scale = TRUE) #standari
env$distance_tertiaryandlink = stand.distance_tertiaryandlink


#########
ppmForm = ~  poly(habit3pc, sbd, temp , clay,elev, degree = 1)
ppmFit = ppmlasso(ppmForm, sp.xy = kolaxy3, env.grid = env, sp.scale = 1)

pred.biasCorrect = predict(ppmFit, newdata=env)
predictions <- cbind(xydatan, pred.biasCorrect)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect")])
plot(xy.rr, las=0)
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = kolaxy4, env.grid = env, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")