library(ppmlasso)
library(spatstat)
library(raster)

# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
load(file="koaladata.RData", .GlobalEnv)
str(koalaData)
backg.env = koalaData$env
# take a random sample
backg.env[is.na(backg.env)] <- 0 # set NA to 0 t0 test
str(backg.env)
koalaxy <- koalaData$koalaxy
xyvariables <- backg.env[c(1,2)] # get x y colomns to be cbind later with predicted valuves.


coordinates(xyvariables) <- ~x+y
coordinates(koalaxy) <- ~x+y
plot(xyvariables)
plot(koalaxy, col="red", add=TRUE)
#standarise observer bias variables
stand.rain = scale.default(backg.env$rain , center = TRUE, scale = TRUE)
stand.temp  = scale.default(backg.env$temp  , center = TRUE, scale = TRUE)



stand.lot_density = scale.default(backg.env$lot_density , center = TRUE, scale = TRUE)
stand.roadk = scale.default(backg.env$roadk , center = TRUE, scale = TRUE) 
 
backg.env$lot_density = stand.lot_density 
backg.env$roadk = stand.roadk
backg.env$rain = stand.rain 
backg.env$temp = stand.temp

# fit a ppm @ res 1k
ppmForm = ~ poly(rain,temp, clay, degree = 2)  
ppmFit =ppmlasso(ppmForm, sp.xy = koalaxy,  env.grid = backg.env , sp.scale = 1)

newEnv = koalaData$vars

# predict using model based control over bias at minimum value at distance_pedestrian 
pred.biasCorrect = predict(ppmFit, newdata=newEnv)
predictions <- cbind(xyvariables, pred.biasCorrect)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = koalaxy,
        env.grid = backg.env, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")
