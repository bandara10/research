library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files( path="wartondata", pattern="\\.tif$", full.names = TRUE) 
myenv.stack <- stack(myenv)

extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)

#changed reolution .5
r <- raster(ncol = 200, nrow = 200)
extent(r) <- extent(myenv.stack[[3]])
distance_tertiaryandlink <- resample(myenv.stack[[1]],r)
elev <- resample(myenv.stack[[2]],r)
temp <- resample(myenv.stack[[3]],r)
uhabit3 <- resample(myenv.stack[[4]],r)
plot(elev)
st <- stack(distance_tertiaryandlink,elev,temp,uhabit3 )

stt <- as.data.frame(st, xy=TRUE, na.rm=T)
colnames(stt)[1] <- "X"
colnames(stt)[2] <- "Y"
# stt[is.na(stt)] <- 0
xydatan <- stt[c(1,2)]
# test resolution on back transformatiotn
# sbd <- rasterFromXYZ(as.data.frame(stt)[, c("X", "Y", "temp")])
# plot(sbd)
scales = c(0.5, 1, 1.2, 2, 4, 8, 16,32)
quad.1 = sample.quad(env.grid =stt , scales, file = "Quad") # this is quadrature points to be use for the analysis.



#A matrix containing locations of species presences in the first two columns and the interpolated
#environmental data in the remaining columns


species.env = env.var(kolaxyT, env.grid = stt, env.scale = .5,
                      file.name = "Sp Env Data") # simialr to extract. isn`t it?

# A matrix dat.ppm with columns representing the latitude and longitude of presence locations and
# quadrature points along with the associated environmental data, as well as a column Pres indicating
# whether either point corresponds to a presence location or a quadrature point, and a column wt of
# observation weights

# determines observation weights
# and sets up the design matrix required for fitting a regularisation path
species.ppm = ppmdat(sp.xy = kolaxyT, back.xy = stt,
                     sp.scale = 1, file.name = "Sp PPM Data")
summary(species.ppm)

#########
#########Pre-standardise observer bias variables
stand.distance_tertiaryandlink=scale.default(stt$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
stt$distance_tertiaryandlink = stand.distance_tertiaryandlink

# To predict using model-based control of observer bias (at min value for distance):
newstt <- stt
newstt$distance_tertiaryandlink = min(stand.distance_tertiaryandlink)

#koala data
kolaxy <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom
kolaxy2 <- subset(kolaxy, X > 442 & X < 540)
kolaxyT <- subset(kolaxy2, Y > 6902 & Y < 7000) # xy within the area only.
#########
ppmForm = ~  poly(temp,elev, degree = 1)
ppmFit = ppmlasso(ppmForm, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):

#newEnv$A.1 = min(stand.A.1)
pred.biasCorrect = predict(ppmFit, newdata=newstt)
#newEnv$A.1 = min(stand.A.1)
#pred.biasCorrect.1 = predict(ppmFit, newdata=sss)
predictions <- cbind(xydatan, pred.biasCorrect)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = kolaxyT, env.grid = stt, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")



#####
