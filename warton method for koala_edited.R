library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files(path="wartondata", pattern="\\.tif$", full.names = TRUE) #path="wartondata", 
myenv.stack <- stack(myenv)
#crs(myenv.stack) <- NA
extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)
plot(myenv.stack)

stt <- as.data.frame(myenv.stack, xy=TRUE, na.rm=T)
#stt <- na.omit(stt)
colnames(stt)[1] <- "X"
colnames(stt)[2] <- "Y"
# stt[is.na(stt)] <- 0
xydatan <- stt[c(1,2)]
stt[] <- lapply(stt, as.integer)
names(stt)

#abc = sample.quad(env.grid =stt , sp.scale= 1, file = "Quad") # this is quadrature points to be use for the analysis.
# 
# 
# 
# #A matrix containing locations of species presences in the first two columns and the interpolated
# #environmental data in the remaining columns
# # 
# species.env = env.var(kolaxyT, env.grid = stt, env.scale = .5,
#                       file.name = "Sp Env Data") # simialr to extract. isn`t it?
# A matrix dat.ppm with columns representing the latitude and longitude of presence locations and
# quadrature points along with the associated environmental data, as well as a column Pres indicating
# whether either point corresponds to a presence location or a quadrature point, and a column wt of
# observation weights.
# determines observation weights and sets up the design matrix required for fitting a regularisation path
# species.ppm = ppmdat(sp.xy = kolaxyT, back.xy = stt,
#                    sp.scale = 1, file.name = "Sp PPM Data")
# summary(species.ppm)

#########
#########Pre-standardise observer bias variables
stt <- stt[c(-3)] # without thses distance varibales.
stand.distance_tertiaryandlink=scale.default(stt$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
stt$distance_tertiaryandlink = stand.distance_tertiaryandlink

stand.dis_visitor=scale.default(stt$dis_visitor, center = TRUE, scale = TRUE) #standarise
stt$dis_visitor = stand.dis_visitor


#koala data
kolaxyT <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom
# coordinates(kolaxyT) <- ~X+Y
# plot(kolaxyT, add=TRUE)
# kolaxy2 <- subset(kolaxy, X > 442 & X < 540)
# kolaxyT <- subset(kolaxy2, Y > 6902 & Y < 7000) # xy within the area only.



######### Step 1.
load(file="datanew.RData", .GlobalEnv)

# Model 1. no distance variables.
ppmForm1 = ~  poly(temp, elev, hpop,lot_density, degree = 2)  

# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, coord = c("X", "Y"), sp.xy = kolaxyT, env.grid = stt, formula = ppmForm1)
## fit the model
ppmFit1 = ppmlasso(ppmForm1, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200)#criterion = "nlgc", alpha= 0.7
### predictions
pred.biasCorrectnot = predict(ppmFit1, newdata=stt)
predictions <- cbind(xydatan, pred.biasCorrectnot)

##### create a raster map.
pred.nct<- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrectnot")])
plot(pred.nct, main=" koala density-warton method/ bias not corrected")
#### residulas model 1.
kenv = envelope(ppmFit1, fun = Kinhom)
resid.plot = diagnose(ppmFit1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias NOT corrected model")


#####    Step 2
# Model 2. bias corrected model with distance variables.
## To predict using model-based control of observer bias (at min value for distance):
# newstt <- stt
# newstt$dis_visitor = min(stand.dis_visitor) # based on minimum distance
# newstt$distance_tertiaryandlink = min(stand.distance_tertiaryandlink)


ppmForm2 = ~  poly(temp,elev,hpop,lot_density, degree = 2) + poly(dis_visitor, distance_tertiaryandlink, degree = 2)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, coord = c("X", "Y"), sp.xy = kolaxyT, env.grid = stt, formula = ppmForm2)

### for the model
ppmFit2 = ppmlasso(ppmForm2, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 100)
pred.biasCorrect = predict(ppmFit2, newdata=newstt)

####  create a raster map
predictions.correct <- cbind(xydatan, pred.biasCorrect)
pred.ct <- rasterFromXYZ(as.data.frame(predictions.correct)[, c("X", "Y", "pred.biasCorrect")])
plot(pred.ct, main=" koala density-warton method/ bias corrected")
### residulas:
kenv = envelope(ppmFit2, fun = Kinhom) # simulated envelop for summary function
resid.plot = diagnose(ppmFit2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
# coordinates(kolaxyT) <- ~X+Y
# plot(kolaxyT, add=TRUE)

##### step 3. 

#species interaction at r km. Provide a avaiability grid is supplied here if some areas are inaccesible.. 
#species.int = point.interactions(dat.ppm, 5)
 load("TestPPM.RData") # load dat.ppm
ai.fit = ppmlasso(ppmForm1, data = dat.ppm, family = "area.inter", r = 1.2)
diagnose(ai.fit, which = "smooth", type = "Pearson")
pred.interaction = predict(ai.fit, newdata=newstt)
pred.inter.action <- cbind(xydatan, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter )
#writeRaster(pred.ct.inter , "pred.ct.inter .tif")

###DO not RUN This is for ArcGIS mapping for better resolution.######
bbr <- raster ("bbr.tif")  # bbr  is a dummy 1k resolution map which isnot in meters.
#use this map and fill with our raster data in to it as seperate new varibales. Becuase this chan be mapped and checked in arcgis.
bbr <- as.data.frame(bbr, xy=TRUE)
bbr <- bbr[c(-3)]
elev <- as.data.frame(myenv.stack[[4]])
ele <- cbind(bbr,elev)
temp <- as.data.frame(myenv.stack[[5]])
uhabit3 <- as.data.frame(myenv.stack[[6]])
stt <- cbind(ele,temp,uhabit3)
stt <- na.omit(stt)
stt[] <- lapply(stt, as.integer)


# ###formal argument "scales" matched by multiple actual arguments
# Becaus the scales arguement is pre-set in the command, we cannot change it directly.
# However, we can alter it with update().
