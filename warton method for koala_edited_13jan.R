library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
library(quickPlot)
# Pre-standardise observer bias variables

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\test")

#full data is in raster stack/newvars
myenv <- list.files(path="wartondata", pattern="\\.tif$", full.names = TRUE) #path="wartondata", 
myenv.stack <- stack(myenv) # these rasters are square areas.
#pb.loc=SpatialPoints(koalafull)
myenv.stack <-  crop(myenv.stack, kolaxyT )
plot(newx)
myenv.stack <- scale(myenv.stack,scale=TRUE,center = TRUE)
# myenv.stack <- scale(myenv.stack) #scalinig results in glm.fit: fitted rates numerically 0 occurred 
# lgalimited.shp <- readShapePoly("LGA10new.shp")
#  myenv.stack<- mask(myenv.stack ,lgalimited.shp )
 
#crs(myenv.stack) <- NA
#stt$Y <- round(stt$Y, digits = -3) # nearest 1000   # make 442329  to 442000 but need it to be above 442300
#stt$X <- round(stt$X, digits = -3) # nearest 1000
extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)
projection(myenv.stack) <- gsub("units=m", "units=km", projection(myenv.stack))
plot(myenv.stack$AnnualMeanTemperature)
# arcGIS dataframe properties, select coordinate system, doubleclick and change to km.
#plot(myenv.stack$fpc.corrected)

stt <- as.data.frame(myenv.stack, xy=TRUE, na.rm=T)
#stt <- na.omit(stt)
colnames(stt)[1] <- "X"
colnames(stt)[2] <- "Y"
# stt[is.na(stt)] <- 0
xydatan <- stt[c(1,2)]
# stt requires xy as integers.
xydata <- as.data.frame(lapply(xydatan, as.integer)) # stt[] <- lapply(stt, as.integer)#this line edited on 04/01# make only xy integer in line with dadta shared with Mark S. 
stt <- cbind(xydata, stt[c(-1,-2)])
names(stt)
head(stt)
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
#stt <- stt[c(-3)] # without thses distance varibales.
stand.distance_tertiaryandlink=scale.default(stt$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
stt$distance_tertiaryandlink = stand.distance_tertiaryandlink
# # 
stand.dis_visitor=scale.default(stt$dis_visitor, center = TRUE, scale = TRUE) #standarise
stt$dis_visitor = stand.dis_visitor

#where is newstt?
newstt <- stt
newstt$dis_visitor = min(stand.dis_visitor)
newstt$distance_tertiaryandlink = min(stand.distance_tertiaryandlink)
#newstt$distance_tertiaryandlink = min(stt$distance_tertiaryandlink)
#koala data
kolaxyT <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom
kolaxyT <- as.data.frame(lapply(kolaxyT, as.integer))
##### Anoher datset below . use only one
koalafull <- read.csv("wartondata\\mydatasighting_cleaned.csv", header = TRUE) # not in decimals.

kolaxyT <- koalafull[c(117,127,128)]
kolaxyT <- subset(kolaxyT, yearnew == 2010)
kolaxyT <- kolaxyT[c(2,3)]
#write.csv()
kolaxyT <- kolaxyT/1000
kolaxy2 <- subset(kolaxyT, X > 443.000 & X < 537.000)
kolaxyT <- subset(kolaxy2, Y > 6902.000 & Y < 6999.000) # xy only within the study area.
kolaxyT <- as.data.frame(lapply(kolaxyT, as.integer))
# plot location of sightings on topof prediction map.
coordinates(kolaxyT) <- ~X+Y
plot(kolaxyT, add=TRUE)

######### Step 1.
#load(file="datanew.RData", .GlobalEnv)

# Model 1. no distance variables.
ppmForm1 = ~  poly(fpcnew,lot_density, AnnualPrecipitation, degree = 2) # AnnualMeanTemperature if included fitted rates numerically 0 occurred   
ppmForm1 = ~  poly(BDW_000_005
                   ,CLY_000_005
                   ,NTO_000_005
                   ,PTO_000_005
                   ,SAWC_000_005
                   ,dem
                   ,habit3percent
                   ,AnnualPrecipitation
                   ,AnnualMeanTemperature, degree = 2) # AnnualMeanTemperature if included fitted rates numerically 0 occurred   

#stt2 <- stt[c(1,2,3,4,12,18)]
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5,1, 2, 4, 8, 16)
findres(scales, coord = c("X", "Y"), sp.xy = kolaxyT, env.grid = stt, formula = ppmForm1)
## fit the model
ppmFit1 = ppmlasso(ppmForm1, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200)#criterion = "nlgc", alpha= 0.7
#print(ppmFit1)
#######
#Block corss validation as a method for choosing LASSO penalty 
#The area-interaction model with radius 2km and LASSO penalty chosen
#by 5-fold block cross-validation in the main text was fitted as follows
#criterion = "blockCV", n.blocks = 5, block.size = 32)
#ppmFit1 = ppmlasso(ppmForm1, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200,criterion = "blockCV", n.blocks = 5, block.size = 32)
### predictions
pred.biasCorrectnot = predict.ppmlasso(ppmFit1, newdata=stt)

predictions <- cbind(xydata, pred.biasCorrectnot) # xydatan was chnaged to xydata.

##### create a raster map.
pred.nct<- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrectnot")])
projection(pred.nct) <- gsub("units=m", "units=km", projection(myenv.stack))
plot(pred.nct, main=" koala density-warton method/ bias not corrected")
#writeRaster(pred.nct,filename = "test")
plot(pred.nct, zlim = c(0.02, 1.2), main=" koala density-warton method/ bias not corrected")
#### residulas model 1.
resid.plot = diagnose(ppmFit1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias NOT corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit1, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit1, which = "y", type = "Pearson", compute.sd = TRUE)
# all four plots together
diagnose.ppmlasso(ppmFit1)

# assessing GOF using 95% simulation envelop
kenv = envelope.ppmlasso(ppmFit1, fun = Kest) # fun= Kinhom
plot(kenv)
####     envelope(ppmFit1,Kinhom, nsim = 10) ## asses the gof via  95% simulation envelope of K (r)


#### model with distance covariates but not bias corrected
####
##### ppmForm1 = ~  poly(AnnualMeanTemperature,fpcnew,lot_density, AnnualPrecipitation, degree = 2) 
ppmForm10 = ~  poly(AnnualMeanTemperature,fpcnew,lot_density, AnnualPrecipitation, dis_visitor, distance_tertiaryandlink, degree = 2)
# To find the resolution (in the range from 0.5 to 16 km):
ppmForm10 = ~  poly(rain_mean_annual,temp_maximum1990_2000, habit2percent,  elevation, degree = 2)  

ppmFit10 = ppmlasso(ppmForm10, sp.xy = kolaxyT, env.grid = stt, sp.scale = 5, n.fits = 100, standardise = TRUE)


pred.distance_vars = predict(ppmFit10, newdata=newstt)

####  create a raster map
predictions.distance_var <- cbind(xydatan, pred.distance_vars)
pred.dv <- rasterFromXYZ(as.data.frame(predictions.distance_var)[, c("X", "Y", "pred.distance_vars")])
plot(pred.dv, main=" koala density-warton method/ with distance ")
plot(pred.dv, zlim = c(0, 0.35), main=" koala density-warton method/corrected")
### residulas:
kenv = envelope(ppmFit10, fun = Kinhom) # simulated envelop for summary function
resid.plot = diagnose(ppmFit10, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit10)

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
ppmFit2 = ppmlasso(ppmForm2, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 100, standardise = TRUE)


pred.biasCorrect = predict(ppmFit2, newdata=newstt)

####  create a raster map
predictions.correct <- cbind(xydatan, pred.biasCorrect)
pred.ct <- rasterFromXYZ(as.data.frame(predictions.correct)[, c("X", "Y", "pred.biasCorrect")])
plot(pred.ct, main=" koala density-warton method/ bias corrected")
plot(pred.ct, zlim = c(0, 0.35), main=" koala density-warton method/corrected")
### residulas:
kenv = envelope(ppmFit2, fun = Kinhom) # simulated envelop for summary function
resid.plot = diagnose(ppmFit2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit2)


#K-envelop
kenv2 = envelope.ppmlasso(ppmFit2, fun = Kest)
plot(kenv2)
# coordinates(kolaxyT) <- ~X+Y
# plot(kolaxyT, add=TRUE)

##### step 3. 

#species interaction at r km. Provide a avaiability grid is supplied here if some areas are inaccesible.. 
#species.int = point.interactions(dat.ppm, 5)
 load("TestPPM.RData") # load dat.ppm
ai.fit = ppmlasso(ppmForm2, data = dat.ppm, family = "area.inter", r = 2)
diagnose(ai.fit, which = "smooth", type = "Pearson")
pred.interaction = predict(ai.fit, newdata=newstt)
pred.inter.action <- cbind(xydatan, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter )

kenv3 = envelope.ppmlasso(ai.fit, fun = Kinhom) # removed fun=Kest
plot(kenv3, main= "99 simulated realisations of fitted Gibbs model")
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
# Plot(myenv.stack)
# Plot(kolaxyT)
#Plot(kolaxyT, addTo = "myenv.stack$hpop")
#https://cran.rstudio.com/web/packages/quickPlot/vignettes/iii-plotting.html

# check demo from spatsta