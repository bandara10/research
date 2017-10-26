#below raster is correct resolution and change raster stack resolution and extenet and write them.
#writeRaster(stack(D), names(D), bylayer=TRUE, format="ascii",overwrite=TRUE)
unzip("vector\\AU_Qld_study_area.zip")
studyareafull.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(studyareafull.shp) <- CRS("+init=epsg:28356") 
windows();plot(studyareafull.shp)
studyareafull <- raster(studyareafull.shp)
res(studyareafull) <- 1000
windows();plot(studyareafull)
----------------------------------------------------------------------------------------------
#changed directory to save files
   #delete this section
  #synchronise raster and write raster
  dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356")
#crop the full study area based on a shape file.
#habitat distance maps are syncronised.
proj4string(habit3.r) <- CRS("+init=epsg:28356")
r2 <- crop(habit3.r, extent(dstudyarea.shp))
r3 <-  mask(r2, dstudyarea.shp)
windows():plot(r3)  
test<-spatial_sync_raster(r3,studyareafull, method = "ngb")
windows(); plot(test)
writeRaster(test,"habit3.TIF", overwrite=TRUE)
----------------------------------------------------------------------------------------  
   
  #synchronise raster and write raster
test<- raster("habit_vector_raster\\habit3.tif")
test<-spatial_sync_raster(test, studyareafull, method = "ngb")
writeRaster(test,"habit_3.TIF", overwrite=TRUE)
windows(); plot(test)
----------------------------------------------------------------------------------------------
test3<- raster("test\\rain_96_20051.tif")
test3<-spatial_sync_raster(test3, studyareafull, method = "ngb")
writeRaster(test3,"AU_Qld_rain_96_20051-MGA56.TIF", overwrite=TRUE)

writeRaster(hpop2,"raster_stack\\hpop2.TIF", overwrite=TRUE)
stack.n <- stack( "aspect2.tif", "awc2.tif","fpc2.tif", "slope2.tif","topo2.tif" ,"hpop2.tif","tpo2.tif", "rain2.tif", "nitro2.tif", "elev2.tif", "temp2.tif", "sbd2.tif")
windows(); hist(stack.n)
windows(); levelplot(stack.n)
windows(); plot(stack.n)
windows(); density(stack.n)
# calculating linne densiy

lines.shp <- readShapeLines("AU_Qld_study_area_mask_roads\\AU_Qld_study_area_mask_lines-MGA56.shp")
windows(); plot(lines.shp)
names(lines.shp@data)
levels(lines.shp@data$highway)
#select colomn secondary and select rivers.
secondary.lines <- lines.shp[lines.shp$highway == "secondary",]
windows(); plot(secondary.lines,lwd=1)

# Convert SpatialLines to psp object using maptools library
pspSl <- as.psp(secondary.lines)
# Pixellate with resolution of 0.5, i.e. 2x2 pixels
px <- pixellate(pspSl, eps=1000)


# This can be converted to raster as desired
rLength <- raster(px)
# read shapefile to get extend to clip
unzip("vector\\AU_Qld_study_area.zip")
studyshp <- readShapePoly("AU_Qld_extended_study_area_mask-MGA56.shp")
proj4string(studyshp) <- CRS("+init=epsg:28356") 
plot(studyshp)

#masking rater
secondary.lines <- crop(x = secondary.lines, y = studyshp, snap = "near")
#
writeRaster(rLength,"rLength3.TIF", overwrite=TRUE)
plot(rLength)
windows(); plot(rLength)

rLength.n<-spatial_sync_raster(rLength, studyshp, method = "ngb")



suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif")
suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif")
suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif")
suitabl_0.r<- raster("habit_vector_raster\\suitable_0.tif")
stack.habit<-stack(suitable_3.r,suitable_2.r,suitable_1.r, suitabl_0.r)
windows(); plot(stack.habit)
library(rasterVis)
#bridle.r<- raster("roads_raster\\bridle.tif")
#cycle.r<- raster("roads_raster\\cycle.tif")
#footway.r<- raster("roads_raster\\footway.tif")
#motorway.r<- raster("roads_raster\\motorway.tif")
#motorwaylink.r<- raster("roads_raster\\motorwaylink.tif")
#path.r<- raster("roads_raster\\path.tif")
#pedestrian.r<- raster("roads_raster\\pedestrian.tif")
#primary.r<- raster("roads_raster\\primary.tif")
#secondry.r<- raster("roads_raster\\secondry.tif")
#secondry_link.r<- raster("roads_raster\\secondry_link.tif")
#tertiary.r<- raster("roads_raster\\tertiary.tif")
#tertiary_link.r<- raster("roads_raster\\tertiary_liink.tif")
#track.r<- raster("roads_raster\\track.tif")
trunck.r<- raster("roads_raster\\trunck.tif")
#unclassified.r<- raster("roads_raster\\unclassified.tif")
#residential.r<- raster("roads_raster\\residential.tif")

lot_density.r<- raster("raster_synchronised\\raster_syn\\lot_density.tif")
windows();levelplot(lot_density.r)
plot3D(lot_density.r)
--------------------------------------------------------------------------------------------------------------------------------------------
library(dichromat)
myTheme <- rasterTheme(region = dichromat(terrain.colors(15)))
windows(); levelplot(lot_density.r, par.setting=myTheme)
meanAug <- cellStats(lot_density.r, mean)
windows();levelplot(lot_density.r - meanAug, par.settings = RdBuTheme)
--------------------------------------------------------------------------------------------------------------------------------------------
stack.road<-stack(cycle.r,footway.r,motorway.r, motorwaylink.r,path.r, pedestrian.r,primary.r, 
               secondry.r,secondry_link.r, tertiary.r, tertiary_liink.r, track.r,trunck.r,unclassified.r, residential.r)

windows(); levelplot(stack.road)
#windows();vectorplot(hpop.r, par.settings=RdBuTheme())
#windows();streamplot(hpop.r)
-----------------------------------------------------
  #start modelling 05/03/2017
dat.ppm000 <- ppm(dat.ppp, trend = ~ 1, covariates = list(cycle = cycle.im,  unclassified = unclassified.im, tertiary=tertiary.im))
summary(dat.ppm000)$coefs.SE.CI
#####----------------------------------
#need glm and pscl to run the ZIP model.Poisson regression coefficients for each of the variables along with standard errors, z-scores, and p-values for the coefficients. A second block follows that corresponds to the inflation model. This includes 
#logit coefficients for predicting excess zeros along with their standard errors, z-scores, and p-values.
library(pscl)
summary(model1<-zeroinfl( koalan ~ hpop  + water+ tpo + rain + nitro + elev +temp +sbd, data = extract.e))
summary(model2<-zeroinfl( koalan ~  water+ tpo  + nitro + elev +temp +sbd, data = extract.e))
summary(model3<-zeroinfl( koalan ~  tpo  + nitro + elev +temp +sbd, data = extract.e))
summary(model4<-zeroinfl( koalan ~  tpo  + nitro + elev +sbd, data = extract.e))
#####ompare with the current model to a null model without predictors using chi-squared test on the difference of log likelihoods.
mnull <- update(model3, . ~ 1)
pchisq(2 * (logLik(model3) - logLik(mnull)), df = 4, lower.tail = FALSE)
#This yields a high significant p-value; thus, our overall model is statistically significant.
#the model output above does not indicate in any way if our zero-inflated model is an improvement over a standard Poisson regression. We can determine this 
#by running the corresponding standard Poisson model and then performing a Vuong test of the two models.
summary(p1 <- glm(koalan ~  tpo  + nitro + elev +temp +sbd, family = poisson,  data = extract.e))
myStr = "koalan ~"
for (i in 1:length(stack.r)){myStr = paste(myStr, "+", stack.r[i])}  
print(myStr)

myglm = glm(koalan ~  tpo  + nitro + elev +temp +sbd, family = poisson,  data = extract.e))
###Vuong test compares the zero-inflated model with an ordinary Poisson regression model. 
vuong(p1, model3)
# we can see that our test statistic is significant, 
#indicating that the zero-inflated model is superior to the standard Poisson model.
summary(model3)
extract.e$res = residuals(p1)
#----------------------------------------------------------------------------------
id <- mydata$yearnew == 2011
table(id)
mydata <- mydata[id,]
dim(mydata)
names(mydata)
id <- mydata$lat > -30 # remove record with wrong coordinate
mydata <- mydata[id,]
mydata.n <- mydata[c(3,4)]
dim(mydata.n)
mydata.n <- unique( mydata.n[,1:2 ] )
dim(mydata.n)
mydata.n = as.data.frame(mydata.n)
mydata.n<-na.omit(mydata.n)
# start distance based selection method
source("Lib_DistEstimatesToItself.r")
myInPtDF = SpatialPointsDataFrame(mydata.n[c("x","y")], mydata.n)

windows();plot(myInPtDF, axes = T)
myInPtDF$disttoitself = Lib_DistEstimatesToItself(myInPtDF$x, myInPtDF$y)
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 500)
windows();plot(myInPtDF2, axes = T)
myFD1 = myInPtDF2[,1:2]
myFD1<-as.data.frame(myFD1)
dim(myFD1)

---------------------------------------------------------------------------------- 
source("owin2sp_source.r")

dat.w <- as(as(dstudyarea.shp, "SpatialPolygons"), "owin")
dat.spw <- owin2SP(dat.w)

# Set projection of the owin as GDA94 / SA Lambert:
#proj4string(dat.spw)
#proj4string(dat.spw) <- CRS("+init=epsg:28356") 
dat.spw <- raster::intersect(x = dstudyarea.shp, y = dat.spw)

# Convert the sp object into a window:
dat.w <- as.owin(dat.spw)
  
   #Make a ppp object:
  dat.ppp <- ppp(x = myFD1[,1], y = myFD1[,2], window = dat.w)

windows(); plot(dat.ppp, axes = TRUE)
 # windows(); plot(studyarea.shp, axes = TRUE)
#points(x = acsel[,1], y = acsel[,2])

# Make a ppp object:

------------------------------------------------------------------------------------------------
  primary.im <- as.im(primary.r)
windows(); plot(primary.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
primary.rho <- rhohat(object = dat.ppp, covariate = primary.im)

windows(); plot(primary.rho,xlab = "primary (units)", main = "")
------------------------------------------------------------------------------------------------
  secondry.r<- raster("roads_raster\\secondry.tif")
secondry.im <- as.im(secondry.r)
windows(); plot(secondry.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
secondry.rho <- rhohat(object = dat.ppp, covariate = secondry.im)

windows(); plot(secondry.rho,xlab = "secondry (units)", main = "")
------------------------------------------------------------------------------------------------
  tertiary.r<- raster("roads_raster\\tertiary.tif")
tertiary.im <- as.im(tertiary.r)
windows(); plot(tertiary.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
tertiary.rho <- rhohat(object = dat.ppp, covariate = tertiary.im)

windows(); plot(tertiary.rho,xlab = "tertiary (units)", main = "")
------------------------------------------------------------------------------------------------
unclassified.r<- raster("roads_raster\\unclassified.tif")
unclassified.im <- as.im(unclassified.r)
windows(); plot(unclassified.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
unclassified.rho <- rhohat(object = dat.ppp, covariate = unclassified.im)

windows(); plot(unclassified.rho,xlab = "unclassified (units)", main = "")
-----------------------------------------------------------------------------------------------
  residential.im <- as.im(residential.r)
windows(); plot(residential.im, axes = TRUE)

# What is the nature of the association between koala sight locations and total phosphorous?
residential.rho <- rhohat(object = dat.ppp, covariate = residential.im)

windows(); plot(residential.rho,xlab = "residential (units)", main = "")
-----------------------------------------------------------------------------------------------
path.r<- raster("roads_raster\\path.tif")
path.im <- as.im(path.r)
windows(); plot(path.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
path.rho <- rhohat(object = dat.ppp, covariate = path.im)

windows(); plot(path.rho,xlab = "path (units)", main = "")
-----------------------------------------------------------------------------------------------
  pedestrian.r<- raster("roads_raster\\pedestrian.tif")
pedestrian.im <- as.im(pedestrian.r)
windows(); plot(pedestrian.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
pedestrian.rho <- rhohat(object = dat.ppp, covariate = pedestrian.im)

windows(); plot(pedestrian.rho,xlab = "pedestrian (units)", main = "")
------------------------------------------------------------------------------------------------
  track.r<- raster("roads_raster\\track.tif")
track.im <- as.im(track.r)
windows(); plot(track.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
track.rho <- rhohat(object = dat.ppp, covariate = track.im)

windows(); plot(track.rho,xlab = "track (units)", main = "")
------------------------------------------------------------------------------------------------
  tertiary_link.r<- raster("roads_raster\\tertiary_liink.tif")
tertiary_link.im <- as.im(tertiary_link.r)
windows(); plot(tertiary_link.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
tertiary_link.rho <- rhohat(object = dat.ppp, covariate = tertiary_link.im)

windows(); plot( tertiary_link.rho,xlab = "tertiary_link(units)", main = "")
------------------------------------------------------------------------------------------------
secondry_link.r<- raster("roads_raster\\secondry_link.tif")
secondry_link.im <- as.im(secondry_link.r)
windows(); plot(secondry_link.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
secondry_link.rho <- rhohat(object = dat.ppp, covariate = secondry_link.im)

windows(); plot(secondry_link.rho,xlab = "secondry_link(units)", main = "")
------------------------------------------------------------------------------------------------
  footway.r<- raster("roads_raster\\footway.tif")
footway.im <- as.im(footway.r)
windows(); plot(footway.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
footway.rho <- rhohat(object = dat.ppp, covariate = footway.im)

windows(); plot(footway.rho,xlab = "footway(units)", main = "")
------------------------------------------------------------------------------------------------
  motorway.r<- raster("roads_raster\\motorway.tif")
motorway.im <- as.im(motorway.r)
windows(); plot(motorway.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
motorway.rho <- rhohat(object = dat.ppp, covariate = motorway.im)

windows(); plot(motorway.rho,xlab = "motorway(units)", main = "")
------------------------------------------------------------------------------------------------
motorwaylink.r<- raster("roads_raster\\motorwaylink.tif")
motorwaylink.im <- as.im(motorwaylink.r)
windows(); plot(motorwaylink.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
motorwaylink.rho <- rhohat(object = dat.ppp, covariate = motorwaylink.im)

windows(); plot(motorwaylink.rho,xlab = "motorwaylink(units)", main = "")
-------------------------------------------------------------------------------------------------
  cycle.r<- raster("roads_raster\\cycle.tif")
cycle.im <- as.im(cycle.r)
windows(); plot(cycle.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
cycle.rho <- rhohat(object = dat.ppp, covariate = cycle.im)

windows(); plot(cycle.rho,xlab = "cycle(units)", main = "")
------------------------------------------------------------------------------------------------  
  trunck.r<- raster("roads_raster\\trunck.tif")
trunck.im <- as.im(trunck.r)
windows(); plot(trunck.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
trunck.rho <- rhohat(object = dat.ppp, covariate = trunck.im)

windows(); plot(trunck.rho,xlab = "trunck(units)", main = "")
------------------------------------------------------------------------------------------------
  bridle.r<- raster("roads_raster\\bridle.tif")
bridle.im <- as.im(bridle.r)
windows(); plot(bridle.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
bridle.rho <- rhohat(object = dat.ppp, covariate = bridle.im)

windows(); plot(bridle.rho,xlab = "bridle(units)", main = "")
-------------------------------------------------------------------------------------------------
 roads_motor.r<- raster("roads_raster\\roads_motor.tif")
roads_motor.im <- as.im(roads_motor.r)
windows(); plot(roads_motor.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
roads_motor.rho <- rhohat(object = dat.ppp, covariate = roads_motor.im)

windows(); plot(roads_motor.rho,xlab = "roads_motor(units)", main = "")

-------------------------------------------------------------------------------------------------
roads_other.r<- raster("roads_raster\\roads_other.tif")

roads_other.im <- as.im(roads_other.r)
windows(); plot(roads_other.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
roads_other.rho <- rhohat(object = dat.ppp, covariate = roads_other.im)

windows(); plot(roads_other.rho,xlab = "roads_other(units)", main = "")
-------------------------------------------------------------------------------------------------

 suitable_3.r<- raster("habit_vector_raster\\suitable_3.tif")
suitable_3.im <- as.im(suitable_3.r)
windows(); plot(suitable_3.im, axes = TRUE)
writeRaster(test,"suitable_3.TIF", overwrite=TRUE)
# Crop the raster to the shape file spatial extent:
suitable_3.r <- crop(x = suitable_3.r, y = r, snap = "near")
windows();plot(suitable_3.r)

# What is the nature of the association between koala sight locations and total phosphorous?
suitable_3.rho <- rhohat(object = dat.ppp, covariate = suitable_3.im)

windows(); plot(suitable_3.rho,xlab = "suitable_3(units)", main = "")
------------------------------------------------------------------------------------------------
suitable_2.r<- raster("habit_vector_raster\\suitable_2.tif")
suitable_2.im <- as.im(suitable_2.r)
windows(); plot(suitable_2.im, axes = TRUE)
# Crop the raster to the shape file spatial extent:
suitable_2.r <- crop(x = suitable_2.r, y = r, snap = "near")
writeRaster(test,"suitable_2.TIF", overwrite=TRUE)
windows();plot(suitable_2.r)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_2.rho <- rhohat(object = dat.ppp, covariate = suitable_2.im)

windows(); plot(suitable_2.rho,xlab = "suitable_2(units)", main = "")
------------------------------------------------------------------------------------------------
suitable_1.r<- raster("habit_vector_raster\\suitable_1.tif")
suitable_1.im <- as.im(suitable_1.r)
windows(); plot(suitable_1.im, axes = TRUE)
writeRaster(test,"suitable_1.TIF", overwrite=TRUE)
# Crop the raster to the shape file spatial extent:
suitable_1.r <- crop(x = suitable_1.r, y = r, snap = "near")
windows();plot(suitable_1.r)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_1.rho <- rhohat(object = dat.ppp, covariate = suitable_1.im)

windows(); plot(suitable_1.rho,xlab = "suitable_1(units)", main = "")
------------------------------------------------------------------------------------------------

suitable_0.r<- raster("habit_vector_raster\\suitable_0.tif")
suitable_0.im <- as.im(suitable_0.r)
windows(); plot(suitable_0.im, axes = TRUE)
writeRaster(test,"suitable_0.TIF", overwrite=TRUE)
# Crop the raster to the shape file spatial extent:
suitable_0.r <- crop(x = suitable_0.r, y = r, snap = "near")
windows();plot(suitable_0.r)
# What is the nature of the association between koala sight locations and total phosphorous?
suitable_0.rho <- rhohat(object = dat.ppp, covariate = suitable_0.im)

windows(); plot(suitable_0.rho,xlab = "suitable_0(units)", main = "")
------------------------------------------------------------------------------------------------
  stream.r<-raster("rivers_raster\\stream.tif")
stream.im <- as.im(stream.r)
windows(); plot(stream.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
stream.rho <- rhohat(object = dat.ppp, covariate = stream.im)

windows(); plot(stream.rho,xlab = "stream(units)", main = "")
----------------------------------------------------------------------------------------------
rivers
rivers.r<-raster("rivers_raster\\rivers.tif")
rivers.im <- as.im(rivers.r)
windows(); plot(rivers.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
rivers.rho <- rhohat(object = dat.ppp, covariate = rivers.im)

windows(); plot(rivers.rho,xlab = "rivers(units)", main = "")
-----------------------------------------------------------------------------------------------
  riverbank.r<-raster("rivers_raster\\riverbank.tif")
riverbank.im <- as.im(riverbank.r)
windows(); plot(riverbank.im, axes = TRUE)
# What is the nature of the association between koala sight locations and total phosphorous?
riverbank.rho <- rhohat(object = dat.ppp, covariate = riverbank.im)
windows(); plot(riverbank.rho,xlab = "riverbank(units)", main = "")
savePlot(filename = "service_dist",type = c( "png"),device = dev.cur()) 
----------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------
  # rhohat analyses - distance to roads:
pedestrian<-raster("raster\\distance_pedestrian.TIF")
pedestrian.im <- as.im(pedestrian)
windows(); plot(pedestrian.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = pedestrian.im)
windows(); plot(roaden.rho, xlab = "pedestrian(distance) ", main = "")
savePlot(filename = "pedestrian_dis",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------
track<-raster("raster\\distance_track.TIF")
track.im <- as.im(track)
windows(); plot(track.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = track.im)
windows(); plot(roaden.rho, xlab = "track (distance m)", main = "")
savePlot(filename = "track_dis",type = c( "png"),device = dev.cur())
--------------------------------------------------------------------------------------------------------------
  trunkandlink<-raster("raster\\distance_trunkandlink.TIF")
trunkandlink.im <- as.im(trunkandlink)
windows(); plot(trunkandlink.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = trunkandlink.im)
windows(); plot(roaden.rho, xlab = "trunkandlink (distance m)", main = "")
savePlot(filename = "trunkandlink_dist",type = c( "png"),device = dev.cur()) 
-----------------------------------------------------------------------------------------------------------------
  motorwayandlink<-raster("raster\\distance_motorwayandlink.TIF")
motorwayandlink.im <- as.im(motorwayandlink)
windows(); plot(motorwayandlink.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = motorwayandlink.im)
windows(); plot(roaden.rho, xlab = "motorwayandlink (distance m)", main = "")
savePlot(filename = "motorwayandlink_dist",type = c( "png"),device = dev.cur()) 
-----------------------------------------------------------------------------------------------------------------
primaryandlink<-raster("raster\\distance_primaryandlink.TIF")
primaryandlink.im <- as.im(primaryandlink)
windows(); plot(primaryandlink.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = primaryandlink.im)
windows(); plot(roaden.rho, xlab = "primaryandlink(distance m)", main = "")
savePlot(filename = "primaryandlink_dist",type = c( "png"),device = dev.cur()) 
-------------------------------------------------------------------------------------------------------------------
secondary<-raster("raster\\distance_secondaryandlink.TIF")
secondary.im <- as.im(secondary)
windows(); plot(secondary.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = secondary.im)
windows(); plot(roaden.rho, xlab = "secondary(distance m)", main = "")
savePlot(filename = "secondary_dist",type = c( "png"),device = dev.cur()) 
-----------------------------------------------------------------------------------------------------------------
  tertiary<-raster("raster\\distance_tertiaryandlink.TIF")
tertiary.im <- as.im(tertiary)
windows(); plot(tertiary.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = tertiary.im)
windows(); plot(roaden.rho, xlab = "tertiary(distance m)", main = "")
savePlot(filename = "tertiary_dist",type = c( "png"),device = dev.cur()) 
-----------------------------------------------------------------------------------------------------------------
residential<-raster("raster\\distance_residentil.TIF")
residential.im <- as.im(residential)
windows(); plot(residential.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = residential.im)
windows(); plot(roaden.rho, xlab = "residential(distance m)", main = "")
savePlot(filename = "residential_dist",type = c( "png"),device = dev.cur())  
-------------------------------------------------------------------------------------------------------------------
footway<-raster("raster\\distance_footway.TIF")
footway.im <- as.im(footway)
windows(); plot(footway.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = footway.im)
windows(); plot(roaden.rho, xlab = "footway(distance m)", main = "")
savePlot(filename = "footway_dist",type = c( "png"),device = dev.cur()) 
-------------------------------------------------------------------------------------------------------------------
  steps<-raster("raster\\distance_steps.TIF")
steps.im <- as.im(steps)
windows(); plot(steps.im, axes = TRUE)
roaden.rho <- rhohat(object = dat.ppp, covariate = steps.im)
windows(); plot(roaden.rho, xlab = "steps(distance m)", main = "")
savePlot(filename = "steps_dist",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
 path<-raster("raster\\distance_path.TIF")
path.im <- as.im(path)
windows(); plot(path.im, axes = TRUE)
path.rho <- rhohat(object = dat.ppp, covariate = path.im)
windows(); plot(path.rho, xlab = "path(distance m)", main = "")
savePlot(filename = "path_dist",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
  cycleway<-raster("raster\\distance_cycleway.TIF")
cycleway.im <- as.im(cycleway)
windows(); plot(cycleway.im, axes = TRUE)
cycleway.rho <- rhohat(object = dat.ppp, covariate = cycleway.im)
windows(); plot(cycleway.rho, xlab = "cycleway(distance m)", main = "")
savePlot(filename = "cycleway_dist",type = c( "png"),device = dev.cur())
-------------------------------------------------------------------------------------------------------------------
 service<-raster("raster\\distance_service.TIF")
service.im <- as.im(service)
windows(); plot(service.im, axes = TRUE)
service.rho <- rhohat(object = dat.ppp, covariate = service.im)
windows(); plot(service.rho, xlab = "service(distance m)", main = "")
savePlot(filename = "service_dist",type = c( "png"),device = dev.cur()) 
-------------------------------------------------------------------------------------------------------------------
  bridleway<-raster("raster\\distance_bridleway.TIF")
bridleway.im <- as.im(bridleway)
windows(); plot(bridleway.im, axes = TRUE)
bridleway.rho <- rhohat(object = dat.ppp, covariate = bridleway.im)
windows(); plot(bridleway.rho, xlab = "bridleway(distance m)", main = "")
savePlot(filename = "bridleway_dist",type = c( "png"),device = dev.cur()) 
------------------------------------------------------------------------------------------------------------------
unclassified<-raster("raster\\distance_unclassified.TIF")
unclassified.im <- as.im(unclassified)
windows(); plot(unclassified.im, axes = TRUE)
unclassified.rho <- rhohat(object = dat.ppp, covariate = unclassified.im)
windows(); plot(unclassified.rho, xlab = "unclassified(distance m)", main = "")
savePlot(filename = "unclassified_dist",type = c( "png"),device = dev.cur()) 
save.image("rastersrack.RData")
-----------------------------------------------------------------------------------------------------------------
