##### load libraries ####
library(spatial.tools)
library(VGAM)
library(mosaic)
library(spatstat)
library(faraway)
library(raster)
library(broom)
library(dismo)
library(randomForest)
library(forestFloor)
library(AUC)
library(rgl)
#library(dplyr)
library(usdm)
library(ROCR)
library(maptools)
library(mapview)
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myfullstack.c <- list.files(pattern="\\.tif$", full.names = TRUE) #select relevant folder to get detection model rasters.
myfullstack <- stack(myfullstack.c)

# get values from raster stack
covariates.m <- getValues(myfullstack)
# convert  valuves as.data.frame 
covariates.df <- as.data.frame(covariates.m)
# read shape file containing fishnet koala locations, whcih then covert to a raster as required for the analysis.


data2011.shp <- readShapePoly("hefley_fishnet_rastermatch2011.shp")
#plot(data2011.shp)

#create a matching raster
mydata.r <- raster("twi.tif")
# Vector to raster conversion:
data2011.r <- rasterize(x = data2011.shp, y = mydata.r, field = "presence")
#now get valuves from
data2011.values <- as.data.frame(data2011.r,xy=TRUE)
#bind two datasets
covariates2011 <-cbind(data2011.values, covariates.df) 

#now omit na rows of layer.
covariates2011 <-  covariates2011[!is.na(covariates2011$layer),]
covariates2011 <- na.omit((covariates2011))
#extract xy co-ordinates from data 
xy <-cbind(covariates2011$x, covariates2011$y) 
xy <- as.data.frame(xy)
#rename colomns to xy
colnames(xy)[2] <- "y"
colnames(xy)[1] <- "x"
#Fit the model with environmental variables 

env_glm <-glm(layer ~ distance_pedestrian + s1_residential_dist + distance_trunkandlink+
                distance_tertiaryandlink,family= "binomial", data=covariates2011) 
residuals <- as.data.frame(env_glm$residuals)
dim(residuals)
summary(env_glm)

rast <-raster(ncol=100, nrow = 100, ymn = 6901098, ymx = 7001098, xmn = 441829, xmx = 541829) 

res(rast) <-1 

#Extract residuals from the model called "env_glm" and map them 

xy_residuals <-cbind(xy, resid(env_glm)) 

rast[cellFromXY(rast, xy_residuals)] <-xy_residuals[,3] 

plot(rast) 
#Calculate residuals autocovariate 

#Focal operations: ngb is neighbourhood size, set to 3 by 3 cells; fun is function, 

#here the mean value within the defined neighbourhood 

focal_rac_rast <-focal(rast, w=matrix(1, ncol=3, nrow=3), fun = mean, na.rm = TRUE)

plot(focal_rac_rast) 
#Extract the values of the focal operation from “focal_rac_rast” rasterfile using the 

#co-ordinates stored in “xy” 
#coordinates(xy) <- ~x+y
plot(xy, add=TRUE)
proj4string(koalan) <- CRS("+init=epsg:28356")
focal_rac_vect <-extract(focal_rac_rast, xy) 

focal_rac_vect <-extract(focal_rac_rast,xy) 
focal_rac_vect <- as.data.frame(focal_rac_vect)
#Add as a column to the data 

covariates2011<-cbind(covariates2011, focal_rac_vect) 

# fit the RAC model using the environmental variables and the residuals autocovariate 
rac_glm <-glm(layer ~ distance_pedestrian + s1_residential_dist + distance_trunkandlink+
                distance_tertiaryandlink + focal_rac_vect, data=covariates2011,family= "binomial" )
summary(rac_glm)
summary(env_glm)
####################################################### 

#Autologistic model 
#Derive the autocovariate from the response variable (ie presence/absence of Snouter) 
#Set up blank rasterfile 
rast_ac <-raster(ncol=100, nrow = 100, ymn = 6901098, ymx = 7001098, xmn = 441829, xmx = 541829) 

res(rast_ac) <-1 
#fill the raster with the values of Snouter presence or absence allocated using “xy” 
rast_ac[cellFromXY(rast_ac, xy)] <-covariates2011[,3] 

#calculate the autocovariate via a focal operation 

focal_response_rast <-focal(rast_ac, w=matrix(1, ncol=3, nrow=3), fun = mean, na.rm = TRUE) 

plot(focal_response_rast) 
#####
focal_rac_vect <-extract(focal_rac_rast,xy) 
focal_rac_vect <- as.data.frame(focal_rac_vect)
#Add as a column to the data 

covariates2011<-cbind(covariates2011, focal_rac_vect) 

#####

ac_vect <-extract(focal_response_rast,xy)
ac_vect <-  as.data.frame(ac_vect)
covariates2011<-cbind(covariates2011, ac_vect) 

#fit the autologistic model 
autolog_glm <- glm(layer ~ distance_pedestrian + s1_residential_dist + distance_trunkandlink+
      distance_tertiaryandlink + ac_vect, data=covariates2011,family= "binomial" )

summary(autolog_glm)
pred <- predict(myfullstack, autolog_glm,type = "response")
myPred1 = predict(myfullstack, Detection.model, type = "response", na.rm=TRUE)
plot(myPred1)
plot(hefleydata.presence, pch=1, add=TRUE)

lambda <- 1
(lambda^9 * 2.71828^-lambda)/(factorial(9))
