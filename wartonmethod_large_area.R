library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(lmtest);library(spatial.tools);library(VGAM);library(mosaic);library(faraway);library(gstat)  #
library(ncf);library(foreign);library(nlme)   ;library(MASS);library(ROCR);library(vcd)
library(RColorBrewer);library(classInt);library(ppmlasso);library(usdm) ; library(ncf); library(epicalc)
library(optiRum)
library(fields)
library(mvtnorm)
library(matrixStats)
#============
library(maxnet); library(maxent);library(dismo);library(rJava);library(maptools)
library(glmnet);library(reshape); library(jsonlite)
##############

setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")

##########Step 1: load koala data frm tow sources. BoalaBASE and Wildnet#####
## step 1 and two are common to all methods.


load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata=mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]

# Get wildnet data

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata_old.RData") # wildnet data.
mydata2 <- wildnetdata[c("X","Y","yearnew")]
# 
# mydata3=rbind(mydata, mydata2)

# Get gold coast data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//koala_gold_coast.RData")
mydata3 <- koala_gold_coast[c("X","Y","yearnew")]
myalldata=rbind(mydata, mydata2,mydata3)

#myalldata[!duplicated(myalldata$X &myalldata$Y), ]

all.locations= subset(myalldata,yearnew >1998, select=X:Y)
pp=SpatialPoints(all.locations)

# Select records based on distance

source("Lib_DistEstimatesToItself.r")## set a minimum distance between koalas
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$X, all.locations$Y)
select.locations = subset(all.locations, all.locations$disttoitself > 500)
selected.locations = select.locations[,1:2]

######Step 2: load rasters and crop to my extent of interest######

myfullstack.a <- list.files(pattern="\\.tif$",full.names=TRUE)
myfullstack.a = scale(stack(myfullstack.a))

#select my extent / myextent=drawExtent()

myextent=c(387900, 553100, 6862400, 7113600) 
myfullstack=crop(myfullstack.a, myextent, snap="near") ###crop to my area of interest
habitat.r<- subset(myfullstack, c(1,2,4,6,13,16,23, 24, 25,26,27,28,29,39,40,41)) # select my variables

##### select uncorrelated variables and save as uncorrelated.r

# myvif=vifstep(habitat.r, th = 3)
# uncorrelated.r=exclude(habitat.r, myvif)
# plot(uncorrelated.r)
# 
# habitat.r=uncorrelated.r

#select koalas from grids

selected.locations=gridSample(selected.locations, habitat.r, n = 100)# disaggregate a raster factor 4 and get one point.
loc=SpatialPoints((selected.locations))
plot(loc, add=TRUE)

#plot the stack and data for visualization.
points=SpatialPoints(select.locations)
fun <- function() {
  plot(points, add = TRUE, col = "red", pch = 1)
}

#plot rasters only##

plot(habitat.rr,  nc = 4, nr =3)

# plot raster and overlay koala locations

plot(habitat.rr,addfun = fun,  nc =4, nr =3)
#--
# ######## All climate rasters .

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = scale(stack(myfullstack.b))
plot(myfullstack.b)
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.b, myextent, snap="near")

habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,34,39,41)) # select my variables


######Step 3. Create quadrature points and species xy data. Remove NA. Keep xy colomns for preditions and rasterize. ######

big.grp <- as.data.frame(habitat.rr, xy=TRUE) # to be used later in Hefley method.

#

bigquad <- as.data.frame(habitat.rr, xy=TRUE,na.rm=TRUE)
xydata <- bigquad[c(1,2)] ### use to get predictions and rasterise.

# ####optional step to create a new varibale. fpc buffer area. 2 km buffer. Then make xy in to km.######
# ##fpc buffer variable.
# bigquadxy=bigquad[c(1,2)]
# fpcnew=myfullstack$fpcnew
# buff=extract(fpcnew,bigquadxy, fun=max,buffer=2,df=TRUE)
# #rename this vvariable as buffer
# colnames(buff)[2] <- "fpcnew_buff"
# buff=buff[c(2)]
# bigquad.t=cbind(buff,bigquad)
# bigquad=bigquad.t[c(2:19,1)]

bigquad =cbind((bigquad[,1:2]/(1000)),bigquad[c(-1,-2)]) # xy convert to km and integers. decimals take days to run..
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'  
xydatan <- as.data.frame(lapply(bigquad[c(1,2)], as.integer)) 
dd=bigquad[c(-1,-2)]
bigquad <- cbind(xydatan,dd )   # Quadrature points XY as integers.

##### Step 4:Create quadrature points with values set to minimum of distance variables 
#             use for bias corretion model
#             model based control of observer bias based on same level of distance

bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 
bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)

####### Step 5.Koala XY data.


#select locations over study area as  dataset is for the full area..

testlayer <- habitat.rr[[1]]   # get a raster

selected.locations.t=cbind(selected.locations,(extract(testlayer,selected.locations))) # locations are NA for outside study area.

# only data from study area. Remove NA rows and keep xy.

selected.locations<- as.data.frame(na.omit(selected.locations.t))

selected.locations=selected.locations[c(1,2)]


###
sp.xy = data.frame(selected.locations)

loc=SpatialPoints(sp.xy)

sp.xy=(sp.xy[,1:2]/(1000))
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

##Check kerneal smoothing maps with spatstat. my.owin=as.owin(my.ppp)
loca=SpatialPoints(sp.xy)
my.ppp=as.ppp(loca)
my.density=density(my.ppp, edge=TRUE) # sigma=1, 
plot(my.density)
plot(loca, add=TRUE)
# plot  LGA 10 map.
LGA.shp <- readShapePoly("LGA10new.shp")
plot(LGA.shp, add=TRUE)


#####step 5: Models: Enviromental only. Test using different env variables.
#                   Enviromental and distance variables
#                   Bias corrected.
#poly() lets you avoid corelation by producing orthogonal polynomials,
#Model 1:  Enviromental varibales only
#ppm.1 = ~ poly(clay,elev,fpcnew,fpcnew_buff, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,degree = 1, raw = TRUE)
ppm.1 = ~ poly(clay
               ,sbd
               ,tpo
               ,AnnualMeanTemperature
               ,awc
               ,habit1decimal
               ,habit2decimal
               ,habit3decimal
               ,fpcnew
               ,AnnualPrecipitation
               ,degree = 2, raw = TRUE)


ppm.1 = ~ poly(Annual_Mean_Temperature
               ,Annual_Precipitation
               ,Max_Temperature_of_Warmest_Month
               ,Min_Temperature_of_Coldest_Month
               ,fpcnew
               ,degree = 2, raw = TRUE)

 #####========================



# p=glmnet(X.des, Pres, family =  "binomial")


scales = c( 0.5, 1, 2, 4, 8,16,32)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.1)#find the spatial resolution best for analysis

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, 
                    criterion = "nlgcv", n.fits = 100)#Fitting a regularisation path of point process models.#a LASSO penalty that optimises non-linear GCV

ppmFit.1$beta  

ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "blockCV", n.blocks = 50,block.size = 50, sp.scale = 1, n.fits = 100)
ppmFit.1$beta

#use interaction terms; family = "area.inter", r = 1    |     interaction = Strauss(r = 6)


ppmFit.1 = ppmlasso(ppm.1, sp.xy = sp.xy, env.grid = bigquad,
                    criterion = "blockCV", n.blocks = 5,
                    block.size = 10, sp.scale = 1, n.fits = 50, family = "area.inter", r = .1)
ppmFit.1$beta

#Check residuals]

diagnose.ppmlasso(ppmFit.1)
resid.plot = diagnose(ppmFit.1, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.1, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.1, which = "y", type = "Pearson", compute.sd = TRUE)

#Check level of interaction

kenv = envelope(ppmFit.1,fun = Kinhom) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

#Predict and plot

pred.fit.e = predict.ppmlasso(ppmFit.1, newdata=bigquad) #get prdictions.
predi.fit.e <- cbind(xydata, pred.fit.e)     #predictions anad xy locations.
pred.model.1<- rasterFromXYZ(as.data.frame(predi.fit.e )[, c("x", "y", "pred.fit.e")])
plot(pred.model.1, main=" koala density-warton method/ env only")
plot(loc, add=TRUE)
(total_koala=cellStats(pred.model.1,sum))        #check total koal in the area


#Model :2 Enviromental and distance covariate which are of original scale (bigquad).

ppm.2 = ~ poly(clay
               ,sbd
               ,tpo
               ,AnnualMeanTemperature
               ,awc
               ,habit1decimal
               ,habit2decimal
               ,habit3decimal
               ,fpcnew
               ,AnnualPrecipitation,degree = 2 
               ,raw = TRUE)+poly(distance_primaryandlink
                                ,distance_motorwayandlink
                                ,degree = 2, raw = TRUE)
### raster data set 2

ppm.2 = ~ poly(Annual_Mean_Temperature
               ,Annual_Precipitation
               ,Max_Temperature_of_Warmest_Month
               ,fpcnew
               ,degree = 2, raw = TRUE)+poly(distance_primaryandlink
                                             ,distance_motorwayandlink
                                             ,degree = 2, raw = TRUE)
# interactions allowed between variables but not between climatic and distance variables.
# linear, quadratic and first order interactions.

ppmFit.2 = ppmlasso(ppm.2, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 100)
diagnose.ppmlasso(ppmFit.2)
#             Check residuals
resid.plot = diagnose(ppmFit.2, which = "smooth", type = "Pearson", main="smoothed pesrson residulas env model")
diagnose.ppmlasso(ppmFit.2, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(ppmFit.2, which = "y", type = "Pearson", compute.sd = TRUE)

#check level of interaction

kenv = envelope(ppmFit.2, fun = Kest) # #check k envelop and see what level of interaction occurs.
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

# Make predictions

pred.fit.dis = predict.ppmlasso(ppmFit.2, newdata=bigquad)
predi.fit.dis <- cbind(xydata, pred.fit.dis) 
pred.model.2<- rasterFromXYZ(as.data.frame(predi.fit.dis)[, c("x", "y", "pred.fit.dis")])
plot(pred.model.2, main=" koala density-warton method/ env_dis")
(total_koala=cellStats(pred.model.2,sum))


# Make predictions to new data. Bias correction step:now predict using quadrature (biguad.2) with minimum distances 

pred.fit.crt = predict.ppmlasso(ppmFit.2, newdata=bigquad.2)
predi.fit.crt <- cbind(xydata, pred.fit.crt) # xydatan was chnaged to xydata.
pred.model.crt  <- rasterFromXYZ(as.data.frame(predi.fit.crt )[, c("x", "y", "pred.fit.crt")])
plot(pred.model.crt , main=" koala density-warton method/ bias corrected")
(total_koala=cellStats(pred.model.crt, sum))

#####step 6.#plot all with same legend######

par(mfrow=c(1,3),oma=c(1,1,1,1))
plot(pred.model.1, zlim = c(0, 7),main="env only")
plot(pred.model.2, zlim = c(0, 7),main="env_dis")
plot(pred.model.crt, zlim = c(0, 7), main="bias corrected")
LGA.shp <- readShapePoly("LGA10new.shp")
plot(LGA.shp, add=TRUE)



##### Detail way of fitting the model:

#A regularisation path of Poisson point process models

quad.1k = sample.quad(bigquad, 1)
ppm.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, criterion = "nlgcv")
diagnose(ppm.fit)

#4.3 Block cross-validation
#block cross-validation as a method for choosing the LASSO penalty
#area interaction model with radius 2k and lasso penalty by5-fold cross validation.

final.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, 
                     criterion = "blockCV", n.blocks = 5, block.size = 10)
#Predict and plot

pred.final.fit = predict.ppmlasso(final.fit, newdata=bigquad.2)

predictions.final.fit <- cbind(xydata, pred.final.fit) # xydatan was chnaged to xydata.
pred.final<- rasterFromXYZ(as.data.frame(predictions.final.fit )[, c("X", "Y", "pred.final.fit")])
plot(pred.final, main=" Koala density-warton method bias corrected blockCV")

#diagnostic plots:residuals and lurking variable plots.

diagnose.ppmlasso(final.fit)
diagnose.ppmlasso(final.fit, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(final.fit, which = "y", type = "Pearson", compute.sd = TRUE)
resid.plot = diagnose(final.fit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas biascorrected interaction model")

########## For large maps.

opar <- par() #make a copy of current settings

#par(opar)          # restore original settings

mypar <- par(mar=c(1.5,1.5,1.5,1.5), oma=c(0.5,0.5,0.5,0.5))


######## =====================================############
#             First run the above code the get objects data to run belwo code.
#
#             IWLR |DWP
#       use same koala and raster dataset prepared above for comparison.
#
########======================================#############

# selected.locations = koala data over the study area.

# habitat.rr = my raster subset for the analysis.

acsel210 = selected.locations   # same as sp.xy in the above model
acsel210$Haskoala <- 1
myfullstack.sub=habitat.rr

selected.loc=SpatialPoints(selected.locations)

#subset raster as df

myfullstack.subdf= as.data.frame(myfullstack.sub) 


#############

# get all presence absences. Put all presences and mark them as 1 and all others zero. Then get 0/1.
# This is the response data in iwlr and dwpr. Generally need presence locations only for other methods but
#this is approximating ppm with logistic regression. So need 0/1.

r <- raster(myfullstack.sub, layer=4)   # get a raster layer

dfull <- as.data.frame(r, xy = TRUE)

dpart = cbind(acsel210,extract(r,acsel210[,1:2]))

dpart <- subset(acsel210, Haskoala==1)

#check minimum valuve of raster r
r

values(r )[values(r) > -1.414443] = 0
rspec <- r

plot(r)
cells <- cellFromXY(rspec, as.matrix(dpart[, c("X", "Y")]))

rspec[cells] <- dpart$Haskoala

plot(rspec)
text(dpart$X, dpart$Y, lab = dpart$Haskoala)

Press.g<- as.data.frame(rspec, xy=TRUE)   

#### Press.g need this for Hefley method for group variable.

#rename variable awc

colnames(Press.g)[3] <- "koala"
ptest=cbind(Press.g,myfullstack.subdf)   ### now remove na
Press= na.omit(ptest)  # same as sp.at=Press ;  Pres=koala

Press.my = Press# save this for future use in comapringlikelihood.

# koala presence locations

Pres <- Press[,3]   # this create a vector, Press[c(3)] create a sublist.

availability.gird=Press[c(1,2,3)]
availability.gird<- rasterFromXYZ(availability.gird[, c("x", "y", "koala")])
plot(availability.gird)
availa.gird=as.matrix(availability.gird)
a.grid=rasterToPolygons(availability.gird)
writeOGR(a.grid, ".", "a.grid", driver="ESRI Shapefile")


available_grdBris=readShapePoly("abc.shp")


av.grid=rasterToPolygons("available_grdBris.shp")


available_grdBris.shp

X.des=Press[,4:10]

# set same distance level to distance covariate as in lasso method.

X.des.2 <- X.des
X.des.2$distance_primaryandlink = min(X.des.2$distance_primaryandlink) 
X.des.2$distance_motorwayandlink = min(X.des.2$distance_motorwayandlink)


X.des.2=as.matrix(X.des.2)  # quadrature points with common level of distance to all locations.
X.des=as.matrix(X.des)  # quadrature points.




####      iwlr      #######==============================

# reall  above ppmlasso use interaction terms. linear, quadatic and first order interactions.         
# QUesiton : how to evaluvate this model? AUS and TSS.
# How to allow interactions among vaiables. to compare for ppmlasso model?
# uncomment to use this section. not too sure of accurasy of it.

X.des.df = as.data.frame(X.des)

# ###
# X.des.poly = polym(X.des.df$Annual_Mean_Temperature, X.des.df$Annual_Precipitation, X.des.df$fpcnew
#       ,X.des.df$Max_Temperature_of_Warmest_Month, X.des.df$Min_Temperature_of_Coldest_Month,degree=2, raw=TRUE)
# X.des.poly2 = polym(X.des.df$distance_motorwayandlink, X.des.df$distance_primaryandlink, degree=2, raw=TRUE)
# X.des = cbind(cbind(X.des.poly,X.des.poly2))

####

X.des = as.matrix(X.des)

###


up.wt = (1.e6)^(1 - Pres)  # up.wt = (10^6)^(1 - Pres) # Pres is a binary vector for presence absence.
iwlr = glm(Pres ~ X.des, family = binomial(), weights = up.wt) #Pres ~ X.des, decide how many quadrature points.

iwlr

dd1 <- as.data.frame(X.des) # get coordinates of the design matrix for predictions. similr to warton method.

pred.iwlr = predict(iwlr, newdata=dd1, response=TRUE) # dd1 =bigquad
#pred.iwlr=logit.prob(pred.iwlr)

# dd2 is the data set for bias correction by setting distance to minimum valuve of distance.
dd2=dd1
dd2$distance_primaryandlink = min(dd1$distance_primaryandlink) 
dd2$distance_motorwayandlink = min(dd1$distance_motorwayandlink)

#r <- raster(myfullstack.sub, layer=2) 

dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only # bias not corrected map

pred.iwlr <- cbind(xydatan, pred.iwlr)

pred.iwlreg <- rasterFromXYZ(as.data.frame(pred.iwlr)[, c("x", "y", "pred.iwlr")])
plot(pred.iwlreg,asp=1)
plot(selected.loc, add=TRUE)

# bias corrected map
# for bias correction new data set`s distance variables are set to a common level as in ppmlasso method.`
# X.de.2 has this common level.

# X.des.df2 = as.data.frame(X.des.2)
# 
# X.des.poly3 = polym(X.des.df2$Annual_Mean_Temperature, X.des.df2$Annual_Precipitation, X.des.df2$fpcnew
#                    ,X.des.df2$Max_Temperature_of_Warmest_Month, X.des.df2$Min_Temperature_of_Coldest_Month,degree=2, raw=TRUE)
# X.des.poly4 = polym(X.des.df$distance_motorwayandlink, X.des.df$distance_primaryandlink, degree=2, raw=TRUE)
# X.des2 = cbind(cbind(X.des.poly3,X.des.poly4))
# X.des2 = as.matrix(X.des2)

#dd2 <- as.data.frame(X.des.2)

pred.iwlr.2 = predict(iwlr, newdata=dd2) # bias corrected check distance variables in  dd2 dataset same as bigquad.2 data.

pred.iwlr.2 <- cbind(xydatan, pred.iwlr.2)
pred.iwlreg.2 <- rasterFromXYZ(as.data.frame(pred.iwlr.2)[, c("x", "y", "pred.iwlr.2")])
plot(pred.iwlreg.2,asp=1)
plot(selected.loc, add=TRUE)




###       DWPR     ######===================================
# QUesiton : how to evaluvate this model? AUC and TSS.
# 

p.wt = rep(1.e-6, length(Pres))
X.des <- X.des[,-6]
p.wt[Pres == 0] = 44415/sum(Pres == 0)
dwpr = glm(Pres/p.wt ~ X.des, family = poisson(), weights = p.wt)
#poly(X.des, degree=3, raw=TRUE)

# check coefficents 
dwpr$coefficients

dd <- as.data.frame(X.des)
pred.dwpr = predict(dwpr, newdata=dd, response=TRUE)
pred.dwpr=logit.prob(pred.dwpr)
dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg,asp=1, main = "dwpr_bias not corrected")

# plot locations over the map
plot(selected.loc, add=TRUE)

# bias correction map
# for bias correction new data set`s distance variables are set to a common level as in ppmlasso method.`
# X.de.2 has this common level.
# dd2 is the data set for bias correction by setting distance to minimum valuve of distance.
dd2=dd
dd2$distance_primaryandlink = min(dd$distance_primaryandlink) 
dd2$distance_motorwayandlink = min(dd$distance_motorwayandlink)


pred.dwpr.2 = predict(dwpr, newdata=dd2) # bias corrected check distance variables in  dd2 dataset.

pred.dwpr.2 <- cbind(xydatan, pred.dwpr.2)
pred.dwpreg.2 <- rasterFromXYZ(as.data.frame(pred.dwpr.2)[, c("x", "y", "pred.dwpr.2")])
plot(pred.dwpreg.2, main = " wepr bias corrected")
plot(selected.loc, add=TRUE)

############## STOP here ###############

# 5.3 Assessing the variability in likelihood for different numbers quadrature points.
# to generate different sie quadrature points.

quad=Press.my[c(-3)]
quad=as.data.frame((quad))
#quad <- as.data.frame(habitrasters, xy=TRUE) 
quad <- na.omit(quad)   ### question
colnames(quad)[1] <- 'X'; colnames(quad)[2] <- 'Y'

#load("Quad100m.RData") #xy and variables.

# another way is to use Pres.my to generate points to generate Quadraatures and use them seperatly for iwlr|dwpr.

n.quad = c(1000, 2000, 5000, 8310) # number of quadrature poiints.
quad.inc = sample(1:dim(quad)[1], 1000)
assign(paste("quad.", n.quad[1], sep = ""), quad[quad.inc[1:n.quad[1]],])
for (i in 2:length(n.quad)){
  quad.inc = c(quad.inc, sample(setdiff(1:dim(quad)[1], quad.inc),
                                (n.quad[i] - n.quad[i - 1])))
  assign(paste("quad.", n.quad[i], sep = ""), quad[quad.inc[1:n.quad[i]],])
}

#3# use them to iwlr and dwpr. 



#compare the likelihood of PPMs fitted using downweighted Poisson regression:
#create species data
#Press has x y koala 0 1  and covariates
# get koala 1 and select variables excluding koala

sp.dat.1=subset(Press.my, koala==1)
sp.dat=sp.dat.1[c(-3)]
colnames(sp.dat)[1] <- 'X'; colnames(sp.dat)[2] <- 'Y'
#to compare likelihoods from iwlr , dwpr and ppmlasso. recall : same as sp.at=Press ;  Pres=koala

# remember to change variable names.

sp.dat$Pres = 1
loglik = rep(NA, length(n.quad))

for (i in 1:length(n.quad)){
  quad = get(paste("quad.", n.quad[i], sep = ""))
  quad$Pres = 1
  all.dat = na.omit(data.frame(rbind(sp.dat, quad)))
  X.des = as.matrix(cbind(poly(  all.dat$AnnualMeanTemperature
                               , all.dat$Annual_Precipitation
                               , all.dat$Max_Temperature_of_Warmest_Month
                               , all.dat$Min_Temperature_of_Coldest_Month  
                               , all.dat$fpcnew
                               , all.dat$distance_motorwayandlink 
                               , all.dat$distance_primaryandlink
                               , degree = 2, raw = TRUE)))   
  p.wt = rep(1.e-8, dim(all.dat)[1])
  p.wt[all.dat$Pres == 0] = 10000/n.quad[i]
  z = all.dat$Pres/p.wt
  dwpr = glm(z ~ X.des, family = poisson(), weights = p.wt)
  mu = dwpr$fitted
  loglik[i] = sum(p.wt*(z*log(mu) - mu))
}
plot(n.quad, loglik, log = "x", type = "o")

#coudn`t get the Renners plot with all lines.


#######   Hefley METHOD                              #########
#                                                     #
#            Use same dataset  ru IWLR |DWPR firt to get data. #
#              Press.g
#               habitat.rr from climatic data folder. ###


#data required are : Koala location as 0 | over the entire study area.
#                  :   raster data for above koala 0 | locations.
# we saved data set "Press.my" has this variables.
#  renamed object 

# Press.g   # x, y, koala
# big.grp # x,y variables
#rename koala colomn to presence:


#### This method is to correct group size so model. same covariate possible.
# create a group variable and make it 1all locations. Otehrwise aggregate koalas in a grid and consider them as one unit or group..

# r is an empty raster creaated above.

group.r=rasterize(selected.locations, r) #  rcomes from line 356.

plot(group.r)

group.r[is.na(group.r[])] <- 0 

group.rr=group.r+r

plot(group.rr)

group=getValues(group.rr)

# remove xy from big.grp  data

myFDc=cbind(Press.g,group,big.grp)
myFDc=myFDc[c(-1,-2)]

myFD1=na.omit(myFDc)
colnames(myFD1)[1] <- "presence"

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 

#####Step 5: select only 1000 absences (monticarlo points as in hefleys method??)####
set.seed(123)                                   # size data comes from absence locations.
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 

##### Step 6: now take a random sample of 80 and assign detected 1 non detected 0.####

train <- sample(seq_len(nrow(ZTGLM.myFD1)), size = 300,replace=FALSE)
detected <- ZTGLM.myFD1[train, ]
notdetected <- ZTGLM.myFD1[-train,] 
#not detected assigned valuve 0
notdetected$presence <- 0

##### Step 7: Create the final data sets for the analysis####

Detection.data= rbind(detected,notdetected) 
str(Detection.data)

IPP.data=rbind(detected, ZTGLM.myFD4) #IPP.data comes from detected data plus no koalas data from myFD1.

ZTGLM.data=(detected)##ZTGLM.data# get the 80 rows selected.


IPP.ignored=glm(presence~Annual_Mean_Temperature
                + Annual_Precipitation
                + Max_Temperature_of_Warmest_Month
                + Min_Temperature_of_Coldest_Month
                + fpcnew
                ,family="binomial",weights=100000^(1-presence),data=IPP.data) # IPP.data2 added.

######poly(Annual_Mean_Temperature, degree=3, raw =TRUE)+ poly(Annual_Precipitation, degree= 3, raw= TRUE)
# plot(IPP.ignored)
# plot(fitted(IPP.ignored),residuals(IPP.ignored))

summary(IPP.ignored)
detecPress.gk=subset(Detection.data, select=c(-2))
Press.gk=as.data.frame(detecPress.gk)

# # ROC curve & GOF metrics
# myPred = prediction(predict(Detection.model, type = "response"), Press.gk$koala)
# perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
# plot(perf, colorize = T)
# myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), Press.gk$koala, 0.5)
# myPredPerfs


#############  
set.seed(125)
#Detection model: steps as in Hefley`s code`
Detection.model=glm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    ,family= "binomial"
                    , data=Detection.data)

unclass(summary(Detection.model))

#Plot models.

myPred0 = predict(habitat.rr, Detection.model, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred0, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main=" detection koala")






###################### Improve the model autocorrelation check and residual approach.#####
#  got to step 4 if want to ignore this section
######## Detection modelimproved by adding 
# 1. Correlogrm prepared, Moran`s` I
# 2. spatial autocorelation corrected using residulas.
# 3. ROC curve 
# 4. LRtest 
###               ###               ###           ### 
Detection.data$res = residuals(Detection.model)

myResCorr <- correlog(Detection.data$x, Detection.data$y, Detection.data$res, na.rm=T, increment=10000, resamp=0, latlon = F)

plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,xlab="distance", ylab="Moran's I")

abline(h=0)



# Crase approach to account for SA
#############################################################
myStack = stack(habitat.rr)
# Plot the predictors
#plot(myStack)
# Creates a mask layer
myMask = myStack$Annual_Precipitation >=-2

plot(myMask)
# Map predictions
AR1 = myMask * 0
plot(AR1)
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

# Crase approach to account for SA
#############################################################
# Build the standard GLM object
source("autoregresive2.r")
Detection.data$AR1 = myLib.AutoRegressiveMean(residuals(Detection.model), Detection.data$x, Detection.data$y, 18000)
AR1

Detection.model.2=glm(presence~ distance_primaryandlink
                      +distance_motorwayandlink
                      + fpcnew
                       +AR1 ,family= "binomial"
                      ,data=Detection.data)
summary(Detection.model.2)


#plot(Detection.model.2)

Detection.data$R1 = residuals(Detection.model.2)
Corr <- correlog(Detection.data$x, Detection.data$y, Detection.data$R1,na.rm=T, increment=18000, resamp=0, latlon = F)              
#Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

#lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram
plot(Corr) 
# Estimates DAIC for each predictor
#############################################################
library(epicalc)
myPredictors=c("distance_motorwayandlink","distance_primaryandlink","fpcnew")
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "presence ~"
  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
  myStrARN = paste(myStrN,"+ AR1")
  Detection.model.2 = glm(as.formula(myStrARN), data = Detection.data, family = "binomial")
  mylrtest = lrtest(Detection.model,Detection.model.2)
  myChi2[i] = mylrtest[[4]]
  myPValue[i] = mylrtest[[6]]
}
myLRTable = data.frame(myChi2,myPValue)
myLRTable = cbind(myPredictors,myLRTable)
myLRTable

# ROC & predictions
#############################################################
# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"


myPred2 = predict(myARStack, Detection.model.2, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred2, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main=" detection koala")

# ROC curve & GOF metrics
myPred = prediction(predict(Detection.model.2, type = "response"), Detection.data$presence)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(Detection.model.2, type = "response"), Detection.data$presence, 0.5)
myPredPerfs

########

#####################
# three  additional varibales, res, AR, R1 created for Detection model.2.  hence, they area not available in IPPdata, ZTGLM.data...
# t=extract(myPred2, ZTGLM.data$x &ZTGLM.data$y)
# To=ilogit(t)
#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(Detection.model,new=ZTGLM.data))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 

# set very small valuves to 0.1
p.det <- ifelse(p.det<0.1,0.1,p.det)

hist(p.det, breaks=70)


# p.det 1 or very low valuves create convergencce issues.
IPP.data$p.det=c(p.det,rep(1,length(ZTGLM.myFD3)))
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
ZTGLM.data$p.det=p.det

######Step 5: - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected=glm(presence~Annual_Mean_Temperature
                + Annual_Precipitation
                + Max_Temperature_of_Warmest_Month
                + Min_Temperature_of_Coldest_Month
                + fpcnew
                ,weights=(1/p.det)*10000^(1-presence),family="binomial",  data=IPP.data)


summary(IPP.corrected)
####Step 6: Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det.####

#use only the significant covariates, tpo +hpop+lot_density+sbd
ZTGLM.corrected=vglm(presence~ distance_primaryandlink
                    +distance_motorwayandlink
                    + fpcnew
                    ,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)

# step 7:  Detecion model Map predictions # 

myPred = predict(habitat.rr, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")

#        IPP ignored model Map predictions

myPred2 = predict(habitat.rr, IPP.ignored, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main=" IPP model |detection ignored- number of koalas")

#        IPP corrected model Map predictions

myPred2.1 = predict(habitat.rr, IPP.corrected, type = "response")
plot(myPred2.1, xlab = "x", ylab= "y",main=" IPP model |detection considered- number of koalas")


#	Calculating mean, standard deviation and 95%, equal-tailed confidence intervals
#	from the empirical distribution. See "Introduction to the Bootstrap" (Efron & Tibshirani 1994) 
#	for more details.
#
#     Estimates of the coefficient for the covariate "x" (0.9520575) of the inhomogeneous Poisson point 
#	process model and  estimate the intercept(0.9832485) and covariate "x" (0.531719) of the zero-truncated 
#	Poisson generalized linear model are close to the true values (1, 1, 0.5).   
#
#	Standard errors when detection bias is accounted for(0.1141398; 0.03696048; 0.0224879), and thus 95%  
#	confidence intervals, and standard errors (0.01185; 0.0099313; 0.0060487) are larger when the variability  
#	in the probability of detection is not accounted for. Note the estimated standard errors(0.01185; 0.0099313; 0.0060487)
#	are obtained from steps 5 & 6 above.	

####### step 8: Two phase algoritham#######################################################################

#### chnage below variable list based on s=identified variables in each model.

set.seed(1234)
tpnbs=function()	{
  bss=resample(1:dim(ZTGLM.data)[1])
  IPP.data.bss=rbind(IPP.data[bss,1:53],IPP.data[which(IPP.data$presence==0),1:53])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  
  Detection.data.bss=resample(Detection.data)
  
  Detection.model=glm(presence ~distance_pedestrian+s1_residential_dist+distance_trunkandlink+
                        distance_tertiaryandlink+scale(group),family="binomial",data=Detection.data.bss)
  
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(ZTGLM.myFD3))) 
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  
  IPP.model= glm(presence~ twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd, 
                 family = "binomial", weights = (1/p.det)*10000^(1-presence), data = IPP.data.bss)
  
  ZTGLM.corrected=vglm(group~twi+tpo+temp+aspect+elev+habit2pc+hpop+lot_density+sbd,
                       weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
 
   c(coef(IPP.model),coef(ZTGLM.corrected))
}

#### step 9: boostrapping ####
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps)*tpnbs() # Ravi changed the moasic::do and tpnbs()
bootstrap.sample <- data.frame(bootstrap.sample)
#bootstrap.sample=as.data.frame (bootstrap.sample) #Ravi created a matrix from the dataframe
######## each step gives mean of coefficients for intercept or covariates
save.image(file="hefleydatapreparation.v3.3.RData")
#### step 10:  Get means and sd ####
colMeans(bootstrap.sample)[1]
sd(bootstrap.sample)[1]
qdata(c=(.025),bootstrap.sample[,1]) # ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[2] 
sd(bootstrap.sample)[2] 
qdata(c=(.025),bootstrap.sample[,2])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[3] 
sd(bootstrap.sample)[3] 
qdata(c=(.025),bootstrap.sample[,3])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[4]
sd(bootstrap.sample)[4]
qdata(c=(.025),bootstrap.sample[,4])#ravi:removed after .975 after .025, add c=

colMeans(bootstrap.sample)[5]
sd(bootstrap.sample)[5]
qdata(c=(.025),bootstrap.sample[,5])#ravi:removed after .975 after .025, add c=




######   Drorazio method  ##############
#This method works with survey data. Small survey datset added.
          #use same data 
####  ####  ####    ###  ####  ####  ####    ###
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_syn\\warton_data")
source("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\Supplimentary materials of papers\\Vira Koshinika paper\\Data\\functions2.r")

## recall rasters being used. "habitat.rr"
## we read rasters directly from the folder.

x.files=c( "AnnualMeanTemperature.tif"
            ,"AnnualPrecipitation.tif"
           
            ,"clay.tif"
            
            ,"habit1decimal.tif"
            )
d <- stack(x.files)

s.occupancy <- scale(d,scale=TRUE,center = TRUE)
plot(s.occupancy)


w.files=c("elev.tif") 

dd <- stack(w.files)

s.detection <- scale(dd, scale=TRUE,center = TRUE)
plot(s.detection)

## we use same datasets and criteria. import data and select in the same manner we did for previous analysis.

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))

#combine with wildnet data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata_old.RData")
mydata2 <- wildnetdata[c("X","Y", "yearnew")] 
table(mydata2$yearnew)
names(mydata2) <- tolower(names(mydata2))

# Get gold coast data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//koala_gold_coast.RData")
mydata3 <- koala_gold_coast[c("X","Y","yearnew")]
names(mydata3) <- tolower(names(mydata3))
all.loc=rbind(mydata, mydata2,mydata3)

(b=SpatialPoints(all.loc))

#crop base on extent of koala locations.
s.occupancy=crop(s.occupancy,b)
s.detection= crop(s.detection,b)

plot(s.occupancy)

# Analyse data 
all.locations= subset(all.loc,yearnew ==2012, select=x:y) #2005 originally
table(all.loc$yearnew)
plot(SpatialPoints(all.locations))
# set a minimum distance between koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$x, all.locations$y)
select.locations = subset(all.locations, all.locations$disttoitself > 400)
selected.locations = select.locations[,1:2]


####################crop raster stack at this stage.
(b=SpatialPoints(all.locations))

##################

pb=selected.locations

pb.loc=SpatialPoints(pb) # over 6629 beyond study area. 
# get locations over rasters only.
pb.occupancy=extract(s.occupancy,pb.loc)  # get covaiates
pb.detection=extract(s.detection,pb.loc)  #get covariates
# only data from study area.
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,] #3663#  distance covariates extract from rasters. only one here
pb.occupancy=pb.occupancy[is.complete.pb,]#3663# env covariates extract from  raster
# upto here basically covariate extraction is done for presence data for occupancy/abunace and detection., .

print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection)) #s.detection is raster stack.

# # bring repeated survey data2010-2013
# so.occupancy <- read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/New folder/dorazio_survey.csv")
# so.occupancy <- so.occupancy[c(1,2)]
# so.occupancy <- extract(s.occupancy, so.occupancy)
# # so.occupancy=read.csv("so_occupancy.csv") # occupancy covariate valuves for survey data
# so.detection=read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/New folder/so.detection.csv") # detection covariate valuves wind sky day 1 and 2 for survey data.
# y.so=read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/New folder/y.so.csv") # survey data or site occupancy, day 1 day 2 presence absence.

########Option 2: ============== 
#bring repeated survey data 1996-1997
survey.data <- read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/Koala Survey Data 1996-2009Cleaned/rpt.survey96_97.csv")
survey.data2 <- read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/Koala Survey Data 2010-2015Cleaned/repeated_surveyDEHP.csv")
survey.data2 <- survey.data2[c(9,10,1:8)]
survey.data <- rbind(survey.data,survey.data2)
# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
survey.data$disttoitself = Lib_DistEstimatesToItself(survey.data$EastingM, survey.data$NorthingsM)
survey.data = subset(survey.data, survey.data$disttoitself > 100)

so.occupancy <- survey.data[c(1,2)]
so.occupancy <- extract(s.occupancy, so.occupancy)

so.detection <- survey.data[c(5:10)]
y.so <- survey.data[c(3,4)] 

#######Option 1: =================================

#  bring repeated survey data2010-2013
survey.data <- read.csv("C:/Users/uqrdissa/ownCloud/koala_survey_DEHP_joerg_96_2015/Koala Survey Data 2010-2015Cleaned/repeated_surveyDEHP.csv")
xysurvey <- survey.data[c(9,10)]
survey.xy <- SpatialPoints(xysurvey)
plot(survey.xy)

so.occupancy <- extract(s.occupancy,xysurvey)
so.detection <- survey.data[c(3:8)]

y.so <- survey.data[c(1,2)] 

###### next step================================
print("removing NA")
is.complete=complete.cases(so.occupancy)&complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*0.04

# print('removing raster files')
# #removing rasters to free the memory
# for (i in c(x.files,w.files)) {
#   do.call(remove,list(i))
# }
# remove(is.complete.pb,po,yb.so,po.loc)


print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection)) #s.detection is raster stack.

# remove all NA values
# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

print("specifying area ")
#area in squared km ??????????????????????? -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell


s.area=area.back*nrow(X.back) #study area


# # adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # v1 and covariate valuves
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection) # v1, probability of detection.

##Preparing Presence-Absence data

# #add a column of ones to the PA covariat
# #y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
so.occupancy <- as.matrix(so.occupancy) # added by me
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
#X.so$v1 <- X.so$`rep(1, nrow(as.matrix(so.occupancy)))
# # # X.so <- X.so[c(-1)]
# # # X.so <- X.so[c(13, 1:12)]
# # #X.so <- as.matrix(X.so)
W.so = array(dim=c(nrow(as.matrix(pb.detection)), J.so, 1))
W.so[,,1] = 1
 #W.so[,,2] = pb.detection# if it changes
# W.so[,,3] = pb.detection2# if it changes



# # Checking whether occupancy and detection rasters have the same resolution -----
# if(sum(res(s.occupancy)!=res(s.detection)))
#   stop("Occupancy and detection raster layers have different resolution")
# 
# if(ncell(s.occupancy)!=ncell(s.detection))
#   stop("Occupancy and detection have different number of cells")
# 
# 
# # Plotting covariates that drive occupancy and detection in PO
# ppi = 300
# png('occupancy-covariates.png', width=9*ppi, height=3*ppi, res=ppi)
# plot(s.occupancy)
# dev.off()
# 
# png('PO-detection -covariates.png', width=9*ppi, height=3*ppi, res=ppi)
# plot(s.detection)
# dev.off()
# adding column of ones - po locations # Is this pb.detection?????
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # pb.occupancy is all presence locations. 
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)

# 2. Analysing the data ========================================================

#Analyzing Presence-Only data
(pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)) 

# W.pb is pb.detection
# X.pb and X.back same covariates. 
# W.back is distance covariate for detection. 
# Analyzing presence-absence data
(so.fit=so.model(X.so,W.so,y.so)) 

# Analyzing presence-only data AND presence-absence data
(poANDso.fit=pbso.integrated(X.pb, W.pb,X.back, W.back,X.so,W.so,y.so))

#estimates are obtained but not se.
# If any of the eigenvalues are zero or negative, there is obviously a problem, 
# perhaps stemming from under identification or boundary conditions
# One thing you might try is to refit the model from its solution, to see 
# if the gradients get smaller, and the NaN's clear up.
# > pb.fit
# $coefs
# Parameter name       Value Standard error
# 1                     beta0 -2.14667766     0.07994286
# 2     AnnualMeanTemperature  0.79053507     0.10671768
# 3                      clay -0.29408512     0.03372959
# 4                      elev -0.77522035     0.13323224
# 5       AnnualPrecipitation  0.35858355     0.03337058
# 6                  habit0pc  0.65510645     0.09293961
# 7                  habit1pc  0.15593408     0.05454184
# 8                  habit2pc  0.11365113     0.05051433
# 9                  habit3pc  0.32729882     0.07020533
# 10                   alpha0  0.14012092     0.18512155
# 11                 dis_city -3.27397922     0.27526465
# 12 distance_tertiaryandlink -0.07361067     0.12179921
# 
# $convergence
# [1] 0
# 
# $value
# [1] 2752.702
# 
# $value
# [1] 2783.46

################ presence ackground IPP model #####
coef <- pb.fit$coefs
coeff <- coef[c(1,2)]

#occupancy rasters

f <- s.occupancy

# get predictions for each raster

f2 <- subset(f,c(1))*coeff$Value[[2]]
f3 <- subset(f,c(2))*coeff$Value[[3]]
f4 <- subset(f,c(3))*coeff$Value[[4]]
f5 <- subset(f,c(4))*coeff$Value[[5]]
# f6 <- subset(f,c(5))*coeff$Value[[6]]
# f7 <- subset(f,c(6))*coeff$Value[[7]]
# f8 <- subset(f,c(1))*coeff$Value[[8]]
# f9 <- subset(f,c(1))*coeff$Value[[9]]
# f10 <- subset(f,c(1))*coeff$Value[[10]]
# f11 <- subset(f,c(1))*coeff$Value[[11]]
# f12 <- subset(f,c(1))*coeff$Value[[12]]

fn <- exp(f2+f3+f4+f5+coeff$Value[[1]])

plot(fn)

plot(fn, main= "intensity 2010 data")


######detection rasters ##### detectoon model

ff <- s.detection

# get predictions for each raster

ff13 <- subset(ff,c(1))*coeff$Value[[6]] 
ff14 <- subset(ff,c(2))*coeff$Value[[9]]
fn2 <-  (ff13)+coeff$Value[[7]]   

fnn <- logit.prob(fn2)

plot(fnn, main =" probability of detection")

########final intensity from two models######

fn3 <- fn*fnn

plot(fn3, main = "bias corrected intensity")


# plot intensity baed on ecological varibales and detection based on observer bias vribales.

par(mfrow=c(2,2))

plot(fn3); plot(fnn) ; plot(fn)

################# Map standarard erros  #######
coef <- pb.fit$coefs
se <- coef[c(1,3)]
se$`Standard error`

#occupancy rasters
s <- s.occupancy

# get predictions for each raster

s2 <- subset(s,c(1))*se$`Standard error`[[2]]
s3 <- subset(s,c(2))*se$`Standard error`[[3]]
s4 <- subset(s,c(3))*se$`Standard error`[[4]]
s5 <- subset(s,c(4))*se$`Standard error`[[5]]
s6 <- subset(s,c(1))*coeff$Value[[6]]
s7 <- subset(s,c(1))*coeff$Value[[7]]
# f8 <- subset(f,c(1))*coeff$Value[[8]]
# f9 <- subset(f,c(1))*coeff$Value[[9]]
# f10 <- subset(f,c(1))*coeff$Value[[10]]
# f11 <- subset(f,c(1))*coeff$Value[[11]]
# f12 <- subset(f,c(1))*coeff$Value[[12]]

sn <- exp(s2+s3+s4+s5+se$`Standard error`[[1]])

plot(sn)

plot(sn, main= "Standard error")

####### integrated model#######

coef <- poANDso.fit$coefs

#occupancy rasters

fi <- s.occupancy

# get predictions for each raster

fi2 <- subset(fi,c(1))*coeff$Value[[2]]
fi3 <- subset(fi,c(2))*coeff$Value[[3]]
fi4 <- subset(fi,c(3))*coeff$Value[[4]]
fi5 <- subset(fi,c(4))*coeff$Value[[5]]
# fi6 <- subset(f,c(5))*coeff$Value[[6]]
# f7 <- subset(f,c(6))*coeff$Value[[7]]
# f8 <- subset(f,c(1))*coeff$Value[[8]]
# f9 <- subset(f,c(1))*coeff$Value[[9]]
# f10 <- subset(f,c(1))*coeff$Value[[10]]
# f11 <- subset(f,c(1))*coeff$Value[[11]]
# f12 <- subset(f,c(1))*coeff$Value[[12]]

fin <- exp(fi2+fi3+fi4+fi5+coeff$Value[[1]])

plot(fin)

plot(fin, main= "intensity 2010 data")


######detection rasters ##### detectoon model

fd <- s.detection

# get predictions for each raster

fd13 <- subset(fd,c(1))*coeff$Value[[7]] 
# fd14 <- subset(fd,c(2))*coeff$Value[[9]]
fdn2 <-  (fd13)+coeff$Value[[6]]   

fnnd <- logit.prob(fdn2)

plot(fnnd, main =" probability of detection")

########final intensity from two models######

fin3 <- fin*fnnd*-0.36960085

plot(fin3, main = "bias corrected intensity")



####### STOP here      ########

f <- s.occupancy
f1 <- subset(f,c(1))*0.6220131 #0.79053507
f2 <- subset(f,c(2))*- 0.3067281#-0.29408512
f3 <- subset(f,c(3))*-0.7443429 #-0.77522035
f4 <- subset(f,c(4))*0.3372941 #0.35858355
f5 <- subset(f,c(5))*0.6024214 # 0.65510645
f6 <- subset(f,c(6))*0.1370172 # 0.15593408
f7 <- subset(f,c(7))*0.1432968 #0.11365113
f8 <- subset(f,c(8))*0.2757422  #0.32729882 
fn <- exp(f1+f2+f3+f4+f5+f6+f7+f8)-1.8802929 #-2.14667766  
plot(fn)
opar <- par() #make a copy of current settings
#par(opar)          # restore original settings
mypar <- par(mar=c(0.5,0.5,0.5,0.5), oma=c(0.5,0.5,0.5,0.5))

ff <- s.detection
f9 <- subset(ff,c(1))*-2.3293703# -0.6404511
f10<- subset(ff,c(2))*-0.6387529# -0.6387529 
fn2 <- (f9+f10)-0.6404511 
library(optiRum)
fnn <- logit.prob(fn2)
plot(fnn)
fn2 <- fn*fnn

plot(fn2)
 
# save model outputs
sink("myoutputs.txt")
print(summary(IPP.model))
print( summary(Detection.model))
print(summary(IPP.corrected))
sink()

# print dpi 300 maps 

jpeg("dorazio.jpeg",width=8,height=4,units="in",res=300) # 8:4 for three mpas. one map 8:10
# png(width=3,height=3,units="in",res=600)

par(mfrow=c(1,3),oma=c(1,1,1,1))
plot(pred.model.1, zlim = c(0, 7),main="env only")
plot(pred.model.2, zlim = c(0, 7),main="env_dis")
plot(pred.model.crt, zlim = c(0, 7), main="bias corrected")

graphics.off()



##### REsults 1.================== parameters, > 2012, distance 50, surveydata96_97.
# $coefs
# Parameter name       Value Standard error
# 1                 beta0  4.01549817     5.25900958
# 2 AnnualMeanTemperature -0.97016671     0.03647652
# 3   AnnualPrecipitation  0.76058893     0.02468097
# 4                  clay  0.47807865     0.03564276
# 5         habit1decimal  0.04782838     0.02754691
# 6                alpha0 -9.27773802     5.25134292
# 7                  elev -2.77694980     0.08155905
# 
# $convergence
# [1] 0
# 
# $value
# [1] 4556.717
# 
# > (so.fit=so.model(X.so,W.so,y.so))
# $coefs
# Parameter name      Value Standard error
# 1                 beta0  1.3098383      1.0900629
# 2 AnnualMeanTemperature -2.7358278      0.6825200
# 3   AnnualPrecipitation  0.1339197      0.8247535
# 4                  clay  0.3944841      0.3751803
# 5         habit1decimal -0.4095235      0.1254893
# 6             alpha0.so  1.0808204      0.1659979
# 
# $convergence
# [1] 0
# 
# $value
# [1] 395.2938
# 
# > (poANDso.fit=pbso.integrated(X.pb, W.pb,X.back, W.back,X.so,W.so,y.so))
# $coefs
# Parameter name       Value Standard error
# 1                 beta0 -0.53605508     0.08592910
# 2 AnnualMeanTemperature -0.98825464     0.03655363
# 3   AnnualPrecipitation  0.78258720     0.02431756
# 4                  clay  0.46210101     0.03584251
# 5         habit1decimal  0.01544906     0.02752850
# 6                alpha0 -4.87370662     0.12206365
# 7                  elev -3.12974096     0.10044113
# 8             alpha0.so  1.01287968     0.17139997
# 
# $convergence
# [1] 0
# 
# $value
# [1] 4973.2


##### Maxent model===============
## we use same datasets and criteria. import data and select in the same manner we did for previous analysis.
## Maxnet is an R package that for fitting Maxent species distribution models using the R package glmnet.

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))

#combine with wildnet data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//wildnetdata_old.RData")
mydata2 <- wildnetdata[c("X","Y", "yearnew")] 
table(mydata2$yearnew)
names(mydata2) <- tolower(names(mydata2))

# Get gold coast data
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//koala_gold_coast.RData")
mydata3 <- koala_gold_coast[c("X","Y","yearnew")]
names(mydata3) <- tolower(names(mydata3))
all.loc=rbind(mydata, mydata2,mydata3)

(b=SpatialPoints(all.loc))

# Analyse data 
all.locations= subset(all.loc,yearnew >2012, select=x:y) #2005 originally
table(all.loc$yearnew)
plot(SpatialPoints(all.locations))
# set a minimum distance between koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$x, all.locations$y)
select.locations = subset(all.locations, all.locations$disttoitself > 400)
selected.loc = select.locations[,1:2]
# loc = SpatialPoints(selected.loc)
# plot(loc, add = TRUE)
# remove duplicates
selected.loc.dups=duplicated(selected.loc[, c("x", "y")])
selected.loc <-selected.loc[!selected.loc.dups, ]


# Get raster data

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = stack(myfullstack.b)
myextent<-extent(min(selected.loc$x)-3000,max(selected.loc$x)+3000,min(selected.loc$y)-3000,max(selected.loc$y)+3000)

habitat.rr=crop(myfullstack.b, myextent, snap="near")

habitat.rr<- subset(habitat.rr, c(1,2,15,18,26,34, 41)) # select my variables
plot(habitat.rr)

habitat.rr=scale(habitat.rr)
plot(habitat.rr)

fold <- kfold(selected.loc, k=5) # add an index that makes five random groups of observations
selected.loctest <- selected.loc[fold == 1, ] # hold out one fifth as test data
selected.loctrain <- selected.loc[fold != 1, ] # the other four fifths are training data

TrainEnv <- extract(habitat.rr,selected.loctrain)

# Background points are for model validation.No partitioning. 
set.seed(0)
backgr <- randomPoints(habitat.rr, 1000)
absvals <- extract(habitat.rr, backgr)

#we make a presence/absence collumn where 1=presence for our occurance points, and 0=absence for our background points
presabs <- c(rep(1, nrow(TrainEnv)), rep(0, nrow(absvals)))

sdmdata <- data.frame(cbind(presabs, rbind(TrainEnv, absvals)))

sdmdata = na.omit(sdmdata)
#we make a subset of that dataset, without the presence and absence values
data <- sdmdata[,-1]
data = na.omit(data)

#we run the maxnet function to fit the SDM
#maxnet fits using glmnet
loc = sdmdata$presabs
koala.model<-maxnet(loc, data,maxnet.formula(loc, data, classes="lq")) 

# Default choices for feature classes and regularization.
#  The second call overrides the default feature classes, using only linear and
# quadratic features for the continuous predictor variables.
####
#fits Maxent models using the same feature classes
#(linear, quadratic, hinge, etc.) and regularization
 #mod <- maxnet(loc, data, maxnet.formula(loc, data, classes="default"))
summary(koala.model)
##three types of response plots
plot(koala.model, type = "exponential")
plot(koala.model, type = "cloglog")
plot(koala.model, type = "logistic")

plot(koala.model, "fpcnew")
plot(koala.model, "Annual_Mean_Temperature")

gg=predict(habitat.rr,koala.model, response=TRUE)
plot(gg)
# ps <- unscale(gg,habitat.rr)
# pp=unscale(gg, center = NULL, scale = NULL)
# plot(ps)

gg=predict(habitat.rr,mod, response=TRUE)
ggg<- logit.prob(gg)
plot(ggg)




###### For Maxent  software coordinates in latitude and longitude.
###### Rasters also in same crs.

myInPtDF = SpatialPointsDataFrame(selected.locations[c("X","Y")],selected.locations)
# Set the projection of the input spatial dataframe
myInProj = CRS("+init=epsg:28356")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+proj=longlat +ellps=WGS84") # UTM zone 45
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
plot(myOutPtDF)
selected.locations$LATITUDE = myOutPtDF@coords[,1]
selected.locations$LONGITUDE = myOutPtDF@coords[,2]

select.maxent <- selected.locations[c(3,4)]
select.maxent $SPECIES <- "koala"
select.maxent <- select.maxent[c(3,2,1)]
write.csv(select.maxent, "select.maxent.csv",row.names=FALSE)

write.csv(selected.locations, "selected.locations.csv")

myfullstack.b <- list.files(path="warton_data_allclimate",pattern="\\.tif$",full.names=TRUE)
myfullstack.b = scale(stack(myfullstack.b))
plot(myfullstack.b)
myextent=c(387900, 553100, 6862400, 7113600) 
habitat.rr=crop(myfullstack.b, myextent, snap="near")

habitat.rr<- subset(myfullstack.b, c(1,2,15,18,26,34,39,41)) # select my variables


habitat.rr
myfullstack.bb <- list.files(pattern="\\.asc$",full.names=TRUE)
myfullstack.bb = scale(stack(myfullstack.bb))


projection(habitat.rr) <- CRS("+init=epsg:28356")
# Reproject
wgs<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
habitat.rrmaxent <- projectRaster(myfullstack.bb, crs = wgs)
writeRaster(habitat.rr, filename = names(habitat.rr), bylayer = TRUE, format= "GTiff")

