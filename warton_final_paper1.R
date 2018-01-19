library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(lmtest)
library(spatial.tools)
library(VGAM)
library(mosaic)
library(faraway)
library(gstat)  #
library(ncf)    #
library(foreign)
library(nlme)   
library(MASS)
library(ROCR)
library(vcd)
library(RColorBrewer) # 
library(classInt)
library(ppmlasso)


##############
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))

# Analyse data for 2011:
all.locations= subset(mydata,yearnew >2000:2015, select=x:y)
plot(all.locations)
# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$x, all.locations$y)
select.locations = subset(all.locations, all.locations$disttoitself > 400)
selected.locations = select.locations[,1:2]

# keep a buffer distance of 2000m as required by this method.
# get  raster dimentions:441829, 541829, 6901098, 7001098  (xmin, xmax, ymin, ymax)

selected.locations2<- subset(selected.locations, x > 443829 & x < 539829)
selected.locations <- subset(selected.locations2, y > 6903098 & y < 6999098) # xy only within the study area.
coordinates(selected.locations) <- ~x+y


#####NOTE: Scaling all rasters gives slightly different predictions.
myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$")
myfullstack = scale(stack(myfullstack.a)) 
#plot rasters 
plot(myfullstack[1:10])
plot(myfullstack, c(1:10))
plot(myfullstack, c(11:20))
plot(myfullstack, c(21:30))
plot(myfullstack, c(31:40))
plot(myfullstack, c(41:45))

##check which distance variables are associated with sightings

Distance_primary <- myfullstack[[17]] 
plot(Distance_primary, main="primary roads"); plot(selected.locations, add=TRUE)
plot(Distance_motorway <- myfullstack[[14]],main="Distance to motorway"); plot(selected.locations, add=TRUE)

# this is the presence data set. extract?.
#get the full raster data set.
#projection(myfullstack) <- gsub("units=m", "units=km", projection(myfullstack))
# myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder
# myfullstack = scale(stack(myfullstack.a)) 
extent(myfullstack) <- extent(c(xmin(myfullstack), xmax(myfullstack), ymin(myfullstack), ymax(myfullstack))/1000)
#quadrature points
habitat.r<- subset(myfullstack, c(1,2,4,5,14, 17,24, 25,27,28,29,32,43,44,45)) # habitat covariates
bigquad <- as.data.frame(habitat.r, xy=TRUE, na.rm=T)

# to predict using model based control of obser bias at minimum distance.
bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 
bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)
#stt <- na.omit(stt)
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'
# stt[is.na(stt)] <- 0
xydatan <- bigquad[c(1,2)]
# stt requires xy as integers.
xydata <- as.data.frame(lapply(xydatan, as.integer)) # stt[] <- lapply(stt, as.integer)#this line edited on 04/01# make only xy integer in line with dadta shared with Mark S. 
bigquad <- cbind(xydata, bigquad[c(-1,-2)])

# get species selected data as a dataframe
selected.locations=as.data.frame(selected.locations)/1000

sp.xy = data.frame(selected.locations)
colnames(sp.xy)[1] <- 'X'; colnames(sp.xy)[2] <- 'Y'
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

ppm.form.e = ~ poly(awc,clay,elev,fpcnew, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi,habit1decimal,habit2decimal,habit3decimal,degree = 2, raw = TRUE)
scales = c( 0.5, 1, 2, 4, 8)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form.e)

ppmFit.e = ppmlasso(ppm.form.e, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 2, n.fits = 200)
#Predict and plot
pred.fit.e = predict.ppmlasso(ppmFit.e, newdata=bigquad)

predictions.fit.e <- cbind(xydata, pred.fit.e) # xydatan was chnaged to xydata.
pred.final0.e<- rasterFromXYZ(as.data.frame(predictions.fit.e )[, c("X", "Y", "pred.fit.e")])
plot(pred.final0.e, main=" koala density-warton method/ env only")
#### Env and distance both

ppm.form = ~ poly(awc,clay,elev,fpcnew, nitro,sbd,AnnualMeanTemperature,AnnualPrecipitation,tpo,twi, habit1decimal,habit2decimal,habit3decimal,degree = 2, raw = TRUE)+poly(distance_primaryandlink,distance_motorwayandlink, degree = 2, raw = TRUE)

scales = c( 0.5, 1, 2, 4, 8)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form)
#4.2 Fitting a regularisation path of point process models
#a LASSO penalty that optimises non-linear GCV
ppmFit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
#Predict and plot
pred.fit = predict.ppmlasso(ppmFit, newdata=bigquad)

predictions.fit <- cbind(xydata, pred.fit) # xydatan was chnaged to xydata.
pred.final0<- rasterFromXYZ(as.data.frame(predictions.fit )[, c("X", "Y", "pred.fit")])
plot(pred.final0, main=" koala density-WM env & dist/ bias not corrected")

# now correct for bias.
pred.fit.correct = predict.ppmlasso(ppmFit, newdata=bigquad.2)

predictions.fit.correct <- cbind(xydata, pred.fit.correct) # xydatan was chnaged to xydata.
pred.final0.correct<- rasterFromXYZ(as.data.frame(predictions.fit.correct )[, c("X", "Y", "pred.fit.correct")])
plot(pred.final0.correct, main=" koala density-warton method/ bias corrected")

### residuals:
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit)
#K-envelop
kenv = envelope(ppmFit, fun = Kinhom, nsim=39) # simulated envelop for summary function 
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

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
###### runs this code without error.

### another way to fit area interaction model. can add: criterion = "blockCV", n.blocks = 5, block.size = 10 
final.fita = ppmlasso(final.fit,sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 5,
                      block.size = 10,sp.scale = 1, n.fits = 100, family = "area.inter", r = .5)
diagnose(final.fita, which = "smooth", type = "Pearson")
pred.interaction = predict(final.fita, newdata=bigquad.2)
pred.inter.action <- cbind(xydata, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter )
diagnose.ppmlasso(final.fita)

#section 2. ############PPM model can be approximated with IWLR / DWPR########################
#Step:1 ###### Number of quadrature points for the analysis. 
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")
#projection(myfullstack) <- gsub("units=m", "units=km", projection(myfullstack))
myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder
myfullstack = scale(stack(myfullstack.a)) 
extent(myfullstack) <- extent(c(xmin(myfullstack), xmax(myfullstack), ymin(myfullstack), ymax(myfullstack))/1000)

habitat.r<- subset(myfullstack, c(1,2,4,5,14, 17,24, 25,32,43,44,45)) # habitat covariates only. removed distanc :14, 17

#now create all background data.
bigquad <- as.data.frame(habitat.r, xy=TRUE,na.rm=T) # if varying size quadrature points are needed chnage bigquad to quad.
#colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'
#quad1 <- quad[c(1,2)]
n.quad = c(50, 100, 200,500, 1000, 1500, 2000,4000, 7000) # number of quadrature poiints.
quad.inc = sample(1:dim(bigquad)[1], 1000)
assign(paste("quad.", n.quad[1], sep = ""), bigquad[quad.inc[1:n.quad[1]],])
for (i in 2:length(n.quad)){
  quad.inc = c(quad.inc, sample(setdiff(1:dim(bigquad)[1], quad.inc),
                                (n.quad[i] - n.quad[i - 1])))
  assign(paste("quad.", n.quad[i], sep = ""), bigquad[quad.inc[1:n.quad[i]],])
}

#compare the likelihood of PPMs fitted using downweighted Poisson regression:
#create species data
#IPP.pre is koala data

IPP.pre2 <- as.data.frame(selected.locations2)/1000 # rasters should be in the same extent.: 
spdata <- cbind(IPP.pre2,(extract(habitat.r, IPP.pre2)))
sp.dat <- as.data.frame(spdata)
sp.dat$Pres = 1
loglik = rep(NA, length(n.quad))

for (i in 1:length(n.quad)){
  quad = get(paste("quad.", n.quad[i], sep = ""))
  quad$Pres = 0
  all.dat = na.omit(data.frame(rbind(sp.dat, quad)))
  X.des = as.matrix(cbind(poly(all.dat$AnnualMeanTemperature, all.dat$AnnualPrecipitation, all.dat$awc
                               ,all.dat$clay, all.dat$elev, all.dat$elev
                               , all.dat$fpcnew,all.dat$nitro,all.dat$sbd,all.dat$tpo
                               ,all.dat$twi ,degree = 2, raw = TRUE)))
  
  p.wt = rep(1.e-8, dim(all.dat)[1])
  p.wt[all.dat$Pres == 0] = 10000/n.quad[i]
  z = all.dat$Pres/p.wt
  dwpr = glm(z ~ X.des, family = poisson(), weights = p.wt)
  mu = dwpr$fitted
  loglik[i] = sum(p.wt*(z*log(mu) - mu))
}
plot(n.quad, loglik, log = "x", type = "o")


###### do a IWLR| DWPR as Renner et al. ##########
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))

# Analyse data for 2011:
all.locations= subset(mydata,yearnew >2000:2015, select=x:y)
plot(all.locations)
# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$x, all.locations$y)
select.locations = subset(all.locations, all.locations$disttoitself > 400)
selected.loc = select.locations[,1:2]

# keep a buffer distance of 2000m as required by this method.
# get  raster dimentions:441829, 541829, 6901098, 7001098  (xmin, xmax, ymin, ymax)

selected.loc2<- subset(selected.loc, x > 443829 & x < 539829)
selected.loc3 <- subset(selected.loc2, y > 6903098 & y < 6999098) # xy only within the study area.
coordinates(selected.loc3) <- ~x+y

selected.loca4=as.data.frame(selected.loc3)/1000
# step 1: preapre data: 
myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder
myfullstack = stack(myfullstack.a) 
habitat.r<- subset(myfullstack, c(1,2,4,5,14, 17,24, 25,32,43,44,45)) # habitat covariates only. distance covariate: 14, 17,
extent(habitat.r) <- extent(c(xmin(habitat.r), xmax(habitat.r), ymin(habitat.r), ymax(habitat.r))/1000)
X.des <- as.data.frame(habitat.r,na.rm=T)
X.des <- as.matrix(X.des)

####need a raster without NA in stack zise.
dummy <- as.data.frame(habitat.r,na.rm=T, xy=TRUE)
r <- rasterFromXYZ(as.data.frame(dummy )[, c("x", "y", "awc")])
plot(dummy.r)


# here comes a problem . varibales length do not match.
# S first use selected.locations and extract valuves, then remove na. get xy to create a vector of 0 1
selected.loc.new <- cbind(selected.loca4,extract(habitat.r,selected.loca4))
#### remove NA at this stage.
selected.loc.new <- as.data.frame(selected.loc.new,na.rm=T)
selected.loc.new <- selected.loc.new[c(1:2)] # this data will be used to create a 0 1 vector.
selected.loc.new$Haskoala <- 1
#xy.select <- selected.loc.new[1,2]

# code to get all presence absences. Put all presences and mark them as 1 and all others zero. Then get 0/1.
# This is the response data in iwlr and dwpr. Generally need presence locations only for other methods but
#this is approximating ppm with logistic regression. So need 0/1.
# r <- raster(habitat.r, layer=4)
dfull <- as.data.frame(r, xy = TRUE,na.rm=T)
dpart = cbind(selected.loc.new,extract(r,selected.loc.new[,1:2]))
dpart <- subset(selected.loc.new, Haskoala==1)
rspec <- raster(r)
rspec[] <- 0
cells <- cellFromXY(rspec, as.matrix(dpart[, c("x", "y")]))
rspec[cells] <- dpart$Haskoala
plot(rspec)
rspec2 <- mask(rspec, r)
plot(rspec2)
text(dpart$x, dpart$y, lab = dpart$Haskoala)
Press<- as.data.frame(rspec2, xy=TRUE,na.rm=T)

Pres <- Press[,3]   
#Pres <- as.data.frame(na.omit(Pres))  # this create a vector, Press[c(3)] create a sublist.

####IWLR#####
up.wt = (10^6)^(1 - Pres) # positive valuves get by using 1.e-6
iwlr = glm(Pres ~ X.des, family = binomial(), weights = up.wt)

dd1 <- as.data.frame(X.des) # get coordinates of the design matrix for predictions. similr to warton method.
pred.iwlr = predict(iwlr, newdata=dd1)
# r <- raster(habitat.r, layer=2)
dummy <- as.data.frame(habitat.r,na.rm=T, xy=TRUE)
r <- rasterFromXYZ(as.data.frame(dummy )[, c("x", "y", "awc")])
dfull <- as.data.frame(r, xy = TRUE,na.rm=T)
xydatan <- dfull[c(1,2)]
# get coordinates only
pred.iwlr <- cbind(xydatan, pred.iwlr)
pred.iwlreg <- rasterFromXYZ(as.data.frame(pred.iwlr)[, c("x", "y", "pred.iwlr")])
plot(pred.iwlreg,asp=1)

####Get area of each cell catogory
r=pred.iwlreg
intervals <- list(c (7,8), c (8 ,9),c (9,14))
sapply(intervals, function(x) { 
  sum(r[] > x[1] & r[] <= x[2])
})

###DWPR#####
p.wt = rep(1.e-6, length(Pres))
p.wt[Pres == 0] = 10000/sum(Pres == 0)
dwpr = glm(Pres/p.wt ~ X.des, family = poisson(), weights = p.wt)

dd <- as.data.frame(X.des)
pred.dwpr = predict(dwpr, newdata=dd)
#r <- raster(habitat.r, layer=2)
dummy <- as.data.frame(habitat.r,na.rm=T, xy=TRUE)
r <- rasterFromXYZ(as.data.frame(dummy )[, c("x", "y", "awc")])
dfull <- as.data.frame(r, xy = TRUE,na.rm=T)
xydatan <- dfull[c(1,2)]
# get coordinates only
pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg,asp=1)









