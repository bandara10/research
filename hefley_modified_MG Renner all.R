library(maptools); library(sp); library(raster); library(rgdal); library(spatstat); library(ggplot2)
library(spatialkernel); library(splancs); library(RColorBrewer); library(dismo); library(spatial.tools)
library(raster) #
library(rgdal)  
library(gstat)  #
library(ncf)    #
library(spatstat)
library(foreign)
library(maptools)
library(nlme)   
library(MASS)
library(ROCR)
library(vcd)
library(epicalc)
library(RColorBrewer) # 
library(classInt)

#load functions code.
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S")

# Load all of the objects in the R workspace up to the point of starting the ppm modelling (starts on line 596):
load(file = "C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\AU_Qld_koala_habitat_v07.RData")

mydata  <- read.csv("mydatasighting.csv", header = TRUE)
mydata <- mydata[c("LAT","LNG","X","Y", "CallerName", "HomePhone", "WorkPhone", "Sighting","SexId", "FateId", "LGA" , "monthnew", "yearnew", "Suburb", "PostCode")]
names(mydata) <- tolower(names(mydata))
# myND <- subset(mydata, yearnew <= 2011, select = c("x","y"))
# mydata$homephone <- ifelse(is.na(mydata$homephone), mydata$workphone, mydata$homephone)
# Analyse data for 2011:
id <- mydata$yearnew == 2011
table(id)
# mydata <- mydata[id,]
# dim(mydata)
# names(mydata)

# Bring in MGA56 study area boundary map:
unzip("vector\\AU_Qld_study_area.zip")
studyarea.shp <- readShapePoly("AU_Qld_study_area-MGA56.shp")
proj4string(studyarea.shp) <- CRS("+init=epsg:28356") 
plot(studyarea.shp)
id<- mydata$lat > -30 # remove record with wrong coordinate
mydata <- mydata[id,]
head(mydata)
# mydata in lat-lon:
coords <- SpatialPoints(mydata[, c("lng", "lat")])
mydata.ll <- SpatialPointsDataFrame(coords, mydata)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))
plot(mydata.mga)
# Add the MGA coordinates to mydata:
mydata$x <- coordinates(mydata.mga)[,1]
mydata$y <- coordinates(mydata.mga)[,2]

windows(); plot(mydata.mga, cex = 0.3, col = 'blue', pch = 15)

# Create a 3 km by 3 km grid:
mydata.r <- raster(studyarea.shp)
res(mydata.r) <- 1000
mydata.r <- extend(mydata.r, extent(mydata.r) + 3000)#extenet should be similar to the selected area not the full dataset
#writeRaster(mydata.r,"testing.asc")
#windows(); plot(mydata.r)
# create the disance matrix
#####Go to 109 for grid based sampling. START Distaance based selection using the function
# Function to estimate the distance to the nearest point of the same vector
source("Lib_DistEstimatesToItself.r")
# start distance based selection method
mydata.n <- mydata.mga[c(3,4)]
mydata.n = as.data.frame(mydata.n)
myInPtDF = SpatialPointsDataFrame(mydata.n[c("x","y")], mydata.n)
windows();plot(myInPtDF, axes = T)
myInPtDF$disttoitself = Lib_DistEstimatesToItself(myInPtDF$x, myInPtDF$y)
myInPtDF2 = subset(myInPtDF, myInPtDF$disttoitself > 500)
windows();plot(myInPtDF2, axes = T)
acsel = myInPtDF[,1:2]
acsel$all = 1
myFD2 = myInPtDF2[,1:2]
myFD2$all = 0
myFD = rbind(acsel,myFD2)
windows();plot(myFD, axes = T,col=c("red","blue")[myFD$all+1])
myFD3<-as.data.frame(myFD) #;acsel <- myFD3 [c(1,2)]
#acsel <- as.data.frame(myFD);acsel<-acsel[c(1,2,3)]; acsel <- subset(acsel, all==1,select=x:y)
#########END of distance based selection#####
# Sample points from within each grid cell:
acsel <- gridSample(mydata.mga, mydata.r, n = 1)
#mydata.p <- rasterToPolygons(mydata.r)
acsel<-as.data.frame(acsel,row.names = NULL)
colnames(acsel) <- c("x","y")# rename colomns to x, y to match negative data colomns
head(acsel)
----
  windows(); plot(mydata.p, border = 'gray')
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(mydata.mga)
# Selected points in red:
points(acsel, cex = 0.2, col = 'red', pch = 15)# go to 103.

dstudyarea.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(dstudyarea.shp) <- CRS("+init=epsg:28356") 

windows(); plot(dstudyarea.shp, axes = TRUE)
plot(studyarea.shp, add = TRUE)
points(acsel, cex = 0.2, col = 'red', pch = 15)
points(x = mydata$x, y = mydata$y,cex = 0.2, col = 'blue', pch = 15)
-----------
  # kdate as a date:
  # dat$kdate <- as.Date(dat$kdate, format = "%d/%m/%Y")
  # dat$kmonth <- format(dat$kdate, format = "%m")
  # 2_Creates a random set of negatives
  #############################################################edit this
  windows();plot(mydata.p, border = 'gray')
# Sets the number of points to generate
n = 1000
#Generates random coordinates
writeRaster(mydata.r,"mydata.r.TIF", overwrite=TRUE)
myMask = raster("mydata.r.TIF") >= 0
x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
# Keep only points according to different conditions
myTDF = as.data.frame(cbind(x,y)) # negative data point
head(myTDF)
acsel2 <- gridSample(myTDF, myMask, n = 1)
points(x = acsel2$x, y = acsel2$y,cex = 0.2, col = 'blue', pch = 15)
#plot(acsel2)
# Merge the two datasets (positives and negatives
#############################################################
names(acsel) = c("x","y")
acsel$Haskoala = 1
acsel2$Haskoala = 0
acsel21 = rbind(acsel,acsel2)
# Remove potential duplicates falling in same pixel
#############################################################
source("Lib_DistEstimatesToItself.r")
acsel21$DistToItself = Lib_DistEstimatesToItself(acsel21$x,acsel21$y)
acsel21 = subset(acsel21, DistToItself > 100)
acsel21 = acsel21[,-4]
# displays the output
plot(mydata.p, border = 'gray')
points(acsel21$x,acsel21$y, pch = 15, col=c("blue","red")[acsel21$Haskoala+1])

# ---------------------------------------------------------------------------------------------------------------------------------
#myPredList<- list.files(pattern="\\.tif$") 
# myfullstack = stack(myfullstack.a)
#myfullstack <- scale(myfullstack)

#myPredictors = c("twi", "tpo", "aspect", "elev", "habit2pc", "hpop", "lot_density", "sbd")
myPredictors = c("s1_residential_dist" , "distance_trunkandlink", "distance_tertiaryandlink")

# add extensions
myPredictorsASC = paste(myPredictors,".TIF", sep = "")
myPredList = as.list(paste(myPredictors,".TIF", sep = ""))  
mystack = stack(myPredList)
#-----
acsel22 = cbind(acsel21,extract(mystack,acsel21[,1:2]))
str(acsel22)
acsel22 = na.omit(acsel22)


  # Build the GLM object
myStr = "Haskoala ~"
for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
print(myStr)
print("##########################################")
# Runs the GLM  
myGLM = glm(as.formula(myStr), data = acsel22, family = "binomial")

  
summary(myGLM)
# Stores the residuals
acsel22$res = residuals(myGLM)
# Create a spatial correlogram of the residuals. coords=acsel22[c(1,2)] ; res=acsel22[c(7)]
myResCorr <- correlog(acsel22$x, acsel22$y, acsel22$res, na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(myResCorr$p<0.05,1,0)  
#plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
#xlab="distance", ylab="Moran's I")
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")
abline(h=0)

# Crase approach to account for SA
#############################################################
myStack = stack(myPredList)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = raster("distance_trunkandlink.TIF") >= 0
plot(myMask)
# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"


# # Stores the residuals
# acsel22$res = residuals(myGLM)
# # Create a spatial correlogram of the residuals
# myResCorr <- correlog(acsel22$x, acsel22$y, acsel22$res,na.rm=T, increment=1000, resamp=0, latlon = F)
# #mySigVec = ifelse(myResCorr$p<0.05,1,0)  
# #plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
# #xlab="distance", ylab="Moran's I")
# plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
#      xlab="distance", ylab="Moran's I")
# abline(h=0)

# Crase approach to account for SA
#############################################################
# Build the standard GLM object

acsel22$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), acsel22$x, acsel22$y, 15000)
# Changes the formular string
myStrAR = paste(myStr,"+ AR1")
# Re-run the model  
myGLMR1 = glm(as.formula(myStrAR), data = acsel22, family = "binomial")
summary(myGLMR1)
plot(myGLMR1)  



# Estimate correlogram of new residuals

acsel22$R1 = residuals(myGLMR1)
Corr <- correlog(acsel22$x, acsel22$y, acsel22$R1,na.rm=T, increment=10000, resamp=0, latlon = F)              
#Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram
plot(Corr)              

# Estimates DAIC for each predictor
#############################################################
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "Haskoala ~"
  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
  myStrARN = paste(myStrN,"+ AR1")
  myGLMR2 = glm(as.formula(myStrARN), data = acsel22, family = "binomial")
  mylrtest = lrtest(myGLMR1,myGLMR2)
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

myPred = predict(myARStack, myGLMR1, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main="koala")

# ROC curve & GOF metrics
myPred = prediction(predict(myGLM, type = "response"), acsel22$Haskoala)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), acsel22$Haskoala, 0.5)
myPredPerfs





# save.image("koalasightingprobability.RData")
# writeRaster(myPred,"koalasightingprobability.TIF")
# koalasightingprob.r<- raster("koalasightingprobability.TIF")
# 
# tpo.crop <- crop(x = koalasightingprob.r, y = r, snap = "near")
# windows(); plot(tpo.crop)
# writeRaster(tpo.crop,"koalasightingprobability2.TIF")
# save.image("koalasightingprobability.RData")
# library(swirl)
# swirl()
# ravi

#now use this detection model to hefley and ccoorect IPP model.
#first create a IPP model. acsel22 is poinit dataset of koala. add habitat rasters create a stack as in hefley method.

myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder.
myfullstack = stack(myfullstack.a)
PP.pro <- acsel22[c(1,2,3)]
IPP.pre <- acsel22[c(1,2)]   # koala location data.
IPP.pre2=as.data.frame(extract(myfullstack,IPP.pre))
IPP.pre3=cbind(IPP.pre2, IPP.pro)

#####Step 4: Estimate the probability of detection for each presence-only location.####
p.det=ilogit(predict(myGLM,new=IPP.pre3))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
hist(p.det, breaks=70)
#)
#IPP.data$p.det=p.det   # IPP.data number of obserarions=1461
IPP.pre3$p.det=p.det


IPP.koala=glm(Haskoala~twi + tpo + aspect + elev+habit2pc+hpop+lot_density+sbd,family="binomial",weights=(1/p.det)*1^(1-Haskoala),data=IPP.pre3) # IPP.data2 added.
summary(IPP.koala)


myPredkoala = predict(myfullstack, IPP.koala, type = "response")
plot(myPredkoala, xlab = "x", ylab= "y",main=" IPP model- number of koalas")
myPred3 = predict(myfullstack, ZTGLM.corrected, type = "response")


###### do a IWLR as Renner et al. 
#######
#########
mydatasighting <- read.csv("mydatasighting.csv", header = TRUE)


# Subset the data:
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
acsel200=mydatasighting_cleaned[c(127,128, 117)]
acsel210 <- subset(acsel200, yearnew==2011, select = c("X","Y"))
acsel21v <- acsel210

# get raaster and subset

myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$",full.names=TRUE) #Mark_s folder.
myfullstack = stack(myfullstack.a)

#subset raster

myfullstack.sub <-subset(myfullstack,c(1,45,44,22,29,32,33))
myfullstack.subdf= as.data.frame(myfullstack.sub) 

# koalas over the study area. create a test lyer, extraxt valuves then remove NA.

testlayer <- myfullstack.sub[[1]]

acsel220=cbind(acsel21v,(extract(testlayer,acsel210)))

# only data from study area. Remove NA rows and keep xy.

acsel220 <- as.data.frame(na.omit(acsel220))

acsel210 <- acsel220[c(1,2)] # get xy

names(acsel210) <- tolower(names(acsel210))

acsel210$Haskoala <- 1

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

values(r )[values(r) > 0] = 0
rspec <- r

plot(r)
cells <- cellFromXY(rspec, as.matrix(dpart[, c("x", "y")]))

rspec[cells] <- dpart$Haskoala

plot(rspec)
text(dpart$x, dpart$y, lab = dpart$Haskoala)
Press<- as.data.frame(rspec, xy=TRUE)   # Has snome NA valuves.

#rename varibael awc

colnames(Press)[3] <- "koala"
ptest=cbind(Press,myfullstack.subdf)   ### now remove na
Press= na.omit(ptest)  # same as sp.at=Press ;  Pres=koala

# save this for future use in comapringlikelihood.
Press.my=Press

# koala presence locations

Pres <- Press[,3]   # this create a vector, Press[c(3)] create a sublist.

### now get design matrix

X.des=Press[,4:10]
X.des=as.matrix(X.des)  # quadrature points.

####iwlr
up.wt = (1.e-6)^(1 - Pres)  # up.wt = (10^6)^(1 - Pres)
iwlr = glm(Pres ~ X.des, family = binomial(), weights = up.wt)

dd1 <- as.data.frame(X.des) # get coordinates of the design matrix for predictions. similr to warton method.
pred.iwlr = predict(iwlr, newdata=dd1)

#r <- raster(myfullstack.sub, layer=2) 

dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only

pred.iwlr <- cbind(xydatan, pred.iwlr)
pred.iwlreg <- rasterFromXYZ(as.data.frame(pred.iwlr)[, c("x", "y", "pred.iwlr")])
plot(pred.iwlreg,asp=1)

###DWPR
p.wt = rep(1.e-6, length(Pres))
p.wt[Pres == 0] = 8310/sum(Pres == 0)
dwpr = glm(Pres/p.wt ~ X.des, family = poisson(), weights = p.wt)

dd <- as.data.frame(X.des)
pred.dwpr = predict(dwpr, newdata=dd)

dfull <- as.data.frame(r, xy = TRUE)

xydatan <- Press[c(1,2)]

# get coordinates only

pred.dwpr <- cbind(xydatan, pred.dwpr)
pred.dwpreg <- rasterFromXYZ(as.data.frame(pred.dwpr)[, c("x", "y", "pred.dwpr")])
plot(pred.dwpreg,asp=1)


##############
#5.3 Assessing the variability in likelihood for different numbers quadrature points.

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
sp.dat$Pres = 1
loglik = rep(NA, length(n.quad))

for (i in 1:length(n.quad)){
  quad = get(paste("quad.", n.quad[i], sep = ""))
  quad$Pres = 1
  all.dat = na.omit(data.frame(rbind(sp.dat, quad)))
  X.des = as.matrix(cbind(poly(all.dat$AnnualMeanTemperature, all.dat$twi, all.dat$tpo,
                               all.dat$distance_trunkandlink, all.dat$habit3decimal,  
                               all.dat$nitro,all.dat$roadk, degree = 2, raw = TRUE)))   
  p.wt = rep(1.e-8, dim(all.dat)[1])
  p.wt[all.dat$Pres == 0] = 10000/n.quad[i]
  z = all.dat$Pres/p.wt
  dwpr = glm(z ~ X.des, family = poisson(), weights = p.wt)
  mu = dwpr$fitted
  loglik[i] = sum(p.wt*(z*log(mu) - mu))
}
plot(n.quad, loglik, log = "x", type = "o")

#coudn`t get the Renners plot with all lines.





#next section starts here.

library(ppmlasso)
                        #4.1 Finding the appropriate spatial resolution for analysis
habitrasters<- subset(myfullstack, c(1,49,48,22,29,32,33)) # habitat covariates
bigquad <- as.data.frame(habitrasters, xy=TRUE, na.rm=T)
#stt <- na.omit(stt)
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'
# bigquad [is.na(bigquad )] <- 0
xydatan <- bigquad[c(1,2)]/1000
# stt requires xy as integers.
xydata <- as.data.frame(lapply(xydatan, as.integer)) # stt[] <- lapply(stt, as.integer)#this line edited on 04/01# make only xy integer in line with dadta shared with Mark S. 
bigquad <- cbind(xydata, bigquad[c(-1,-2)])

#IPP.pre2 <- acsel210 /1000
sp.xy = data.frame(IPP.pre2)
colnames(sp.xy)[1] <- 'X'; colnames(sp.xy)[2] <- 'Y'
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))
ppm.form = ~ poly(aspect, tpo, twi, elev, lot_density,habit2pc,hpop,lot_density, degree = 2, raw = TRUE)#+ poly(sqrt(habit2pc), sqrt(hpop), degree = 2, raw = TRUE)
scales = c( 0.5, 1, 2, 4, 8)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form)
                        #4.2 Fitting a regularisation path of point process models
#a LASSO penalty that optimises non-linear GCV
ppmFit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
 
#A regularisation path of Poisson point process models
quad.1k = sample.quad(bigquad, 1)
ppm.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, criterion = "nlgcv")

                      #4.3 Block cross-validation
#block cross-validation as a method for choosing the LASSO penalty
#area interaction model with radius 2k and lasso penalty by5-fold cross validation.
final.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, 
                     criterion = "blockCV", n.blocks = 5, block.size = 10)
#Predict and plot
pred.final.fit = predict.ppmlasso(final.fit, newdata=bigquad)

predictions.final.fit <- cbind(xydata, pred.final.fit) # xydatan was chnaged to xydata.
pred.final<- rasterFromXYZ(as.data.frame(predictions.final.fit)[, c("X", "Y", "pred.final.fit")])

plot(pred.final, main=" koala density-warton method/ bias not corrected")
###### runs this code without error.





####################### Generate quadrature points of different sizes. check R enviroment for genetted data.
library(ppmlasso)
load("Quad100m.RData")# this is from Renners example. add my quadratures.
# quad1 <- quad[c(1,2)]
# generte varying zise of quDRATURE POINTS BELOW SECTION
n.quad = c(1000, 2000, 5000, 10000, 20000)
quad.inc = sample(1:dim(quad)[1], 1000)
assign(paste("quad.", n.quad[1], sep = ""), quad[quad.inc[1:n.quad[1]],])
for (i in 2:length(n.quad)){
  quad.inc = c(quad.inc, sample(setdiff(1:dim(quad)[1], quad.inc),
                                (n.quad[i] - n.quad[i - 1])))
  assign(paste("quad.", n.quad[i], sep = ""), quad[quad.inc[1:n.quad[i]],])
}
#################################

# method in Renner et al 2015 Point process models for presence data. 

mydata = read.csv("Eucalyptus sparsifolia.csv")
XY <- mydata[c(3,4)]
XY=XY/1000 
save(XY, file = "Eucalyptus sparsifolia2.RData")
load("Eucalyptus sparsifolia2.RData")

#