setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack\\wartondata")
source("functions2.r")
library(raster)
x.files=c( "AnnualMeanTemperature.tif"
           ,"clay.tif"
           ,"elev.tif"
           ,"AnnualPrecipitation.tif"
           ,"habit0pc.tif"
           ,"habit1pc.tif"
           ,"habit2pc.tif"
           ,"habit3pc.tif"
          )
d <- stack(x.files)
s.occupancy <- scale(d,scale=TRUE,center = TRUE)

# distance_tertiaryandlink <- raster("distance_tertiaryandlink.tif")
# distance_tertiaryandlink<- sqrt (distance_tertiaryandlink)
# writeRaster(distance_tertiaryandlink, "distance_tertiaryandlink2.tif",overwrite=TRUE)
# dis_city <- raster("dis_city.tif")
# dis_city <- dis_city*.001
# writeRaster(dis_city,"dis_city2.tif",overwrite=TRUE)



w.files=c("dis_city2.tif","distance_tertiaryandlink2.tif")# get distance in km.
#w.files=c("hpop.tif")
dd <- stack(w.files)

s.detection <- scale(dd,scale=TRUE, center = TRUE)
#"distance_to_roads_100_no",
# for (i in c(x.files, w.files)) {
#   do.call("=", list(i, raster(i)))
# }
# s.occupancy= raster(get(x.files[1])) # no values associated with this RasterLayer
# 
# for (i in x.files) {
#   temp=get(i)
#   names(temp) = i
#   s.occupancy = addLayer(s.occupancy, temp)
# } 
# 
# s.detection= raster(get(w.files[1]))
# for (i in w.files) {
#   temp=get(i)
#   names(temp) = i
#   s.detection = addLayer(s.detection, temp)
# }
#plot(s.detection) #0-10
# #plotting occupancy and detection covariates distributions in the study area
# print('plotting occupancy and detection rasters')
# ppi = 300
# png('occupancy-covariates.png', width=9*ppi, height=round(length(x.files)/2+2)*ppi, res=ppi)
# plot(s.occupancy)
# dev.off()
# 
# png('pb-detection-covariates.png', width=9*ppi, height=round(length(w.files)/2+2)*ppi, res=ppi)
# plot(s.detection)
# dev.off()
load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
pb=mydatasighting_cleaned[c(127,128)]
# or
pb=read.csv("kolaxy.csv")
pb=read.csv("kolaxyT2010.csv")

# kolaxy2 <- subset(kolaxyT, X > 442000 & X < 540000)
# pb <- subset(kolaxy2, Y > 6902000 & Y < 7000000) # xy within the area only.
# pb <- as.matrix(pb)
pb.loc=SpatialPoints(pb) # over 6629 beyond study area. 
# get locations over rasters only.
pb.occupancy=extract(s.occupancy,pb.loc) 
pb.detection=extract(s.detection,pb.loc) 
# only data from study area.
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,] #3663#  distance covariates extract from rasters. only one here
pb.occupancy=pb.occupancy[is.complete.pb,]#3663# env covariates extract from  raster
# upto here basically covariate extraction is done for presence data for occupancy/abunace and detection., .

print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection)) #s.detection is raster stack.

# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

print("specifying area ")
#area in squared km ??????????????????????? -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell


s.area=area.back*nrow(X.back) #study area


# # adding column of ones - po locations
# X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # v1 and covariate valuves
# W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection) # v1, probability of detection.

# #add a column of ones to the PA covariat
# #y.so # matrix of presences and absences (when the species was and wasn't present)
# J.so=ncol(y.so)
# so.occupancy <- as.matrix(so.occupancy) # added by me
# X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
# #X.so$v1 <- X.so$`rep(1, nrow(as.matrix(so.occupancy)))
# # X.so <- X.so[c(-1)]
# # X.so <- X.so[c(13, 1:12)]
# #X.so <- as.matrix(X.so)
# W.so = array(dim=c(nrow(as.matrix(pb.detection)), J.so, 1))
# W.so[,,1] = 1
# W.so[,,2] = pb.detection# if it changes
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


# 2. Analising the data ========================================================

#Analyzing Presence-Only data
pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back) 
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

################
coef <- pb.fit$coefs
coeff <- coef[c(1,2)]
f <- s.occupancy
f1 <- subset(f,c(1))*coeff$Value[[2]]
f2 <- subset(f,c(1))*coeff$Value[[3]]
f3 <- subset(f,c(1))*coeff$Value[[4]]
f4 <- subset(f,c(1))*coeff$Value[[5]]
f5 <- subset(f,c(1))*coeff$Value[[6]]
f6 <- subset(f,c(1))*coeff$Value[[7]]
f7 <- subset(f,c(1))*coeff$Value[[8]]
f8 <- subset(f,c(1))*coeff$Value[[9]]
fn <-  exp(f1+f2+f3+f4+f5+f6+f7+f8)+coeff$Value[[1]]
plot(fn, main= "intensity 2010 data")

ff <- s.detection
ff1 <- subset(ff,c(1))*coeff$Value[[11]] 
ff2 <- subset(ff,c(2))*coeff$Value[[12]]
fn2 <-  (f1+f2)+coeff$Value[[10]]   

library(optiRum)
fnn <- logit.prob(fn2)
plot(fnn, main =" probability of detection")
fn3 <- fn*fnn

plot(fn3, main = "bias corrected intensity")
#################
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


ff <- s.detection
f9 <- subset(ff,c(1))*-2.3293703# -0.6404511
f10<- subset(ff,c(2))*-0.6387529# -0.6387529 
fn2 <- (f9+f10)-0.6404511 
library(optiRum)
fnn <- logit.prob(fn2)
plot(fnn)
fn2 <- fn*fnn
plot(fn2)

lgalimited.shp <- readShapePoly("LGA10.shp")
lgal <- crop(fn2,lgalimited.shp)
d <- mask(lgal,lgalimited.shp )
plot(d, axes=FALSE, box=FALSE)
plot(lgalimited.shp, add=TRUE)

plot(fn2, zlim = c(1), main=" koala density-warton method/ bias not corrected")

plot(pb.loc, add=TRUE)
# fn3 <- as.data.frame(fn2, xy=TRUE)
# fn4 <- fn3$layer
# fn5 <- Logit(fn4,min = 0, max = 1) #library(DescTools)
# plot(fn5)
# xx <- cbind(fn3, fn4)
d <- as.data.frame(s.detection, xy=TRUE)   # wrong approach.
b <- d[c(-3,-4)]
f <- d[c(-1,-2)]
f1 <- f$dis_city2
f2 <- f$distance_tertiaryandlink2
a <- exp(f1*-2.3293703+f2*-0.6387529-0.6404511)
pp <- cbind(b,a)
pred.ct <- rasterFromXYZ(as.data.frame(pp)[, c("x", "y", "a")])
plot(pred.ct)

plot(pb.loc, add=TRUE)


fn3 <- fn*fn2
plot(fn3)
plot(fn3, zlim=c(-3, 5))