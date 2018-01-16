###################################################################################################

###################################################################################################
######### Data required for the function ##########################################################
###################################################################################################
## All the data required for this function is stored in a file data.rda
## s.occupancy - raster with background covatiates that effect occupancy
## s.detection - raster with background covariates that effect detection
## pb.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
## pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
##
## y.so - matrix of detection/non detection of the PA surveys
## so.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
## so.detection - matrix with covariates that effect detection in the locations of PA survey sites in each survey
###################################################################################################

library(raster)
library(fields)
library(mvtnorm)
library(matrixStats)

source("functions2.r")

#  C:\Users\uqrdissa\ownCloud\Covariates_analysis\Mark_S\Supplimentary materials of papers\Vira Koshinika paper\Data

# Data Preparation -----------------------------------------------------------------------------
#### read spatial covariates #plot(visible_sky_sept2012_100_no.grd)
x.files=c( "Dem75mInteger_100_no.grd"
           ,"MinTemp_Jul_100_no.grd"
           ,"MaxTemp_Jan_100_no.grd"
           ,"wetness_index_saga_sept2012_100_no.grd"
           ,"Evaporation_Jul_100_no.grd"
           ,"Evaporation_Jan_100_no.grd"
           ,"RainDays_Jul_100_no.gri"
           ,"Rainfall_Jan_100_no.grd"
           ,"visible_sky_sept2012_100_no.grd")
#,"Rainfall_Jul_100_no.grd"
#,"RainDays_Jan_100_no.gri"
w.files=c("log_vertical_distance_major_streams_sept2012_100_no.grd", "terrain_ruggedness_index_sept2012_100_no.grd")
#"distance_to_roads_100_no",
for (i in c(x.files, w.files)) {
  do.call("=", list(i, raster(i)))
}
s.occupancy= raster(get(x.files[1])) # no values associated with this RasterLayer

for (i in x.files) {
  temp=get(i)
  names(temp) = i
  s.occupancy = addLayer(s.occupancy, temp)
} 
#plot(s.occupancy) # this is raster stack.

s.detection= raster(get(w.files[1]))
for (i in w.files) {
  temp=get(i)
  names(temp) = i
  s.detection = addLayer(s.detection, temp)
}
#plot(s.detection) #0-10
#plotting occupancy and detection covariates distributions in the study area
print('plotting occupancy and detection rasters')
ppi = 300
png('occupancy-covariates.png', width=9*ppi, height=round(length(x.files)/2+2)*ppi, res=ppi)
plot(s.occupancy)
dev.off()

png('pb-detection-covariates.png', width=9*ppi, height=round(length(w.files)/2+2)*ppi, res=ppi)
plot(s.detection)
dev.off()

# read pb and PA survey data
#List of locations of sightings of yellow-bellied glider in opportunistic surveys
pb=read.csv("yb-pb-location.csv") # xy locations over 6000.
#s <- raster("MinTemp_Jul_100_no.grd")
#plot(pb.loc, add=TRUE)

pb.loc=SpatialPoints(pb) # over 6629 beyond study area. 
# get locations over rasters only.
pb.occupancy=extract(s.occupancy,pb.loc) 
pb.detection=extract(s.detection,pb.loc) 
# only data from study area.
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,] #3663#  distance covariates extract from rasters. only one here
pb.occupancy=pb.occupancy[is.complete.pb,]#3663# env covariates extract from  raster
# upto here basically covariate extraction is done for presence data for occupancy/abunace and detection., .
#######################################
#######################################
so.occupancy=read.csv("so_occupancy.csv") # occupancy covariate valuves for survey data
so.detection=read.csv("so_detection.csv") # detection covariate valuves wind sky day 1 and 2 for survey data.
y.so=read.csv("yb-so.csv") # survey data or site occupancy, day 1 day 2 presence absence.

print("removing NA")
is.complete=complete.cases(so.occupancy)&complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*0.04

print('removing raster files')
#removing rasters to free the memory
for (i in c(x.files,w.files)) {
  do.call(remove,list(i))
}
remove(is.complete.pb,po,yb.so,po.loc)


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


# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # v1 and covariate valuves
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection) # v1, probability of detection.



#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
so.occupancy <- as.matrix(so.occupancy) # added by me
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
#X.so$v1 <- X.so$`rep(1, nrow(as.matrix(so.occupancy)))
# X.so <- X.so[c(-1)]
# X.so <- X.so[c(13, 1:12)]
#X.so <- as.matrix(X.so)
W.so = array(dim=c(nrow(as.matrix(pb.detection)), J.so, 3))
W.so[,,1] = 1
W.so[,,2] = pb.detection# if it changes
W.so[,,3] = pb.detection2# if it changes


######################didnt work up to tis level.


#Checking whether occupancy and detection rasters have the same resolution -----
if(sum(res(s.occupancy)!=res(s.detection)))
  stop("Occupancy and detection raster layers have different resolution")

if(ncell(s.occupancy)!=ncell(s.detection))
  stop("Occupancy and detection have different number of cells")


# Plotting covariates that drive occupancy and detection in PO
ppi = 300
png('occupancy-covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.occupancy)
dev.off()

png('PO-detection -covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.detection)
dev.off()

# 1. Preparing the data ========================================================
#Preparing Presence-Only data ------------------------------------------------repeat???

#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
colnames(X.back)=c("",names(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))
colnames(W.back)=c("",names(s.detection))
# remove all NA values
tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell
s.area=area.back*nrow(X.back) #study area

# adding column of ones - po locations # Is this pb.detection?????
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # pb.occupancy is all presence locations. 
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)



#Preparing Presence-Absence data# seems to be a repat ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)

W.so = array(dim=c(nrow(as.matrix(so.detection)), J.so, 1))
W.so[,,1] = 1
W.so[,,2] = so.detection[,1:2]# if it changes
W.so[,,3] = so.detection[,3:4]# if it changes
# 2. Analising the data ========================================================

#Analyzing Presence-Only data
pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)  
# W.pb is pb.detection
# X.pb and X.back same covariates. 
# W.back is distance covariate for detection. 
# Analyzing presence-absence data
so.fit=so.model(X.so,W.so,y.so) 

# Analyzing presence-only data AND presence-absence data
poANDso.fit=pbso.integrated(X.pb, W.pb,X.back, W.back,X.so,W.so,y.so)


# pb.fit
#$coefs
# Parameter name      Value Standard error
# 1                                                    beta0 -5.2043587     0.13225605
# 2                                 Dem75mInteger_100_no.grd -1.0188289     0.09002529
# 3                                   MinTemp_Jul_100_no.grd -0.7245087     0.05534356
# 4                                   MaxTemp_Jan_100_no.grd -0.4412986     0.08397245
# 5                   wetness_index_saga_sept2012_100_no.grd -0.2431656     0.02429644
# 6                               Evaporation_Jul_100_no.grd  0.6241968     0.15377184
# 7                               Evaporation_Jan_100_no.grd  0.2788767     0.08262562
# 8                                  RainDays_Jul_100_no.gri  0.6352531     0.07357142
# 9                                  Rainfall_Jan_100_no.grd  1.6488527     0.03502068
# 10                         visible_sky_sept2012_100_no.grd -0.1527510     0.01963205
# 11                                                  alpha0  2.2552370     0.37952045
# 12 log_vertical_distance_major_streams_sept2012_100_no.grd  0.3345543     0.11314596
# 13            terrain_ruggedness_index_sept2012_100_no.grd -1.0139156     0.08707117
# 
# $convergence
# [1] 0
# 
# $value
# [1] 13028.05



# > so.fit
# $coefs
# Parameter name        Value Standard error
# 1                                                beta0 -35.32066783     13.8974792
# 2                                 Dem75mInteger_100_no  -3.85970335      2.1075431
# 3                                   MinTemp_Jul_100_no   7.03926349      3.1381702
# 4                                   MaxTemp_Jan_100_no   2.97951068      2.8663198
# 5  log_vertical_distance_major_streams_sept2012_100_no   1.49008050      0.7158995
# 6                   wetness_index_saga_sept2012_100_no   0.11506930      0.2674998
# 7                           Evaporation_Jul_100_no.grd -10.10350950      4.9318864
# 8                           Evaporation_Jan_100_no.grd -23.97880500     10.0095411
# 9                              RainDays_Jul_100_no.gri -11.04918569      2.9825634
# 10                             Rainfall_Jul_100_no.grd   1.83391730      0.8050790
# 11                             RainDays_Jan_100_no.gri  -2.87105727      2.6243461
# 12                             Rainfall_Jan_100_no.grd   2.61632927      1.4959734
# 13                         visible_sky_sept2012_100_no  -0.25346487      0.1857578
# 14                                           alpha0.so   0.02252261      0.3022473
# 
# $convergence
# [1] 0
# 
# $value
# [1] 114.0533