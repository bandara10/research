setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack\\wartondata")
source("functions2.r")
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
s.occupancy <- scale(d)
#w.files=c("dis_city.tif", "distance_tertiaryandlink.tif")
w.files=c("hpop.tif")
dd <- stack(w.files)
s.detection <- scale(dd)
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

s.detection= raster(get(w.files[1]))
for (i in w.files) {
  temp=get(i)
  names(temp) = i
  s.detection = addLayer(s.detection, temp)
}
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


# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # v1 and covariate valuves
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection) # v1, probability of detection.

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



# #Checking whether occupancy and detection rasters have the same resolution -----
# if(sum(res(s.occupancy)!=res(s.detection)))
#   stop("Occupancy and detection raster layers have different resolution")
# 
# if(ncell(s.occupancy)!=ncell(s.detection))
#   stop("Occupancy and detection have different number of cells")


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


################
#Belwo code not run.



myenv <- list.files(path="wartondata", pattern="\\.tif$", full.names = TRUE) #path="wartondata",
myenv.stack <- stack(myenv)
#turning rasters into tables - background, adding column of ones
s.occupancy <- dropLayer(myenv.stack,c(4,5,6))
s.detection <- dropLayer(myenv.stack,c(1,2,3,7,8,9,10,4))


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

#this is presence locations and covariates do by extract.
#Use X.back and extract valuves.
kolaxyT <- read.csv("koalaxy.csv", header = TRUE)
kolaxyT <- kolaxyT *1000
#kolaxyT <- kolaxyT[c(-1)]
kolaxy2 <- subset(kolaxyT, X > 442000 & X < 540000)
kolaxyT <- subset(kolaxy2, Y > 6902000 & Y < 7000000) # xy within the area only.
pb.occupancy = extract(s.occupancy, kolaxyT)

#Use X.back and extract vlauve.

pb.detection = extract(s.detection, kolaxyT)
# adding column of ones - po locations # Is this pb.detection?????
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # pb.occupancy is all presence locations. 
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)

pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)
