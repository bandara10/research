
min_iter <- 5 # Minimum number of iterations
max_iter <- 10 # Maximum number of iterations
p_train <- 0.025 # Subsampling of the initial data


res_rf <- dissever(
  coarse = r, # stack of fine resolution covariates
  fine = pc, # coarse resolution raster
  method = "rf", # regression method used for disseveration
  p = p_train, # proportion of pixels sampled for training regression model
  min_iter = min_iter, # minimum iterations
  max_iter = max_iter # maximum iterations
 )

rrt=raster.downscale(r, pc, scatter=TRUE)


elev <- raster::getData('alt', country='SWZ', mask=TRUE)
tmax <- raster::getData('worldclim', var='tmax', res=10, lon=8.25, lat=46.8)
tmax <- crop(tmax[[1]], extent(elev))

tmax.ds <- raster.downscale(mytp, temp.crop, scatter=TRUE)
par(mfrow=c(2,2))
plot(tmax, main="Temp max")
plot(elev, main="elevation")
plot(tmax.ds$downscale, main="Downscaled Temp max")


library(raster)
library(spatialEco)

#step :1

temp=raster("meanann.txt")# load temperature australia map wgs 84

# step:2

myt=habitat.rr$Annual_Mean_Temperature #  my raster temp

mytp<- projectRaster(from = myt, crs = CRS("+init=epsg:4326")) # reproject to wgs.

temp.crop=crop(temp, mytp)   # crop to myt area.

# step 3: downscale

tmax.ds <- raster.downscale(mytp, temp.crop, scatter=TRUE)

# plot
par(mfrow=c(2,2))
plot(temp.crop)
plot(tmax.ds$downscale, main="Downscaled Temp max")
plot(myt)


