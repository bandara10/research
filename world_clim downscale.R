library(spatialEco)
#####raster downscaling####
#doesnt work
min_iter <- 5 # Minimum number of iterations
max_iter <- 10 # Maximum number of iterations
p_train <- 0.025 # Subsampling of the initial data


res_rf <- dissever::dissever(
  coarse =temp.crop, # stack of fine resolution covariates
  fine = mytp, # coarse resolution raster
  method = "rf", # regression method used for disseveration
  p = p_train, # proportion of pixels sampled for training regression model
  min_iter = min_iter, # minimum iterations
  max_iter = max_iter # maximum iterations
 )



########raster.downscale#######

elev <- raster::getData('alt', country='SWZ', mask=TRUE)
tmax2<- raster::getData('worldclim', var='tmax', res=10, lon=8.25, lat=46.8)
tmax <- crop(tmax[[1]], extent(elev))

tmax.ds <- raster.downscale(mytp, temp.crop, scatter=TRUE)
par(mfrow=c(2,2))
plot(tmax, main="Temp max")
plot(elev, main="elevation")
plot(tmax.ds$downscale, main="Downscaled Temp max")


library(raster)
library(spatialEco)

#step :1


elev <- raster::getData('alt', country='AUS', mask=TRUE)
tmax2<- raster::getData('worldclim', var='tmax', res=10, lon=8.25, lat=46.8)
tmax <- crop(tmax[[1]], extent(elev))



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


##### get global datasets/admin boundaries.######
elev <- raster::getData('alt', country='AUS', mask=TRUE)
australia0 <- raster::getData('GADM', country="AUS", level=1)

CRS("+init=epsg:28356")



plot(australia0, add=TRUE)

# get australia
austria0 <- getData('GADM', country='AUS', level=0)
# Select Dataset: The first argument specifies the dataset. 'GADM' returns the global administrative boundaries.
# Select country: The second argument provides the country name of the boundaries by using its ISO A3 country code (more info here)
# Specify level: The third argument specifies the level of of administrative subdivision (0=country, 1=first level subdivision).


climate <- getData('worldclim', var='bio', res=2.5)
# Select Dataset: The first argument specifies the dataset. 'worldclim' returns the World Climate Data.
# Select variable: The second argument specifies the variable: 'tmin', 'tmax', 'prec' and 'bio' (more info here).
# Specify resolution:  0.5, 2.5, 5, and 10 (minutes of a degree). In the case of res=0.5, you must also provide a lon and lat argument for a tile



####### SRTM 90 Elevation#######
# Last but not least, lets have a look at the SRTM 90 Data. We will use the getData() function one last time:
#   
  srtm <- getData('SRTM', lon=16, lat=48)
  
  
# Select Dataset: The first argument specifies the dataset. 'SRTM' returns the SRTM 90 elevation data.
# Specify Lon: The second argument specifies the lon of the SRTM tile.
# Specify Lat:  The second argument specifies the lat of the SRTM tile.
  
  
  #Download two more tiles
  srtm2 <- getData('SRTM', lon=13, lat=48)
  srtm3 <- getData('SRTM', lon=9, lat=48)
  
  #Mosaic/merge srtm tiles
  srtmmosaic <- mosaic(srtm, srtm2, srtm3, fun=mean)
  
    
    plot(srtmmosaic, main="Elevation (SRTM)")
  plot(austria0, add=TRUE)
  
  #http://www.gis-blog.com/r-raster-data-acquisition/
