
#vector projection.


# mydata in lat-lon:
coords <- SpatialPoints(mydata[, c("lng", "lat")])
mydata.ll <- SpatialPointsDataFrame(coords, mydata)
proj4string(mydata.ll) <- CRS("+init=epsg:4326") 

# Transform to MGA 56:
mydata.mga <- spTransform(mydata.ll, CRS("+init=epsg:28356"))


#raster projection

projection(tpo.r) <- CRS("+init=epsg:4326")
windows(); plot(tpo.r)

# Reproject to MGA-56:
tpo.r <- projectRaster(from = tpo.r, crs = CRS("+init=epsg:28356"))