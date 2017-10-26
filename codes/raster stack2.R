#below raster is correct resolution and change raster stack resolution and extenet and write them.

#writeRaster(stack(D), names(D), bylayer=TRUE, format="ascii",overwrite=TRUE)

unzip("vector\\AU_Qld_study_area.zip")
studyareafull.shp <- readShapePoly("AU_Qld_detail_studyarea_outline-MGA56.shp")
proj4string(studyareafull.shp) <- CRS("+init=epsg:28356") 
windows();plot(studyareafull.shp)
studyareafull <- raster(studyareafull.shp)
res(studyareafull) <- 1000

hpop2<-spatial_sync_raster(hpop.r, studyareafull, method = "ngb")
writeRaster(hpop2,"raster_stack\\hpop2.TIF", overwrite=TRUE)

stack.n <- stack( "aspect2.tif", "awc2.tif","fpc2.tif", "slope2.tif","topo2.tif" ,"hpop2.tif","tpo2.tif", "rain2.tif", "nitro2.tif", "elev2.tif", "temp2.tif", "sbd2.tif")

windows(); hist(stack.n)
windows(); levelplot(stack.n)
windows(); plot(stack.n)
windows(); density(stack.n)
# calculating linne densiy
# Create some sample lines
l1 <- Lines(Line(cbind(c(0,1),c(.25,0.25))), ID="a")
l2 <- Lines(Line(cbind(c(0.25,0.25),c(0,1))), ID="b")
sl <- SpatialLines(list(l1,l2))
pspSl <- as.psp(sl)

lines.shp <- readShapeLines("AU_Qld_study_area_mask_roads\\AU_Qld_study_area_mask_lines-MGA56.shp")
proj4string(lines.shp) <- CRS("+init=epsg:28356") 
river.lines <- lines.shp[lines.shp$waterway == "river",]



intersect(lines.shp, studyareafull.shp)
windows(); plot(studyareafull.shp)
names(lines.shp@data)
#select colomn waterway and  select catogory rivers.
river.lines <- lines.shp[lines.shp$waterway == "river",]
windows(); plot(river.lines,lwd=1)




# Convert SpatialLines to psp object using maptools library
pspSl <- as.psp(lines.shp)
# Pixellate with resolution of 0.5, i.e. 2x2 pixels
px <- pixellate(pspSl, eps=1000)
# This can be converted to raster as desired
rLength <- raster(px)
rLength <- crop(x = rLength, y = dstudyarea.shp)
writeRaster(rLength,"rLength2.TIF", overwrite=TRUE)
plot(rLength)
plot(dstudyarea.shp)
rLength.n<-spatial_sync_raster(rLength, studyareafull, method = "ngb")
