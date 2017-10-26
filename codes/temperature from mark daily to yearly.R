library(spatial.tools)
library(raster)

setwd("C:/Users/uqrdissa/Downloads/BOM_temps/BOM_temps/maximum_temperature_2001-2010/maximum_temperature_2001-2010")

grids <- list.files("C:/Users/uqrdissa/Downloads/BOM_temps/BOM_temps/maximum_temperature_2001-2010/maximum_temperature_2001-2010/New folder/", pattern = "*.txt")

for (i in 1:length(grids)){
  r <- raster(paste0("C:/Users/uqrdissa/Downloads/BOM_temps/BOM_temps/maximum_temperature_2001-2010/maximum_temperature_2001-2010/New folder/", grids[i]))  
  #perform the reclassifcation
  #write each reclass to a new file 
  writeRaster(r, filename = paste0 ("C:/Users/uqrdissa/Downloads/BOM_temps/BOM_temps/maximum_temperature_2001-2010/maximum_temperature_2001-2010/New folder/", grids[i]),
              format="GTiff", overwrite=TRUE)
}


grids.r <- list.files("C:/Users/uqrdissa/Downloads/BOM_temps/BOM_temps/maximum_temperature_2001-2010/maximum_temperature_2001-2010/New folder/",
                      pattern = "*.tif$",full.names = TRUE)
s <- stack(grids.r)
#####-----------------------------This section doesn`t work. ignor for the momment`
# assign a coordinate reference system then transform
proj4string(s) <- CRS("+init=epsg:28356")
mask <- raster("BDW_000_005.tif")
ss <- crop(x = s, y = mask, snap="near")

rc2 <- crop(r, extent(r, 5, 10, 7, 15))

extent.new <-  extent(150, 153, -30, -27)
s.new <- crop(x = s, y = extent.newE)
writeRaster(s.new, filename ="s.new.tif", overwrite=TRUE)
-------------------------
plot(s)
indices<-rep(1:8, each=8) # 8 is the number of layers
s.sum<-stackApply(s, indices, fun = sum)
s.sum2<-s.sum/8
plot(s.sum2)
writeRaster(s.sum2, filename ="s.sum2.tif")
# Get average yearly temperature .txt to a folder
# convert all to rasters
# stack rasters
#use the function to divide staack by number of layers. For a year 365 or the number of lraster layers.
# Crop the area required as this is whole australia.