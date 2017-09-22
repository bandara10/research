#First,example raster with values of [1,2] representing tow habitat types. and generate a random point sample
#if we create random negatives, same can be calculated for negative locations. This then becocme a new variables if necessasary.
landcover.prop <- list()
for( i in 1:3) {
  landcover.prop[[i]] <- extract(rc, hefleydata.presence, buffer=2000, small=TRUE, fun=function(hefleydata.presence, p = i) 
  { prop.table( ifelse(hefleydata.presence == p, 1, 0) ) } )
}

str(landcover.prop)
land <- data.frame(landcover.prop)
 

#The elements in the list object are vectors that are ordered the same as your points
#so, you can just add them to the SpatialPointsDataFrame object.
##lets reclassify a raster for this 
m <- c(0, 10, 1, 10, 150, 2,  150, 190, 3)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(fpc, rclmat)
plot(rc)
plot(hefleydata.presence, add=TRUE)
