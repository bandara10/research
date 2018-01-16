
#example http://rpubs.com/cyclemumner/294656
r <- raster(matrix(1:30, 5, 6))

## this is the full data set in data frame form
dfull <- as.data.frame(r, xy = TRUE)

## this is the partial data set, only the points with a valid value
## row-order doesn't matter, but we keep it for illustration
set.seed(10)
dpart <- dfull[sort(sample(seq_len(nrow(dfull)), 22)), ]

rspec <- raster(r)  ## this drops the data, keeps the structure

rspec[] <- NA_real_
## these names were nominated above, and might be different for a different
## input
cells <- cellFromXY(rspec, as.matrix(dpart[, c("x", "y")]))

rspec[cells] <- dpart$layer

plot(rspec)
text(dpart$x, dpart$y, lab = dpart$layer)

####real koala data
#for Renners IWLR and DWPR.


r <- raster(mystack, layer=2)
dfull <- as.data.frame(r, xy = TRUE)
dpart = cbind(acsel21,extract(r,acsel21[,1:2]))
dpart <- subset(acsel21, Haskoala==1)
rspec <- raster(r)
rspec[] <- 0


cells <- cellFromXY(rspec, as.matrix(dpart[, c("x", "y")]))
rspec[cells] <- dpart$Haskoala

plot(rspec)
text(dpart$x, dpart$y, lab = dpart$Haskoala)
f <- as.data.frame(rspec, xy=TRUE)
