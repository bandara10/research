r.val <- getValues(myfullstack)# all raster valuves
is.na(r.val) <- sapply(r.val, is.infinite) # convert infinite to na
r.val = na.omit(r.val) # omint na valuves
pca <- princomp(r.val, cor=T)
summary(pca)
plot(pca)
biplot(pca)
biplot(pca, choices = 2:3)


# for the raster stack first tow PCA.
r.val <- getValues(myfullstack)# all raster valuves
is.na(r.val) <- sapply(r.val, is.infinite) # convert infinite to na
r.val = na.omit(r.val) # omint na valuves
idx <- which(!is.na(r.val))
ncomp <- 2 # first two principle components
r.pca <- myfullstack[[1:ncomp]]
for(i in 1:ncomp) { r.pca[[i]][idx] <- pca$scores[,i] } 
plot(r.pca, asp=3)

pc1 = psych::principal(r.val, nfactors = length(1), rotate="none")
plot(pc1$values, type = 'b') 