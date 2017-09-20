#
#https://www.linkedin.com/pulse/naive-principal-component-analysis-r-pablo-bernabeur.val <- getValues(myfullstack)# all raster valuves
dim(myfullstack)
is.na(r.val) <- sapply(r.val, is.infinite) # convert infinite to na
r.val = na.omit(r.val) # omint na valuves
pca <- princomp(r.val, cor=T)
summary(pca)
plot(pca)
biplot(pca)
biplot(pca, choices = 2:3)


# for the raster stack first two PCA.
r.val <- getValues(myfullstack)# all raster valuves
is.na(r.val) <- sapply(r.val, is.infinite) # convert infinite to na
r.val = na.omit(r.val) # omint na valuves
idx <- which(!is.na(r.val))
ncomp <- 2 # first two principle components
r.pca <- myfullstack[[1:ncomp]]
for(i in 1:ncomp) { r.pca[[i]][idx] <- pca$scores[,i] } 
plot(r.pca, asp=3)
#####


#Two sources will determine the number of components to select for the next stage:
  
#Kaiser's criterion: components with SS loadings > 1. In our example, only PC1( 9.87).
#A more lenient alternative is Joliffe's criterion, SS loadings > .7.
#Scree plot: the number of points after point of inflexion. For this plot, call:
pc1 = psych::principal(r.val, nfactors = length(1), rotate="none")
plot(pc1$values, type = 'b') 

pc2 = psych::principal(r.val, nfactors=2, 
                       rotate = "varimax", scores = TRUE)
