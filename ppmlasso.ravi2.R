library(ppmlasso)
library(spatstat)
library(raster)

# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

start_y <- 441829 #100
start_x <- 6901098 #200
griddf <- expand.grid(X = seq(from = start_y, by = 1, l = 100),
                      Y = seq(from = start_x, by = 1, l = 100))

A <- rnorm(10000, 55, 6)
B <- rnorm(10000, 28, 20)

envs<- cbind(griddf,A,B)

xy.rr <- rasterFromXYZ(as.data.frame(envs)[, c("X", "Y", "A")])

plot(xy.rr)
xy.rr2 <- disaggregate(xy.rr, fact=100)


#coordinates(griddf) <- ~X+Y
#plot(griddf)

#SPP locations



D <- griddf[sample(nrow(griddf), 500), ] 


#coordinates(E) <- ~X+Y
#plot(E, add=TRUE, col="red")
  

ppmForm = ~ poly(A,B, degree = 1) 
ppmFit = ppmlasso (ppmForm, sp.xy = D,
                  env.grid = envs, sp.scale =1)


pred.biasCorrect.1 = predict(ppmFit, newdata=envs)
xydatan <- envs[c(1,2)]

predictions <- cbind(xydatan, pred.biasCorrect.1)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = D,
        env.grid = envs, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")









