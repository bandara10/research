library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files( path="wartondata", pattern="\\.tif$", full.names = TRUE) 
myenv.stack <- stack(myenv)

crs(myenv.stack) <- NA



 
extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)


plot(myenv.stack[[3]])
####### create a new raster of the same extent as stack and assign values 
start_y <- 441829 #100
start_x <- 6901098 #200
griddf <- expand.grid(X = seq(from = start_y, by = 1, l = 100),
                      Y = seq(from = start_x, by = 1, l = 100))

A <- rnorm(10000, 55, 6)
B <- rnorm(10000, 28, 20)

envs<- cbind(griddf,A,B)

xy.rr <- rasterFromXYZ(as.data.frame(envs)[, c("X", "Y", "A")])

plot(xy.rr)


######
ab.r <- myenv.stack[[3]]
tem <- getValues(ab.r) # give this valuve to new raster created in simulated dataset.
temp <- as.vector (tem)
ab.rnew<- setValues(xy.rr, temp)
plot(ab.rnew)
ab <- as.data.frame(ab.rnew, xy=TRUE)

#######
abc.r <- myenv.stack[[4]]
rain <- getValues(abc.r)
rain <- as.vector(rain)
abc.rnew <- setValues(xy.rr, rain)
abc<- as.data.frame(abc.rnew, xy=TRUE)


xydatan <- abc[c(1,2)]


######
abcd <- cbind(ab,abc)
env.v <- abcd[c(1:3,6)]
env.v[is.na(env.v)] <- 0
colnames(env.v)[1] <- "X"
colnames(env.v)[2] <- "Y"
#plot(env.v, pch = ".", add = T)
str(env.v)
kolaxy <- read.csv("wartondata\\kolaxy.csv", header = TRUE)
kolaxy <- kolaxy[,1:ncol(kolaxy)]/1000 # selects every row and 2nd to last columns

coordinates(kolaxy)= ~X+Y
plot(kolaxy, add=TRUE)

points(kolaxy$X, kolaxy$Y, pch = 15, col = "red")
mydata <- list(koala=kolaxy ,  env.var=env.v)
save(mydata, file="mydata.RData")
load(file="mydata.RData")
mydata$koala

ppmForm = ~  poly(A, A.1, degree = 2)
ppmFit = ppmlasso(ppmForms, sp.xy = kolaxy, env.grid = env.v, sp.scale = 1.5)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = BlueMountains$env
head(newEnv)
#newEnv$D_MAIN_RDS = min(stand.D_MAIN_RDS)
pred.biasCorrect.1 = predict(ppmFits, newdata=env.v)
predictions <- cbind(xydatan, pred.biasCorrect.1)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("x", "y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = kolaxy,
        env.grid = env.v, formula = ppmForms)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFits, fun = Kinhom)
resid.plot = diagnose(ppmFits, which = "smooth", type = "Pearson")



#
#lapply outputs a list of rasters
# m is data matrix from getvalues stack.
#r is an empry raster
ll <- getValues(myenv.stack)
lll <- as.data.frame(ll)

#create an empty raster
start_y <- 441829 #100
start_x <- 6901098 #200
griddf <- expand.grid(X = seq(from = start_y, by = 1, l = 100),
                      Y = seq(from = start_x, by = 1, l = 100))

A <- rnorm(10000, 55, 6)
B <- rnorm(10000, 28, 20)

envs<- cbind(griddf,A,B)

xy.rr <- rasterFromXYZ(as.data.frame(envs)[, c("X", "Y", "A")])

####
l <- lapply(1:ncol(lll), function(i) {
 setValues(xy.rr,lll[ , i])
} )

# stack the list
s <- stack(l)
ss <- as.data.frame(s, xy=TRUE)
ss[is.na(ss)] <- 0
colnames(ss)[1] <- "X"
colnames(ss)[2] <- "Y"


## Pre-standardise observer bias variables
stand.A.1=scale.default(ss$A.1, center = TRUE, scale = TRUE) #standarise
stand.A.2=scale.default(ss$A.2, center = TRUE, scale = TRUE) #standarise
ss$A.1 = stand.A.1
ss$A.2 = stand.A.2

#########+ poly (A.3, A.4, degree = 2)
ppmForm = ~  poly(A.3, A.4, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = kolaxy, env.grid = ss, sp.scale = 1)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = ss
newEnv$A.1 = min(stand.A.1)
pred.biasCorrect = predict(ppmFit, newdata=newEnv)
#newEnv$D_MAIN_RDS = min(stand.D_MAIN_RDS)
pred.biasCorrect.1 = predict(ppmFit, newdata=ss)
predictions <- cbind(xydatan, pred.biasCorrect.1)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("x", "y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = kolaxy,
        env.grid = env.v, formula = ppmForms)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFits, fun = Kinhom)
resid.plot = diagnose(ppmFits, which = "smooth", type = "Pearson")

