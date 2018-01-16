setwd("C:/Users/uqrdissa/ownCloud/Covariates_analysis/Mark_S/raster_syn")#/renamed_vars

###################################################################################################

###################################################################################################

#likelihood functions
# Utility functions
logit = function(pp) { log(pp) - log(1-pp) }
expit = function(eta) {1/(1+exp(-eta))}

# function that calculates a probability of occupancy for each location in the dataset X
predict=function(mymodel, X){
  beta=mymodel$coefs[1:ncol(X),'Value']
  
  lambda=exp(X %*% beta)
  # psi =1- exp(-lambda*area)
  return(lambda)
}


#Function that fits IPP model
pb.ipp=function(X.po, W.po,X.back, W.back){
  
  beta.names=colnames(X.back)
  beta.names[1]='beta0'
  
  alpha.names=colnames(W.back)
  alpha.names[1]='alpha0'
  
  par.names=c(beta.names,	alpha.names)
  
  
  
  minrecipCondNum = 1e-6
  paramGuess = c(rep(.1, ncol(X.po)), rep(-.1, ncol(W.po)))
  
  
  fit.po = optim(par=paramGuess, fn=negLL.po, method='BFGS', hessian=FALSE
                 , X.po=X.po, W.po=W.po,X.back=X.back,W.back=W.back ) # params for likelyhood function
  
  
  # calculating se with Hessian matrix
  recipCondNum.po = NA
  se.po = rep(NA, length(fit.po$par))
  if (fit.po$convergence==0) {
    hess = ObsInfo.po(fit.po$par, X.po, W.po,X.back, W.back)
    ev = eigen(hess)$values
    recipCondNum.po = ev[length(ev)]/ev[1]
    if (recipCondNum.po>minrecipCondNum) {
      vcv = chol2inv(chol(hess))
      se.po = sqrt(diag(vcv))
    }
  }
  
  #printing PO results
  tmp= data.frame(par.names,fit.po$par,se.po)
  names(tmp)=c('Parameter name', 'Value', 'Standard error')
  p=NULL
  p$coefs=tmp
  p$convergence=fit.po$convergence
  p$optim_message=fit.po$message
  p$value=fit.po$value
  # print("Estimated parameters beta and alpha", quote=FALSE)
  # print(p)
  return(p)
}
# negative loglikelihood function for Poisson point process
negLL.pp = function(param) {
  
  beta = param[1:dim(X.pp)[2]]
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  
  logL.pp = sum(X.pp %*% beta) - sum(mu)
  
  (-1)*sum(logL.pp)
}

# negative loglikelihood function for thinned Poisson point process
negLL.po = function(param, X.po, W.po,X.back, W.back) {
  
  beta = param[1:dim(X.po)[2]]
  alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
  # dim(X.back)
  # length(beta)
  # length(area.back)
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  logL.po = sum(X.po %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
  
  (-1)*sum(logL.po)
}

# Observed hessian matrix of negative loglikelihood function for thinned Poisson point process
ObsInfo.po = function(param, X.po,W.po,X.back, W.back) {
  
  beta = param[1:dim(X.back)[2]]
  alpha = param[(dim(X.back)[2]+1):(dim(X.back)[2]+dim(W.back)[2])]
  
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  p.po = expit(W.po %*% alpha)
  
  nxcovs = length(beta)
  nwcovs = length(alpha)
  
  nparams = nxcovs + nwcovs
  Hmat = matrix(nrow=nparams, ncol=nparams)
  
  #  beta partials
  for (i in 1:nxcovs) {
    for (j in 1:i) {
      Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
      Hmat[j,i] = Hmat[i,j]
    }
  }
  
  # alpha partials
  for (i in 1:nwcovs) {
    for (j in 1:i) {
      Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.po[,i] * W.po[,j] * p.po * (1-p.po))
      Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
    }
  }
  
  # alpha-beta partials
  for (i in 1:nwcovs) {
    for (j in 1:nxcovs) {
      Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
      Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
    }
  }
  
  Hmat
}

# Expected hessian matrix of negative loglikelihood function for thinned Poisson point process
FisherInfo.po = function(param) {
  
  beta = param[1:dim(X.back)[2]]
  alpha = param[(dim(X.back)[2]+1):(dim(X.back)[2]+dim(W.back)[2])]
  
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  
  nxcovs = length(beta)
  nwcovs = length(alpha)
  
  nparams = nxcovs + nwcovs
  Hmat = matrix(nrow=nparams, ncol=nparams)
  
  #  beta partials
  for (i in 1:nxcovs) {
    for (j in 1:i) {
      Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
      Hmat[j,i] = Hmat[i,j]
    }
  }
  
  # alpha partials
  for (i in 1:nwcovs) {
    for (j in 1:i) {
      Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.back[,i] * W.back[,j] * p * (1-p) * mu * p)
      Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
    }
  }
  
  # alpha-beta partials
  for (i in 1:nwcovs) {
    for (j in 1:nxcovs) {
      Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
      Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
    }
  }
  
  Hmat
}

########
library(raster)
library(maptools)
x.files = c(

"AnnualPrecipitation.tif"
,"AnnualMeanTemperature.tif"
, "habit0percent.tif"
,"habit1percent.tif"
,"habit2percent.tif"
,"habit3percent.tif"
,"dem.tif"
,"CLY_100_200.tif"
)
d <- stack(x.files)
lgalimited.shp <- readShapePoly("LGA10new.shp") # uncomment to use 10 LGA shape file.
d <- crop(d ,lgalimited.shp ) #uncomment to use 10 LGA 
s.occupancy <- scale(d,scale=TRUE,center = TRUE)

w.files=c(
"distance_tertiaryandlink.tif"
)
dd <- stack(w.files)
lgalimited.shp <- readShapePoly("LGA10new.shp")
dd <- crop(dd ,lgalimited.shp)##mask =crop/uncomment to use 10 LGA 
s.detection <- scale(dd,scale=TRUE, center = TRUE)

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
pb=mydatasighting_cleaned[c(127,128)]
## or
pb=read.csv("kolaxy.csv")
## or
pb=read.csv("kolaxyT2010.csv")

# kolaxy2 <- subset(kolaxyT, X > 442000 & X < 540000)
# pb <- subset(kolaxy2, Y > 6902000 & Y < 7000000) # xy within the area only.
# pb <- as.matrix(pb)
pb.loc=SpatialPoints(pb) # over 6629 beyond study area. 
# get locations over rasters only.
pb.occupancy=extract(s.occupancy,pb.loc) 
pb.detection=extract(s.detection,pb.loc) 
# only data from study area.
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,] #3663#  distance covariates extract from rasters. only one here
pb.occupancy=pb.occupancy[is.complete.pb,]#3663# env covariates extract from  raster
# upto here basically covariate extraction is done for presence data for occupancy/abunace and detection., .

print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection)) #s.detection is raster stack.

# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

print("specifying area ")
#area in squared km ??????????????????????? -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell


s.area=area.back*nrow(X.back) #study area


# # adding column of ones - po locations
# X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # v1 and covariate valuves
# W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection) # v1, probability of detection.

# #add a column of ones to the PA covariat
# #y.so # matrix of presences and absences (when the species was and wasn't present)
# J.so=ncol(y.so)
# so.occupancy <- as.matrix(so.occupancy) # added by me
# X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
# #X.so$v1 <- X.so$`rep(1, nrow(as.matrix(so.occupancy)))
# # X.so <- X.so[c(-1)]
# # X.so <- X.so[c(13, 1:12)]
# #X.so <- as.matrix(X.so)
# W.so = array(dim=c(nrow(as.matrix(pb.detection)), J.so, 1))
# W.so[,,1] = 1
# W.so[,,2] = pb.detection# if it changes
# W.so[,,3] = pb.detection2# if it changes



# # Checking whether occupancy and detection rasters have the same resolution -----
# if(sum(res(s.occupancy)!=res(s.detection)))
#   stop("Occupancy and detection raster layers have different resolution")
# 
# if(ncell(s.occupancy)!=ncell(s.detection))
#   stop("Occupancy and detection have different number of cells")
# 
# 
# # Plotting covariates that drive occupancy and detection in PO
# ppi = 300
# png('occupancy-covariates.png', width=9*ppi, height=3*ppi, res=ppi)
# plot(s.occupancy)
# dev.off()
# 
# png('PO-detection -covariates.png', width=9*ppi, height=3*ppi, res=ppi)
# plot(s.detection)
# dev.off()
# adding column of ones - po locations # Is this pb.detection?????
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy) # pb.occupancy is all presence locations. 
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)


# 2. Analising the data ========================================================

#Analyzing Presence-Only data
(pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)) 
#estimates are obtained but not se.
# If any of the eigenvalues are zero or negative, there is obviously a problem, 
# perhaps stemming from under identification or boundary conditions
# One thing you might try is to refit the model from its solution, to see 
# if the gradients get smaller, and the NaN's clear up.
# > pb.fit
# $coefs
# Parameter name       Value Standard error
# 1                     beta0 -2.14667766     0.07994286
# 2     AnnualMeanTemperature  0.79053507     0.10671768
# 3                      clay -0.29408512     0.03372959
# 4                      elev -0.77522035     0.13323224
# 5       AnnualPrecipitation  0.35858355     0.03337058
# 6                  habit0pc  0.65510645     0.09293961
# 7                  habit1pc  0.15593408     0.05454184
# 8                  habit2pc  0.11365113     0.05051433
# 9                  habit3pc  0.32729882     0.07020533
# 10                   alpha0  0.14012092     0.18512155
# 11                 dis_city -3.27397922     0.27526465
# 12 distance_tertiaryandlink -0.07361067     0.12179921
# 
# $convergence
# [1] 0
# 
# $value
# [1] 2752.702
# 
# $value
# [1] 2783.46

coef <- pb.fit$coefs
coeff <- coef[c(1,2)]
f <- s.occupancy
f1 <- subset(f,c(1))*coeff$Value[[2]]
f2 <- subset(f,c(1))*coeff$Value[[3]]
f3 <- subset(f,c(1))*coeff$Value[[4]]
f4 <- subset(f,c(1))*coeff$Value[[5]]
f5 <- subset(f,c(1))*coeff$Value[[6]]
fn <-  exp(f1+f2+f3+f4+f5)+coeff$Value[[1]]
fn <- mask(fn ,lgalimited.shp)
plot(fn)

ff <- s.detection
f1 <- subset(ff,c(1))*coeff$Value[[7]] 
f2 <- subset(ff,c(2))*coeff$Value[[8]]
fn2 <-  (f1+f2)+coeff$Value[[6]]   

library(optiRum)
fnn <- logit.prob(fn2)
fnn <- mask(fnn ,lgalimited.shp)
plot(fnn)
fn3 <- fn*fnn
fn3 <- mask(fn3 ,lgalimited.shp)##mask
#plot(fn3,col=topo.colors(10,1),axes=FALSE, box=FALSE)
plot(fn3)
plot(fn3, zlim = c(.8, 15), main=" koala density-warton method/corrected")
#####
lgalimited.shp <- readShapePoly("LGA10new.shp")
lgal <- crop(fn3,lgalimited.shp)
d <- mask(lgal,lgalimited.shp )
plot(d,col=rainbow(12),axes=FALSE, box=FALSE)
plot(lgalimited.shp,border="black", lwd=.2, add=TRUE)

plot(lgalimited.shp, border="black", lwd=.2, add=TRUE)

