library(raster)
library(fields)
library(mvtnorm)

# Utility functions
logit = function(pp) { log(pp) - log(1-pp) }
expit = function(eta) {1/(1+exp(-eta))}



# negative loglikelihood function for Poisson point process

negLL.pp = function(param) {
  
  beta = param[1:dim(X.pp)[2]]
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  
  logL.pp = sum(X.pp %*% beta) - sum(mu)
  
  (-1)*sum(logL.pp)
}




# negative loglikelihood function for thinned Poisson point process

negLL.po = function(param) {
  
  beta = param[1:dim(X.po)[2]]
  alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
  
  lambda = exp(X.back %*% beta)
  mu = lambda * area.back
  p = expit(W.back %*% alpha)
  
  logL.po = sum(X.po %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)
  
  (-1)*sum(logL.po)
}



# Observed hessian matrix of negative loglikelihood function for thinned Poisson point process

ObsInfo.po = function(param) {
  
  beta = param[1:dim(X.po)[2]]
  alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
  
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



#  negative of log-likelihood function for N-mixture model of point counts parameterized for spatial resolution

negLL.pc = function(param) {
  beta = param[1:dim(X.pc)[2]]
  lambda.pc = exp(X.pc %*% beta)
  alpha = param[(dim(X.pc)[2]+1):(dim(X.pc)[2]+dim(W.pc)[3])]
  p.pc = matrix(nrow=n.pc, ncol=J.pc)
  for (j in 1:J.pc) {
    p.pc[, j] = expit(as.matrix(W.pc[,j,], nrow=n.pc) %*% alpha)
  }
  
  Ninfty = 100
  logL.pc = rep(NA,n.pc)
  for (i in 1:n.pc) {
    
    siteSum = 0
    Nmin = ymax.pc[i]
    if (Nmin==0) {
      siteSum = exp(-lambda.pc[i]*area.pc[i])
      Nmin = 1
    }
    N = Nmin:Ninfty
    logSum = log(dpois(N, lambda=lambda.pc[i]*area.pc[i]))
    
    yvec = y.pc[i, jind.pc[i,]]
    pvec = p.pc[i, jind.pc[i,]]
    
    for (j in 1:length(yvec)) {
      logSum = logSum + dbinom(yvec[j], size=N, prob=pvec[j], log=TRUE)
    }
    
    siteSum = siteSum + sum(exp(logSum))
    logL.pc[i] = log(siteSum)
  }
  
  (-1)*sum(logL.pc)
}


negLL.1pc = function(param) {
  beta = param[1:dim(X.pc)[2]]
  lambda.pc = exp(X.pc %*% beta)
  alpha = param[(dim(X.pc)[2]+1):(dim(X.pc)[2]+dim(W.pc)[3])]
  p.pc = expit(Wmat.pc %*% alpha)
  
  logL.pc = dpois(y.pc, lambda=lambda.pc*area.pc*p.pc, log=TRUE)
  
  (-1)*sum(logL.pc)
}



negLL.poANDpc = function(param)  {
  
  param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
  param.pc = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pc)[3]))]
  negLL.po(param.po) + negLL.pc(param.pc)
}


negLL.poAND1pc = function(param)  {
  
  param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
  param.pc = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pc)[3]))]
  negLL.po(param.po) + negLL.1pc(param.pc)
}




# Simulate a population of individuals distributed over rectangular region S
# and analyze samples of this population obtained via:
# 1. presence-only observations
# 2. point-count surveys



# Define a rectangular region S
s.xmin = -1
s.xmax =  1
s.ymin = -1
s.ymax =  1

s.area =  (s.xmax-s.xmin)*(s.ymax-s.ymin)

npixels.x = 1000
npixels.y = 1000

s = raster(ncol=npixels.x, nrow=npixels.y, xmn=s.xmin, xmx=s.xmax, ymn=s.ymin, ymx=s.ymax)
s.loc = xyFromCell(s, 1:ncell(s))



# Simulate covariate values over discretization of region S
mu1.x = s.xmin + 0.75*(s.xmax-s.xmin)
mu1.y = s.ymin + 0.40*(s.ymax-s.ymin)
sigma1.x = 0.25*abs(s.xmax-s.xmin)
sigma1.y = 0.50*abs(s.ymax-s.ymin)
rho1.xy = 0.5
mu1 = c(mu1.x, mu1.y)
Sigma1 = matrix(c(sigma1.x^2, rep(rho1.xy*sigma1.x*sigma1.y, 2), sigma1.y^2), ncol=2)

mu2.x = s.xmin + 0.15*(s.xmax-s.xmin)
mu2.y = s.ymin + 0.80*(s.ymax-s.ymin)
sigma2.x = 0.50*abs(s.xmax-s.xmin)
sigma2.y = 0.25*abs(s.ymax-s.ymin)
rho2.xy = -0.4
mu2 = c(mu2.x, mu2.y)
Sigma2 = matrix(c(sigma2.x^2, rep(rho2.xy*sigma2.x*sigma2.y, 2), sigma2.y^2), ncol=2)

xcov = 0.4 * dmvnorm(s.loc, mean=mu1, sigma=Sigma1) + 0.6 * dmvnorm(s.loc, mean=mu2, sigma=Sigma2)
xcov = (xcov - mean(xcov))/sd(xcov)

values(s) = xcov
names(s) = 'x'

mu3.x = s.xmin + 0.25*(s.xmax-s.xmin)
mu3.y = s.ymin + 0.65*(s.ymax-s.ymin)
sigma3.x = 0.25*abs(s.xmax-s.xmin)
sigma3.y = 0.50*abs(s.ymax-s.ymin)
rho3.xy = 0.1
mu3 = c(mu3.x, mu3.y)
Sigma3 = matrix(c(sigma3.x^2, rep(rho3.xy*sigma3.x*sigma3.y, 2), sigma3.y^2), ncol=2)

wcov = dmvnorm(s.loc, mean=mu3, sigma=Sigma3)
wcov = (wcov - mean(wcov))/sd(wcov)


temp = raster(s)
values(temp) = wcov
## values(temp) = xcov   #   NOTE:  UNCOMMENT THIS LINE TO MAKE COVARIATES x and w IDENTICAL
names(temp) = 'w'
s = addLayer(s, temp)



# Assign values of model parameters and design matrices

beta.param = c(log(8000), 0.5)
X = cbind(rep(1, ncell(s)), values(s)[,'x'])

alpha.param = c(-1.0, -1.0)
W = cbind(rep(1, ncell(s)), values(s)[,'w'])


# ... compute expected limiting density of individuals and true detection probabilities over discretization of region S
temp = raster(s)
values(temp) = exp(X %*% beta.param)
names(temp) = 'lambda'
s = addLayer(s, temp)

temp = raster(s)
values(temp) = expit(W %*% alpha.param)
names(temp) = 'pTrue'
s = addLayer(s, temp)

maxlambda = max(values(s)[,'lambda'])

plot(temp)


# Begin loop to simulate the population and to analyze samples of the population

nsims = 10

CR = '\n'

header.po  = c('n','conv','recip','beta0','beta1','alpha0','alpha1','se.beta0','se.beta1','se.alpha0','se.alpha1')

cat(header.po, sep=',', file='output-po.csv')
cat(CR, file='output-po.csv', append=TRUE)


header.pc  = c('n.pc', 'conv','recip','beta0','beta1','alpha0','alpha1','se.beta0','se.beta1','se.alpha0','se.alpha1')

cat(header.pc, sep=',', file='output-pc.csv')
cat(CR, file='output-pc.csv', append=TRUE)


header.poANDpc  = c('n.pc','conv','recip','beta0','beta1','alpha0','alpha1','alpha0.pc','alpha1.pc','se.beta0','se.beta1','se.alpha0','se.alpha1','se.alpha0.pc','se.alpha1.pc')

cat(header.poANDpc, sep=',', file='output-poANDpc.csv')
cat(CR, file='output-poANDpc.csv', append=TRUE)


minrecipCondNum = 1e-6

for (sim in 1:nsims) {
  
  
  # ... simulate point pattern of individuals over discretization of S
  N.hpp = rpois(1, maxlambda*s.area)
  
  ind.hpp = sample(1:ncell(s), size=N.hpp, replace=FALSE)   #  sampling w/o replacement ensures only 1 indiv per pixel
  loc.hpp = s.loc[ind.hpp, ]
  lambda.hpp = values(s)[,'lambda'][ind.hpp]
  
  ind.ipp = runif(N.hpp, 0,1) <= lambda.hpp/maxlambda
  N.ipp = sum(ind.ipp)
  loc.ipp = loc.hpp[ind.ipp, ]
  
  
  
  # ... simulate presence-only data (= detections of individuals as a thinned point process)
  
  pTrue.ipp = values(s)[,'pTrue'][ind.hpp][ind.ipp]
  y.ipp = rbinom(length(pTrue.ipp), size=1, prob=pTrue.ipp)
  
  ind.po = (y.ipp == 1)
  loc.po = loc.ipp[ind.po, ]
  
  
  # ... analyze simulated presence-only data
  
  # ... first establish discrete grid for integrating over region S;
  # ... then compute area, average value of covariates, and number of individuals in each element of grid
  
  sgrid = aggregate(s, fact=2, fun=mean)
  
  sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))
  sgrid.xcov = values(sgrid)[,'x']
  sgrid.wcov = values(sgrid)[,'w']
  
  
  area.back = rep(res(sgrid)[1]*res(sgrid)[2], ncell(sgrid))
  X.back = cbind(rep(1, ncell(sgrid)), sgrid.xcov)
  W.back = cbind(rep(1, ncell(sgrid)), sgrid.wcov)
  
  
  X.po = X[ind.hpp[ind.ipp][ind.po], ]  # thinned observations
  W.po = W[ind.hpp[ind.ipp][ind.po], ]  # thinned observations
  
  betaGuess = rep(0, dim(X.po)[2])
  alphaGuess = rep(0, dim(W.po)[2])
  paramGuess = c(betaGuess, alphaGuess)
  
  fit.po = optim(par=paramGuess, fn=negLL.po, method='BFGS', hessian=FALSE)
  
  recipCondNum.po = NA
  se.po = rep(NA, length(fit.po$par))
  if (fit.po$convergence==0) {
    hess = ObsInfo.po(fit.po$par)
    ev = eigen(hess)$values
    recipCondNum.po = ev[length(ev)]/ev[1]
    if (recipCondNum.po>minrecipCondNum) {
      vcv = chol2inv(chol(hess))
      se.po = sqrt(diag(vcv))
    }
  }
  
  cat(c(sum(ind.po), fit.po$convergence, recipCondNum.po, fit.po$par, se.po),  sep=',', file='output-po.csv', append=TRUE)
  cat(CR, file='output-po.csv', append=TRUE)
  
  
  
  # .... simulate point-count data
  
  # .... first compute abundances within a sample frame
  spop = dropLayer(s, c('lambda','pTrue'))
  temp = raster(spop)
  z = rep(0, ncell(spop))
  z[ ind.hpp[ind.ipp] ] = 1
  values(temp) = z
  names(temp) = 'presence'
  spop = addLayer(spop, temp)
  
  # .... form sample frame (= grid of sample units) by aggregating every no. cols and every no. rows in spop
  gridfact = c(10,10)
  sgrid = aggregate(spop, fact=gridfact, fun=mean)
  abund = aggregate(subset(spop, 'presence'), fact=gridfact, fun=sum)
  names(abund) = 'N'
  sgrid = addLayer(sgrid, abund)
  rm(abund)
  
  
  sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))
  sgrid.xcov = values(sgrid)[,'x']
  sgrid.wcov = values(sgrid)[,'w']
  sgrid.N = values(sgrid)[,'N']
  
  
  # .... select sample units randomly;
  # .... then compute observed count(s) for each sample unit
  
  sampleFraction.pc = c(0.005, 0.01, 0.02, 0.04, 0.08)
  
  for (samp in 1:length(sampleFraction.pc)) {
    
    n.pc = floor(ncell(sgrid)*sampleFraction.pc[samp])  
    ind.pc = sample(1:ncell(sgrid), size=n.pc, replace=FALSE)  # simple random sample w/o replacement
    N.pc = values(sgrid)[,'N'][ind.pc]
    area.pc = rep(res(sgrid)[1]*res(sgrid)[2], length(N.pc))
    xcov.pc = values(sgrid)[,'x'][ind.pc]
    wcov.pc = values(sgrid)[,'w'][ind.pc]
    
    X.pc = cbind(rep(1,n.pc), xcov.pc)
    
    # .... compute observed counts
    J.pc = 4
    alpha.pc = c(0, -1.0)
    W.pc = array(dim=c(n.pc, J.pc, 2))
    W.pc[,,1] = 1
    W.pc[,,2] = matrix(rep(wcov.pc,J.pc), ncol=J.pc)
    Wmat.pc = cbind(rep(1,n.pc), wcov.pc)
    
    pTrue.pc = matrix(nrow=n.pc, ncol=J.pc)
    y.pc = matrix(nrow=n.pc, ncol=J.pc)
    for (j in 1:J.pc) {
      pTrue.pc[, j] = expit(as.matrix(W.pc[,j,], nrow=n.pc) %*% alpha.pc)
      y.pc[, j] = rbinom(n.pc, size=N.pc, prob=pTrue.pc[, j])
    }
    
    
    # .... analyze simulated point count data
    
    jind.pc = !is.na(y.pc)  # index for non-missing observations of y
    ymax.pc = apply(y.pc,1,max, na.rm=TRUE)
    betaGuess = rep(0, dim(X.pc)[2])
    alphaGuess.pc = rep(0, dim(W.pc)[3])
    paramGuess = c(betaGuess, alphaGuess.pc)
    fit.pc = NA
    if (J.pc>1) {
      fit.pc = optim(par=paramGuess, fn=negLL.pc, method='BFGS', hessian=TRUE)
    }
    else {
      fit.pc =  optim(par=paramGuess, fn=negLL.1pc, method='BFGS', hessian=TRUE)
    }
    
    
    recipCondNum.pc = NA
    se.pc = rep(NA, length(fit.pc$par))
    if (fit.pc$convergence==0) {
      hess = fit.pc$hessian
      ev = eigen(hess)$values
      recipCondNum.pc = ev[length(ev)]/ev[1]
      if (recipCondNum.pc>minrecipCondNum) {
        vcv = chol2inv(chol(hess))
        se.pc = sqrt(diag(vcv))
      }
    }
    
    cat(c(n.pc, fit.pc$convergence, recipCondNum.pc, fit.pc$par, se.pc),  sep=',', file='output-pc.csv', append=TRUE)
    cat(CR, file='output-pc.csv', append=TRUE)
    
    
    # .... analyze simulated presence-only data AND point count data
    
    paramGuess = c(betaGuess, alphaGuess, alphaGuess.pc)
    
    fit.poANDpc = NA
    if (J.pc>1) {
      fit.poANDpc = optim(par=paramGuess, fn=negLL.poANDpc, method='BFGS', hessian=TRUE)
    }
    else {
      fit.poANDpc = optim(par=paramGuess, fn=negLL.poAND1pc, method='BFGS', hessian=TRUE)
    }
    
    
    recipCondNum.poANDpc = NA
    se.poANDpc = rep(NA, length(fit.poANDpc$par))
    if (fit.poANDpc$convergence==0) {
      hess = fit.poANDpc$hessian
      ev = eigen(hess)$values
      recipCondNum.poANDpc = ev[length(ev)]/ev[1]
      if (recipCondNum.poANDpc>minrecipCondNum) {
        vcv = chol2inv(chol(hess))
        se.poANDpc = sqrt(diag(vcv))
      }
    }
    
    cat(n.pc, c(fit.poANDpc$convergence, recipCondNum.poANDpc, fit.poANDpc$par, se.poANDpc),  sep=',', file='output-poANDpc.csv', append=TRUE)
    cat(CR, file='output-poANDpc.csv', append=TRUE)
    
    
    
    
    
  }  # end of sampleFraction loop
  
  cat('Completed analysis of data set #', sim, " ", CR)
  
  
} # end of sim loop

