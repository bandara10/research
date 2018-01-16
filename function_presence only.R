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

