#Spatial model

model{
  # process model
  #random-effect for a - loop through years
  
  for (i in 1:4)
  {
    a[i] ~ dlnorm(amean,tau1)
  }
  # loop through sites
  
  for (i in 1:NSites)
  {
    #intercept random-effect
    
    intmean[i] <- sum(Y[i,1,] * alpha)
    int[i] ~ dnorm(intmean[i],tau2)
    
    # sample density for first time step
    Lambda[i,1] <- exp(int[i] + sum(X[i,1,] * beta))
    # D[i,j] is density
    # mean of dgamma(a, b) is a/b, so if Lambda[i] = a/b, we sample a and
    # calculate b deterministically using Lambda
    b[i,1] <- a[trunc((Year[i,1] - 1995)/5.1) + 1]/Lambda[i,1]
    D[i,1] ~ dgamma(a[trunc((Year[i,1] - 1995)/5.1) + 1],b[i,1])
    # simulate for model adequacy checks
    #D.sim[i,1] ~ dgamma(a[trunc((Year[i,1] - 1995)/5.1) + 1],b[i,1])
    # loop through remaining time steps
    for (j in 2:NSteps)
    {
      # sample density for remaining time step
      Lambda[i,j] <- exp(int[i] + sum(X[i,j,] * beta))
      # D[i,j] is density
      # mean of dgamma(a, b) is a/b, so if Lambda[i] = a/b, we sample a and
      
      # calculate b deterministically using Lambda
      b[i,j] <- a[trunc((Year[i,j] - 1995)/5.1) + 1]/Lambda[i,j]
      D[i,j] ~ dgamma(a[trunc((Year[i,j] - 1995)/5.1) + 1],b[i,j])
      # simulate for model adequacy checks (uncomment to use)
      #D.sim[i,j] ~ dgamma(a[trunc((Year[i,j] - 1995)/5.1) + 1],b[i,j])
    }
  }
  # observation model - strip transects
  # count data
  for (i in 1:NStrips)
  {
    #sample observed counts
    SNTrue[i] <- ifelse((round(D[SIndex[i,1],SIndex[i,2]] * Area[i]) - (D[SIndex[i,1],SIndex[i,2]] *
                                                                          Area[i])) <= 0.5,round(D[SIndex[i,1],SIndex[i,2]] * Area[i]),trunc(D[SIndex[i,1],SIndex[i,2]] *
                                                                                                                                               Area[i]))
    NAll[i] ~ dbin(p,SNTrue[i])
    #goodness-of-fit & model adequacy (uncomment to use)
    #compute discrepancy statistics for observed data
    #expect[i] <- p * Lambda[SIndex[i,1],SIndex[i,2]] * Area[i]
    #res.obs[i] <- NAll[i] - expect[i]
    #R.obs[i] <- pow((NAll[i] - expect[i]),2) / expect[i]
    #compute discrepancy statistics for simulated data
    #SNT.sim[i] <- ifelse((round(D.sim[SIndex[i,1],SIndex[i,2]] * Area[i]) -
    (D.sim[SIndex[i,1],SIndex[i,2]] * Area[i])) <= 0.5,round(D.sim[SIndex[i,1],SIndex[i,2]] *
                                                               #Area[i]),trunc(D.sim[SIndex[i,1],SIndex[i,2]] * Area[i]))
                                                               #y.sim[i] ~ dbin(p,SNT.sim[i])
                                                               #res.sim[i] <- y.sim[i] - expect[i]
                                                               #R.sim[i] <- pow((y.sim[i] - expect[i]),2) / expect[i]
  }
  # observation model - line transects
  # perpendicular distance data
  # Distance Sampling (estimate f0) - this is separate from the density estimate
  # because this should not be applied to any records with 0 observations, whereas
  # the density estimate is applied to all line transects
  for (i in 1:NDists)
  {
    #likelihood function for half-normal
    zeros[i] ~ dpois(phi[i])
    phi[i] <- -(log(2 * DistLambda / 3.141593) / 2 - DistLambda * pow(PDist[i],2) / 2)
  }
  f0 <- sqrt(2 * DistLambda / 3.141593)
  # count data
  for (i in 1:NLines)
  {
    #sample observed counts
    NAll[NStrips + i] ~ dpois((D[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000))
    #goodness-of-fit
    #compute discrepancy statistics for observed data
    #expect[NStrips + i] <- (Lambda[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000)
    #res.obs[NStrips + i] <- NAll[NStrips + i] - expect[NStrips + i]
    #R.obs[NStrips + i] <- pow((NAll[NStrips + i] - expect[NStrips + i]),2) / expect[NStrips + i]
    #compute discrepancy statistics for simulated data
    #y.sim[NStrips + i] ~ dpois((D.sim[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000))
    #res.sim[NStrips + i] <- y.sim[NStrips + i] - expect[NStrips + i]
    #R.sim[NStrips + i] <- pow((y.sim[NStrips + i] - expect[NStrips + i]),2) / expect[NStrips + i]
  }
  #calculate goodness-of-fit statistics (uncomment to use)
  #fit.obs <- sum(R.obs[])
  #fit.sim <- sum(R.sim[])
  #fit.test <- fit.obs - fit.sim
  # priors
  amean ~ dnorm(0,0.001)
  tau1 <- sig1^-2
  sig1 ~ dunif(0,10)
  tau2 <- sig2^-2
  sig2 ~ dunif(0,10)
  for (i in 1:Ny)
  {
    alpha[i] ~ dnorm(0,0.001)
  }
  for (i in 1:Nx)
  {
    beta[i] ~ dnorm(0,0.001)
  }
  p ~ dunif(0,1)
  DistLambda ~ dgamma(0.001,0.001)
}
Trend model
model {
  # process model
  #random-effect for a - loop through years
  for (i in 1:4)
  {
    a[i] ~ dlnorm(amean,tau1)
  }
  # loop through sites
  for (i in 1:NSites)
  {
    #intercept for random-effect on initial population size
    intmean[i] <- sum(Y[i,1] * alpha)
    int[i] ~ dnorm(intmean[i],tau2)
    # sample density for first time step
    Lambda[i,1] <- exp(int[i])
    # D[i,j] is density
    # mean of dgamma(a, b) is a/b, so if Lambda[i] = a/b, we sample a and
    # calculate b deterministically using Lambda
    b[i,1] <- a[trunc((Year[i,1] - 1995)/5.1) + 1]/Lambda[i,1]
    D[i,1] ~ dgamma(a[trunc((Year[i,1] - 1995)/5.1) + 1],b[i,1])
    # simulate for model adequacy checks
    #D.sim[i,1] ~ dgamma(a[trunc((Year[i,1] - 1995)/5.1) + 1],b[i,1])
    #intercept for growth rate
    intr[i] <- sum(X[i,1,] * beta)
    # loop through remaining time steps
    for (j in 2:NSteps)
    {
      # sample density for remaining time step
      r[i,j - 1] <- exp(intr[i] + sum(Z[i,j - 1,] * epsilon))
      Lambda[i,j] <- D[i,j - 1] * r[i,j - 1]
      # D[i,j] is density
      # mean of dgamma(a, b) is a/b, so if Lambda[i] = a/b, we sample a and
      # calculate b deterministically using Lambda
      b[i,j] <- a[trunc((Year[i,j] - 1995)/5.1) + 1]/Lambda[i,j]
      D[i,j] ~ dgamma(a[trunc((Year[i,j] - 1995)/5.1) + 1],b[i,j])
      # simulate for model adequacy checks (uncomment to use)
      #D.sim[i,j] ~ dgamma(a[trunc((Year[i,j] - 1995)/5.1) + 1],b[i,j])
    }
  }
  # observation model - strip transects
  # count data
  for (i in 1:NStrips)
  {
    #sample observed counts
    SNTrue[i] <- ifelse((round(D[SIndex[i,1],SIndex[i,2]] * Area[i]) - (D[SIndex[i,1],SIndex[i,2]] *
                                                                          Area[i])) <= 0.5,round(D[SIndex[i,1],SIndex[i,2]] * Area[i]),trunc(D[SIndex[i,1],SIndex[i,2]] * Area[i]))
    NAll[i] ~ dbin(p,SNTrue[i])
    #goodness-of-fit & model adequacy (uncomment to use)
    #compute discrepancy statistics for observed data
    #expect[i] <- p * Lambda[SIndex[i,1],SIndex[i,2]] * Area[i]
    #res.obs[i] <- NAll[i] - expect[i]
    #R.obs[i] <- pow((NAll[i] - expect[i]),2) / expect[i]
    #compute discrepancy statistics for simulated data
    #SNT.sim[i] <- ifelse((round(D.sim[SIndex[i,1],SIndex[i,2]] * Area[i]) -
    (D.sim[SIndex[i,1],SIndex[i,2]] * Area[i])) <= 0.5,round(D.sim[SIndex[i,1],SIndex[i,2]] *
                                                               #Area[i]),trunc(D.sim[SIndex[i,1],SIndex[i,2]] * Area[i]))
                                                               #y.sim[i] ~ dbin(p,SNT.sim[i])
                                                               #res.sim[i] <- y.sim[i] - expect[i]
                                                               #R.sim[i] <- pow((y.sim[i] - expect[i]),2) / expect[i]
  }
  # observation model - line transects
  # perpendicular distance data
  # Distance Sampling (estimate f0) - this is separate from the density estimate
  # because this should not be applied to any records with 0 observations, whereas
  # the density estimate is applied to all line transects
  for (i in 1:NDists)
  {
    #likelihood function for half-normal
    zeros[i] ~ dpois(phi[i])
    phi[i] <- -(log(2 * DistLambda / 3.141593) / 2 - DistLambda * pow(PDist[i],2) / 2)
  }
  f0 <- sqrt(2 * DistLambda / 3.141593)
  # count data
  for (i in 1:NLines)
  {
    #sample observed counts
    NAll[NStrips + i] ~ dpois((D[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000))
    #goodness-of-fit & model adequacy (uncomment to use)
    #compute discrepancy statistics for observed data
    #expect[NStrips + i] <- (Lambda[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000)
    #res.obs[NStrips + i] <- NAll[NStrips + i] - expect[NStrips + i]
    #R.obs[NStrips + i] <- pow((NAll[NStrips + i] - expect[NStrips + i]),2) / expect[NStrips + i]
    #compute discrepancy statistics for simulated data
    #y.sim[NStrips + i] ~ dpois((D.sim[LIndex[i,1],LIndex[i,2]] * 2 * LL[i]) / (f0 * 10000))
    #res.sim[NStrips + i] <- y.sim[NStrips + i] - expect[NStrips + i]
    #R.sim[NStrips + i] <- pow((y.sim[NStrips + i] - expect[NStrips + i]),2) / expect[NStrips + i]
  }
  #calculate goodness-of-fit statistics (uncomment to use)
  #fit.obs <- sum(R.obs[])
  #fit.sim <- sum(R.sim[])
  #fit.test <- fit.obs - fit.sim
  # priors
  amean ~ dnorm(0,0.001)
  tau1 <- sig1^-2
  sig1 ~ dunif(0,10)
  tau2 <- sig2^-2
  sig2 ~ dunif(0,10)
  for (i in 1:Ny)
  {
    alpha[i] ~ dnorm(0,0.001)
  }
  for (i in 1:Nx)
  {
    beta[i] ~ dnorm(0,0.001)
  }
  for (i in 1:Nz)
  {
    epsilon[i] ~ dnorm(0,0.001)
  }
  p ~ dunif(0,1)
  DistLambda ~ dgamma(0.001,0.001)
}