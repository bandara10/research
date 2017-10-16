
# =========================================================================
#
#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#
#   Marc K???ry & J. Andy Royle
#
#   *** This is the text file with all R and BUGS code from the book ***
#
#   Created 2 Dec 2015 based on draft from 21 Oct 2015
#
# =========================================================================
### Last change: 19 May 2017 by Mike Meredith
# Incorporated errata up to 19 May 2017
# Code updated to implement new names and changes to functions in AHMbook package 0.1.4
# Use built-in data sets instead of .csv files
# In Chapter 11, replaced 'Y', 'Ysum' and 'Yaug' with lower case 'y', 'ysum' and 'yaug'
#  to match the code in the printed book.

# =========================================================================
#
#  1 Distribution, abundance and species richness in ecology
#
# =========================================================================



# 1.1 Point processes, distribution, abundance and species richness
# ------------------------------------------------------------------------


sim.fn(quad.size = 10, cell.size = 1, intensity = 1)


set.seed(82)
tmp <- sim.fn(quad.size = 16, cell.size = 2, intensity = 0.5)


# Effect of grain size of study on abundance and occupancy (intensity constant)
tmp <- sim.fn(quad.size = 10, cell.size = 1, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 2, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 5, intensity = 0.5)
tmp <- sim.fn(quad.size = 10, cell.size = 10, intensity = 0.5)


# Effect of intensity of point pattern (intensity) on abundance and occupancy
tmp <- sim.fn(intensity = 0.1) # choose default quad.size = 10, cell.size = 1
tmp <- sim.fn(intensity = 1)
tmp <- sim.fn(intensity = 5)
tmp <- sim.fn(intensity = 10)


simrep <- 100                 # Run 100 simulation reps ## See errata
grain <- c(0.1,0.2,0.25,0.5,1,2) # values will be fed into 'cell.size' argument
int <- seq(0.1, 3,,6)         # values will be fed into 'lambda' argument
n.levels <- length(grain)     # number of factor levels in simulation
results <- array(NA, dim = c(n.levels, n.levels, 2, simrep)) # 4-D array !
for(i in 1:n.levels){         # Loop over levels of factor grain
  for(j in 1:n.levels){       # Loop over levels of factor intensity
    for(k in 1:simrep){
      cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
      tmp <- sim.fn(cell.size = grain[i], intensity = int[j], show.plot = F)
      results[i,j,1:2,k] <- c(mean(tmp$N), tmp$psi)
    }
  }
}


# Plot these two prediction matrices (NOT IN BOOK)
par(mfrow = c(2, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
# Plot mean abundance in sampled quadrats
z1 <- apply(results[,,1,], c(1,2), mean)   # mean abundance
image(x=grain, y=int, z=z1, col = mapPalette(100), axes = T, xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot mean occupancy in sampled quadrats
z2 <- apply(results[,,2,], c(1,2), mean)   # mean occupancy
image(x=grain, y=int, z=z2, col = mapPalette(100), axes = T, xlab = "Grain size (cell.size)", ylab = "Intensity of PPP")
contour(x=grain, y=int, z=z2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
# Plot relationship between occupancy and abundance for whole range of abundance
plot(results[,,1,], results[,,2,], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
lines(smooth.spline(results[,,2,] ~ results[,,1,], df = 4), lwd = 3, col = "blue")
#abline(0, 1, lwd = 3)
# ... and only for very small abundance
keep <- results[,,1,] < 0.25
plot(results[,,1,][keep], results[,,2,][keep], xlab = "Mean abundance", ylab = "Occupancy", frame = FALSE)
abline(0, 1, lwd = 3)
lines(smooth.spline(results[,,2,][keep] ~ results[,,1,][keep], df = 4), lwd = 3, col = "blue")


# 1.2 Meta-population designs
# ------------------------------------------------------------------------


# 1.3 State and rate parameters
# ------------------------------------------------------------------------


# 1.4 Measurement error models in ecology
# ------------------------------------------------------------------------


# 1.5 Hierarchical models for distribution, abundance, and species richness
# ------------------------------------------------------------------------


# 1.6 Summary and outlook
# ------------------------------------------------------------------------


# 1.7 Exercises
# ------------------------------------------------------------------------










# =========================================================================
#
# 2. What are hierarchical models and how do we analyze them ?
#
# =========================================================================





# 2.1 Introduction
# ------------------------------------------------------------------------



# 2.2 Random variables, probability density functions, statistical models, probability, and statistical inference
# ------------------------------------------------------------------------


dbinom(0:5, size = 5, prob = 0.2)

pnorm(200, mean = 190, sd = 10)

pnorm(200, mean = 190, sd = 10) - pnorm(180, mean = 190, sd = 10)

f <- function(x, mu, sigma){
  (1 / sqrt(2*pi*sigma^2)) * exp( -((x-mu)^2)/(2*sigma^2))
}

integrate(f, lower = 180, upper = 200, mu = 190, sigma = 10)


# 2.2.1 Statistical models
# ------------------------------------------------------------------------

# 2.2.2 Joint, marginal, and conditional distributions
# ------------------------------------------------------------------------
Y <- 0:5   # Possible values of Y (# surveys with peregrine sightings)
X <- 0:5   # Possible values of X (# fledged young)
p <- plogis(-1.2 + 2*X) # p as function of X
round(p, 2)


# Joint distribution [Y, X]
lambda <- 0.4
joint <- matrix(NA, length(Y), length(X))
rownames(joint) <- paste("y=", Y, sep="")
colnames(joint) <- paste("x=", X, sep="")
for(i in 1:length(Y)) {
  joint[,i] <- dbinom(Y, 5, p[i]) * dpois(X[i], lambda)
}

round(joint, 3)

margX <- colSums(joint)
round(margX, 4)

margY <- rowSums(joint)
round(margY, 4)

YgivenX <- joint / matrix(margX, nrow(joint), ncol(joint), byrow=TRUE)
round(YgivenX, 2)


# 2.2.3 Statistical inference
# ------------------------------------------------------------------------

# 2.3 Hierarchical Models
# ------------------------------------------------------------------------

# 2.3.1 Two canonical hierarchical models in ecology
# ------------------------------------------------------------------------

# 2.3.1.1 The occupancy model for species distribution
# ------------------------------------------------------------------------

# 2.3.1.2 The N-mixture model for abundance
# ------------------------------------------------------------------------

# 2.3.2 The process of hierarchical modeling
# ------------------------------------------------------------------------

# 2.3.3 Inference for hierarchical models
# ------------------------------------------------------------------------

# 2.4 Classical inference based on likelihood
# ------------------------------------------------------------------------

# 2.4.1 The frequentist interpretation
# ------------------------------------------------------------------------

# 2.4.2 Properties of MLEs
# ------------------------------------------------------------------------

# 2.4.3 Delta approximation
# ------------------------------------------------------------------------

# 2.4.4 Example: Classical inference for logistic regression
# ------------------------------------------------------------------------
# Simulate a covariate called vegHt for 100 sites
set.seed(2014)  # Set seed so we all get the same values of vegHt
M <- 100        # Number of sites surveyed
vegHt <- runif(M, 1, 3) # uniform from 1 to 3

# Suppose that occupancy probability increases with vegHt
# The relationship is described by an intercept of -3 and
#    a slope parameter of 2 on the logit scale
beta0 <- -3
beta1 <- 2
psi <- plogis(beta0 + beta1*vegHt) # apply inverse logit

# Now we go to 100 sites and observe presence or absence
z <- rbinom(M, 1, psi)

# Definition of negative log-likelihood.
negLogLike <- function(beta, y, x) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  psi <- plogis(beta0 + beta1*x)
  likelihood <- psi^y * (1-psi)^(1-y) # same as next line:
  #   likelihood <- dbinom(y, 1, psi)
  return(-sum(log(likelihood)))
}

# Look at (negative) log-likelihood for 2 parameter sets
negLogLike(c(0,0), y=z, x=vegHt)
negLogLike(c(-3,2), y=z, x=vegHt) # Lower is better!

# Let's minimize it formally by function minimisation
starting.values <- c(beta0=0, beta1=0)
opt.out <- optim(starting.values, negLogLike, y=z, x=vegHt, hessian=TRUE)
(mles <- opt.out$par)     # MLEs are pretty close to truth

# Alternative 1: Brute-force grid search for MLEs
mat <- as.matrix(expand.grid(seq(-10,10,0.1), seq(-10,10,0.1)))
# above: Can vary resolution
nll <- array(NA, dim = nrow(mat))
for (i in 1:nrow(mat)){
  nll[i] <- negLogLike(mat[i,], y = z, x = vegHt)
}
which(nll == min(nll))
mat[which(nll == min(nll)),]

# Produce a likelihood surface, shown in Fig. 2-2.
library(raster)
r <- rasterFromXYZ(data.frame(x = mat[,1], y = mat[,2], z = nll))
mapPalette <- colorRampPalette(rev(c("grey", "yellow", "red")))
plot(r, col = mapPalette(100), main = "Negative log-likelihood",
     xlab = "Intercept (beta0)", ylab = "Slope (beta1)")
contour(r, add = TRUE, levels = seq(50, 2000, 100))

# Alternative 2: Use canned R function glm as a shortcut
(fm <- glm(z ~ vegHt, family = binomial)$coef)

# Add 3 sets of MLEs into plot
# 1. Add MLE from function minimisation
points(mles[1], mles[2], pch = 1, lwd = 2)
abline(mles[2],0)  # Put a line through the Slope value
lines(c(mles[1],mles[1]),c(-10,10))
# 2. Add MLE from grid search
points(mat[which(nll == min(nll)),1], mat[which(nll == min(nll)),2],
       pch = 1, lwd = 2)
# 3. Add MLE from glm function
points(fm[1], fm[2], pch = 1, lwd = 2)

# Note they are essentially all the same (as they should be)

Vc <- solve(opt.out$hessian)         # Get variance-covariance matrix
ASE <- sqrt(diag(Vc))                # Extract asymptotic SEs
print(ASE)

# Compare to SEs reported by glm() function (output thinned)
summary(glm(z ~ vegHt, family = binomial))

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        ASE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$ASE
mle.table$upper <- mle.table$Est + 1.96*mle.table$ASE
mle.table

# Plot the actual and estimated response curves
plot(vegHt, z, xlab="Vegetation height", ylab="Occurrence probability")
plot(function(x) plogis(beta0 + beta1 * x), 1.1, 3, add=TRUE, lwd=2)
plot(function(x) plogis(mles[1] + mles[2] * x), 1.1, 3, add=TRUE,
     lwd=2, col="blue")
legend(1.1, 0.9, c("Actual", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=2)


# 2.4.5 Bootstrapping
# ------------------------------------------------------------------------
nboot <- 1000   # Obtain 1000 bootstrap samples
boot.out <- matrix(NA, nrow=nboot, ncol=3)
dimnames(boot.out) <- list(NULL, c("beta0", "beta1", "psi.bar"))

for(i in 1:1000){
  # Simulate data
  psi <- plogis(mles[1] + mles[2] * vegHt)
  z <- rbinom(M, 1, psi)
  
  # Fit model
  tmp <- optim(mles, negLogLike, y=z, x=vegHt, hessian=TRUE)$par
  psi.mean <- plogis(tmp[1] + tmp[2] * mean(vegHt))
  boot.out[i,] <- c(tmp, psi.mean)
}

SE.boot <- sqrt(apply(boot.out, 2, var))  # Get bootstrap SE
names(SE.boot) <- c("beta0", "beta1", "psi.bar")

# 95% bootstrapped confidence intervals
apply(boot.out,2,quantile,c(0.025,0.975))

# Boostrap SEs
SE.boot

# Compare these with the ASEs for regression parameters
mle.table


# 2.4.6 Likelihood analysis of hierarchical models
# ------------------------------------------------------------------------

# 2.4.6.1 Discrete random variable
# ------------------------------------------------------------------------
set.seed(2014)
M <- 100                             # number of sites
vegHt <- runif(M, 1, 3)              # uniform from 1 to 3
psi <- plogis(beta0 + beta1 * vegHt) # occupancy probability
z <- rbinom(M, 1, psi)               # realised presence/absence
p <- 0.6                             # detection probability
J <- 3                               # sample each site 3 times
y <-rbinom(M, J, p*z)                # observed detection frequency

# Define negative log-likelihood.
negLogLikeocc <- function(beta, y, x, J) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  p <- plogis(beta[3])
  psi <- plogis(beta0 + beta1*x)
  marg.likelihood <- dbinom(y, J, p)*psi + ifelse(y==0,1,0)*(1-psi)
  return(-sum(log(marg.likelihood)))
}

starting.values <- c(beta0=0, beta1=0,logitp=0)
(opt.out <- optim(starting.values, negLogLikeocc, y=y, x=vegHt,J=J,
                  hessian=TRUE))

sqrt(diag(solve(opt.out$hessian)))


# 2.4.6.2 A continuous latent variable
# ------------------------------------------------------------------------
marg <- rep(NA, J+1)
for(j in 0:J){
  marg[j] <- integrate(
    function(x){
      dbinom(j, J, plogis(x)) * dnorm(x, mu, sigma)},
    lower=-Inf,upper=Inf)$value
}
}

# nx = encounter frequencies, number inds. encountered 1, 2, ..., 14 times
nx <- c(34, 16, 10, 4, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0)
nind <- sum(nx)      # Number of individuals observed
J <- 14              # Number of sample occasions

# Model Mh likelihood
Mhlik <- function(parms){
  mu <- parms[1]
  sigma <- exp(parms[2])
  # n0 = number of UNobserved individuals: N = nind + n0
  n0 <- exp(parms[3])
  
  # Compute the marginal probabilities for each possible value j=0,..,14
  marg <- rep(NA,J+1)
  for(j in 0:J){
    marg[j+1] <- integrate(
      function(x){dbinom(j, J, plogis(x)) * dnorm(x, mu, sigma)},
      lower=-Inf,upper=Inf)$value
  }
  
  # The negative log likelihood involves combinatorial terms computed
  # using lgamma()
  -1*(lgamma(n0 + nind + 1) - lgamma(n0 + 1) + sum(c(n0, nx) * log(marg)))
}
(tmp <- nlm(Mhlik, c(-1, 0, log(10)), hessian=TRUE))


(SE <- sqrt( (exp(tmp$estimate[3])^2)* diag(solve(tmp$hessian))[3] ) )


# 2.4.7 The R package ???unmarked???
# ------------------------------------------------------------------------

# 2.5 Bayesian Inference
# ------------------------------------------------------------------------

# 2.5.1 Bayes??? rule
# ------------------------------------------------------------------------

# 2.5.2 Principles of Bayesian inference
# ------------------------------------------------------------------------

# 2.5.3 Prior distributions
# ------------------------------------------------------------------------

# 2.5.4 Computing posterior distributions
# ------------------------------------------------------------------------

# 2.6 Basic Markov chain Monte Carlo (MCMC)
# ------------------------------------------------------------------------

# 2.6.1 Metropolis-Hastings algorithm
# ------------------------------------------------------------------------

# 2.6.2 Illustration: using MH for a binomial model
# ------------------------------------------------------------------------
# Simulate data
set.seed(2016)
y <- rbinom(2, size=10, p = 0.5)

# Define the joint distribution (= likelihood) which we will maximize
jointdis <- function(data, K, p){
  prod(dbinom(data, size=K, p=p))
}

# Posterior is proportional to likelihood times prior
posterior <- function(p, data, K, a, b){
  prod(dbinom(data, size=K, p=p)) * dbeta(p, a, b)
}

# Do 100,000 MCMC iterations using Metropolis algorithm
# Assume vague prior which is beta(1,1) = Unif(0,1)
mcmc.iters <- 100000
out <- rep(NA, mcmc.iters)

# Starting value
p <- 0.2

# Begin the MCMC loop
for(i in 1:mcmc.iters){
  
  # Use a uniform candidate generator (not efficient)
  p.cand <- runif(1, 0, 1)
  
  # Alternative: random walk proposal
  # p.cand <- rnorm(1, p, .05)  # Need to reject if > 1 or < 0
  # if(p.cand < 0 | p.cand > 1 ) next
  
  r <- posterior(p=p.cand, y, K=10, a=1, b=1) / posterior(p=p, y, K=10, a=1, b=1)
  # Generate a uniform r.v. and compare with "r", this imposes the
  #    correct probability of acceptance
  if(runif(1) < r)
    p <- p.cand
  
  # Save the current value of p
  out[i] <- p
}


mean(out)

sd(out)

quantile(out, c(0.025, 0.975))

# Evaluate likelihood for a grid of values of p
p.grid <- seq(0.1, 0.9, , 200)
likelihood <- rep(NA, 200)

for(i in 1:200){
  likelihood[i] <- jointdis(y, K=10, p=p.grid[i])
}

par(mfrow=c(2,1), mar = c(5,5,3,2))
plot(p.grid, likelihood, xlab="", ylab="Likelihood", xlim=c(0,1), ty = "l", main = "Likelihood function")
p.hat <- p.grid[likelihood == max(likelihood)]
abline(v = p.hat)
text(p.hat, 0.005, paste("MLE = ", round(p.hat, 3), sep= ""))

plot(density(out), xlim=c(0,1), main = "Posterior distribution", xlab = "p", ylab = "Posterior")
p.mean <- mean(out)
abline(v = p.mean)
text(p.mean, 0.5, paste("Post. mean = ", round(p.mean, 3),sep=" "))


# 2.6.3 Metropolis algorithm for multi-parameter models
# ----------------------------------------------------------------------------------------
log.posterior <- function(beta0, beta1, z, vegHt){
  # Note: "z" and "vegHt" must be input
  loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)
  logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)
  return(loglike + logprior[1] + logprior[2])
}

niter <- 50000
out <- matrix(NA, niter, 2, dimnames = list(NULL, c("beta0", "beta1")))

# Initialize parameters
beta0 <- rnorm(1)
beta1 <- rnorm(1)

# Current value of the log(posterior)
logpost.curr <- log.posterior(beta0, beta1, z, vegHt)

# Run MCMC algorithm
for(i in 1:niter){
  if(i %% 1000 == 0)                     # report progress
    cat("iter", i, "\n")
  # Update intercept (beta0)
  # Propose candidate values of beta
  # If the proposal was not symmetric, would be Metrop-*Hastings*
  beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
  # Evaluate the log(posterior)
  logpost.cand <- log.posterior(beta0.cand, beta1, z, vegHt)
  # Compute Metropolis acceptance probability, r
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta0 <- beta0.cand
    logpost.curr <- logpost.cand
  }
  
  # Update slope (beta1)
  beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
  # Evaluate the log(posterior)
  logpost.cand <- log.posterior(beta0, beta1.cand, z, vegHt)
  
  # Compute Metropolis acceptance probability
  r <- exp(logpost.cand - logpost.curr)
  # Keep candidate if it meets criterion (u < r)
  if(runif(1) < r){
    beta1 <- beta1.cand
    logpost.curr <- logpost.cand
  }
  out[i,] <- c(beta0, beta1) # Save samples for iteration
}

# Plot
layout(rbind(c(1,1),
             c(2,2),
             c(3,4)), respect = T) # <- play with these settings

par(oma=c(0,0,0,0),mar=c(5,4,1,1))
plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
plot(density(out[,1]), xlab="beta0", main="")
plot(density(out[,2]), xlab="beta1", main="")


# 2.6.4 Why we need to use MCMC for logistic regression
# ------------------------------------------------------------------------

# 2.6.5 Gibbs sampling
# ------------------------------------------------------------------------

# 2.6.6 Convergence and mixing of Markov chains
# ------------------------------------------------------------------------

# 2.6.6.1 Slow mixing and thinning of Markov chains
# ------------------------------------------------------------------------

# 2.6.6.2 Effective sample size and Monte Carlo error
# ------------------------------------------------------------------------

# 2.6.7 Bayesian analysis of hierarchical models
# ------------------------------------------------------------------------

# The occupancy model
# ------------------------------------------------------------------------

# Simulate the data set
set.seed(2014)
M <- 100                        # number of sites
vegHt <- runif(M, 1, 3)         # uniform from 1 to 3
psi <- plogis(-3 + 2*vegHt)     # occupancy probability
z <- rbinom(M, 1, psi)          # realised presence/absence
p <- 0.6                        # detection probability
J <- 3                          # sample each site 3 times
y <-rbinom(M, J, p*z)           # observed detection frequency

# Number of MCMC iterations to to
niter <- 50000

# Matrix to hold the simulated values
out <- matrix(NA, niter, 3, dimnames = list(NULL, c("beta0", "beta1", "p")))

# Initialize parameters, likelihood, and priors
starting.values <- c(beta0=0, beta1=0)
beta0 <- starting.values[1]
beta1 <- starting.values[2]
z <-ifelse(y>0, 1, 0)
p <- 0.2

# NOTE: using logistic reg. likelihood function here!
loglike <- -1*negLogLike(c(beta0, beta1), z, vegHt)
logprior <- dnorm(c(beta0, beta1), 0, 10, log=TRUE)

# Run MCMC algorithm
for(i in 1:niter) {
  if(i %% 1000 ==0 ) # report progress
    cat("iter", i, "\n")
  
  # PART 1 of algorithm -- same as before
  # Update intercept (beta0)
  # propose candidate values of beta
  beta0.cand <- rnorm(1, beta0, 0.3) # 0.3 is tuning parameter
  # evaluate likelihood and priors for candidates
  loglike.cand <- -1*negLogLike(c(beta0.cand, beta1), z, vegHt)
  logprior.cand <- dnorm(beta0.cand, 0, 10, log=TRUE)
  # Compute Metropolis acceptance probability (r)
  r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[1]))
  # Keep candidate if it meets the criterion
  if(runif(1) < r){
    beta0 <- beta0.cand
    loglike <- loglike.cand
    logprior[1] <- logprior.cand
  }
  # Update slope (beta1)
  beta1.cand <- rnorm(1, beta1, 0.3) # 0.3 is tuning parameter
  # evaluate likelihood and priors for candidates
  loglike.cand <- -1*negLogLike(c(beta0,beta1.cand), z, vegHt)
  logprior.cand <- dnorm(beta1.cand, 0, 10, log=TRUE)
  # Compute Metropolis acceptance probability r
  r <- exp((loglike.cand+logprior.cand) - (loglike + logprior[2]))
  # Keep the candidates if they meet the criterion
  if(runif(1) < r) {
    beta1 <- beta1.cand
    loglike <- loglike.cand
    logprior[2] <- logprior.cand
  }
  
  # Part 2 of the algorithm
  # update z. Note we only need to update z if y=0.
  # The full conditional has known form
  
  psi <- plogis(beta0 + beta1 * vegHt)
  psi.cond <- dbinom(0,J,p) * psi /(dbinom(0, J, p) * psi + (1-psi))
  z[y==0] <- rbinom(sum(y==0), 1, psi.cond[y==0])
  loglike <- -1 * negLogLike(c(beta0, beta1), z, vegHt)
  
  # Part 3: update p
  ## The commented code will update p using Metropolis
  ## loglike.p <- sum(log(dbinom(y[z==1],J,p)))
  ## p.cand <- runif(1, 0, 1)
  ## loglike.p.cand <- sum(log(dbinom(y[z==1], J, p.cand)))
  ## if(runif(1) < exp(loglike.p.cand-loglike.p))
  ##    p<-p.cand
  ## This bit draws p directly from its full conditional
  p <- rbeta(1, 1+ sum(y), sum(z)*J +1 - sum(y) )
  
  # Save MCMC samples
  out[i,] <- c(beta0,beta1,p)
}

# Plot bivariate representation of joint posterior
pairs(out)

# Trace/history plots for each parameter (Fig. 2.8)
op <- par(mfrow=c(2,1))
plot(out[,1], type="l", xlab="Iteration", ylab="beta0")
abline(h = mean(out[,1]), col = "blue", lwd = 2)
abline(h = -3, col = "red", lwd = 2)
plot(out[,2], type="l", xlab="Iteration", ylab="beta1")
abline(h = mean(out[,2]), col = "blue", lwd = 2)
abline(h = 2, col = "red", lwd = 2)



# 2.7 Model selection and averaging
# ------------------------------------------------------------------------

# 2.7.1 Model selection by AIC
# ------------------------------------------------------------------------

# 2.7.2 Model selection by DIC
# ------------------------------------------------------------------------

# 2.7.3 Bayesian model averaging with indicator variables
# ------------------------------------------------------------------------



# 2.8 Assessment of model fit
# ------------------------------------------------------------------------


# 2.8.1 Parametric bootstrapping example
# ------------------------------------------------------------------------

sim.data <- function(beta0 = -3, beta1 = 2, p = 0.6, x=NULL){
  # Function allows input of covariate "x", or simulates new
  
  M <- 100
  if(is.null(x))
    vegHt <- runif(M, 1, 3) # uniform from 1 to 3
  
  # Suppose that occupancy probability increases with vegHt
  # The relationship is described (default) by an intercept of -3 and
  #    a slope parameter of 2 on the logit scale
  # plogis is the inverse-logit (constrains us back to the [0-1] scale)
  psi <- plogis(beta0 + beta1*vegHt)
  
  # Now we simulated true presence/absence for 100 sites
  z <- rbinom(M, 1, psi)
  
  # Now generate observations
  J <- 3 # sample each site 3 times
  y <- rbinom(M,J,p*z)
  
  list(y=y, J=J, vegHt=vegHt)
}

# This is the negative log-likelihood based on the marginal distribution
# of y. It is the pmf of a zero-inflated binomial random variable.
#
negLogLikeocc <- function(beta, y, x, J) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  p<- plogis(beta[3])
  psi <- plogis(beta0 + beta1*x)
  marg.likelihood <- dbinom(y, J, p) * psi + ifelse(y==0, 1, 0) * (1-psi)
  return(-sum(log(marg.likelihood)))
}

data <- sim.data()        # Generate a data set

# Let's minimize it
starting.values <- c(beta0=0, beta1=0, logitp=0)
opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x=data$vegHt,J=data$J, hessian=TRUE)
(mles <- opt.out$par)

# Make a table with estimates, SEs, and 95% CI
mle.table <- data.frame(Est=mles,
                        SE = sqrt(diag(solve(opt.out$hessian))))
mle.table$lower <- mle.table$Est - 1.96*mle.table$SE
mle.table$upper <- mle.table$Est + 1.96*mle.table$SE
mle.table


# Define a fit statistic
fitstat <- function(y, Ey){
  sum((sqrt(y) - sqrt(Ey)))
}
# Compute it for the observed data
T.obs <- fitstat(y, J*plogis(mles[1] + mles[2]*vegHt)*plogis(mles[3]))

# Get bootstrap distribution of fit statistic
T.boot <- rep(NA, 100)
for(i in 1:100){
  # Simulate a new data set and extract the elements. Note we use
  # the previously simulated "vegHt" covariate
  data <- sim.data(beta0=mles[1],beta1=mles[2],p=plogis(mles[3]),x=vegHt)
  # Next we fit the model
  starting.values <- c(0,0,0)
  opt.out <- optim(starting.values, negLogLikeocc, y=data$y, x= data$vegHt, J=data$J, hessian=TRUE)
  (parms <- opt.out$par)
  # Obtain the fit statistic
  T.boot[i]<- fitstat(y, J*plogis(parms[1] + parms[2]*vegHt)*plogis(parms[3]) )
}

(T.obs)

summary(T.boot)


# 2.8.2 Bayesian p-value
# ------------------------------------------------------------------------

2.9 Summary and outlook
# ------------------------------------------------------------------------









# ============================================================================
#
# 3. Linear models, generalised linear models (GLMs) and random effects models:
#    the components of hierarchical models
#
# ============================================================================







# 3.1 Introduction
# ------------------------------------------------------------------------


# Define data
pop <- factor(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)), levels = c("Navarra", "Aragon", "Catalonia"))         # Population
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9) # Wing span
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2) # Body length
sex <- factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F"))
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6)      # Number of ectoparasites
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57) # Color intensity
damage <- c(0,2,0,0,4,2,1,0,1)                 # Number of wings damaged

cbind(pop, sex, wing, body, mites, color, damage) # Look at data

str(pop)

par(mfrow = c(1, 3), cex = 1.2)
colorM <- c("red", "red", "blue", "green", "green")  # Pop color code males
colorF <- c("red", "blue", "blue", "green", "green") # Pop color code females
plot(body[sex == "M"], wing[sex == "M"], col =colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col =colorF, pch = 16)
text(6.8, 13.8, "A", cex = 1.5)
plot(body[sex == "M"], mites[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 10), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Parasite load")
points(body[sex == "F"], mites[sex == "F"], col = colorF, pch = 16)
text(6.8, 9.7, "B", cex = 1.5)
plot(body[sex == "M"], damage[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(0, 4), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Damaged wings")
points(body[sex == "F"], damage[sex == "F"], col = colorF, pch = 16)
text(6.8, 3.9, "C", cex = 1.5)



# 3.2 Linear models
# ------------------------------------------------------------------------


# 3.2.1 Linear models with main effects of one factor and one continuous covariate
# ------------------------------------------------------------------------
summary(fm1 <- lm(wing ~ pop + body))

summary(fm2 <- lm(wing ~ pop-1 + body))

cbind(model.matrix(~pop+body) %*% fm1$coef, predict(fm1)) # Compare two solutions

model.matrix(~ pop + body) # Effects parameterisation

model.matrix(~ pop-1 + body) # Means parameterization

par(mfrow = c(1, 3), mar = c(5,4,2,2), cex = 1.2, cex.main = 1)
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm2)[1], coef(fm2)[4], col = "red", lwd = 2)
abline(coef(fm2)[2], coef(fm2)[4], col = "blue", lwd = 2)
abline(coef(fm2)[3], coef(fm2)[4], col = "green", lwd = 2)
text(6.8, 14, "A", cex = 1.5)


# 3.2.2 Linear models with interaction between one factor and one continuous covariate
# ------------------------------------------------------------------------

model.matrix(~ pop*body)  # Effects parameterisation

model.matrix(~ pop*body-1-body)  # Means parameterisation

summary(fm3 <- lm(wing ~ pop*body-1-body))

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm3)[1], coef(fm3)[4], col = "red", lwd = 2)
abline(coef(fm3)[2], coef(fm3)[5], col = "blue", lwd = 2)
abline(coef(fm3)[3], coef(fm3)[6], col = "green", lwd = 2)
text(6.8, 14, "B", cex = 1.5)

# Create new design matrix
(DM0 <- model.matrix(~ pop*body-1-body)) # Original DM for means param
DM0[7:9,5] <- DM0[7:9,6]                 # Combine slopes for Ar and Cat
(DM1 <- DM0[, -6])                       # Delete former slope column for Cat

# Fit model with partial interaction
summary(fm4 <- lm(wing ~ DM1-1))

# Do significance test
anova(fm3, fm4)             # F test between two models

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5), ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17, xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm4)[1], coef(fm4)[4], col = "red", lwd = 2)
abline(coef(fm4)[2], coef(fm4)[5], col = "blue", lwd = 2)
abline(coef(fm4)[3], coef(fm4)[5], col = "green", lwd = 2)
text(6.8, 14, "C", cex = 1.5)


# 3.2.3 Linear models with two factors
# ------------------------------------------------------------------------
model.matrix(~ pop+sex)  # Design matrix of main-effects 2-way ANOVA

# Fit linear model with that design matrix
summary(fm5 <- lm(wing ~ pop + sex))

model.matrix(~ pop+sex-1)  # Design matrix of the main-effects 2-way ANOVA

# Fit linear model with that design matrix
summary(fm6 <- lm(wing ~ pop + sex-1))

# Variant 1: Effects parameterisation (R default)
model.matrix(~ pop*sex)
#model.matrix(~ pop + sex + pop:sex)     # Same

# Variant 2: Means param. for pop, effects param. for sex
model.matrix(~ pop*sex-1)

# Variant 3 (output slightly trimmed): full means parameterisation
model.matrix(~ pop:sex-1)


# 3.2.4 Linear models with two continuous covariates and including polynomials
# ------------------------------------------------------------------------
model.matrix(~ body + color)  # main-effects of covariates

summary(fm7 <- lm(wing ~ body + color))  # Fit that model

model.matrix(~ body*color)  # Interaction between two covariates

summary(fm8 <- lm(wing ~ body*color))  # Fit that model

# Cubic polynomial of body in R
body2 <- body^2           # Squared body length
body3 <- body^3           # Cubed body length
model.matrix(~ body + body2 + body3)

summary(fm9 <- lm(wing ~ body + body2 + body3))  # Fit that model
# summary(fm9 <- lm(wing ~ body + I(body^2) + I(body^3))) # same



# 3.3 Generalised linear models (GLMs)
# ------------------------------------------------------------------------


# 3.3.1 Poisson generalised linear model (GLM) for unbounded counts
# ------------------------------------------------------------------------
summary(fm10 <- glm(mites ~ pop-1 + body, family = poisson))


# 3.3.2 Offsets in the Poisson GLM
# ------------------------------------------------------------------------
summary(fm10a <- glm(mites ~ pop-1 + wing, offset = log(body), family = poisson))
# summary(fm10a <- glm(mites ~ offset(log(body)) + pop-1 + wing, family = poisson))     # same


# 3.3.3 Overdispersion and underdispersion
# ------------------------------------------------------------------------

# 3.3.4 Zero-inflation
# ------------------------------------------------------------------------

# 3.3.5 Bernoulli GLM: logistic regression for a binary response
# ------------------------------------------------------------------------
presence <- ifelse(mites > 0, 1, 0)  # convert abundance to presence/absence
summary(fm11 <- glm(presence ~ pop-1 + body, family = binomial))

# 3.3.6 Modeling a Poisson process from presence/absence data using a Bernoulli GLM with cloglog link
# ------------------------------------------------------------------------
summary(fm11a <- glm(presence ~ pop-1 + body, family = binomial(link = "cloglog")))

summary(fm10)

# 3.3.7 Binomial GLM: logistic regression for bounded counts
# -------------------------------------------------------------------------
summary(fm12 <- glm(cbind(damage, 4-damage) ~ pop + body -1, family = binomial))


# 3.3.8 The GLM as the quintessential statistical model
# ------------------------------------------------------------------------



# 3.4 Random effects (mixed) models
# ------------------------------------------------------------------------


# 3.4.1 Random effects for a normal data distribution: normal-normal generalised linear mixed model (GLMM)
# ------------------------------------------------------------------------
# Plot data without distinguishing sex
plot(body, wing, col = rep(c("red", "blue", "green"), each = 3), xlim = c(6.5, 9.5), ylim = c(10, 14), cex = 1.5, lwd = 2, frame.plot = FALSE, las = 1, pch = 16, xlab = "Body length", ylab = "Wing span")

summary(lm <- lm(wing ~ pop-1 + body))     # Same as fm2

library(lme4)
summary(lmm1 <- lmer(wing ~ (1|pop) + body))  # Fit the model
ranef(lmm1)                                   # Print random effects

alpha_j <- fixef(lmm1)[1]+ranef(lmm1)$pop[,1]
cbind(fixed = coef(lm)[1:3], random = alpha_j)

par(lwd = 3)
abline(lm$coef[1], lm$coef[4], col = "red", lty = 2)
abline(lm$coef[2], lm$coef[4], col = "blue", lty = 2)
abline(lm$coef[3], lm$coef[4], col = "green", lty = 2)
abline(alpha_j[1], fixef(lmm1)[2], col = "red")
abline(alpha_j[2], fixef(lmm1)[2], col = "blue")
abline(alpha_j[3], fixef(lmm1)[2], col = "green")
abline(fixef(lmm1), col = "black")
legend(6.5, 14, c("Catalonia", "Aragon", "Navarra"), col=c("blue", "green", "red"), lty = 1, pch = 16, bty = "n", cex = 1.5)

summary(lmm2 <- lmer(wing ~ body + (1|pop) + (0+body|pop)))


# 3.4.2 Random effects for a Poisson data distribution: normal-Poisson generalised linear mixed model (GLMM)
# ------------------------------------------------------------------------
summary(glmm <- glmer(mites ~ body + (1|pop), family = poisson))



# 3.5 Summary and outlook
# ------------------------------------------------------------------------


# 3.6 Exercises
# ------------------------------------------------------------------------
# Define and plot data (10 times larger data set than the toy data set)
clone.size <- 10               # clone size
pop <- factor(rep(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)), levels = c("Navarra", "Aragon", "Catalonia"), clone.size))
wing <- rep(c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9), clone.size)
body <- rep(c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2), clone.size)
sex <- rep(factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F")), clone.size)
mites <- rep(c(0, 3, 2, 1, 0, 7, 0, 9, 6), clone.size)
color <- rep(c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57), clone.size)
damage <- rep(c(0,2,0,0,4,2,1,0,1), clone.size)









# =========================================================================
#
# 4. Introduction to data simulation
#
# =========================================================================





# 4.1 What do we mean by data simulation and why is it so tremendously useful ?
# ------------------------------------------------------------------------




# 4.2 Generation of a typical point count data set
# ------------------------------------------------------------------------


# 4.2.1 Initial steps: sample size and covariate values
# ------------------------------------------------------------------------
M <- 267         # Number of spatial replicates (sites)
J <- 3           # Number of temporal replicates (counts)

set.seed(24)     # Can choose seed of your choice

elev <- runif(n = M, -1, 1)             # Scaled elevation of a site
forest <- runif(n = M, -1, 1)           # Scaled forest cover at each site
wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed


# 4.2.2 Simulating the ecological process and its outcome: great tit abundance
# ------------------------------------------------------------------------
mean.lambda <- 2               # Mean expected abundance of great tits
beta0 <- log(mean.lambda)      # Same on log scale (= log-scale intercept)
beta1 <- -2                    # Effect (slope) of elevation
beta2 <- 2                     # Effect (slope) of forest cover
beta3 <- 1                     # Interaction effect (slope) of elev and forest

log.lambda <- beta0 + beta1 * elev + beta2 * forest + beta3 * elev * forest
lambda <- exp(log.lambda)      # Inverse link transformation

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
curve(exp(beta0 + beta1*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 18), xlab = "Elevation", ylab = "lambda", lwd = 2)
text(-0.9, 17, "A", cex = 1.5)
plot(elev, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Elevation", ylab = "")
text(-0.9, 36, "B", cex = 1.5)
curve(exp(beta0 + beta2*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 18), xlab = "Forest cover", ylab = "lambda", lwd = 2)
text(-0.9, 17, "C", cex = 1.5)
plot(forest, lambda, frame.plot = FALSE, ylim = c(0, 38), xlab = "Forest cover", ylab = "")
text(-0.9, 36, "D", cex = 1.5)


# Compute expected abundance for a grid of elevation and forest cover
cov1 <- seq(-1, 1,,100)                       # Values for elevation
cov2 <- seq(-1,1,,100)                        # Values for forest cover
lambda.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix, for every combination of values of elevation and forest cover
for(i in 1:100){
  for(j in 1:100){
    lambda.matrix[i, j] <- exp(beta0 + beta1 * cov1[i] + beta2 * cov2[j] + beta3 * cov1[i] * cov2[j])
  }
}

par(mfrow = c(1, 2), mar = c(5,4,3,2), cex.main = 1.6)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x = cov1, y = cov2, z = lambda.matrix, col = mapPalette(100), xlab = "Elevation", ylab = "Forest cover", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = lambda.matrix, add = TRUE, lwd = 1)
title(main = "A")
matpoints(elev, forest, pch="+", cex=0.8)


N <- rpois(n = M, lambda = lambda)   # Realised abundance
sum(N)                               # Total population size at M sites
table(N)                             # Frequency distribution of tit abundance

# 4.2.3 Simulating the observation process and its outcome: point counts of great tits
# ------------------------------------------------------------------------
mean.detection <- 0.3            # Mean expected detection
alpha0 <- qlogis(mean.detection) # same on logit scale (intercept)
alpha1 <- 1                      # Effect (slope) of elevation
alpha2 <- -3                     # Effect (slope) of wind speed
alpha3 <- 0                      # Interaction effect (slope) of elev and wind

logit.p <- alpha0 + alpha1 * elev + alpha2 * wind + alpha3 * elev * wind
p <- plogis(logit.p)             # Inverse link transform
mean(p)                          # average per-individual p is 0.38

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1.1), xlab = "Elevation", ylab = "p", lwd = 2)
text(-0.9, 1.05, "A", cex = 1.5)
matplot(elev, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1), xlab = "Elevation", ylab = "")
text(-0.9, 1.05, "B", cex = 1.5)
curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", frame.plot = FALSE, ylim = c(0, 1.1), xlab = "Wind speed", ylab = "p", lwd = 2)
text(-0.9, 1.05, "C", cex = 1.5)
matplot(wind, p, pch = "*", frame.plot = FALSE, ylim = c(0, 1.1), xlab = "Wind speed", ylab = "p")
text(-0.9, 1.05, "D", cex = 1.5)


# Compute expected detection probability for a grid of elevation and wind speed
cov1 <- seq(-1, 1,,100)                  # Values of elevation
cov2 <- seq(-1,1,,100)                   # Values of wind speed
p.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix which combines every value in cov 1 with every other in cov2
for(i in 1:100){
  for(j in 1:100){
    p.matrix[i, j] <- plogis(alpha0 + alpha1 * cov1[i] + alpha2 * cov2[j] + alpha3 * cov1[i] * cov2[j])
  }
}
image(x = cov1, y = cov2, z = p.matrix, col = mapPalette(100), xlab = "Elevation", ylab = "Wind speed", cex.lab = 1.2)
contour(x = cov1, y = cov2, z = p.matrix, add = TRUE, lwd = 1)
title(main = "B")
matpoints(elev, wind, pch="+", cex=0.7, col = "black")

C <- matrix(NA, nrow = M, ncol = J)      # Prepare array for counts
for (i in 1:J){                          # Generate counts
  C[,i] <- rbinom(n = M, size = N, prob = p[,i])
}

head(cbind("True N" = N, "1st count" = C[,1], "2nd count" = C[,2], "3rd count" = C[,3]), 10)                    # First 10 rows (= sites)

table(C)

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
matplot(elev, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38), xlab = "Elevation", ylab = "Count (C)")
text(-0.9, 36, "A", cex = 1.5)
matplot(forest, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38), xlab = "Forest cover", ylab = "Count (C)")
text(-0.9, 36, "B", cex = 1.5)
matplot(wind, C, pch = "*", frame.plot = FALSE, ylim = c(0, 38), xlab = "Wind speed", ylab = "Count (C)")
text(-0.9, 36, "C", cex = 1.5)
hist(C, breaks = 50, col = "grey", ylim = c(0, 460), main = "", xlab = "Count (C)")
text(3, 450, "D", cex = 1.5)

sum(N)                   # True total abundance (all sites)
sum(apply(C, 1, max))    # 'Observed' total abundance (all sites)

sum(N>0)                 # True number of occupied sites
sum(apply(C, 1, max)>0)  # 'Observed' number of occupied sites



# 4.3 Packaging everything in a function
# ------------------------------------------------------------------------
# Function definition with set of default values
data.fn <- function(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0, show.plot = TRUE){
  #
  # Function to simulate point counts replicated at M sites during J occasions.
  # Population closure is assumed for each site.
  # Expected abundance may be affected by elevation (elev),
  # forest cover (forest) and their interaction.
  # Expected detection probability may be affected by elevation,
  # wind speed (wind) and their interaction.
  # Function arguments:
  #     M: Number of spatial replicates (sites)
  #     J: Number of temporal replicates (occasions)
  #     mean.lambda: Mean abundance at value 0 of abundance covariates
  #     beta1: Main effect of elevation on abundance
  #     beta2: Main effect of forest cover on abundance
  #     beta3: Interaction effect on abundance of elevation and forest cover
  #     mean.detection: Mean detection prob. at value 0 of detection covariates
  #     alpha1: Main effect of elevation on detection probability
  #     alpha2: Main effect of wind speed on detection probability
  #     alpha3: Interaction effect on detection of elevation and wind speed
  #     show.plot: if TRUE, plots of the data will be displayed;
  #        set to FALSE if you are running simulations.
  
  # Create covariates
  elev <- runif(n = M, -1, 1)                         # Scaled elevation
  forest <- runif(n = M, -1, 1)                       # Scaled forest cover
  wind <- array(runif(n = M*J, -1, 1), dim = c(M, J)) # Scaled wind speed
  
  # Model for abundance
  beta0 <- log(mean.lambda)               # Mean abundance on link scale
  lambda <- exp(beta0 + beta1*elev + beta2*forest + beta3*elev*forest)
  N <- rpois(n = M, lambda = lambda)      # Realised abundance
  Ntotal <- sum(N)                        # Total abundance (all sites)
  psi.true <- mean(N>0)                   # True occupancy in sample
  
  # Plots
  if(show.plot){
    par(mfrow = c(2, 2), cex.main = 1)
    devAskNewPage(ask = TRUE)
    curve(exp(beta0 + beta1*x), -1, 1, col = "red", main = "Relationship lambda-elevation \nat average forest cover", frame.plot = F, xlab = "Scaled elevation")
    plot(elev, lambda, xlab = "Scaled elevation", main = "Relationship lambda-elevation \nat observed forest cover", frame.plot = F)
    curve(exp(beta0 + beta2*x), -1, 1, col = "red", main = "Relationship lambda-forest \ncover at average elevation", xlab = "Scaled forest cover", frame.plot = F)
    plot(forest, lambda, xlab = "Scaled forest cover", main = "Relationship lambda-forest cover \nat observed elevation", frame.plot = F)
  }
  
  # Model for observations
  alpha0 <- qlogis(mean.detection)        # mean detection on link scale
  p <- plogis(alpha0 + alpha1*elev + alpha2*wind + alpha3*elev*wind)
  C <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C[,i] <- rbinom(n = M, size = N, prob = p[,i])
  }
  summaxC <- sum(apply(C,1,max))          # Sum of max counts (all sites)
  psi.obs <- mean(apply(C,1,max)>0)       # Observed occupancy in sample
  
  # More plots
  if(show.plot){
    par(mfrow = c(2, 2))
    curve(plogis(alpha0 + alpha1*x), -1, 1, col = "red", main = "Relationship p-elevation \nat average wind speed", xlab = "Scaled elevation", frame.plot = F)
    matplot(elev, p, xlab = "Scaled elevation", main = "Relationship p-elevation\n at observed wind speed", pch = "*", frame.plot = F)
    curve(plogis(alpha0 + alpha2*x), -1, 1, col = "red", main = "Relationship p-wind speed \n at average elevation", xlab = "Scaled wind speed", frame.plot = F)
    matplot(wind, p, xlab = "Scaled wind speed", main = "Relationship p-wind speed \nat observed elevation", pch = "*", frame.plot = F)
    
    matplot(elev, C, xlab = "Scaled elevation", main = "Relationship counts and elevation", pch = "*", frame.plot = F)
    matplot(forest, C, xlab = "Scaled forest cover", main = "Relationship counts and forest cover", pch = "*", frame.plot = F)
    matplot(wind, C, xlab = "Scaled wind speed", main = "Relationship counts and wind speed", pch = "*", frame.plot = F)
    desc <- paste('Counts at', M, 'sites during', J, 'surveys')
    hist(C, main = desc, breaks = 50, col = "grey")
  }
  
  # Output
  return(list(M = M, J = J, mean.lambda = mean.lambda, beta0 = beta0, beta1 = beta1, beta2 = beta2, beta3 = beta3, mean.detection = mean.detection, alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, elev = elev, forest = forest, wind = wind, lambda = lambda, N = N, p = p, C = C, Ntotal = Ntotal, psi.true = psi.true, summaxC = summaxC, psi.obs = psi.obs))
}


data.fn()                  # Execute function with default arguments
data.fn(show.plot = FALSE) # same, without plots
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, alpha1 = 1, alpha2 = -3, alpha3 = 0) # Explicit defaults
data <- data.fn()          # Assign results to an object called 'data'

simrep <- 10000
NTOTAL <- SUMMAXC <- numeric(simrep)
for(i in 1:simrep){
  data <- data.fn(show.plot = FALSE)
  NTOTAL[i] <- data$Ntotal
  SUMMAXC[i] <- data$summaxC
}
plot(sort(NTOTAL), ylim = c(min(SUMMAXC), max(NTOTAL)), ylab = "", xlab = "Simulation", col = "red", frame = FALSE)
points(SUMMAXC[order(NTOTAL)], col = "blue")

data.fn(J = 2)              # Only 2 surveys
data.fn(J = 1)              # No temporal replicate
data.fn(M = 1, J = 100)     # No spatial replicates, but 100 counts
data.fn(alpha3 = 1)          # With interaction elev-wind on p ## see errata
data.fn(M = 267, J = 3, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 1)         # No obs. process (i.e., p = 1, perfect detection)
data.fn(mean.lambda = 50)   # Really common species
data.fn(mean.lambda = 0.05) # Really rare species

set.seed(24)
data <- data.fn()    # Default arguments
str(data)            # Look at the object

attach(data)         # Make objects inside of 'data' accessible directly

detach(data)         # Make clean up



# 4.4 Summary and outlook
# ------------------------------------------------------------------------



# 4.5 Exercises
# ------------------------------------------------------------------------


# Generate data set and fit Poisson and Bernoulli models
data <- data.fn(M = 267, J = 1, mean.lambda = 2, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 1, show.plot = FALSE)
summary(fmPois <- glm(C ~ elev*forest, family = poisson, data = data))
presence <- ifelse(data$C > 0, 1, 0)
summary(fmBern <- glm(presence ~ elev*forest, family = binomial(link = "cloglog"), data = data))

# Compare Poisson and Bernoulli estimates with truth
print(cbind(Truth = c(data$beta0, data$beta1, data$beta2, data$beta3), summary(fmPois)$coef[,1:2], summary(fmBern)$coef[,1:2]),3)










# =========================================================================
#
# 5. Fitting models using the Bayesian modeling software BUGS and JAGS
#
# =========================================================================





# 5.1 Introduction
# ------------------------------------------------------------------------



# 5.2 Introduction to BUGS software: WinBUGS, OpenBUGS, and JAGS
# ------------------------------------------------------------------------



# 5.3 Linear model with normal response (normal GLM): multiple linear regression
# ------------------------------------------------------------------------


# Generate data with data.fn from chapter 4
set.seed(24)
data <- data.fn()
str(data)
attach(data)

# Summarize data by taking mean at each site and plot
Cmean <- apply(C, 1, mean)
par(mfrow = c(1,3))
hist(Cmean, 50)               # Very skewed
plot(elev, Cmean)
plot(forest, Cmean)

# Package the data needed in a bundle
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, forest = forest)
str(win.data)                    # Check what???s in win.data


# Write text file with model description in BUGS language
cat(file = "multiple_linear_regression_model.txt",
    "   # --- Code in BUGS language starts with this quotation mark ---
    model {
    
    # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)                   # Precision tau = 1/(sd^2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # dispersion tau is precision (1/variance)
    mu[i] <- alpha0 + alpha1*elev[i] + alpha2*forest[i] + alpha3*elev[i]*forest[i]
    }
    
    # Derived quantities
    for (i in 1:M){
    resi[i] <- Cmean[i] - mu[i]
    }
    }"#  --- Code in BUGS language ends on this line ---
)


# Initial values (have to give for at least some estimands)
inits <- function() list(alpha0 = rnorm(1,0,10), alpha1 = rnorm(1,0,10), alpha2 = rnorm(1,0,10), alpha3 = rnorm(1,0,10))


# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "resi")


# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
# ni <- 10   ;   nt <- 1   ;   nb <- 0   ;  nc <- 8 # not run


# Call WinBUGS from R (approximate run time (ART) <1 min)
library(R2WinBUGS)
bugs.dir <- "C:/WinBUGS14/"          # Place where your WinBUGS installed
out1B <- bugs(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())


# Overview of the object created by bugs
names(out1B)
str(out1B, 1)


# Call OpenBUGS from R (ART <1 min)
library(R2OpenBUGS)
out1OB <- bugs(data=win.data, inits=inits, parameters.to.save = params, model.file = "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, working.directory = getwd())
detach("package:R2OpenBUGS", unload=T)  # Otherwise R2WinBUGS is 'masked'


# Call JAGS from R (ART <1 min)
library(jagsUI)
?jags                 # Look at main function
out1J <- jags(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

out1J <- jags(win.data, inits, params, parallel = TRUE, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


par(mfrow = c(3,2))
traceplot(out1J, param = c('alpha1', 'alpha2', 'resi[c(1,3, 5:6)]')) # Subset
# traceplot(out1J)                  # All params

# Overview of object created by jags()
names(out1J)

# Summarize posteriors from WinBUGS run
print(out1B, 2)

# Summarize posteriors from JAGS run
print(out1J, 2)

n <- 10                   # maximum lag
par(mfrow = c(2, 3), mar = c(5,4,2,2), cex.main = 1)
for(k in 1:5){
  matplot(0:n, autocorr.diag(as.mcmc(out1B$sims.array[,,k]), lags = 0:n),
          type = "l", lty = 1, xlab = "lag (n)", ylab = "autocorrelation",
          main = colnames(out1B$sims.matrix)[k], lwd = 2)
  abline(h = 0)
}

plot(out1B)              # For WinBUGS analysis from R2WinBUGS
plot(out1J)              # For JAGS analysis from jagsUI

par(mfrow = c(1, 2), mar = c(5,4,2,2), cex.main = 1)
whiskerplot(out1J, param = c('alpha0', 'alpha1', 'alpha2', 'alpha3', 'sd', 'resi[c(1,3, 5:7)]'))    # For JAGS analysis from jagsUI
library(denstrip)      # Similar, but more beautiful, with package denstrip
plot(out1J$alpha0, xlim=c(-4, 4), ylim=c(1, 5), xlab="", ylab="", type="n", axes = F, main = "Density strip plots")
axis(1)
axis(2, at = 1:5, labels = c('alpha0','alpha1','alpha2','alpha3','sd'), las = 1)
abline(v = c(-4,-2,2,4), col = "grey")  ;  abline(v = 0)
for(k in 1:5){
  denstrip(unlist(out1J$sims.list[k]), at = k, ticks = out1J$summary[k, c(3,5,7)])
}


(fm <- summary(lm(Cmean ~ elev*forest)))


print(cbind(out1B$summary[1:5, 1:2], out1J$summary[1:5, 1:2], rbind(fm$coef[,1:2], c(fm$sigma, NA))), 4)

plot(lm(Cmean ~ elev*forest))

mu <- out1B$mean$alpha0 + out1B$mean$alpha1 * elev + out1B$mean$alpha2 * forest + out1B$mean$alpha3 * elev * forest       # Compute the posterior mean of mu

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
plot(1:M, out1B$summary[6:272, 1], xlab = "Order of values", ylab = "Residual", frame.plot = F, ylim = c(-10, 15))
abline(h = 0, col = "red", lwd = 2)
segments(1:267, out1B$summary[6:272, 3], 1:267, out1B$summary[6:272, 7], col = "grey")
text(10, 14, "A", cex = 1.5)
hist(out1B$summary[6:272, 1], xlab = "Residual", main = "", breaks = 50, col = "grey", xlim = c(-10, 15))
abline(v = 0, col = "red", lwd = 2)
text(-9, 48, "B", cex = 1.5)
qq <- qnorm(seq(0,0.9999,,data$M), mean = 0, sd = out1B$summary[5, 1])
plot(sort(qq), sort(out1B$summary[6:272, 1]), xlab = "Theoretical quantile", ylab = "Residual", frame.plot = F, ylim = c(-10, 15)) # could also use qqnorm()
abline(0, 1, col = "red", lwd = 2)
text(-4.5, 14, "C", cex = 1.5)
plot(mu, out1B$summary[6:272, 1], xlab = "Predicted values", ylab = "Residual", frame.plot = F, ylim = c(-10, 15))
abline(h = 0, col = "red", lwd = 2)
segments(mu, out1B$summary[6:272, 3], mu, out1B$summary[6:272, 7], col = "grey")
text(-1, 14, "D", cex = 1.5)

fm

confint(lm(Cmean ~ elev*forest))

par(mfrow = c(2, 2), mar = c(5,4,2,2), cex.main = 1)
hist(out1B$sims.list$alpha1, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1B$sims.list$alpha1, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(-2.4, 1.8, "A", cex = 1.5)
hist(out1B$sims.list$alpha2, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1B$sims.list$alpha2, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(1.7, 2, "B", cex = 1.5)
hist(out1B$sims.list$alpha3, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1B$sims.list$alpha3, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(-2.2, 1.2, "C", cex = 1.5)
hist(out1B$sims.list$sd, main = "", breaks = 100, col = "grey", freq=F)
abline(v = quantile(out1B$sims.list$sd, prob = c(0.025, 0.975)), col = "red", lwd = 2)
text(1.6, 4.9, "D", cex = 1.5)


HPDinterval(as.mcmc(out1B$sims.list$sd), prob = 0.95)  # HPDI
quantile(out1B$sims.list$sd, prob = c(0.025, 0.975))   # Percentile-based CRI


cbind(confint(lm(Cmean ~ elev*forest))[2:4,], out1B$summary[2:4, c(3,7)])

mean(out1B$sims.list$alpha1 < -1.6)
mean(out1B$sims.list$alpha1 < -1.6 & out1B$sims.list$alpha1 > -1.8)

plot(out1B$sims.list$alpha1, out1B$sims.list$alpha2)
abline(h = c(2.5, 2.8), col = "red", lwd = 2)
abline(v = c(-1.9, -1.6), col = "red", lwd = 2)


mean(out1B$sims.list$alpha1 < -1.6 & out1B$sims.list$alpha1 > -1.9 & out1B$sims.list$alpha2 > 2.5 & out1B$sims.list$alpha2 < 2.8)

crazy.ratio <- out1B$sims.list$alpha2 / abs(out1B$sims.list$alpha1)
hist(crazy.ratio, main = "", breaks = 100, col = "grey", freq = F)
abline(v = quantile(crazy.ratio, prob = c(0.025, 0.975)), col = "red", lwd = 3)


mean(abs(out1B$sims.list$alpha2 / out1B$sims.list$alpha1) > 1)


# Compute expected abundance for a grid of elevation and forest cover
elev.pred <- seq(-1, 1,,100)                       # Values of elevation
forest.pred <- seq(-1,1,,100)                      # Values of forest cover
pred.matrix <- array(NA, dim = c(100, 100)) # Prediction matrix
for(i in 1:100){
  for(j in 1:100){
    pred.matrix[i, j] <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred[i] + out1J$mean$alpha2 * forest.pred[j] + out1J$mean$alpha3 * elev.pred[i] * forest.pred[j]
  }
}

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.main = 1.6, cex.axis = 1.5, cex.lab = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=elev.pred, y= forest.pred, z=pred.matrix, col = mapPalette(100), xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z=pred.matrix, add = TRUE, lwd = 1, cex = 1.5)
title(main = "A")
matpoints(elev, forest, pch="+", cex=1.5)
abline(h = c(-1, -0.5, 0, 0.5, 1))

# Predictions for elev. at specific values of forest cover (-1,-0.5,0,0.5,1)
pred1 <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred + out1J$mean$alpha2 * (-1) + out1J$mean$alpha3 * elev.pred * (-1)
pred2 <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred + out1J$mean$alpha2 * (-0.5) + out1J$mean$alpha3 * elev.pred * (-0.5)
pred3 <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred + out1J$mean$alpha2 * 0 + out1J$mean$alpha3 * elev.pred * 0
# pred3b <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred   # same
pred4 <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred + out1J$mean$alpha2 * 0.5 + out1J$mean$alpha3 * elev.pred * 0.5
pred5 <- out1J$mean$alpha0 + out1J$mean$alpha1 * elev.pred + out1J$mean$alpha2 * 1 + out1J$mean$alpha3 * elev.pred * 1
matplot(seq(-1, 1,,100), cbind(pred1, pred2, pred3, pred4, pred5), type = "l", lty= 1, col = "blue", ylab = "Prediction of mean count", xlab = "Elevation", ylim = c(-1.5, 7), lwd = 2)
title(main = "B")


pred.mat <- array(dim = c(length(elev.pred), length(out1J$sims.list$alpha0)))
for(j in 1:length(out1J$sims.list$alpha0)){
  pred.mat[,j] <- out1J$sims.list$alpha0[j] + out1J$sims.list$alpha1[j] * elev.pred + out1J$sims.list$alpha2[j] * 0.5 + out1J$sims.list$alpha3[j] * elev.pred * 0.5
}

CL <- apply(pred.mat, 1, function(x){quantile(x, prob = c(0.025, 0.975))})
plot(seq(-1, 1,,100), pred4, type = "l", lty= 1, col = "blue", ylab = "Prediction of mean count at forest = -0.5", xlab = "Elevation", las =1, ylim = c(-1.5, 7), lwd = 3)
matlines(seq(-1, 1,,100), t(CL), lty = 1, col = "blue", lwd = 2)
title(main = "C")

pred <- predict(lm(Cmean ~ elev*forest), newdata = data.frame(elev = seq(-1, 1,,100), forest = 0.5), se.fit = TRUE, interval = "confidence")
lines(seq(-1, 1,,100), pred$fit[,1], lty= 2, col = "red", lwd = 3)
matlines(seq(-1, 1,,100), pred$fit[,2:3], lty = 2, col = "red", lwd = 2)



# 5.4 The R package rjags
# ------------------------------------------------------------------------


library(rjags)
load.module("glm")     # Careful with that package, see JAGS discussion list
load.module("dic")

# Have to explicitly list the deviance if want samples
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "deviance")

# Adaptative phase to maximize MCMC efficiency
model <- jags.model(file = "multiple_linear_regression_model.txt", data = win.data, inits = inits, n.chains = nc, n.adapt = 1000)

# Burnin
update(model, nb)

# Generate posterior samples
samples <- coda.samples(model = model, variable.names = params, n.iter = ni - nb, thin = nt)

# Get the summary statistics for the posterior samples
summfit <- summary(samples)
print(summfit, 2)

# Traceplots and posterior densities
plot(samples[,1:4])

# Compute the Brooks-Gelman-Rubin statistic (R-hat)
gelman.diag(samples)

# Compute the effective sample size
effectiveSize(samples)


# Secondary burnin can be applied (e.g. another 500 samples tossed out)
#samples <- window(samples, start = nb + 500 + 1, end = ni)

# More samples can be drawn (starting where the old chains stopped, not starting from 0)
newsamples <- coda.samples(model = model, variable.names = params, n.iter = 1500, thin = nt)

# Combine the new samples with the old ones (ugly but works)
mc1 <- as.mcmc(rbind(samples[[1]], newsamples[[1]]))
mc2 <- as.mcmc(rbind(samples[[2]], newsamples[[2]]))
mc3 <- as.mcmc(rbind(samples[[3]], newsamples[[3]]))
allsamples <- as.mcmc.list(list(mc1, mc2, mc3))

# Mean deviance
Dbar <- summfit$statistics["deviance","Mean"]

# Variance of the deviance
varD <- summfit$statistics["deviance","SD"]^2

# Compute pD and DIC (according to A. Gelman, implemented in R2jags)
pD <- varD/2
DIC <- Dbar + pD

# Another DIC computation (according to M. Plummer). DIC = Penalized deviance
(dic.pD <- dic.samples(model, 2000, "pD"))



# 5.5 Missing values in a Bayesian analysis
# ------------------------------------------------------------------------


# 5.5.1 Some responses missing
# ------------------------------------------------------------------------
# Copy mean counts and turn first 10 into NAs
Cm <- Cmean         # Copy Cmean into Cm
Cm[1:10] <- NA      # turn first 10 into missing

# Bundle data (inside BUGS use Cm for Cmean)
win.data <- list(Cmean = Cm, M = length(Cm), elev = elev, forest = forest)

# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean", "mu")

# ??? or this to get a subset of the parameters
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "Cmean[1:10]", "mu[1:10]")

# Call WinBUGS or JAGS from R (ART <1 min) and summarize posteriors
out1.1 <- bugs(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd(), DIC = FALSE)

out1.1 <- jags(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1.1, 2)

print(cbind(Truth=Cmean[1:10], out1.1$summary[6:15,c(1:3,7)], out1.1$summary[16:25, c(1:3,7)]),3)


# 5.5.2 All responses missing
# ------------------------------------------------------------------------
# Bundle data: simply drop the response from list
win.data <- list(M = length(Cm), elev = elev, forest = forest)

# Alternatively, add all-NA data vector
win.data$Cmean <- as.numeric(rep(NA, length(Cmean)))
str(win.data)  # Cmean is numeric

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "mu[1:2]")

# Call WinBUGS from R (ART <1 min) and summarize posteriors
out1.2 <- bugs(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd(), DIC = FALSE)

print(out1.2, 2)

par(mfrow = c(3, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out1.2$sims.list$alpha0, breaks = 100, col = "grey", main = "")
text(-3500, 600, "A", cex = 2)
hist(out1.2$sims.list$alpha1, breaks = 100, col = "grey", main = "", ylab = "")
text(-3500, 580, "B", cex = 2)
hist(out1.2$sims.list$alpha2, breaks = 100, col = "grey", main = "")
text(-3700, 580, "C", cex = 2)
hist(out1.2$sims.list$sd, breaks = 100, col = "grey", main = "", ylim = c(0, 230), ylab = "")
text(30, 220, "D", cex = 2)
hist(out1.2$sims.list$mu[,1], breaks = 100, col = "grey", main = "")
text(-4000, 480, "E", cex = 2)
hist(out1.2$sims.list$mu[,2], breaks = 100, col = "grey", main = "", ylab = "")
text(-4000, 480, "F", cex = 2)



# 5.5.3 Missing values in a covariate
# ------------------------------------------------------------------------
# Shoot 'holes' in the covariate data
ele <- elev          # copy of elevation covariate
ele[1:10] <- NA      # missing values in covariate elevation

# Bundle data: feed new 'ele' into 'elev' covariate inside of BUGS model
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = ele, forest = forest)

# Specify model in BUGS language
cat(file = "missing_cov_imputation_model_1.txt","
    model {
    
    # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
    }
    
    # Model for missing covariates
    for (i in 1:M){
    elev[i] ~ dnorm(0, 0.001)
    }
    }")

# Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10), alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS from R (ART <1 min)
out1.3 <- bugs(win.data, inits, params, "missing_cov_imputation_model_1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())


# Specify model in BUGS language
cat(file = "missing_cov_imputation_model_2.txt","
    model {
    
    # Priors
    alpha0 ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1 ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2 ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3 ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)      # precision tau = 1 / variance
    mu[i] <- alpha0 + alpha1 * elev[i] + alpha2 * forest[i] + alpha3 * elev[i] * forest[i]
    }
    
    # Covariate mean as a model for missing covariates
    for (i in 1:M){
    elev[i] ~ dnorm(mu.elev, tau.elev)    # Assume elevation normally distributed
    }
    mu.elev ~ dnorm(0, 0.0001)
    tau.elev <- pow(sd.elev, -2)
    sd.elev ~ dunif(0, 100)
    }")

# Initial values
inits <- function() list(alpha0 = rnorm(1,,10), alpha1 = rnorm(1,,10), alpha2 = rnorm(1,,10), alpha3 = rnorm(1,,10))

# Parameters monitored
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd", "elev", "mu.elev", "sd.elev")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS from R (ART <1 min)
out1.4 <- bugs(win.data, inits, params, "missing_cov_imputation_model_2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())


par(cex = 1.5, lwd = 2)
plot(elev[1:10]-0.01, out1.3$summary[6:15,1], ylim = c(-10, 10), col = "red", xlab = "True value of covariate", ylab = "Estimate (with 95% CRI)", frame.plot =F)
segments(elev[1:10]-0.01, out1.3$summary[6:15,3], elev[1:10]-0.01, out1.3$summary[6:15,7], col = "red")
points(elev[1:10]+0.01, out1.4$summary[6:15,1], ylim = c(-3, 3), col = "blue")
segments(elev[1:10]+0.01, out1.4$summary[6:15,3], elev[1:10]+0.01, out1.4$summary[6:15,7], col = "blue")
abline(0,1)



# 5.6 Linear model with normal response (normal GLM): analysis of covariance (ANCOVA)
# -----------------------------------------------------------------------------------


# Generate factor and plot raw data in boxplot as function of factor A
facFor <- as.numeric(forest < -0.5)         # Factor level 1
facFor[forest < 0 & forest > -0.5] <- 2     # Factor level 2
facFor[forest < 0.5 & forest > 0] <- 3      # Factor level 3
facFor[forest > 0.5] <- 4                   # Factor level 4
table(facFor)                               # every site assigned a level OK

par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Cmean ~ factor(facFor), col = c("red", "blue", "green", "grey"), xlab = "Forest cover class", ylab = "Mean count of great tits", frame.plot = F, ylim = c(0,20))
text(0.8, 20, "A", cex=1.6)


# Bundle data
win.data <- list(Cmean = Cmean, M = length(Cmean), elev = elev, facFor = facFor)

# Specify model in BUGS language in effects parameterisation
cat(file = "ANCOVA1.txt","
    model {
    
    # Priors
    alpha ~ dnorm(0, 1.0E-06)            # Prior for intercept = effect of level 1 of forest factor
    beta2 ~ dnorm(0, 1.0E-06)            # Prior for slope = effect of elevation for level 1 of forest factor
    beta1[1] <- 0                        # Set to zero effect of first level of facFor
    beta3[1] <- 0                        # Set to zero effect of first level of facFor of elevation
    for(k in 2:4){
    beta1[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
    beta3[k] ~ dnorm(0, 1.0E-06)       # Prior for effects of factor facFor
    }
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha + beta1[facFor[i]] + beta2 * elev[i] + beta3[facFor[i]] * elev[i]
    }
    }
    ")

# Initial values
inits <- function() list(alpha = rnorm(1,,10), beta1 = c(NA, rnorm(3,,10)), beta2 = rnorm(1,,10), beta3 = c(NA, rnorm(3,,10)))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "sd")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min)
out3 <- bugs(win.data, inits, params, "ANCOVA1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

out3J <- jags(win.data, inits, params, "ANCOVA1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
traceplot(out3J)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev)))

# Summarize posteriors
print(out3, 3)


# Specify model in BUGS language
cat(file = "ANCOVA2.txt","
    model {
    
    # Priors
    for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)       # Priors for intercepts
    beta[k] ~ dnorm(0, 1.0E-06)        # Priors for slopes
    }
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu[i], tau)          # precision tau = 1 / variance
    mu[i] <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    }
    
    # Derived quantities: comparison of slopes (now you can forget the delta rule !)
    for(k in 1:4){
    diff.vs1[k] <- beta[k] - beta[1]    # Differences relative to beta[1]
    diff.vs2[k] <- beta[k] - beta[2]    # ... relative to beta[2]
    diff.vs3[k] <- beta[k] - beta[3]    # ... relative to beta[3]
    diff.vs4[k] <- beta[k] - beta[4]    # ... relative to beta[4]
    }
    }
    ")

# Initial values
inits <- function() list(alpha = rnorm(4,,10), beta = rnorm(4,,10))

# Parameters monitored
params <- c("alpha", "beta", "sd", "diff.vs1", "diff.vs2", "diff.vs3", "diff.vs4")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min) and summarize posteriors
out4 <- bugs(win.data, inits, params, "ANCOVA2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

system.time(out4J <- jags(win.data, inits, params, "ANCOVA2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb))
traceplot(out4J)

print(out4, 2)

# Fit model using least-squares (produces MLEs)
(fm <- summary(lm(Cmean ~ as.factor(facFor)*elev-1-elev)))


plot(elev[facFor==1], Cmean[facFor==1], col = "red", ylim = c(0, 20), xlab = "Elevation", ylab = "", frame.plot = F)
points(elev[facFor==2], Cmean[facFor==2], col = "blue")
points(elev[facFor==3], Cmean[facFor==3], col = "green")
points(elev[facFor==4], Cmean[facFor==4], col = "black")
abline(fm$coef[1,1], fm$coef[5,1], col = "red")
abline(fm$coef[2,1], fm$coef[6,1], col = "blue")
abline(fm$coef[3,1], fm$coef[7,1], col = "green")
abline(fm$coef[4,1], fm$coef[8,1], col = "black")
text(-0.8, 20, "B", cex=1.6)


attach.bugs(out4)     # Allows to directly address the sims.list
str(diff.vs3)
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(diff.vs3[,1], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.2, 0.8, "A", cex = 2)
hist(diff.vs3[,2], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-1.4, 0.8, "B", cex = 2)
hist(diff.vs3[,4], col = "grey", breaks = 100, main = "", freq=F, ylim = c(0, 0.8))
abline(v = 1, lwd = 3, col = "red")
text(-2.2, 0.8, "C", cex = 2)

# Prob. difference greater than 1
mean(diff.vs3[,1] > 1)
mean(diff.vs3[,2] > 1)
mean(diff.vs3[,4] > 1)



# 5.7 Proportion of variance explained (R2)
# ------------------------------------------------------------------------


cat(file = "Model0.txt","
    model {
    # Priors
    mu ~ dnorm(0, 1.0E-06)
    tau <- pow(sd, -2)
    sd ~ dunif(0, 1000)
    # Likelihood
    for (i in 1:M){
    Cmean[i] ~ dnorm(mu, tau)
    }
    }
    ")
inits <- function() list(mu = rnorm(1))
params <- c("mu", "sd")
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3
out0 <- jags(win.data, inits, params, "Model0.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

print(out0)

# Compute R2 from BUGS analysis
(total.var <- mean(out0$sims.list$sd^2))       # Total variance around the mean
(unexplained.var <- mean(out3$sims.list$sd^2)) # Not explained by the ANCOVA
(prop.explained <- (total.var - unexplained.var)/total.var)



# 5.8 Fitting a model with non-standard likelihood using the zeros or the ones tricks
# ------------------------------------------------------------------------------------



# Package the data needed in a bundle
win.data <- list(Cmean1 = Cmean, Cmean2 = Cmean, zeros = rep(0, M), ones = rep(1, M), M = length(Cmean), elev = elev, forest = forest) # note 2 copies of response

# Write text file with model description in BUGS language
cat(file = "multiple_linear_regression_model.txt",
    "model {
    
    # Priors
    for(k in 1:3){ # Loop over three ways to specify likelihood
    alpha0[k] ~ dnorm(0, 1.0E-06)           # Prior for intercept
    alpha1[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of elev
    alpha2[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of forest
    alpha3[k] ~ dnorm(0, 1.0E-06)           # Prior for slope of interaction
    sd[k] ~ dunif(0, 1000)                  # Prior for dispersion on sd scale
    }
    var1 <- pow(sd[1], 2)                      # Variance in zeros trick
    var2 <- pow(sd[2], 2)                      # Variance in ones trick
    tau <- pow(sd[3], -2)                      # Precision tau = 1/(sd^2)
    
    C1 <- 10000 # zeros trick: make large enough to ensure lam >= 0
    C2 <- 10000 # ones trick: make large enough to ensure p <= 1
    pi <- 3.1415926
    
    # Three variants of specification of the likelihood
    for (i in 1:M){
    # 'Zeros trick' for normal likelihood
    zeros[i] ~ dpois(phi[i])  # likelihood contribution is exp(-phi)
    #   negLL[i] <- log(sd[1]) + 0.5 * pow((Cmean1[i] - mu1[i]) / sd[1],2 )
    negLL[i] <- -log(sqrt(1/(2*pi*var1))) + pow(Cmean1[i]-mu1[i],2)/(2*var1)
    phi[i] <- negLL[i] + C1
    mu1[i] <- alpha0[1] + alpha1[1]*elev[i] + alpha2[1]*forest[i] + alpha3[1]*elev[i]*forest[i]
    
    # 'Ones trick' for normal likelihood
    ones[i] ~ dbern(p[i])  # likelihood contribution is p directly
    L[i] <- sqrt(1/(2*pi*var2)) * exp(-pow(Cmean1[i]-mu2[i],2)/(2*var2))
    p[i] <- L[i] / C2
    mu2[i] <- alpha0[2] + alpha1[2]*elev[i] + alpha2[2]*forest[i] + alpha3[2]*elev[i]*forest[i]
    
    # Standard distribution function for the normal
    Cmean2[i] ~ dnorm(mu3[i], tau)
    mu3[i] <- alpha0[3] + alpha1[3]*elev[i] + alpha2[3]*forest[i] + alpha3[3]*elev[i]*forest[i]
    }
    }"
)

# Initial values
inits <- function() list(alpha0 = rnorm(3, 0, 10), alpha1 = rnorm(3,0,10), alpha2 = rnorm(3,0,10), alpha3 = rnorm(3,0,10))

# Parameters monitored (i.e., for which estimates are saved)
params <- c("alpha0", "alpha1", "alpha2", "alpha3", "sd")

# MCMC settings
ni <- 1200   ;   nt <- 1   ;   nb <- 200   ;  nc <- 3    # For JAGS

# Call JAGS
library(jagsUI)
outX <- jags(win.data, inits, params, "multiple_linear_regression_model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(outX)


# Define negative log-likelihood function
neglogLike <- function(param) {
  alpha0 = param[1]
  alpha1 = param[2]
  alpha2 = param[3]
  alpha3 = param[4]
  sigma = exp(param[5])    # Estimate sigma on log-scale
  mu = alpha0 + alpha1*elev + alpha2*forest + alpha3*elev*forest
  #   -sum(dnorm(Cmean, mean=mu, sd=sigma, log=TRUE))  # cheap quick way
  sum(-log(sqrt(1/(2*3.1415926*sigma^2))) + (Cmean-mu)^2/(2*sigma^2))
}

# Find parameter values that minimize function value
(fit <- optim(par = rep(0, 5), fn = neglogLike, method = "BFGS"))

exp(fit$par[5])            # Backtransform to get sigma



# 5.9 Poisson generalized linear model (Poisson GLM)
# --------------------------------------------------


# Summarize data by taking max at each site
Cmax <- apply(C, 1, max)
table(Cmax)

# Bundle data
win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001)

# Specify model in BUGS language
cat(file = "Poisson_GLM.txt","
    model {
    
    # Priors
    for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)       # Prior for intercepts
    beta[k] ~ dnorm(0, 1.0E-06)        # Prior for slopes
    }
    
    # Likelihood
    for (i in 1:M){
    Cmax[i] ~ dpois(lambda[i])         # note no variance parameter
    log(lambda[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    resi[i] <- (Cmax[i]-lambda[i]) / (sqrt(lambda[i])+e)   # Pearson resi
    }
    }
    ")


# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("alpha", "beta", "lambda", "resi")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R and summarize posteriors
out5 <- bugs(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

system.time(out5J <- jags(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb))
par(mfrow = c(4,2))    ;    traceplot(out5J, c("alpha[1:4]", "beta[1:4]"))
print(out5J, 3)

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex = 1.3, cex.lab = 1.5, cex.axis = 1.5)
hist(out5$summary[276:542, 1], xlab = "Pearson residuals", col = "grey", breaks = 50, main = "", freq = F, xlim = c(-5, 5), ylim = c(0, 0.57))
abline(v = 0, col = "red", lwd = 2)
text(-4.7, 0.54, "A", cex = 1.5)

plot(1:267, out5$summary[276:542, 1], main = "", xlab = "Order of data", ylab = "Pearson residual", frame.plot = F)
abline(h = 0, col = "red", lwd = 2)
text(8, 4, "B", cex = 1.5)

plot(out5$summary[9:275, 1],out5$summary[276:542, 1], main = "", xlab = "Predicted values", ylab = "Pearson residual", frame.plot = F, xlim = c(-1, 14))
abline(h = 0, col = "red", lwd = 2)
text(-0.5, 4, "C", cex = 1.5)


summary(glm(Cmax ~ factor(facFor)*elev-1-elev, family = poisson))


lambda2 <- array(dim = c(15000, 267))
for(j in 1:267){                            # Loop over sites
  lambda2[,j] <- exp(out5$sims.list$alpha[,facFor[j]] + out5$sims.list$beta[,facFor[j]] * elev[j]) # linear regression/backtransform
}
plot(out5$sims.list$lambda ~ lambda2, pch = ".")  # Check the two are identical
lm(c(out5$sims.list$lambda) ~ c(lambda2))

sorted.ele1 <- sort(elev[facFor == 1])
sorted.y1 <- out5$summary[9:275,][facFor == 1,][order(elev[facFor == 1]),]

# Plot A
par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(elev[facFor == 1], jitter(Cmax[facFor ==1]), ylab = "Maximum count", xlab = "Elevation (scaled)", frame.plot=F), ylim = c(0, 6))
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2) # Post. mean
lines(sorted.ele1, sorted.y1[,3], col = "grey", lwd = 2) # Lower 95% CL
lines(sorted.ele1, sorted.y1[,7], col = "grey", lwd = 2) # Upper 95% CL
text(-0.8, 6, "A", cex = 2)

# Plot B
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
polygon(c(sorted.ele1, rev(sorted.ele1)), c(sorted.y1[,3], rev(sorted.y1[,7])), col='grey', border=NA)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "B", cex = 2)

# Plot C
elev.pred <- seq(-1,1, length.out = 200)  # Cov. for which to predict lambda
n.pred <- 50                             # Number of prediction profiles
pred.matrix <- array(NA, dim = c(length(elev.pred), n.pred))
for(j in 1:n.pred){
  sel <- sample(1:length(out5$sims.list$alpha[,1]),1) # Choose one post. draw
  pred.matrix[,j] <- exp(out5$sims.list$alpha[sel,1] + out5$sims.list$beta[sel,1] * elev.pred)
}
plot(sorted.ele1, sorted.y1[,1], type='n', xlab = "Elevation (scaled)", ylab = "", frame.plot = F, ylim = c(0, 6))
matlines(elev.pred, pred.matrix, col = "grey", lty = 1, lwd = 1)
lines(sorted.ele1, sorted.y1[,1], col = "blue", lwd = 2)
text(-0.8, 6, "C", cex = 2)



# 5.10 Goodness of fit assessment: posterior predictive checks and the parametric bootstrap
# -----------------------------------------------------------------------------------------


# Bundle data
win.data <- list(Cmax = Cmax, M = length(Cmax), elev = elev, facFor = facFor, e = 0.0001)

# Specify model in BUGS language
cat(file = "Poisson_GLM.txt","
    model {
    # Priors
    for(k in 1:4){
    alpha[k] ~ dnorm(0, 1.0E-06)
    beta[k] ~ dnorm(0, 1.0E-06)
    }
    
    # Likelihood and computations for posterior predictive check
    for (i in 1:M){
    Cmax[i] ~ dpois(lambda[i])
    log(lambda[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    
    # Fit assessments: Chi-squared test statistic and posterior predictive check
    chi2[i] <- pow((Cmax[i]-lambda[i]),2) / (sqrt(lambda[i])+e)         # obs.
    Cmax.new[i] ~ dpois(lambda[i])      # Replicate (new) data set
    chi2.new[i] <- pow((Cmax.new[i]-lambda[i]),2) / (sqrt(lambda[i])+e) # exp.
    }
    # Add up discrepancy measures for entire data set
    fit <- sum(chi2[])                     # Omnibus test statistic actual data
    fit.new <- sum(chi2.new[])             # Omnibus test statistic replicate data
    
    # range of data as a second discrepancy measure
    obs.range <- max(Cmax[]) - min(Cmax[])
    exp.range <- max(Cmax.new[]) - min(Cmax.new[])
    }
    ")

# Initial values
inits <- function() list(alpha = rnorm(4,,3), beta = rnorm(4,,3))

# Parameters monitored
params <- c("chi2", "fit", "fit.new", "obs.range", "exp.range")
params <- c("Cmax.new", "chi2.new", "chi2", "fit", "fit.new", "obs.range", "exp.range")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R and summarize posteriors
out5.1 <- jags(win.data, inits, params, "Poisson_GLM.txt", n.chains = nc,
               n.thin = nt, n.iter = ni, n.burnin = nb, modules = 'bugs')

print(out5.1, 2)

par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(out5.1$sims.list$fit, out5.1$sims.list$fit.new, xlim = c(200, 1000), ylim = c(200, 1000), main = "", xlab = "Discrepancy observed data", ylab = "Discrepancy expected data", frame.plot = F, cex = 1.5)
abline(0,1, lwd = 2)
text(240, 1000, "A", cex = 2)



(bpv <- mean(out5.1$sims.list$fit.new > out5.1$sims.list$fit))

plot(Cmax, out5.1$mean$chi2, xlab = "Observed data", ylab = "Chi2 contribution", frame.plot = F, ylim = c(0, 70), cex = 1.5)
lines(0:30, sqrt(0:30), lwd = 2)
text(2, 70, "B", cex = 2)


# Fit model to actual data and compute two fit statistics
fm <- glm(Cmax ~ as.factor(facFor)*elev, family = poisson) # Fit model
observed <- Cmax
expected <- predict(fm, type = "response")
plot(observed, expected)
abline(0,1)
chi2.obs <- sum((observed - expected)^2 / (expected + 0.0001)) # fit stat 1
range.obs <- diff(range(observed))                             # fit stat 2

# Generate reference distribution of fit stat under a fitting model
simrep <- 100000      # Might want to try 1000 first
chi2vec <- rangevec <- numeric(simrep)  # vectors for chi2 and maximum
for(i in 1:simrep){
  cat(paste(i, "\n"))
  Cmaxrep <- rpois(n = 267, lambda = expected)   # Generate replicate data set
  fmrep <- glm(Cmaxrep ~ as.factor(facFor)*elev, family = poisson) # Refit model
  expectednew <- predict(fmrep, type = "response")
  chi2vec[i] <- sum((Cmaxrep - expectednew)^2 / (expectednew + 0.0001))
  rangevec[i] <- diff(range(Cmaxrep))
}

# Summarize bootstrap results and compare with posterior predictive dist.
par(mfrow = c(2, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
hist(out5.1$sims.list$fit.new, col = "grey", main = "", breaks = 100, xlim = c(180, 900), freq = F, ylim = c(0, 0.01))
abline(v = mean(out5.1$sims.list$fit), col = "red", lwd = 2)
text(200, 0.009, "A", cex = 2)
hist(out5.1$sims.list$exp.range, col = "grey", main = "", breaks = 50, xlim = c(10, 32), freq = F, ylim = c(0, 0.40))
abline(v = mean(out5.1$sims.list$obs.range), col = "red", lwd = 2)
text(11, 0.38, "B", cex = 2)
hist(chi2vec, col = "grey", main = "", breaks = 100, xlim = c(180, 900), freq = F, ylim = c(0, 0.02))
abline(v = chi2.obs, col = "red", lwd = 2)
text(200, 0.018, "C", cex = 2)
hist(rangevec, col = "grey", main = "", breaks = 50, xlim = c(10, 32), freq = F)
abline(v = range.obs, col = "red", lwd = 2)
text(11, 0.38, "D", cex = 2)

# Lack of fit ratio in PPD and parboot
mean(out5.1$sims.list$fit/out5.1$sims.list$fit.new) # ppc
mean(chi2.obs/chi2vec)                              # parboot

(pval1 <- 1-rank(c(chi2vec, chi2.obs))[simrep+1]/(simrep+1))
(pval2 <- 1-rank(c(rangevec, range.obs))[simrep+1]/(simrep+1))



# 5.11 Binomial generalised linear model (binomial GLM, logistic regression)
# --------------------------------------------------------------------------


# Quantize counts from first survey and describe
y1 <- as.numeric(C[,1] > 0)  # Gets 1 if first count greater than zero
table(y1)

mean(N > 0)          # True occupancy
mean(y1)             # Observed occupancy after first survey

# Bundle data
win.data <- list(y1 = y1, M = length(y1), elev = elev, facFor = facFor)

# Specify model in BUGS language
cat(file = "Bernoulli_GLM.txt","
    model {
    
    # Priors
    for(k in 1:4){
    alpha[k] <- logit(mean.psi[k])     # intercepts
    mean.psi[k] ~ dunif(0,1)
    beta[k] ~ dnorm(0, 1.0E-06)        # slopes
    }
    
    # Likelihood
    for (i in 1:M){
    y1[i] ~ dbern(theta[i])
    logit(theta[i]) <- alpha[facFor[i]] + beta[facFor[i]] * elev[i]
    }
    }
    ")

# Initial values
inits <- function() list(mean.psi = runif(4), beta = rnorm(4,,3))   # Priors 2

# Parameters monitored
params <- c("mean.psi", "alpha", "beta", "theta")

# MCMC settings
ni <- 6000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call WinBUGS or JAGS from R (ART <1 min)
out6 <- bugs(win.data, inits, params, "Bernoulli_GLM.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

out6J <- jags(win.data, inits, params, "Bernoulli_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(4,2))    ;    traceplot(out6J, c("alpha[1:4]", "beta[1:4]"))

print(out6, 2)

# Compare with MLEs
summary(glm(y1 ~ factor(facFor)*elev-1-elev, family = binomial))


# Plot of observed response vs. two covariates
par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
F1 <- facFor == 1 ; F2 <- facFor == 2 ; F3 <- facFor == 3 ; F4 <- facFor == 4
plot(jitter(y1,,0.05) ~ facFor, xlab = "Forest factor", ylab = "Observed occupancy probability", frame.plot = F, ylim = c(0, 1.15))
lines(1:4, out6$summary[1:4,1], lwd = 2)
segments(1:4, out6$summary[1:4,3], 1:4, out6$summary[1:4,7])
text(1.15, 1.1, "A", cex=1.6)

plot(elev[F1], jitter(y1,,0.1)[F1], xlab = "Elevation", ylab = "", col = "red", frame.plot = F)
points(elev[F2], jitter(y1,,0.05)[F2], col = "blue")
points(elev[F3], jitter(y1,,0.05)[F3], col = "green")
points(elev[F4], jitter(y1,,0.05)[F4], col = "grey")
lines(sort(elev[F1]), out6$mean$theta[F1][order(elev[F1])], col="red", lwd=2)
lines(sort(elev[F2]), out6$mean$theta[F2][order(elev[F2])], col="blue", lwd=2)
lines(sort(elev[F3]), out6$mean$theta[F3][order(elev[F3])], col="green", lwd=2)
lines(sort(elev[F4]), out6$mean$theta[F4][order(elev[F4])], col="grey", lwd=2)
text(-0.9, 1.1, "B", cex=1.6)



# 5.12 Moment-matching in a binomial GLM to accommodate underdispersion
# ---------------------------------------------------------------------


# Bundle data
win.data <- list(y = data$C, M = nrow(data$C), J = ncol(data$C), elev = elev, N = 32)

# Specify model in BUGS language
cat(file = "squeezed_count_GLM.txt","
    model {
    
    # Priors
    alpha ~ dnorm(0, 1.0E-06)
    beta ~ dnorm(0, 1.0E-06)
    
    # Likelihood
    for (i in 1:M){
    mu[i] <- alpha + beta * elev[i] # linear model for expected response
    logit(p[i]) <- logit(mu[i] / N) # express param as function of first moment
    for(j in 1:J){
    y[i,j] ~ dbin(p[i], N)
    }
    }
    }
    ")

# Initial values
inits <- function() list(alpha = 1.7, beta = -1.2)          # works always
inits <- function() list(alpha = runif(1), beta = runif(1)) # works sometimes

# Parameters monitored
params <- c("alpha", "beta", "mu")

# MCMC settings
ni <- 3000   ;   nt <- 1   ;   nb <- 1000   ;  nc <- 3

# Call JAGS from R (ART <1 min)
library(jagsUI)
out7 <- jags(win.data, inits, params, "squeezed_count_GLM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(4,2))   ;   traceplot(out7)
print(out7, 3)



# 5.13 Random-effects Poisson GLM (Poisson GLMM)
# ----------------------------------------------


# Bundle data
win.data <- list(C = C, M = nrow(C), J = ncol(C), elev = elev, forest = forest, elev.forest = elev * forest, wind = wind)

# Specify model in BUGS language
cat(file = "RE.Poisson.txt","
    model {
    
    # Priors
    mu.alpha ~ dnorm(0, 0.001)                # Mean hyperparam
    tau.alpha <- pow(sd.alpha, -2)
    sd.alpha ~ dunif(0, 10)                   # sd hyperparam
    for(k in 1:4){
    alpha[k] ~ dunif(-10, 10)              # Regression params
    }
    
    # Likelihood
    for (i in 1:M){
    alpha0[i] ~ dnorm(mu.alpha, tau.alpha) # Random effects and hyperparams
    re0[i] <- alpha0[i] - mu.alpha         # zero-centered random effects
    for(j in 1:J){
    C[i,j] ~ dpois(lambda[i,j])
    log(lambda[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + alpha[3] * elev.forest[i] + alpha[4] * wind[i,j]
    }
    }
    }")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha = rnorm(4)) # Inits
params <- c("mu.alpha", "sd.alpha", "alpha0", "alpha", "re0") # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call WinBUGS or JAGS from R (ART 6-7 min) and summarize posteriors
out8 <- bugs(win.data, inits, params, "RE.Poisson.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

out8 <- jags(win.data, inits, params, "RE.Poisson.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(3,2))  ;  traceplot(out8, c("mu.alpha", "sd.alpha", "alpha[1:3]"))

print(out8, 3)

Cvec <- as.vector(C)            # Vector of M*J counts
elev.vec <- rep(elev, J)        # Vectorized elevation covariate
forest.vec <- rep(forest, J)    # Vectorized forest covariate
wind.vec <- as.vector(wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:M, J)) # Site indicator (factor)
cbind(Cvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood (NOTE: glmer uses ML instead of REML)
library(lme4)
summary(fm <- glmer(Cvec ~ elev.vec*forest.vec + wind.vec + (1| fac.site), family = poisson))              # Fit model
ranef(fm)                       # Print zero-centered random effects


# Compare fixed-effects estimates (in spite of the confusing naming in glmer output), Bayesian post. means and sd left, frequentist MLEs and SEs right
print(cbind(out8$summary[c(1:2, 270:273), 1:2], rbind(summary(fm)$coef[1,1:2], c(sqrt(summary(fm)$varcor$fac.site), NA), summary(fm)$coef[c(2,3,5,4),1:2])), 3)

# Compare graphically non-Bayesian and Bayesian random effects estimates
Freq.re <- ranef(fm)$fac.site[,1]         # Non-Bayesian estimates (MLEs)
Bayes.re <- out8$summary[274:540,]        # Bayesian estimates

par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(Freq.re, Bayes.re[,1], xlab = "Non-Bayesian (glmer)", ylab = "Bayesian (BUGS)", xlim = c(-0.4, 0.42), ylim = c(-2, 2), frame.plot = F, type = "n")
segments(Freq.re, Bayes.re[,3], Freq.re, Bayes.re[,7], col = "grey", lwd = 0.5)
abline(0, 1, lwd = 2)
points(Freq.re, Bayes.re[,1])
text(-0.38, 2, "A", cex=1.6)

wind.pred <- seq(-1, 1, , 1000)     # Covariate values for prediction
pred <- array(NA, dim = c(1000, 267))
for(i in 1:267){
  pred[,i]<- exp(out8$mean$alpha0[i] + out8$mean$alpha[1] * 0 + out8$mean$alpha[2] * 0 + out8$mean$alpha[3] * 0 + out8$mean$alpha[4] * wind.pred)    # Predictions for each site
}

matplot(wind.pred, pred, type = "l", lty = 1, col = "grey", xlab = "Wind speed", ylab = "Expected count", frame.plot = F, ylim = c(0, 4))
lines(wind.pred, exp(out8$mean$mu.alpha + out8$mean$alpha[4] * wind.pred), col = "black", lwd = 3)
text(-0.9, 4, "B", cex=1.6)



# 5.14 Random-effects binomial GLM (binomial GLMM)
# ------------------------------------------------


# Get detection/nondetection response
y <- C
y[y > 0] <- 1

# Bundle data
win.data <- list(y = y, M = nrow(y), J = ncol(y), elev = elev, forest = forest, elev.forest = elev * forest, wind = wind)
str(win.data)

# Specify model in BUGS language
cat(file = "RE.Bernoulli.txt","
    model {
    
    # Priors
    mu.alpha0 <- logit(mean.theta)              # Random intercepts
    mean.theta ~ dunif(0,1)
    tau.alpha0 <- pow(sd.alpha0, -2)
    sd.alpha0 ~ dunif(0, 10)
    mu.alpha4 ~ dnorm(0, 0.001)                 # Random slope on wind
    tau.alpha4 <- pow(sd.alpha4, -2)
    sd.alpha4 ~ dunif(0, 10)
    for(k in 1:3){
    alpha[k] ~ dnorm(0, 0.001)               # Slopes
    }
    
    # Likelihood
    for (i in 1:M){
    alpha0[i] ~ dnorm(mu.alpha0, tau.alpha0) # Intercept random effects
    re00[i] <- alpha0[i] - mu.alpha0         # same zero-centered
    alpha4[i] ~ dnorm(mu.alpha4, tau.alpha4) # Slope random effects
    re04[i] <- alpha4[i] - mu.alpha4         # same zero-centered
    for(j in 1:J){
    y[i,j] ~ dbern(theta[i,j])
    logit(theta[i,j]) <- alpha0[i] + alpha[1] * elev[i] + alpha[2] * forest[i] + alpha[3] * elev.forest[i] + alpha4[i] * wind[i,j]
    }
    }
    }")

# Other model run preparations
inits <- function() list(alpha0 = rnorm(M), alpha4 = rnorm(M))# Inits
params <- c("mu.alpha0", "sd.alpha0", "alpha0", "alpha", "mu.alpha4", "sd.alpha4", "alpha4", "re00", "re04")                        # Params
ni <- 30000 ; nt <- 25 ; nb <- 5000 ; nc <- 3                 # MCMC settings

# Call WinBUGS from R .... and crash !
out9 <- bugs(win.data, inits, params, "RE.Bernoulli.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Call JAGS from R (ART 2.5 min)
out9 <- jags(win.data, inits, params, "RE.Bernoulli.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
par(mfrow = c(2,2))
traceplot(out9, c("mu.alpha0", "sd.alpha0", "alpha[1:3]", "mu.alpha4", "sd.alpha4"))
print(out9, 3)

yvec <- as.vector(y)            # Vector of M*J counts
elev.vec <- rep(elev, J)        # Vectorized elevation covariate
forest.vec <- rep(forest, J)    # Vectorized forest covariate
wind.vec <- as.vector(wind)     # Vectorized wind covariate
fac.site <- factor(rep(1:M, J)) # Site indicator (factor)
cbind(yvec, fac.site, elev.vec, forest.vec, wind.vec) # Look at data

# Fit same model using maximum likelihood
library(lme4)                   # Load package
summary(frem <- glmer(yvec ~ elev.vec*forest.vec + wind.vec + (wind.vec || fac.site), family = binomial))              # Fit model


# Compare Bayesian and non-Bayesian estimates
print(out9$summary[c(1:2, 270:274),c(1:3,7:9)], 4)

(re <- ranef(frem))                 # Print zero-centered random effects

par(mfrow = c(1, 3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
pop.mean.int <- summary(frem)$coef[1,1]
pop.mean.slope <- summary(frem)$coef[4,1]
plot(sort(wind.vec), plogis(pop.mean.int + sort(wind.vec) * pop.mean.slope), type = "l", xlab = "Wind speed", ylab = "Prob. to count >0 great tits", lwd = 3, frame.plot = F)
for(i in 1:267){
  lines(sort(wind.vec), plogis(pop.mean.int + re$fac.site[i,1] + sort(wind.vec) * (pop.mean.slope + re$fac.site[i,2])), lwd = 1, col = i)
}
title(main = "A", cex.main = 2)

# Compute expected detection/nondetection probability for a grid of elevation and forest cover, at wind-speed = 0 (covariate average) and for hypermean of intercepts alpha0
n.sims <- length(out9$sims.list$mu.alpha0)
elev.pred <- seq(-1, 1,,100)                       # Values of elevation
forest.pred <- seq(-1,1,,100)                      # Values of forest cover
pred.array <- array(NA, dim = c(100, 100, n.sims)) # Prediction array
for(i in 1:100){
  for(j in 1:100){
    pred.array[i,j,] <- plogis(out9$sims.list$mu.alpha0 + out9$sims.list$alpha[,1] * elev.pred[i] + out9$sims.list$alpha[,2] * forest.pred[j] + out9$sims.list$alpha[,3] * elev.pred[i] * forest.pred[j])
  }
}
pm.pred.array <- apply(pred.array, c(1,2), mean)   # Get posterior mean
psd.pred.array <- apply(pred.array, c(1,2), sd)    # Get posterior sd

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
image(x=elev.pred, y= forest.pred, z=pm.pred.array, col = mapPalette(100), xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z= pm.pred.array, add = TRUE, lwd = 1)
title(main = "B", cex.main = 2)

image(x=elev.pred, y= forest.pred, z=psd.pred.array, col = mapPalette(100), xlab = "Elevation", ylab = "Forest cover")
contour(x=elev.pred, y=forest.pred, z= psd.pred.array, add = TRUE, lwd = 1)
title(main = "C", cex.main = 2)



# 5.15 General strategy of model building with BUGS
# --------------------------------------------------

# 5.16 Summary and outlook
# --------------------------------------------------










# =========================================================================
#
# 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
#
# =========================================================================






# 6.1. Introduction to the modeling of abundance
# ------------------------------------------------------------------------



# 6.2 An exercise in hierarchical modeling: Derivation of binomial N-mixture models from first principles
# ------------------------------------------------------------------------



# 6.3. Simulation and analysis of the simplest possible N-mixture model
# ------------------------------------------------------------------------


# Choose sample sizes and prepare observed data array C
set.seed(24)                # So we all get same data set
M <- 150                    # Number of sites
J <- 2                      # Number of abu. measurements per site (rep. counts)
C <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data

# Parameter values
lambda <- 2.5               # Expected abundance
p <- 0.4                    # Probability of detection (per individual)

# Generate local abundance data (the truth)
N <- rpois(n = M, lambda = lambda)

# Conduct repeated measurements (generate replicated counts)
for(j in 1:J){
  C[,j] <- rbinom(n = M, size = N, prob = p)
}

# Look at data
# The truth ....
table(N)                    # True abundance distribution
sum(N)                      # True total population size at M sites
sum(N>0)                    # True number of occupied sites
mean(N)                     # True mean abundance (estimate of lambda)

# ... and the observations
table(apply(C, 1, max))     # Observed abundance distribution (max count)
sum(apply(C, 1, max))       # Observed total population size at M sites
sum(apply(C, 1, max)>0)     # Observed number of occupied sites
mean(apply(C, 1, max))      # Observed mean "relative abundance"

head(cbind(N=N, count1=C[,1], count2=C[,2])) # First 6 sites

cor(C)[1,2]

library(unmarked)                  # Load package
umf <- unmarkedFramePCount(y = C)  # Create um data frame
summary(umf)                       # Summarize
(fm1 <- pcount(~1 ~1, data = umf)) # Fit model: get estimates on link scale
backTransform(fm1, "state")        # Get estimates on natural scale
backTransform(fm1, "det")


# Bundle and summarize data set
win.data <- list(C = C, M = nrow(C), J = ncol(C))
str(win.data)                      # Look at data

# Specify model in BUGS language
sink("model1.txt")
cat("
    model {
    # Priors
    lambda ~ dgamma(0.001, 0.001)
    p ~ dunif(0, 1)
    # Likelihood
    for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
    C[i,j] ~ dbin(p, N[i]) # Observation model
    }
    }
    }
    ",fill = TRUE)
sink()

# Initial values
Nst <- apply(C, 1, max)       # Avoid data/model/inits conflict
inits <- function(){list(N = Nst)}

# Parameters monitored
params <- c("lambda", "p")

# MCMC settings
ni <- 25000   ;   nt <- 20   ;   nb <- 5000   ;   nc <- 3

# Call JAGS (ART 1 min) and summarize posteriors
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model1.txt", n.chains = nc,
            n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)



# 6.4 A slightly more complex N-mixture model with covariates
# ------------------------------------------------------------------------



# Choose sample sizes and prepare obs. data array y
set.seed(1)                   # So we all get same data set
M <- 100                      # Number of sites
J <- 3                        # Number of repeated abundance measurements
C <- matrix(NA, nrow = M, ncol = J) # to contain the observed data

# Create a covariate called vegHt
vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience

# Choose parameter values for abundance model and compute lambda
beta0 <- 0                    # Log-scale intercept
beta1 <- 2                    # Log-scale slope for vegHt
lambda <- exp(beta0 + beta1 * vegHt) # Expected abundance
plot(vegHt, lambda, type = "l", lwd = 3)  # Expected abundance

# Draw local abundance and look at data so far
N <- rpois(M, lambda)
points(vegHt, N)              # Add realized abundance to plot
table(N)

# Plot the true system state (Fig. 6???2, left)
par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
plot(vegHt, N, xlab="Vegetation height", ylab="True abundance (N)", frame = F, cex = 1.5)
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")


# Create a covariate called wind
wind <- array(runif(M * J, -1, 1), dim = c(M, J))

# Choose parameter values for measurement error model and compute detectability
alpha0 <- -2                        # Logit-scale intercept
alpha1 <- -3                        # Logit-scale slope for wind
p <- plogis(alpha0 + alpha1 * wind) # Detection probability
#plot(p ~ wind, ylim = c(0,1))       # Look at relationship

# Take J = 3 abundance measurements at each site
for(j in 1:J) {
  C[,j] <- rbinom(M, N, p[,j])
}

# Plot observed data and effect of wind on det. probability (Fig. 6???2, middle)
plot(wind, C/max(C), xlab="Wind", ylab="Scaled counts: C/max(C)", frame = F, cex = 1.5)
lines(seq(-1,1,,100), plogis(alpha0 + alpha1*seq(-1,1,,100)), lwd=3, col="red")


# Expected (lambda) and realized abundance (N) and measurements (C)
cbind(lambda=round(lambda,2), N=N, C1=C[,1], C2=C[,2], C3=C[,3])


# Create factors
time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # assumes M = 100


# Load unmarked, format data in unmarked data frame and summarize
library(unmarked)
umf <- unmarkedFramePCount(
  y=C,                                            # Counts matrix
  siteCovs= data.frame(vegHt = vegHt, hab = hab), # Site covariates
  obsCovs = list(time = time, wind = wind))       # Observation covs
summary(umf)


# Fit model and extract estimates
# linear model for p follows first tilde, then comes linear model for lambda
summary(fm.Nmix1 <- pcount(~wind ~vegHt, data=umf, control=list(trace=T, REPORT=1)))


fm.Nmix2 <- pcount(~wind ~vegHt, data=umf, mixture="NB", control=list(trace=TRUE, REPORT=5))
fm.Nmix3 <- pcount(~wind ~vegHt, data=umf, mixture="ZIP", control=list(trace=TRUE, REPORT=5))
cbind(AIC.P=fm.Nmix1@AIC, AIC.NB=fm.Nmix2@AIC, AIC.ZIP=fm.Nmix3@AIC)


# Predictions of lambda for specified values of vegHt, say 0.2 and 2.1
newdat <- data.frame(vegHt=c(0.2, 1.1))
predict(fm.Nmix1, type="state", newdata=newdat, append = T)


# ... or of p for values of wind of -1 to 1
newdat <- data.frame(wind=seq(-1, 1, , 5))
predict(fm.Nmix1, type="det", newdata=newdat, append = T)


# Predict lambda and detection for actual data set
(lambda.hat <- predict(fm.Nmix1, type="state"))     # lambda at every site
(p.hat <- predict(fm.Nmix1, type="det"))            # p during every survey


# Predict lambda and detection as function of covs
newdat <- data.frame(vegHt=seq(-1, 1, 0.01))
pred.lam <- predict(fm.Nmix1, type="state", newdata=newdat)
newdat <- data.frame(wind=seq(-1, 1, 0.1))
pred.det <- predict(fm.Nmix1, type="det", newdata=newdat)


# Fit detection-naive GLM to counts and plot comparison (Fig. 6???2, right)
summary(fm.glm <- glm(c(C) ~ rep(vegHt, 3), family=poisson)) # p-naive  model
matplot(vegHt, C, xlab="Vegetation height", ylab="Counts", frame = F, cex = 1.5, pch = 1, col = "black")
lines(seq(-1,1,,100), exp(beta0 + beta1* seq(-1,1,,100)), lwd=3, col = "red")
curve(exp(coef(fm.glm)[1]+coef(fm.glm)[2]*x), -1, 1, type ="l", lwd=3, add=TRUE)
lines(vegHt, predict(fm.Nmix1, type="state")[,1], col = "blue", lwd = 3)
legend(-1, 7, c("Truth", "'Poisson GLM' with p", "Poisson GLM without p"), col=c("red", "blue", "black"), lty = 1, lwd=3, cex = 1.2)


ranef(fm.Nmix1)

# calculate lambda.hat: exp(a0 + a1*vegHt)
lambda.hat <- predict(fm.Nmix1, type="state")[,1]

# calculate p.hat: plogis(b0 + b1*wind)
p.hat <- matrix(predict(fm.Nmix1, type="det")[,1], ncol=ncol(C), byrow=TRUE)

Ngrid <- 0:(100+max(umf@y, na.rm = TRUE))
posterior <- matrix(NA, nrow=nrow(C), ncol=length(Ngrid))
bup2 <- array(NA, dim = M)

for(i in 1:nrow(C)){      # Loop over sites
  # Compute prior using MLE
  gN <- dpois(Ngrid, lambda.hat[i])
  gN <- gN/sum(gN)
  
  # Compute likelihood for each possible value of N
  fy <- rep(NA, length(Ngrid))
  for(j in 1:length(Ngrid)){
    fy[j]<- prod(dbinom(C[i,], Ngrid[j], p.hat[i,]))
  }
  
  # Compute marginal of y. for denominator of Bayes rule
  qy <- sum(fy * gN)
  
  # Posterior
  posterior[i,] <- fy * gN / qy
  
  # N can't be less than max(C)
  if(max(C[i,] > 0))
    posterior[i,0:max(C[i,])]<- 0
  
  # Compute posterior mean (BUP)
  bup2[i] <- sum(posterior[i,] * Ngrid)
}

# Compare BUPS with true N and counts for first and last 5 sites
(bup1 <- bup(ranef(fm.Nmix1)))
cbind(N=N,count1=C[,1],count2=C[,2],count3=C[,3],BUP1=bup1,BUP2=bup2)[c(1:5, 91:95),]


plot(ranef(fm.Nmix1), xlim = c(0,12))[sort(sample(1:100, 12))]


# Main-effects ANCOVA: additive effects of factor and covariate
summary(fm.Nmix2 <- pcount(~ wind+time-1 ~ vegHt+hab-1, data=umf))

# Interaction-effects ANCOVA: multiplicative effects of factor and covariate
summary(fm.Nmix3 <- pcount(~ wind*time-1-wind ~ vegHt*hab-1-vegHt, data=umf))

# Get predictions for factor levels at average values of covariates
newdat <- data.frame(vegHt=0, hab = c("A", "B", "C"))
predict(fm.Nmix2, type="state", newdata=newdat, appendData = T)


newdat <- data.frame(time = c("1", "2", "3"), wind = 0)
predict(fm.Nmix3, type="det", newdata=newdat, appendData = T)


newdat <- data.frame(vegHt=seq(0, 2 ,by = 0.1), hab = factor("A", levels = c("A", "B", "C")))
predict(fm.Nmix2, type="state", newdata=newdat, appendData = T)


LRT(fm.Nmix3, fm.Nmix1)


# Bundle data
win.data <- list(C = C, M = nrow(C), J = ncol(C), wind = wind, vegHt = vegHt, hab = as.numeric(factor(hab)), XvegHt = seq(-1, 1,, 100), Xwind = seq(-1, 1,,100) )
str(win.data)

# Specify model in BUGS language
cat(file = "model2.txt", "
    model {
    # Priors
    for(k in 1:3){                # Loop over 3 levels of hab or time factors
    alpha0[k] ~ dunif(-10, 10) # Detection intercepts
    alpha1[k] ~ dunif(-10, 10) # Detection slopes
    beta0[k] ~ dunif(-10, 10)  # Abundance intercepts
    beta1[k] ~ dunif(-10, 10)  # Abundance slopes
    }
    
    # Likelihood
    # Ecological model for true abundance
    for (i in 1:M){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0[hab[i]] + beta1[hab[i]] * vegHt[i]
    # Some intermediate derived quantities
    critical[i] <- step(2-N[i])# yields 1 whenever N is 2 or less
    z[i] <- step(N[i]-0.5)     # Indicator for occupied site
    # Observation model for replicated counts
    for (j in 1:J){
    C[i,j] ~ dbin(p[i,j], N[i])
    logit(p[i,j]) <- alpha0[j] + alpha1[j] * wind[i,j]
    }
    }
    
    # Derived quantities
    Nocc <- sum(z[])         # Number of occupied sites among sample of M
    Ntotal <- sum(N[])       # Total population size at M sites combined
    Nhab[1] <- sum(N[1:33])  # Total abundance for sites in hab A
    Nhab[2] <- sum(N[34:66]) # Total abundance for sites in hab B
    Nhab[3] <- sum(N[67:100])# Total abundance for sites in hab C
    for(k in 1:100){         # Predictions of lambda and p ...
    for(level in 1:3){    #    ... for each level of hab and time factors
    lam.pred[k, level] <- exp(beta0[level] + beta1[level] * XvegHt[k])
    logit(p.pred[k, level]) <- alpha0[level] + alpha1[level] * Xwind[k]
    }
    }
    N.critical <- sum(critical[]) # Number of populations with critical size
    }")

# Initial values
Nst <- apply(C, 1, max)+1   # Important to give good inits for latent N
inits <- function() list(N = Nst, alpha0 = rnorm(3), alpha1 = rnorm(3), beta0 = rnorm(3), beta1 = rnorm(3))

# Parameters monitored
params <- c("alpha0", "alpha1", "beta0", "beta1", "Nocc", "Ntotal", "Nhab", "N.critical", "lam.pred", "p.pred") # could also estimate N, bayesian counterpart to BUPs before: simply add "N" to the list

# MCMC settings
nc <- 3   ;   ni <- 22000   ;   nb <- 2000   ;   nt <- 10

# Call JAGS, time run (ART 1 min) and summarize posteriors
system.time(out <- jags(win.data, inits, params, "model2.txt", n.chains = nc,
                        n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
)
traceplot(out, param = c('alpha0', 'alpha1', 'beta0', 'beta1', 'Nocc', 'Ntotal', 'Nhab', 'N.critical'))
print(out, 2)


plot(table(out$sims.list$N.critical), xlab="Number of populations with critical size", ylab="Frequency", frame = F)
abline(v = 74.5, col = "red", lwd = 3)


(metapop.extinction.risk <- mean(out$sims.list$N.critical > 74))


par(mfrow = c(1,2), mar = c(5,5,3,2), cex.axis = 1.5, cex.lab = 1.5)
X <- seq(-1, 1,, 100)
plot(X, out$summary[219:318,1], xlab = "Vegetation Height", ylab = "Expected abundance (lambda)", ylim = c(0, 11), frame = F, type = "l")
polygon(c(X, rev(X)), c(out$summary[219:318,3], rev(out$summary[219:318,7])), col = "gray", border = F)
lines(X, out$summary[219:318,1], lty = 1, lwd = 3, col = "blue")
plot(X, out$summary[519:618,1], xlab = "Wind speed", ylab = "Detection probability (p)", ylim = c(0, 1), frame = F, type = "l")
polygon(c(X, rev(X)), c(out$summary[519:618,3], rev(out$summary[519:618,7])), col = "gray", border = F)
lines(X, out$summary[519:618,1], lty = 1, lwd = 3, col = "blue")



# 6.5 A very general data simulation function for N-mixture models: simNmix
# -------------------------------------------------------------------------


# Execute function and inspect results
data <- simNmix()                   # Default arguments
data <- simNmix(show.plot = FALSE)  # Default args, no plots
set.seed(24)
str(data <- simNmix(nsite = 267, nvisit = 3, mean.theta = 1, mean.lam = 2, mean.p = 0.6, area = FALSE, beta1.theta = 0, beta2.theta = 0, beta3.theta = 0, beta2.lam = 0, beta3.lam = 0, beta4.lam = 0, beta3.p = 0, beta5.p = 0, beta6.p = 0, beta.p.survey = 0, beta.p.N = 0, sigma.lam = 0, dispersion = 10, sigma.p.site = 0, sigma.p.visit = 0, sigma.p.survey = 0, sigma.p.ind = 0, Neg.Bin = FALSE, open.N = FALSE, show.plot = TRUE))  # All default args explicit


str(data <- simNmix())                  # Null data-generating model
str(data <- simNmix(mean.theta = 0.60)) # ZIP with 40% structural zeroes
str(data <- simNmix(sigma.lam = 1))     # Poisson-lognormal (PLN) mixture
str(data <- simNmix(Neg.Bin = TRUE))    # Negative-binomial mixture
str(data <- simNmix(mean.theta = 0.6, sigma.lam = 1))  # Zero-inflated PLN
str(data <- simNmix(mean.theta = 0.6, Neg.Bin = TRUE)) # Zero-infl. NegBin
str(data <- simNmix(mean.p = 1))        # Perfect detection (p = 1)
str(data <- simNmix(mean.theta = 0.6, mean.p = 1))     # ZIP with p = 1
str(data <- simNmix(sigma.lam = 1, mean.p = 1))        # PLN with p = 1


areas <- runif(267, 1, 2)                # Generate vector with site area
str(data <- simNmix(nsite = 267, area = areas)) # Sites with variable area
str(data <- simNmix(nvisit = 1))         # Only one visit
str(data <- simNmix(sigma.p.site = 1))   # Random site effects in p
str(data <- simNmix(sigma.p.visit = 1))  # Random visit (= time) effects in p
str(data <- simNmix(sigma.p.survey = 1)) # Random site-by-visit effects in p
str(data <- simNmix(sigma.p.ind = 1))    # Random individual effects in p
str(data <- simNmix(mean.theta = 0.5, beta1.theta = 1)) # Site cov 1 in suit.
str(data <- simNmix(beta2.lam = 1))      # Site covariate 2 in abundance process
str(data <- simNmix(beta3.p = 1))        # Site covariate 3 in detection process
str(data <- simNmix(beta.p.N = 1))       # Positive density-dep. in p
str(data <- simNmix(beta.p.N = -1))      # Negative density-dep. in p
# Same covariate in suitab. and abund. (see Phillips & Elith, Ecology, 2014 !)
str(data <- simNmix(mean.theta = 0.5, beta2.theta = 1, beta2.lam = -1))
# Same covariate in abundance and detection (see K???ry, Auk, 2008)
str(data <- simNmix(beta3.lam = 1, beta3.p = -1))
# Same covariate in all three levels of model (ouch !)
str(data <- simNmix(mean.theta = 0.5, beta3.theta = 1, beta3.lam = 1, beta3.p = -1))


# Use unmarked to fit some models to these data sets
cov <- data$site.cov
summary(umf <- unmarkedFramePCount(
  y=data$C, siteCovs= data.frame(cov1=cov[,1], cov2=cov[,2], cov3=cov[,3],
                                 cov4=cov[,4], cov5=cov[,5], cov6=cov[,6], area = data$area),
  obsCovs = list(survey.cov = data$survey.cov)))
summary(fm <- pcount(~1 ~1, umf))
summary(fm <- pcount(~1 ~1, umf, mixture = "ZIP"))
summary(fm <- pcount(~1 ~1, umf, mixture = "NB"))
summary(fm <- pcount(~cov1+cov2+cov3 ~ cov1+cov2+cov3, umf))



# 6.6 Study design, and bias and precision of the N-mixture estimator
# -------------------------------------------------------------------


# Define simulation settings and arrays for sim results
simreps <- 1000                    # Simulate and analyse 1000 data sets
nsites <- c(20, 120, 250)          # Levels for nsites factor
nreps <- c(2, 5, 10)               # Levels of nrep factor
estimates <- array(NA, dim = c(2, simreps, 3, 3))

# Fill p with random numbers between 0.01 and 0.99
p <- array(runif(n=simreps*3*3, 0.01, 0.99), dim = c(simreps, 3, 3))

# Launch simulation (takes about 6.3 hours)
for(s in 1:3){                     # Loop over levels of nsites factor
  for(r in 1:3){                   # Loop over levels of nreps factor
    for(i in 1:simreps){           # Simulate and analyse 1000 data sets
      cat("*** Simrep number", i, "***\n")
      data <- simNmix(nsite=nsites[s], nvisit=nreps[r], mean.lam = 5,
                      mean.p=p[i,s,r], show.plot = F)        # Generate data set
      umf <- unmarkedFramePCount(y = data$C) # Bundle data for unmarked
      fm <- pcount(~1 ~1, umf)               # Fit model
      estimates[,i,s,r] <- coef(fm)          # Save estimates
    }
  }
}

# Visualisation
par(mfrow = c(3,3), mar = c(4.5,4.5,2,2), cex.lab = 1.5, cex.axis = 1.3)
for(s in 1:3){                     # Loop over nsites
  for(r in 1:3){                   # Loop over nreps
    plot(p[,s,r], exp(estimates[1,,s,r]), xlab = "Detection probability",
         ylab = "lambda_hat", main = "", ylim = c(0, 75), frame = F)
    text(0.75, 60, paste("M = ", nsites[s], ", J = ", nreps[r], sep = ""),
         cex = 1.5)
    abline(h = 5, col = "red", lwd = 2)
    lines(smooth.spline(exp(estimates[1,,s,r])~p[,s,r]), col="blue", lwd=2)
  }
}



# 6.7 Study of some assumption violations using function simNmix
# --------------------------------------------------------------


simreps <- 1000                   # Number of data sets created/analysed
MLE <- array(dim = c(5, simreps)) # Array to hold MLEs

for(i in 1:simreps){              # Create and analyse 1000 data sets
  cat("*** Simrep number", i, "***\n")
  # Create data set with some extra (here: open populations)
  data <- simNmix(mean.lam=exp(1), mean.p=0.5, beta2.lam=1,
                  beta3.p=1, beta.p.survey=1, open.N=TRUE, show.plot=F)
  # Analyse data set with standard model (here: assuming closure)
  umf <- unmarkedFramePCount(y=data$C, siteCovs =
                               data.frame(cov2=data$site.cov[,2], cov3=data$site.cov[,3]),
                             obsCovs = list(survey.cov = data$survey.cov))
  fm <- pcount(~cov3+survey.cov ~cov2, umf, se = F)
  # Save MLEs
  MLE[,i] <- coef(fm)
}



# 6.8 Goodness of fit
# -------------------


# Case 1: Test GoF of correct model
library(AICcmodavg)
par(mfrow = c(3,3))
for(i in 1:9){
  data <- simNmix(show.plot = F)          # Create data set
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100) # 100 bootstrap reps
}


# Case 2: Simulate data with zero inflation and analyse without
val.range <- seq(0.1, 1,,9)            # Much to no zero-inflation
for(i in 1:9){
  data <- simNmix(mean.theta = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}

# Case 3: Extra-Poisson dispersion in lambda
val.range <- seq(1, 0,,9)            # Some to no extra-Poisson dispersion
for(i in 1:9){
  data <- simNmix(sigma.lam = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}

# Case 4: Site covariate in lambda
val.range <- seq(3, 0,,9)            # Strong to no effect of covariate
for(i in 1:9){
  data <- simNmix(beta3.lam = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}

# Case 5: Extra-binomial dispersion in p (survey random effect)
val.range <- seq(1, 0,,9)            # Strong to no effect extra-dispersion
for(i in 1:9){
  data <- simNmix(sigma.p.survey = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}

# Case 6: Site covariate in p
val.range <- seq(3, 0,,9)            # Strong to no covariate effect
for(i in 1:9){
  data <- simNmix(beta3.p = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}

# Case 7: Observational covariate in p
val.range <- seq(3, 0,,9)            # Strong to no covariate effect
for(i in 1:9){
  data <- simNmix(beta.p.survey = val.range[i], show.plot = F)
  fm <- pcount(~1 ~1, unmarkedFramePCount(y = data$C)) # Fit model
  pb.gof <- Nmix.gof.test(fm, nsim = 100)
}


# Bundle and summarize data set
str( win.data <- list(C = data$C, M = nrow(data$C), J = ncol(data$C), e = 0.001))

# Specify model in BUGS language
sink("model.txt")
cat("
    model {
    # Priors
    lambda ~ dgamma(0.001, 0.001)
    p ~ dunif(0, 1)
    # Likelihood
    for (i in 1:M) {
    N[i] ~ dpois(lambda)      # State model
    for (j in 1:J) {
    C[i,j] ~ dbin(p, N[i]) # Observation model
    }
    }
    # Posterior predictive distributions of chi2 discrepancy
    for (i in 1:M) {
    for (j in 1:J) {
    C.sim[i,j] ~ dbin(p, N[i]) # Create new data set under model
    e.count[i,j] <- N[i] * p   # Expected datum
    # Chi-square discrepancy for the actual data
    chi2.actual[i,j] <- pow((C[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
    # Chi-square discrepancy for the simulated ('perfect') data
    chi2.sim[i,j] <- pow((C.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
    # Add small value e to denominator to avoid division by zero
    }
    }
    # Add up individual chi2 values for overall fit statistic
    fit.actual <- sum(chi2.actual[,])  # Fit statistic for actual data set
    fit.sim <- sum(chi2.sim[,])        # Fit statistic for a fitting model
    c.hat <- fit.actual / fit.sim      # c-hat estimate
    bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
    }
    ",fill = TRUE)
sink()

# Do other preps and run model with JAGS
inits <- function(){list(N = apply(data$C, 1, max)+1)}
params <- c("lambda", "p", "fit.actual", "fit.sim", "c.hat", "bpv")
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3
fm <- jags(win.data, inits, params, "model.txt", n.chains = nc,
           n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm, dig = 3)

ppc.plot(fm)               # Produces Fig. 6???7



# 6.9 Abundance mapping of Swiss Great tits with unmarked
# ------------------------------------------------------------------------


# 6.9.1 Set up of the analysis
# ------------------------------------------------------------------------
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
?SwissTits
str(SwissTits)
SwissTits$species  # Available species

# Select Great tit and covariate data from 2013 and
#   drop 4 sites not surveyed in 2013
y0 <- SwissTits$counts[, , '2013', 'Great tit']
( NA.sites <- which(rowSums(is.na(y0)) == 3) ) # Unsurveyed sites
y <- y0[-NA.sites, ]                 # Drop them from the count data
tits <- SwissTits$sites[-NA.sites, ] # Also drop from the site covariates
str(y)  # Check the matrix of count data
# Get date and duration data for 2013, without the NA.sites rows:
date <- SwissTits$date[-NA.sites, , '2013']
dur <- SwissTits$dur[-NA.sites, , '2013']

# Plot observed data: counts vs survey date (Fig. 6-9)
matplot(t(date), t(y), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey data (1 = April 1)", ylab = "Count of Great Tits")

# Load unmarked, create unmarked data frame and inspect result
library(unmarked)
time <- matrix(rep(as.character(1:3), nrow(y)), ncol = 3, byrow = TRUE)
umf <- unmarkedFramePCount(y = y,
                           siteCovs=data.frame(elev=scale(tits[,"elev"]), forest=scale(tits[,"forest"]), iLength=1/tits[,"rlength"]),
                           obsCovs=list(time = time, date = scale(date), dur = scale(dur)))
summary(umf)                            # Summarize unmarked data frame
summary(apply(y, 1, max, na.rm = TRUE)) # Summarize max counts


# 6.9.2 Model fitting
# ------------------------------------------------------------------------
fm1 <- pcount(~ (elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
              ~ (elev+I(elev^2)) * (forest+I(forest^2)) + iLength,
              umf, control=list(trace=TRUE, REPORT=5))
summary(fm1)   ;   fm1@AIC


fm1.K500 <- pcount(fm1@formula, umf, control=list(trace=T, REPORT=5), K = 500)
summary(fm1.K500)   ;   fm1.K500@AIC


fm2 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
              - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
              - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
              - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
              ~ (elev+I(elev^2)) * (forest+I(forest^2))
              + iLength, starts = coef(fm1)[1:31],
              umf, control=list(trace=TRUE, REPORT=5))
summary(fm2)                      # AIC = 3695.792


fm3 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
              - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
              - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
              - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
              - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
              ~ (elev+I(elev^2)) * (forest+I(forest^2))
              + iLength, starts = coef(fm2)[-c(23, 27, 31)],
              umf, control=list(trace=TRUE, REPORT=5))
summary(fm3)                      # AIC = 3691.184


fm4 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
              - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
              - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
              - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
              - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
              - elev:I(date^2) - I(date^2):dur
              ~ (elev+I(elev^2)) * (forest+I(forest^2))
              + iLength, starts = coef(fm3)[-c(21, 28)],
              umf, control=list(trace=TRUE, REPORT=5))
summary(fm4)                      # AIC = 3687.565


fm5 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
              - elev:date:dur - elev:date:I(dur^2) - elev:I(date^2):dur
              - elev:I(date^2):I(dur^2) - I(elev^2):date:dur - I(elev^2):date:I(dur^2)
              - I(elev^2):I(date^2):dur - I(elev^2):I(date^2):I(dur^2)
              - I(elev^2):I(date^2) - I(elev^2):I(dur^2) - I(date^2):I(dur^2)
              - elev:I(date^2) - I(date^2):dur
              ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength
              - I(elev^2):forest - I(elev^2):I(forest^2),
              starts = coef(fm4)[-c(9:10)],
              umf, control=list(trace=TRUE, REPORT=5))
summary(fm5)                      # AIC = 3686.094


# Negative binomial (NB) mixture
fm5NB <- pcount(fm5@formula, starts = c(coef(fm5),0),
                umf, control=list(trace=TRUE, REPORT=5), mixture = "NB")
summary(fm5NB)                      # AIC = 3181.046

# Zero-inflated Poisson (ZIP) mixture
fm5ZIP <- pcount(fm5@formula, starts = c(coef(fm5),0),
                 umf, control=list(trace=TRUE, REPORT=5), mixture = "ZIP")
summary(fm5ZIP)                      # AIC = 3636.058


cbind(rbind("Poisson" = exp(coef(fm5)[1]), "NegBin" = exp(coef(fm5NB)[1]), "ZIP" = exp(coef(fm5ZIP)[1])), rbind(plogis(coef(fm5)[15:17]), plogis(coef(fm5NB)[15:17]), plogis(coef(fm5ZIP)[15:17])))



# 6.9.3 Model criticism and goodness of fit
# ------------------------------------------------------------------------
library(AICcmodavg)
system.time(gof.P <- Nmix.gof.test(fm5, nsim=100))      # 65 min
system.time(gof.NB <- Nmix.gof.test(fm5NB, nsim=100))   # 131 min
system.time(gof.ZIP <- Nmix.gof.test(fm5ZIP, nsim=100)) # 69 min
gof.P   ;   gof.NB   ;   gof.ZIP                        # print results


# Look at data, fitted values and residuals and produce plots
print(cbind(y, fitted(fm5), residuals(fm5)), 2)  # For Poisson model
plot_Nmix_resi(fm5, fm5NB, fm5ZIP)               # Produces Fig. 6???10 ## function renamed

plot_Nmix_resi <- function(fmP, fmNB, fmZIP){
  # Function does diagnostic plots for one Nmix model fitted with all three
  #   mixture distributions currently availabe in unmarked:
  #   Poisson, negative binomial and zero-inflated Poisson
  # For each, fitted values vs. observed data and
  #   residuals vs. fitted values are plotted.
  library(unmarked)
  
  # Plot fitted vs. observed data
  par(mfrow = c(2,3), mar = c(4,4,2,2), cex = 1.2)
  tmp1 <- range(c(fitted(fmP), fitted(fmNB), fitted(fmZIP)), na.rm = T)
  limits1 = round(c(tmp1[1], tmp1[2]))
  tmp2 <- range(c(residuals(fmP), residuals(fmNB), residuals(fmZIP)), na.rm = T)
  limits2 = round(c(tmp2[1], tmp2[2]))
  
  plot(fitted(fmP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (P)", frame = F, ylim = limits1)
  abline(0,1, lwd = 3 )
  abline(lm(c(fitted(fmP))~ c(fmP@data@y)), col = "blue", lwd = 3)
  plot(fitted(fmNB)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (NB)", frame = F, ylim = limits1)
  abline(0,1, lwd = 3)
  abline(lm(c(fitted(fmNB))~ c(fmP@data@y)), col = "blue", lwd = 3)
  plot(fitted(fmZIP)~ fmP@data@y, xlab = "Observed data", ylab = "Fitted values (ZIP)", frame = F, ylim = limits1)
  abline(0,1, lwd = 3)
  abline(lm(c(fitted(fmZIP)) ~ c(fmP@data@y)), col = "blue", lwd = 3)
  
  # Plot residuals vs. fitted values
  plot(residuals(fmP)~ fitted(fmP), xlab = "Fitted values (P)", ylab = "Residuals", frame = F, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmP)) ~ c(fitted(fmP))), col = "blue", lwd = 3)
  plot(residuals(fmNB)~ fitted(fmNB), xlab = "Fitted values (NB)", ylab = "Residuals", frame = F, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmNB)) ~ c(fitted(fmNB))), col = "blue", lwd = 3)
  plot(residuals(fmZIP)~ fitted(fmZIP), xlab = "Fitted values (ZIP)", ylab = "Residuals", frame = F, xlim = limits1, ylim = limits2)
  abline(h = 0, lwd = 2)
  abline(lm(c(residuals(fmZIP)) ~ c(fitted(fmZIP))), col = "blue", lwd = 3)
}


# Compute RMSE for all three models
(RMSEP <- sqrt(mean((y - fitted(fm5))^2, na.rm = TRUE)))      # Poisson
(RMSENB <- sqrt(mean((y - fitted(fm5NB))^2, na.rm = TRUE)))   # NB
(RMSEZIP <- sqrt(mean((y - fitted(fm5ZIP))^2, na.rm = TRUE))) # ZIP


map.Nmix.resi <- function(fm, x = tits$coordx, y = tits$coordy){
  # Function produces a map of the mean residuals from an N-mixture model
  #    object named fm, which was fit by function pcount in unmarked
  # Function arguments are the fitted model object and the x and y coordinates
  #    of every site
  library(sp)
  mean.resi <- apply(residuals(fm), 1, mean, na.rm = TRUE)
  mean.resi[mean.resi == "NaN"] <- mean(mean.resi, na.rm = TRUE)
  spdata <- data.frame(residuals = mean.resi, x = x, y = y)
  coordinates(spdata) <- c("x", "y")
  plot(bubble(spdata, "residuals", col = c("blue", "red"), main = paste("Average residuals of fitted N-mixture model")))
}

map.Nmix.resi(fm5, x = tits$coordx, y = tits$coordy)    # Map of average residuals for Poisson model ## function defaults changed in AHMbook 0.1.3
map.Nmix.resi(fm5NB, x = tits$coordx, y = tits$coordy)  # Map of average residuals for NB model
map.Nmix.resi(fm5ZIP, x = tits$coordx, y = tits$coordy) # Map of average residuals for ZIP model


binFittedP <- fitted(fm5) %/% 2.5                     # Bin fitted values
mMeanP <- tapply(fitted(fm5)^2, binFittedP, mean, na.rm = T)   # Mean mean
mVarP <- tapply(residuals(fm5)^2, binFittedP, mean, na.rm = T) # Mean variance
nsampleP <- table(binFittedP)                         # Sample size
binFittedNB <- fitted(fm5NB) %/% 2.5
mMeanNB <- tapply(fitted(fm5NB)^2, binFittedNB, mean, na.rm = T)
mVarNB <- tapply(residuals(fm5NB)^2, binFittedNB, mean, na.rm = T)
nsampleNB <- table(binFittedNB)
plot(mMeanP, mVarP, xlab = "Binned mean fitted response", ylab = "Mean variance of response", frame = F, cex = log(nsampleP), col = "grey", pch = 16, cex.lab = 1.5)
points(mMeanNB, mVarNB, cex = log(nsampleNB), col = "green", pch = 16)


# 6.9.4 Analysis of results
# ------------------------------------------------------------------------
# Two new data sets for prediction: for lambda and for p (200 data points each)
lamNewData <- data.frame(elev = (seq(200, 2250,,200) - mean(tits[,"elev"]))/ sd(tits[,"elev"]), forest = 0, iLength = 1/5.1)
pNewData <- data.frame(elev = 0, time = factor("2", levels = c("1", "2", "3")), dur = 0, date = (seq(1,90,,200) - mean(date, na.rm = T))/sd(date, na.rm = T))


# Predictions for lambda and p, with SE and CIs, but no overdispersion
predict(fm5, type="state", newdata=lamNewData)    # Poisson model
predict(fm5, type="det", newdata=pNewData)
predict(fm5NB, type="state", newdata=lamNewData)  # NegBin model
predict(fm5NB, type="det", newdata=pNewData)
predict(fm5ZIP, type="state", newdata=lamNewData) # ZIP
predict(fm5ZIP, type="det", newdata=pNewData)


# Predictions for lambda only, incl. SE and CIs, with overdispersion
predictSE(fm5ZIP, newdata=lamNewData, print.matrix = TRUE, type="response", parm.type = "lambda", c.hat = 2.47)


# Predictions for lambda and p, incl. SE and with overdispersion, but no CIs
# Poisson model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 3.82)
modavgPred(cand.set = list(fm5), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 3.82)    # Could be used to get 95% CIs
modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 3.82)
modavgPred(cand.set = list(fm5), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 3.82)      # Could be used to get 95% CIs

# NegBin model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 1.79)
modavgPred(cand.set = list(fm5NB), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 1.79)    # Could be used to get 95% CIs
modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 1.79)
modavgPred(cand.set = list(fm5NB), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 1.79)         # Could be used to get 95% CIs

# ZIP model, for natural and for link scale of both lambda and p
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda", type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=lamNewData, parm.type = "lambda", type = "link", c.hat = 2.47)     # Not yet implemented (May 2015)
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect", type = "response", c.hat = 2.47)
modavgPred(cand.set = list(fm5ZIP), newdata=pNewData, parm.type = "detect", type = "link", c.hat = 2.47)         # this works, so we could get 95% CIs


rlength <- seq(1, 30, 0.01)         # Vary route length from 1 to 30 kms
newData <- data.frame(elev=0, forest=0, iLength=1/rlength)
pred <- predictSE(fm5ZIP, parm.type="lambda", newdata=newData, c.hat = 2.47)
par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.3)
plot(rlength, pred[[1]], type = "l", lwd = 3, col = "blue", frame = F, xlab = "Transect length (km)", ylab = "Exposed population (lambda)", ylim = c(0, 16), axes = F)
axis(1, at = seq(2,30,2))       ;      axis(2)
abline(v = c(1.2, 5.135, 9.4), lwd = 2)
matlines(rlength, cbind(pred[[1]]-pred[[2]], pred[[1]]+pred[[2]]), type = "l", lty = 1, lwd = 2, col = "gray")


sat.pred <- predictSE(fm5ZIP, parm.type="lambda", newdata= data.frame(elev=0, forest=0, iLength=0), c.hat = 2.47)
abline(h = sat.pred$fit, lwd = 2, lty = 2)


# Inspect the numbers
print(cbind("Route length" = rlength, "Exp. pop" = pred[[1]], "Rel. exp. pop" = pred[[1]] / sat.pred$fit), 3)


# Create covariate vectors for prediction and standardise as in analysis
ep.orig <- seq(200, 2250, length.out=100)  # Elevation between 200 and 2250 m
(elev.mean <- mean(tits[,"elev"]))
(elev.sd <- sd(tits[,"elev"]))
ep <- (ep.orig - elev.mean) / elev.sd      # Standardized for prediction
fp.orig <- seq(0, 100, length.out=100)     # Forest cover between 0 and 100%
(forest.mean <- mean(tits[,"forest"]))
(forest.sd <- sd(tits[,"forest"]))
fp <- (fp.orig - forest.mean) / forest.sd  # Standardised for prediction
date.orig <- seq(1, 90, length.out=100)  # Survey date from 1 April - 1 July
(date.mean <- mean(date, na.rm = TRUE))
(date.sd <- sd(date, na.rm = TRUE))
datep <- (date.orig - date.mean) / date.sd # Standardised for prediction
dur.orig <- seq(90, 420, length.out=100)   # Survey duration from 90 - 420 min
(dur.mean <- mean(dur, na.rm = TRUE))
(dur.sd <- sd(dur, na.rm = TRUE))
durp <- (dur.orig - dur.mean) / dur.sd     # Standardised for prediction

# Do predictions along single covariate gradient
newData1 <- data.frame(elev=ep, forest=0, iLength=0, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))
pred1 <- predictSE(fm5ZIP, newdata=newData1, c.hat = 2.47)
pred2 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData1, parm.type = "detect", type = "response", c.hat = 2.47)
newData3 <- data.frame(elev=0, forest=fp, iLength=0, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))
pred3 <- predictSE(fm5ZIP, newdata=newData3, c.hat = 2.47)
newData4 <- data.frame(elev=0, forest=0, iLength=0, date=datep, dur=0, time = factor("2", levels = c("1", "2", "3")))
pred4 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData4, parm.type = "detect", type = "response", c.hat = 2.47)
newData5 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=durp, time = factor("2", levels = c("1", "2", "3")))
pred5 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData5, parm.type = "detect", type = "response", c.hat = 2.47)
newData6 <- data.frame(elev=0, forest=0, iLength=0, date=0, dur=0,
                       time = c("1", "2", "3"))
pred6 <- modavgPred(cand.set = list(fm5ZIP), newdata=newData6, parm.type = "detect", type = "response", c.hat = 2.47)

# Plot these predictions along single covariate gradient
par(mfrow = c(3,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
plot(ep.orig, pred1[[1]], type = "l", lwd = 2, col = "blue", xlab = "Elevation (m)", ylab = "Expected abundance", las = 1, ylim = c(0,50), frame = F)
matlines(ep.orig, cbind(pred1[[1]]-pred1[[2]], pred1[[1]]+pred1[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")
plot(fp.orig, pred3[[1]], type = "l", lwd = 2, col = "blue", xlab = "Forest cover (%)", ylab = "Expected abundance", las = 1, ylim = c(0, 18), frame = F)
matlines(fp.orig, cbind(pred3[[1]]-pred3[[2]], pred3[[1]]+pred3[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")
plot(ep.orig, pred2[[1]], type = "l", lwd = 2, col = "blue", xlab = "Elevation (m)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F)
matlines(ep.orig, cbind(pred2[[1]]-pred2[[2]], pred2[[1]]+pred2[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")
plot(date.orig, pred4[[1]], type = "l", lwd = 2, col = "blue", xlab = "Survey date (1 = April 1)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F)
matlines(date.orig, cbind(pred4[[1]]-pred4[[2]], pred4[[1]]+pred4[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")
plot(dur.orig, pred5[[1]], type = "l", lwd = 2, col = "blue", xlab = "Survey duration (min)", ylab = "Expected detection", las = 1, ylim = c(0,1), frame = F)
matlines(dur.orig, cbind(pred5[[1]]-pred5[[2]], pred5[[1]]+pred5[[2]]), type = "l", lty = 1, lwd = 1, col = "gray")
barplot(pred6[[1]], names.arg = c("1", "2", "3"), ylim = c(0,1),
        ylab = "Expected detection", xlab = "Survey Number")
segments(c(0.7,1.9, 3.1), pred6[[1]]-pred6[[2]], c(0.7,1.9, 3.1), pred6[[1]]+pred6[[2]], lwd = 2)


# Make predictions along two covariate gradients
# (1) Expected abundance (lambda) for forest and elevation
pred.matrix1 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(x=0, y=0, elev=ep[i], forest=fp[j], iLength=0)
    pred.matrix1[i,j] <- predict(fm5ZIP,type="state", newdata=newData)[1,1]
  }
}

# (2) Expected detection (p) for elevation and survey date
pred.matrix2 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], date=datep[j], dur=0,
                          time = factor("2", levels = c("1", "2", "3")))
    pred.matrix2[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1]
  }
}

# (3) Expected detection (p) for elevation and survey duration
pred.matrix3 <- array(NA, dim = c(100, 100))
for(i in 1:100){
  for(j in 1:100){
    newData <- data.frame(elev=ep[i], date=0, dur=durp[j],
                          time = factor("2", levels = c("1", "2", "3")))
    pred.matrix3[i,j] <- predict(fm5ZIP, type="det", newdata=newData)[1,1]
  }
}

# Plot these prediction matrices
par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

image(x=ep.orig, y=fp.orig, z=pred.matrix1, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Forest cover (%)")
contour(x=ep.orig, y=fp.orig, z=pred.matrix1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(0, 100, by = 10))
box()
title(main = "(A)", font.main = 2)
points(tits$elev, tits$forest, pch="+", cex=1.5)

image(x=ep.orig, y=date.orig, z=pred.matrix2, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Date (1 = April 1)")
contour(x=ep.orig, y=date.orig, z=pred.matrix2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(10, 120, by = 10))
box()
title(main = "(B)", font.main =2)
matpoints(tits$elev, date, pch="+", cex=1.5)

image(x=ep.orig, y=dur.orig, z=pred.matrix3, col = mapPalette(100), axes = F, xlab = "Elevation (m)", ylab = "Duration (min)")
contour(x=ep.orig, y=dur.orig, z=pred.matrix3, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
axis(1, at = seq(min(ep.orig), max(ep.orig), by = 250))
axis(2, at = seq(90, 420, by = 20))
box()
title(main = "(C)", font.main = 2)
matpoints(tits$elev, dur, pch="+", cex=1.5)


data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland


# Predictions for lambda, with overdispersion
newData <- data.frame(elev = (CH$elev-elev.mean)/elev.sd, forest = (CH$forest-forest.mean)/forest.sd, iLength = 0, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))
predCH <- predictSE(fm5ZIP, newdata=newData, print.matrix = TRUE, type="response", parm.type = "lambda", c.hat = 2.47)


max(predCH[,1])                 # Look at the max prediction  --- 43.8
sum(predCH[,1] > 60)            # How many are > 60 ?  --- none
plot(CH$elev, predCH[,1])       # Relationship with elevation
predCH[1][predCH[1] > 60] <- 60 # Censor freak predicions (not req'd here)

# Prepare Swiss coordinates and produce map
library(raster)
library(rgdal)

# Define a new dataframe with coordinates and outcome to be plotted
PARAM1 <- data.frame(x = CH$x, y = CH$y, z = predCH[,1])

# Convert the DataFrame into a raster object
r1 <- rasterFromXYZ(PARAM1)

# Create mask for elevation (mask areas > 2250 m)
elev <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)

# Create custom color palette
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Map expected abundance of great tits in Switzerland in 2013
par(mfrow = c(1,2), mar = c(1,1,2,4))
plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


# Prepare raster with prediction SEs
r2 <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = predCH[,2]))
elev <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elev[elev > 2250] <- NA
r2 <- mask(r2, elev)

# Map prediction SEs of expected abundance
plot(r2, col = mapPalette(100), axes = F, box = FALSE, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


# Predictions for p with overdispersion
newData <- data.frame(elev = (CH$elev-elev.mean)/elev.sd, date=0, dur=0, time = factor("2", levels = c("1", "2", "3")))
predCHp <- modavgPred(cand.set = list(fm5ZIP), newdata = newData, parm.type = "detect", type = "response", c.hat = 2.47)


(N <- sum(predCH[-which(CH$elev > 2250),1]))    # Nat'l population size

out <- which(CH$water > 50 | CH$elev > 2250)
(N <- sum(predCH[-out,1]))  # 'Terrestrial' population size < 2250 m elevation


# Remind ourselves of the relevant coefficients in the model
cbind(coef(fm5ZIP)[25])     # Zero-inflation model
cbind(coef(fm5ZIP)[1:8])    # Abundance model


pelev <- (CH$elev - elev.mean)/elev.sd
pforest <- (CH$forest - forest.mean)/forest.sd


Nhat <- function(fm = fm5ZIP, iLength = 0, area = 1) {
  betavec <- coef(fm)[1:8]
  DM <- cbind(rep(1,length(pelev)), pelev, pelev^2, pforest,
              pforest^2, rep(iLength,length(pelev)), pelev*pforest, pelev*pforest^2)
  pred <- exp(DM %*% betavec) * (1-plogis(coef(fm)[25])) * (1/area)
  pred2 <- pred
  N1 <- sum(pred[-which(CH$elev > 2250),])# Drop quads > 2250 m
  pred2[pred2 > 100] <- 100      # Censor freak-high estimates
  N2 <- sum(pred2[-which(CH$elev > 2250),]) # Estimate with freaks censored
  out <- which(CH$water > 50 | CH$elev > 2250)
  N3 <- sum(pred2[-out,])  # Estimate excluding water bodies
  return(c(N1 = N1, N2 = N2, N3 = N3))
}


(Nest <- Nhat(fm5ZIP))            # For default saturation density
(Nest <- Nhat(fm5ZIP, 1/10))      # For 10 km routes

# Launch the bootstrap (takes about 30h)
system.time(pb.N <- parboot(fm5ZIP, Nhat, nsim = 2500, report=5))


bs <- pb.N@t.star[,3]             # Extract the bootstrapped vals of N3

# Sample statistics
summary(bs)
quantile(bs, prob = c(0.025, 0.975))  # Get 95% CI

# Plot
hist(bs, breaks = 100, col = "grey", xlab = "National population size of great tits", cex = 1.5, main = "", xlim = c(400000, 1000000))
abline(v = Nest[3], lwd = 3)
abline(v = mean(bs), lwd = 3, col = "blue")
abline(v = quantile(bs, prob = c(0.025, 0.975)), lwd = 3, col = "blue", lty = 2)



(CI1 <- quantile(bs, prob = c(0.025, 0.975)))     # Percentile-based CI
(CI2 <- c(mean(bs-2*sd(bs)), mean(bs+2*sd(bs))))  # Normal approx.-based CI


c.hat <- gof.ZIP$c.hat.est
(CI <- c(mean(bs-2*sqrt(c.hat)*sd(bs)), mean(bs+2*sqrt(c.hat)*sd(bs))))


# 6.9.5 Conclusions on the analysis with unmarked
# ------------------------------------------------------------------------



# 6.10 The issue of space, or: what is your effective sample area ?
# ------------------------------------------------------------------------


AHR <- seq(0.001, 100, 0.1)         # Area home range (in km2): 0.1-100
RHR <- sqrt(AHR/pi)                 # Translate to home range radius (in km)
ESA <- (1 + 2*RHR)^2                # Eff. sample area for 1km2 nominal area
par(mfrow = c(1,2), mar = c(5,5,3,2), cex.lab = 1.3, cex.axis = 1.3)
plot(AHR, ESA, xlab = "Home range size (km2)", ylab = "Effective sample area (km2)", type = "l", lwd = 3, frame = F)
abline(h = 1, col = "red", lwd = 3) # Nominal sample area
abline(h = 0, col = "grey", lwd = 3)


GT.HR <- seq(0.003, 0.1,,1000)     # Great tit home ranges in km2
GT.rad <- sqrt(GT.HR/pi)           # Great tit home range radius (in km)
ESA.GT <- (1 + 2*GT.rad)^2         # Effective sample area for Great tit
NTOT <- numeric(length(ESA.GT))    # Adjusted national total population
for(i in 1:length(NTOT)){
  NTOT[i] <- Nhat(fm5ZIP, 0, ESA.GT[i])[3]
}
plot(100*GT.HR, NTOT, xlab = "Great tit home-range (ha)", ylab = "Adjusted national total (pairs)", type = "l", lwd = 3, frame = F)



# 6.11 Bayesian modeling of Swiss Great tits with BUGS
# ------------------------------------------------------------------------


# 6.11.1 Bayesian fitting of the basic ZIP N-mixture model
# ------------------------------------------------------------------------
# Prepare data for BUGS data bundle
elev <- umf@siteCovs$elev   ;   elev2 <- elev^2
forest <- umf@siteCovs$forest   ;   forest2 <- forest^2
date <- matrix(umf@obsCovs$date, ncol = 3, byrow = TRUE)
dur <- matrix(umf@obsCovs$dur, ncol = 3, byrow = TRUE)
date[is.na(date)] <- 0   ;   date2 <- date^2
dur[is.na(dur)] <- 0   ;   dur2 <- dur^2
iRoute <- umf@siteCovs$iLength

# Design matrix for abundance model (no intercept)
lamDM <- model.matrix(~ elev + elev2 + forest + forest2 + elev:forest + elev:forest2 + iRoute)[,-1]


# Specify model in BUGS language
sink("ZIPNmix.txt")
cat("
    model {
    
    # Specify priors
    # zero-inflation/suitability
    phi ~ dunif(0,1)          # proportion of suitable sites
    theta <- 1-phi            # zero-inflation (proportion of unsuitable)
    ltheta <- logit(theta)
    
    # abundance
    beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
    for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
    }
    tau.lam <- pow(sd.lam, -2)
    sd.lam ~ dunif(0, 2)      # site heterogeneity in lambda
    
    # detection
    for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1)# p intercept for occasions 1-3
    }
    for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
    }
    tau.p.site <- pow(sd.p.site, -2)
    sd.p.site ~ dunif(0, 2)   # site heterogeneity in p
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 2) # site-survey heterogeneity in p
    
    # ZIP model for abundance
    for (i in 1:nsite){
    a[i] ~ dbern(phi)
    eps.lam[i] ~ dnorm(0, tau.lam)       # Random site effects in log(abundance)
    loglam[i] <- beta0 + inprod(beta[], lamDM[i,]) + eps.lam[i] * hlam.on
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # ???Stabilize??? log
    lam[i] <- exp(loglam.lim[i])
    mu.poisson[i] <- a[i] * lam[i]
    N[i] ~ dpois(mu.poisson[i])
    }
    
    # Measurement error model
    for (i in 1:nsite){
    eps.p.site[i] ~ dnorm(0, tau.p.site) # Random site effects in logit(p)
    for (j in 1:nrep){
    y[i,j] ~ dbin(p[i,j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
    lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # ???Stabilize??? logit
    lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
    alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
    alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
    alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
    alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
    alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
    alpha[13] * date[i,j] * dur2[i,j] +
    eps.p.site[i] * hp.site.on + eps.p.survey[i,j] * hp.survey.on
    eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
    }
    # Posterior predictive distributions of chi2 discrepancy
    for (i in 1:nsite) {
    for (j in 1:nrep) {
    y.sim[i,j] ~ dbin(p[i,j], N[i]) # Create new data set under model
    e.count[i,j] <- N[i] * p[i,j]   # Expected datum
    # Chi-square discrepancy for the actual data
    chi2.actual[i,j] <- pow((y[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
    # Chi-square discrepancy for the simulated ('perfect') data
    chi2.sim[i,j] <- pow((y.sim[i,j]-e.count[i,j]),2) / (e.count[i,j]+e)
    # Add small value e to denominator to avoid division by zero
    }
    }
    # Add up individual chi2 values for overall fit statistic
    fit.actual <- sum(chi2.actual[,])  # Fit statistic for actual data set
    fit.sim <- sum(chi2.sim[,])        # Fit statistic for a fitting model
    bpv <- step(fit.sim-fit.actual)    # Bayesian p-value
    c.hat <- fit.actual/fit.sim        # c-hat estimate
    
    # Derived parameters: Total abundance at 263 sampled sites
    Ntotal263 <- sum(N[])
    }
    ",fill = TRUE)
sink()


# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, beta0 = 0, mean.p = rep(0.5,3), beta = runif(7, 0,0), alpha = runif(13, 0,0))}

# Parameters monitored
params <- c("theta", "ltheta", "phi", "beta0", "beta", "sd.lam", "alpha0", "mean.p", "alpha", "sd.p.site", "sd.p.survey", "fit.actual", "fit.sim", "bpv", "c.hat", "Ntotal263")


# Bundle data and choose to fit simple ZIP model (model 1)
win.data1 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 0,
                  hp.survey.on = 0)

# MCMC settings
ni <- 50000    ;    nt <- 4    ;    nb <- 10000    ;    nc <- 3

# Call WinBUGS from R (ART 93 min) and summarize posteriors
out1 <- bugs(win.data1, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
print(out1, dig = 3)


# Compare MLEs and Bayesian posterior means (order first): table and graph
tmp <- summary(fm5ZIP)
ord.MLE <- rbind(tmp$psi[,1:2], tmp$state[,1:2], tmp$det[c(7:9, 1:6, 10:16),1:2])
ord.Bayes <- out1$summary[-c(1,3,12,16:18,32:39), 1:2]
cbind(ord.MLE, ord.Bayes)
par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(ord.MLE[,1], ylim = c(-3,3), pch = 16, col = "black", main = "", frame = F, xlab = "Parameters (zero-inflation, abundance, detection)", ylab = "Parameter estimate (+/- 1 SE)", cex = 1.5)
segments(1:25, ord.MLE[,1]-ord.MLE[,2], 1:25, ord.MLE[,1]+ord.MLE[,2], lwd = 2)
abline(h = 0)
abline(v = c(1.5, 9.5), col = "grey")
points((1:25)+0.3, ord.Bayes[,1], pch = 16, col = "blue", cex = 1.5)
segments(1:25+0.3, ord.Bayes[,1]-ord.Bayes[,2], 1:25+0.3, ord.Bayes[,1]+ord.Bayes[,2], col = "blue", lwd = 2)


library(unmarked)
data(Switzerland)             # Load Swiss landscape data in unmarked
CH <- Switzerland


ELEV <- (CH$elev-elev.mean)/elev.sd
ELEV2 <- ELEV^2
FOREST <- (CH$forest-forest.mean)/forest.sd
FOREST2 <- FOREST^2
CHdata <- cbind(elev = ELEV, elev2 = ELEV2, forest = FOREST, forest2 = FOREST^2, iRoute = rep(0, length(CH$elev)), elev.forest = ELEV * FOREST, elev.forest2 = ELEV * FOREST2)
str(CHdata)                  # This is a design matrix


MCMCout <- out1              # Choose results output from model 1
(nsamp <- length(MCMCout$sims.list$theta))  # how many MCMC samples do we have ?

# Subsample
sub.sample.size <- 3000      # choose sample of 3000
selection <- sort(sample(1:nsamp, sub.sample.size))

# Array to hold predictions for every Swiss 1km2 quadrat
lamPred <- array(NA, dim =c(length(CH[,1]), sub.sample.size))

# Fill the array
for(i in 1:sub.sample.size){
  MCMCstep <- selection[i]
  lamPred[,i] <- (1-MCMCout$sims.list$theta[MCMCstep]) *
    exp(MCMCout$sims.list$beta0[MCMCstep] +
          CHdata %*% MCMCout$sims.list$beta[MCMCstep,1:7])
}


# Get posterior means for every quadrat and check if sensible
meanlam <- apply(lamPred, 1, mean)   # Get posterior mean
max(meanlam)                         # Check maximum
sum(meanlam > 100)                   # Are any predictions >100 ?
plot(CH$elev, meanlam, ylim = c(0, max(meanlam))) ; abline(v=2250, col="red")
plot(CH$elev, meanlam, ylim = c(0,100)) ; abline(v = 2250, col = "red", lwd = 2)
meanlam[meanlam > 100] <- 100        # censor all predictions
lamPred[lamPred > 100] <- 100

# Produce map of posterior mean of lambda
library(raster)
library(rgdal)
r <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = meanlam))
elevation <- rasterFromXYZ(cbind(CH$x, CH$y,CH$elevation))
elevation[elevation > 2250] <- NA
r <- mask(r, elevation)
mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
par(mfrow = c(1,2), mar = c(5,5,1,5))
plot(r, col = mapPalette(100), axes = F, box = F, main ="")
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


elev.class <- 100*(CH$elev %/% 100 + 1)      # elevation class of each km2
tmp <- aggregate(lamPred, by = list(elev.class), FUN = sum)
N.elev <- as.matrix(tmp[,-1])              # Posterior sample of Ntotal per band
band <- tmp[,1]                            # elevation band (in m)
meanN <- apply(N.elev, 1, mean)
barplot(meanN, col = "grey", horiz = T, xlab = "Number of Great tit territories", ylab = "Elevation band (100m)", xlim = c(0, 200000))
axis(2, at = 1:length(band), labels = band)


# Posterior distribution of total number of great tit territories in 2013
keep <- which((CH$water < 50) & (CH$elev < 2251))
Ntot <- apply(lamPred[keep,], 2, sum)
hist(Ntot, breaks = 100, col = "grey", main = "Posterior of national population size")

# Point estimate and 95% CRI
mean(Ntot)
quantile(Ntot, prob = c(0.025, 0.975))


# 6.11.2 Adding random effects in BUGS
# ------------------------------------------------------------------------

# 6.11.2.1 Accounting for overdispersion at multiple scales
# ------------------------------------------------------------------------
# MCMC settings
ni <- 10^6    ;    nt <- 80    ;    nb <- 200000    ;    nc <- 3

# Bundle data and select model 2
win.data2 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
                  hp.survey.on = 0)

# Call WinBUGS from R (ART 4050 min) and summarize posteriors
out2 <- bugs(win.data2, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out2, dig = 3)


# Bundle data and select model 3
win.data3 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 1,
                  hp.survey.on = 0)

# Call WinBUGS from R (ART 4200 min) and summarize posteriors
out3 <- bugs(win.data3, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out3, dig = 3)


# Bundle data and select model 4
win.data4 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 0, hp.site.on = 0,
                  hp.survey.on = 1)

# Call WinBUGS from R (ART 4020 min) and summarize posteriors
out4 <- bugs(win.data4, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out4, dig = 3)


# Bundle data and select model 5
win.data5 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
                  hp.survey.on = 0)

# Call WinBUGS from R (ART 4250 min) and summarize posteriors
out5 <- bugs(win.data5, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out5, dig = 3)


# Bundle data and select model 6
win.data6 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 0,
                  hp.survey.on = 1)

# Call WinBUGS from R (ART 4230 min) and summarize posteriors
out6 <- bugs(win.data6, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out6, dig = 3)


# Bundle data and select model 7
win.data7 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                  lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                  date2 = date2, dur2 = dur2, e = 1e-06, hlam.on = 1, hp.site.on = 1,
                  hp.survey.on = 1)

# Call WinBUGS from R (ART 4625 min) and summarize posteriors
out7 <- bugs(win.data7, inits, params, "ZIPNmix.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = F, bugs.directory = bugs.dir)
print(out7, dig = 3)

# Look at posteriors for random effects
MCMCout <- out7       # Choose which model you want to plot random effects
par(mfrow = c(1,3))
hist(MCMCout$sims.list$sd.lam, breaks = 60, col = "grey")
hist(MCMCout$sims.list$sd.p.site, breaks = 60, col = "grey")
hist(MCMCout$sims.list$sd.p.survey, breaks = 60, col = "grey")


# 6.11.2.2 Linear modeling of a variance in the N-mixture model
# ------------------------------------------------------------------------
# Bundle and summarize data set
str( win.data8 <- list(y = y, nsite = nrow(y), nrep = ncol(y),
                       lamDM = lamDM, elev = elev, date = date, dur = dur, elev2 = elev2,
                       date2 = date2, dur2 = dur2) )

# Specify model in BUGS language
sink("Nmix.special.txt")
cat("
    model {
    
    # Specify priors
    # abundance
    beta0 ~ dnorm(0, 0.1)     # log(lambda) intercept
    for(k in 1:7){            # Regression params in lambda
    beta[k] ~ dnorm(0, 1)
    }
    # Model for unexplained variance in lambda among sites
    for (i in 1:nsite){
    tau.lam[i] <- 1/var.lam[i]
    log(var.lam[i]) <- alpha.var.lam + beta.var.lam * elev[i]
    }
    # Priors for intercept and slope of linear model for variance
    alpha.var.lam ~ dunif(-1, 1)
    beta.var.lam ~ dunif(0, 3)
    
    # detection
    for(j in 1:3){
    alpha0[j] <- logit(mean.p[j])
    mean.p[j] ~ dunif(0, 1)# p intercept for occasions 1-3
    }
    for(k in 1:13){           # Regression params in p
    alpha[k] ~ dnorm(0, 1)
    }
    tau.p.survey <- pow(sd.p.survey, -2)
    sd.p.survey ~ dunif(0, 1) # site-survey heterogeneity in p
    
    # Poisson-lognormal model for abundance
    for (i in 1:nsite){
    eps.lam[i] ~ dnorm(0, tau.lam[i]) # Random site effects in log(abundance)
    loglam[i] <- beta0 + inprod(beta[], lamDM[i,]) + eps.lam[i]
    loglam.lim[i] <- min(250, max(-250, loglam[i]))  # ???Stabilize??? log
    mu.poisson[i] <- exp(loglam.lim[i])
    N[i] ~ dpois(mu.poisson[i])
    }
    
    # Binomial measurement error model with extra-binomial dispersion
    for (i in 1:nsite){
    for (j in 1:nrep){
    y[i,j] ~ dbin(p[i,j], N[i])
    p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
    lp.lim[i,j] <- min(250, max(-250, lp[i,j]))  # ???Stabilize??? logit
    lp[i,j] <- alpha0[j] + alpha[1] * elev[i] + alpha[2] * elev2[i] +
    alpha[3] * date[i,j] + alpha[4] * date2[i,j] +
    alpha[5] * dur[i,j] + alpha[6] * dur2[i,j] +
    alpha[7] * elev[i] * date[i,j] + alpha[8] * elev2[i] * date[i,j] +
    alpha[9] * elev[i] * dur[i,j] + alpha[10] * elev[i] * dur2[i,j] +
    alpha[11] * elev2[i] * dur[i,j] + alpha[12] * date[i,j] * dur[i,j] +
    alpha[13] * date[i,j] * dur2[i,j] + eps.p.survey[i,j]
    eps.p.survey[i,j] ~ dnorm(0, tau.p.survey) # Random site-survey effects
    }
    }
    }
    ",fill = TRUE)
sink()


# Initial values
Nst <- apply(y, 1, max, na.rm = T) + 1
Nst[is.na(Nst)] <- round(mean(y, na.rm = TRUE))
Nst[Nst == "-Inf"] <- round(mean(y, na.rm = TRUE))
inits <- function(){ list(N = Nst, beta0 = 0, mean.p = rep(0.5,3), beta = runif(7, 0,0), alpha = runif(13, 0,0), alpha.var.lam = 0, beta.var.lam = 1.5, sd.p.survey = 0.3)}

# Parameters monitored
params <- c("beta0", "beta", "alpha.var.lam", "beta.var.lam", "alpha0", "mean.p", "alpha", "sd.p.survey")

# MCMC settings
ni <- 180000    ;    nt <- 100    ;    nb <- 10000    ;    nc <- 3

# Call WinBUGS from R (ART 374 min) and summarize posteriors
out8 <- bugs(win.data8, inits, params, "Nmix.special.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
print(out8, dig = 3)


print(cbind(out6$summary[c(4:11,13:31,33),1:2], out8$summary[c(1:8,11:30),1:2]),2)


# Predict site variance as a function of elevation
# Get posterior distribution of predictions first
orig.elev.pred <- seq(200, 2250, 50)
elev.pred <- (orig.elev.pred - mean(tits$elev)) / sd(tits$elev)
(n.mcmc <- length(out8$sims.list$alpha.var.lam))  # how many MCMC samples ?
post.sd.lam <- array(NA, dim = c(length(elev.pred), n.mcmc))
for(i in 1:length(elev.pred)){
  post.sd.lam[i,] <- sqrt(exp(out8$sims.list$alpha.var.lam + out8$sims.list$beta.var.lam * elev.pred[i]))
}

# Plot posterior mean and a sample of 500 regression lines from posterior
show <- sample(1:n.mcmc, 500)
matplot(orig.elev.pred, post.sd.lam[,show], xlab = "Elevation (m)", ylab = " sd.lam", type = "l", lty = 1, lwd = 1, col = "grey", frame = F, ylim = c(0, 6))
lines(orig.elev.pred, apply(post.sd.lam, 1, mean), lwd = 3, col = "blue")



# 6.12 Time-for-space substitution
# ------------------------------------------------------------------------


simpleNmix(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1)

# Define function to simulate such data
simpleNmix <- function(nyear = 12, nrep = 4, beta0 = 2, beta1 = 0.1, alpha0 = 0.5, alpha1 = -0.1, alpha2 = 1){
  # Simple function simulates data under binomial N-mixture model where you have
  # a single site that is survyed over 'nyear' primary sampling periods
  # ('seasons', 'years'), within which there are 'nrep' secondary samples each
  # alpha0, alpha1 are the logit-linear coefficients of detection (p) on Time
  #    and on a survey-specific covariate such as temperature (temp).
  # beta0 and beta1 are the log-linear coefficients of expected abundance
  #   (lambda) on Time.
  
  Time <- 1:nyear
  temp <- matrix(runif(nyear*nrep, -2, 2), ncol = nrep)
  N <- rpois(n = nyear, lambda = exp(beta0 + beta1 * Time))
  C <- array(NA, dim = c(nyear, nrep))
  p <- plogis(alpha0 + alpha1*Time + alpha2*temp)
  for(j in 1:nrep){
    C[,j] <- rbinom(n = nyear, size = N, prob =p[,j])
  }
  par(mfrow = c(3, 2))
  curve(exp(beta0 + beta1 * x), 1, nyear, main = "Expected abundance (lambda) over time", frame = F, lwd = 2, ylab = "lambda", xlab = "Time")
  plot(Time, N, main = "Realized abundance (N) over time", frame = F)
  curve(plogis(alpha0 +alpha1 * x), 1, nyear, main = "p over time", frame = F, lwd = 2, xlab = "Time", ylab = "p (at averate temp)")
  matplot(Time, C, main = "Counts (C) over time", frame = F)
  curve(plogis(alpha0 + alpha2 * x), -2, 2, main = "p vs. Temperature", frame = F, lwd = 2, xlab = "Temperature", ylab = "p (at start of study)")
  matplot(temp, C, main = "Counts (C) over time", frame = F)
  
  return(list(nyear=nyear, nrep=nrep, beta0=beta0, beta1=beta1, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, N=N, C=C, Time=Time, temp = temp, p = p))
}


library(unmarked)
simrep <- 2500                  # Number of simreps
results <- array(NA, dim = c(simrep, 8)) # Array for results
for(i in 1:simrep){
  cat("Simrep", i, "\n")
  data <- simpleNmix(nyear = 12, nrep = 4) # Simulate a data set
  umf <- unmarkedFramePCount(y = data$C,
                             siteCovs = data.frame(Time = data$Time), obsCov = list(temp = data$temp))
  fm1 <- pcount(~Time+temp ~Time, data = umf)
  fm2 <- glm(c(data$C)~rep(data$Time,data$nrep)+c(data$temp),family='poisson')
  results[i, 1:5] <- coef(fm1)
  results[i, 6:8] <- coef(fm2)
}
colnames(results) <- c(names(coef(fm1)), names(coef(fm2)))


par(mfrow = c(2,2), mar = c(4,4,3,2))
hist(results[,1], breaks = 100, col = "grey", main = "lambda(Int) Nmix", xlim = c(0, 4))
abline(v = data$beta0, col = "red", lwd = 3)
abline(v = mean(results[,1]), col = "blue", lwd = 3)
hist(results[,2], breaks = 100, col = "grey", main = "lambda(Slope Time) Nmix", xlim = c(-0.1, 0.3))
abline(v = data$beta1, col = "red", lwd = 3)
abline(v = mean(results[,2]), col = "blue", lwd = 3)
hist(results[,6], breaks = 100, col = "grey", main = "lambda(Int) GLM", xlim = c(0, 4))
abline(v = data$beta0, col = "red", lwd = 3)
abline(v = mean(results[,6]), col = "blue", lwd = 3)
hist(results[,7], breaks = 100, col = "grey", main = "lambda(Slope Time) GLM", , xlim = c(-0.1, 0.3))
abline(v = data$beta1, col = "red", lwd = 3)
abline(v = mean(results[,7]), col = "blue", lwd = 3)



# 6.13 The Royle-Nichols model and other non-standard N-mixture models
# ------------------------------------------------------------------------


# 6.13.1 The Royle-Nichols or Poisson/Bernoulli N-mixture model
# ------------------------------------------------------------------------
playRN <- function(M = 267, J = 3, mean.abundance = 1, mean.detection = 0.3){
  # Function generates replicated count data under the Nmix model of Royle (2004),
  #   then 'degrades' the data to detection/nondetection and fits the RN model
  #   (Royle & Nichols 2003) using unmarked and estimates site-specific abundance.
  #   Requires function simNmix and package unmarked.
  #
  devAskNewPage(ask = FALSE)
  #
  # Simulate Nmix data under a range of abundance levels
  data <- simNmix(nsite = M, nvisit = J, mean.lam = mean.abundance, mean.p = mean.detection, beta2.lam = 1, beta3.p = -1, beta.p.survey = -1, show.plot = FALSE)
  # Turn counts into detection/nondetection data
  y <- data$C          # Copy counts C into y
  y[y>0] <- 1          # Turn counts >0 into 1
  # Load unmarked, format data and summarize
  library(unmarked)
  umf <- unmarkedFrameOccu(y=y, siteCovs= data.frame(cov2 = data$site.cov[,2], cov3 = data$site.cov[,3]), obsCovs = list(obscov = data$survey.cov))
  # Fit data-generating model
  fm <- occuRN(~cov3+obscov ~cov2, data=umf)
  # Estimate local abundance N and plot against true N (known in simulation)
  Nest <- bup(ranef(fm, K = ), "mean")
  par(mfrow = c(1,1))
  plot(data$N, Nest, xlab = "True local abundance", ylab = "Estimated local abundance", frame = F)
  abline(0,1, lwd = 3)                              # 1:1 line
  abline(lm(Nest ~ data$N), col = "blue", lwd = 3)  # Regression
  slope <- coef(lm(Nest ~ data$N))[2]               # Is 1 if model perfect
  return(list(nsite = M, nvisit = J, coef = coef(fm), slope = slope))
}

# Execute the function using various settings
playRN(M = 100, J = 3, mean.abundance = 0.1)  # Increasing abundance
playRN(M = 100, J = 3, mean.abundance = 1)
playRN(M = 100, J = 3, mean.abundance = 5)
playRN(M = 100, J = 3, mean.detection = 0.3)  # Increasing detection
playRN(M = 100, J = 3, mean.detection = 0.5)
playRN(M = 100, J = 3, mean.detection = 0.7)
playRN(M = 100, J = 20)                       # More visits
playRN(M = 1000, J = 3)                       # More sites


# Run simulation
lam <- c(0.1, 0.5, 1, 2.5, 5, 10)    # 6 levels of mean abundance
simrep <- 100
results <- array(NA, dim = c(length(lam), simrep, 6))
for(i in 1:6){
  for(j in 1:simrep){
    cat(paste("\n *** lambda level", lam[i], ", simrep", j, "***\n"))
    tmp <- playRN(mean.abundance = lam[i])
    results[i,j,1:5] <- tmp$coef   # Coefficients of RN model
    results[i,j,6] <- tmp$slope    # Slope of regression of Nest on Ntrue
  }
}

# Summary of results for abundance (Fig. 6-23)
par(mfrow = c(1, 2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
boxplot(t(exp(results[,,1])), names = as.character(lam), outline = F, frame = F, col = "grey", xlab = "Mean abundance (lambda)", ylab = "lambda intercept", ylim = c(0,10))
points(1:6, lam, pch = "*", col = "red", cex = 3)
boxplot(t(results[,,2]), names = as.character(lam), outline = F, frame = F, col = "grey", xlab = "Mean abundance (lambda)", ylab = "lambda slope", ylim = c(0.5, 1.5))
abline(h = 1, col = "red", lwd = 3)


# Load data on Swiss tits
## Code modified to use the SwissTits data set included in the AHMbook package
data(SwissTits)
?SwissTits
str(SwissTits)

# Select the count data for 2013 (all species)
y0 <- SwissTits$count[, , '2013', ]
str(y0)
# We keep the sites with count data, remove those with 3 NAs
# See which sites have counts in 2013 for (say) Great tits:
keep <- which(rowSums(is.na(y0[, , "Great tit"])) != 3)
length(keep)
y3D <- y0[keep, , ]

# Get covariate data (site and observational) and drop unsurveyed sites
elev <- SwissTits$sites$ele[keep]
route <- SwissTits$sites$rlength[keep]
forest <- SwissTits$sites$forest[keep]
date <- SwissTits$date[keep, , '2013']  # Survey date
dur <- SwissTits$dur[keep, , '2013']    # Survey duration

# 'Degrade' counts to mere detection/nondetection data
y3DRN <- y ; y3DRN[y3DRN > 1] <- 1  # Overwrite any count >1 with 1 (for RN model)

( spec.names <- paste0(SwissTits$species$name, "s") )

library(unmarked)

# Loop over 6 species of tits
par(mfrow = c(2,2), mar = c(5,4,3,1))
for(k in 1:6){
  cat("\n*** Analysis for ", spec.names[k], "***\n")
  # Plot observed data: counts vs survey date
  matplot(t(date), t(y3D[,,k]), type = "l", lwd = 3, lty = 1, frame = F, xlab = "Survey date (1 = April 1)", ylab = "Observed counts", main = paste("Counts of", spec.names[k], "as a function of survey date"))
  
  # Fit standard Nmix model (Nmix1)
  time <- matrix(rep(as.character(1:3), 263), ncol = 3, byrow = T)
  summary(umf1 <- unmarkedFramePCount(y = y3D[,,k], siteCovs=data.frame(elev=scale(elev), forest=scale(forest), iLength=1/route), obsCovs=list(time = time, date = scale(date), dur = scale(dur))) )
  Nmix1 <- pcount(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
                  ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength,
                  umf1, control=list(trace=TRUE, REPORT=5, maxit = 250))
  (tmp1 <- summary(Nmix1))
  
  # Fit RN model (Nmix2)
  summary(umf2 <- unmarkedFrameOccu(y = y3DRN[,,k], siteCovs=data.frame(elev=scale(elev), forest=scale(forest), iLength=1/route), obsCovs=list(time = time, date = scale(date), dur = scale(dur))))
  # Use solutions from Nmix1 as inits for Nmix2
  Nmix2 <- occuRN(~(elev+I(elev^2)) * (date+I(date^2)) * (dur+I(dur^2)) + time-1
                  ~ (elev+I(elev^2)) * (forest+I(forest^2))+ iLength,
                  umf2, control=list(trace=TRUE, REPORT=5), starts = coef(Nmix1))
  (tmp2 <- summary(Nmix2))
  
  # Compare estimates under both models
  # Table with MLEs and SEs
  print(cbind(rbind(tmp1$state[,1:2], tmp1$det[,1:2]), rbind(tmp2$state[,1:2], tmp2$det[,1:2])))
  
  # Plot of all RN estimates versus all Nmix estimates
  plot(coef(Nmix1), coef(Nmix2), xlab = "Coefficients Nmix", ylab = "Coefficients RN", main = spec.names[k])
  abline(0,1, lwd = 2)
  abline(h = 0, lwd = 1, col = "grey")
  abline(v = 0, lwd = 1, col = "grey")
  abline(lm(coef(Nmix2) ~ coef(Nmix1)), lwd = 2, col = "blue")
  browser()
  
  # Overall discrepancy measure (for state model only): slope of regression and r2
  print(slope <- coef(lm(coef(Nmix2)[1:10] ~ coef(Nmix1)[1:10]))[2]) # Slope
  print(r <- cor(coef(Nmix2)[1:10], coef(Nmix1)[1:10]))       # Correlation
}


# 6.13.2 The Poisson/Poisson N-mixture model
# ------------------------------------------------------------------------



# 6.14 Multi-scale N-mixture models
# ------------------------------------------------------------------------



# 6.15 Summary and outlook
# ------------------------------------------------------------------------



6.16 Exercises
# ------------------------------------------------------------------------


# Solution A:
range(mhbdata[,12:14], na.rm = TRUE)
day.mean <- mean(as.matrix(mhbdata[,12:14]), na.rm = TRUE)
day.sd <- sd(c(as.matrix(mhbdata[,12:14])), na.rm = TRUE)
original.pred.day <- 15:110
pred.day <- (original.pred.day - day.mean) / day.sd
new<- data.frame(day=pred.day)
pred<-predict(fm31,type="det",newdata=new,appendData=TRUE)
head(pred)

plot(Predicted ~ original.pred.day, pred,type="l",xlab="Date (1 = 1 April)", ylab="Expected detection prob",ylim=c(0,1), lwd = 2)
lines(lower ~ original.pred.day, pred,type="l",col="red", lwd = 2)
lines(upper ~ original.pred.day, pred,type="l",col="red", lwd = 2)










# =========================================================================
#
# 7. Modeling abundance using multinomial N-mixture models
#
# =========================================================================





# 7.1 Introduction
# ------------------------------------------------------------------------



# 7.2 Multinomial N-Mixture Models in ecology
# ------------------------------------------------------------------------


# 7.2.1 Covariate models
# ------------------------------------------------------------------------


# 7.2.2 Types of multinomial models
# ------------------------------------------------------------------------



# 7.3 Simulating multinomial observations in R
# ------------------------------------------------------------------------



rmultinom(10, 5, c(0.1, 0.2, 0.7))

(rmn <- sample(1:3, 5, replace=TRUE, prob = c(0.1, 0.2, 0.7)))

set.seed(2015)                         # Initialize RNG

# Simulate covariate values and local population size for each point
x <- rnorm(100)
N <- rpois(100, lambda=exp(-1 + 1*x) ) # Intercept and slope equal to 1
table(N)                               # Summarize

# Define detection probabilities (p) for both observers
p1 <- 0.8
p2 <- 0.6

# Construct the multinomial cell probabilities (pi)
cellprobs <- c(p1*p2, p1*(1-p2), (1-p1)*p2, (1-p1)*(1-p2))

# Create a matrix to hold the data
y <- matrix(NA, nrow=100, ncol=4)
dimnames(y) <- list(1:100, c("11", "10", "01", "00"))

# Loop over sites and generate data with function rmultinom()
for(i in 1:100){
  y[i,] <- rmultinom(1, N[i], cellprobs)
}

# Remove 4th column ("not detected") and summarize results
y <- y[,-4]
apply(y, 2, sum)


# Generate specific pseudo-random data set
set.seed(2014)
data <- simNmix(mean.lam = exp(1), beta3.lam = 1, mean.p = plogis(0),
                sigma.p.visit = 1, show.plot=FALSE)
str(data$DH)

# View detection histories for site with max abundance (here, N = 30)
t(data$DH[min(which(data$N == max(dim(data$DH)[3]))),,])


# Get detection history frequencies for each site (for exactly 3 surveys)
dhfreq <- array(NA, dim = c(data$nsite, 7),
                dimnames = list(NULL, c("100", "010", "001", "110", "101", "011", "111")))
for(i in 1:data$nsite){
  dhfreq[i,] <- table(factor(paste(data$DH[i,1,], data$DH[i,2,],
                                   data$DH[i,3,], sep = ""),
                             levels = c("100", "010", "001", "110", "101", "011", "111")))
}

head(dhfreq)      # Data for first 6 sites

# Get occasions with first detection of each individual
f <- apply(data$DH, c(1,3), function(x) min(which(x != 0)))
head(f)   ;   str(f)   ;   table(f)    # Inspect result

# Produce removal counts
y <- array(NA, dim = c(data$nsite, data$nvisit), dimnames = list(NULL,
                                                                 as.factor(1:data$nvisit)))
for(i in 1:data$nsite){
  y[i,] <- table(factor(f[i,], levels = as.character(1:data$nvisit)))
}
head(y)      # Data for first 6 sites

set.seed(24)
data <- simNmix(mean.lam = exp(0), beta2.lam = 1, beta3.lam =-1,
                beta4.lam = 0.2, dispersion = 1, mean.p = plogis(0),
                beta3.p = 1, sigma.p.visit = 1, Neg.Bin = TRUE)



# 7.4 Likelihood inference for multinomial N-mixture models
# ------------------------------------------------------------------------



# 7.5 Example 1: Bird point counts based on removal sampling
# ------------------------------------------------------------------------


data(ovendata)
ovendata.list$data[11:20,]   # Look at a snippet of data set

apply(ovendata.list$data,2,sum)  # Removals in occasion 1-4


# 7.5.1 Setting up the data for analysis
# ------------------------------------------------------------------------
library(unmarked)
data(ovendata)
ovenFrame <- unmarkedFrameMPois(y = ovendata.list$data,
                                siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
                                type = "removal")


# 7.5.2 Fitting models using function multinomPois
# ------------------------------------------------------------------------
# Fit models: multinomPois order of formulas: detection, abundance
fm0 <- multinomPois(~ 1 ~ 1, ovenFrame)
fm1 <- multinomPois(~ 1 ~ ufc, ovenFrame)
fm2 <- multinomPois(~ 1 ~ trba, ovenFrame)
fm3 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)
fm4 <- multinomPois(~ 1 ~ ufc + trba + ufc:trba, ovenFrame)
fm5 <- multinomPois(~ ufc ~ ufc + trba, ovenFrame)
fm6 <- multinomPois(~ ufc ~ ufc + trba + ufc:trba, ovenFrame)

# Rank models by AIC
ms <- fitList(
  "lam(.)p(.)"                                = fm0,
  "lam(ufc)p(.)"                              = fm1,
  "lam(trba)p(.)"                             = fm2,
  "lam(ufc+trba)p(.)"                         = fm3,
  "lam(ufc+trba+ufc:trba)p(.)"                = fm4,
  "lam(ufc+trba)p(ufc)"                       = fm5,
  "lam(ufc+trba+ufc:trba)p(ufc)"              = fm6)

(ms1 <- modSel(ms))

# Table with everything you could possibly need
coef(ms1)[,1:4]  # Only first 4 columns shown

output <- as(ms1, "data.frame")



# 7.5.3 Fitting models using function gmultmix
# ------------------------------------------------------------------------
ovenFrame <- unmarkedFrameGMM(ovendata.list$data,
                              siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
                              numPrimary=1,type = "removal")


fm0 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                data=ovenFrame)


# Fit Poisson models
fm1 <- gmultmix(~ ufc, ~ 1, ~  1, data = ovenFrame)
fm2 <- gmultmix(~ trba, ~ 1, ~ 1, data = ovenFrame)
fm3 <- gmultmix(~ ufc + trba, ~ 1, ~ 1, data = ovenFrame)
fm4 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, data = ovenFrame)
# Maybe p also depends on understory foliage?
fm5 <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, data = ovenFrame)
fm6 <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, data = ovenFrame)

# Fit analogous NegBin models
fm0nb <- gmultmix(~ 1, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm1nb <- gmultmix(~ ufc, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm2nb <- gmultmix(~ trba, ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm3nb <- gmultmix(~ ufc + trba , ~ 1, ~ 1, mixture = "NB", data = ovenFrame)
fm4nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ 1, mixture = "NB",
                  data = ovenFrame)
# maybe p also depends on understory foliage?
fm5nb <- gmultmix(~ ufc + trba, ~ 1, ~ ufc, mixture = "NB",
                  data = ovenFrame)
fm6nb <- gmultmix(~ ufc + trba + ufc:trba, ~ 1, ~ ufc, mixture = "NB",
                  data = ovenFrame)

# Rank models by AIC
gms <- fitList(
  "lam(.)p(.)"                                = fm0,
  "lam(ufc)p(.)"                              = fm1,
  "lam(trba)p(.)"                             = fm2,
  "lam(ufc+trba)p(.)"                         = fm3,
  "lam(ufc+trba+ufc:trba)p(.)"                = fm4,
  "lam(ufc+trba)p(ufc)"                       = fm5,
  "lam(ufc+trba+ufc:trba)p(ufc)"              = fm6,
  "NB,lam(.)p(.)"                             = fm0nb,
  "NB,lam(ufc)p(.)"                           = fm1nb,
  "NB,lam(trba)p(.)"                          = fm2nb,
  "NB,lam(ufc+trba)p(.)"                      = fm3nb,
  "NB,lam(ufc+trba+ufc:trba)p(.)"             = fm4nb,
  "NB,lam(ufc+trba)p(ufc)"                    = fm5nb,
  "NB,lam(ufc+trba+ufc:trba)p(ufc)"           = fm6nb)

(gms1 <- modSel(gms))

# Table with everything you could possibly need
output <- as(gms1, "data.frame")

# Summary results
gms1

fm2nb


print(coef(gms1), digits = 2)


# 7.5.4 Assessing model fit in unmarked
# ------------------------------------------------------------------------
set.seed(2015)
(gof <- parboot(fm2, fitstats, nsim = 1000, report = 1))



#  7.6 Bayesian Analysis in BUGS using the Conditional Multinomial (3-part model)
# ------------------------------------------------------------------------


# Set-up data with a missing value for element "not captured"
y[i,] <- c(y[i,],NA)

# Then, in BUGS, do this:
y[i,] ~ dmulti(probs[i,], N[i])
N[i] ~ dpois(lambda[i])


# Harvest the data and bundle it up for sending to BUGS
y <- as.matrix(getY(ovenFrame))
ncap <- apply(y, 1, sum)   # number of individuals removed per point
data <- list(y = y, M = nrow(y), n = ncap, X=as.matrix(siteCovs(ovenFrame)))
str(data)                  # Good practice to always inspect your BUGS data

# Write BUGS model
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    for(i in 1:M){ # Loop over sites
    # Conditional multinomial cell probabilities
    pi[i,1] <- p[i]
    pi[i,2] <- p[i]*(1-p[i])
    pi[i,3] <- p[i]*(1-p[i])*(1-p[i])
    pi[i,4] <- p[i]*(1-p[i])*(1-p[i])*(1-p[i])
    pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3] + pi[i,4])
    pcap[i] <- 1 - pi0[i]
    for(j in 1:4){
    pic[i,j] <- pi[i,j] / pcap[i]
    }
    
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    
    # Model specification, three parts:
    y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
    #    cell probabilities
    n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
    #    observed sample size
    N[i] ~ dpois(lambda[i])             # component 3 is the process model
    
    # log-linear model for abundance: UFC + TRBA + UFC:TRBA
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
    }
    }
    ",fill=TRUE, file="model.txt")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2)
}

# Parameters monitored

parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3")

# MCMC settings
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1

# Experience the power of BUGS and print posterior summary
library(R2WinBUGS)
bd <- "c:/Program Files/WinBUGS14/"
out <- bugs(data, inits, parameters, "model.txt", n.thin = nt,
            n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE,
            bugs.directory = bd)
print(out, 3)


# To make trace plots using CODA functionality
plot(as.mcmc.list(out), ask = TRUE)

library(jagsUI)
out <- jags(data, inits, parameters, "model.txt", n.thin=nt,
            n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, 3)

set.factory("bugs::Conjugate", FALSE, type="sampler")


# 7.6.1 Goodness-of-fit using Bayesian p-values
# ------------------------------------------------------------------------
parameters <- c("N", "p0", "beta0", "beta1", "beta2", "beta3",
                "fit1.data", "fit1.pred", "fit2.data", "fit2.pred")
out <- bugs (data, inits, parameters, "model.txt",n.thin=nt, n.chains=nc,
             n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)

mean(out$sims.list$fit1.pred > out$sims.list$fit1.data)

mean(out$sims.list$fit2.pred > out$sims.list$fit2.data)


# 7.6.2 Model selection in BUGS
# ------------------------------------------------------------------------
# Write BUGS model
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- log(p0/(1-p0))
    alpha1 ~ dnorm(0, 0.1)
    beta0 ~ dnorm(0, 0.1)
    beta1 ~ dnorm(0, 0.1)
    beta2 ~ dnorm(0, 0.1)
    beta3 ~ dnorm(0, 0.1)
    w1 ~ dbern(0.5)
    w2 ~ dbern(0.5)
    w3 ~ dbern(0.5)
    
    for(i in 1:M){
    # Conditional multinomial cell probabilities
    pi[i,1] <- p[i]
    pi[i,2] <- p[i] * (1-p[i])
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i])
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i])
    pi0[i] <- 1 - (pi[i,1] + pi[i,2] + pi[i,3]+ pi[i,4])
    pcap[i] <- 1 - pi0[i]
    for(j in 1:4){
    pic[i,j] <- pi[i,j] / pcap[i]
    }
    
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    
    # Model specification, 3 parts:
    y[i,1:4] ~ dmulti(pic[i,1:4], n[i]) # component 1 uses the conditional
    #    cell probabilities
    n[i] ~ dbin(pcap[i], N[i])          # component 2 is a model for the
    #    observed sample size
    N[i] ~ dpois(lambda[i])             # component 3 is the process model
    
    # log-linear model for abundance: effects of UFC and TRBA, with weights
    log(lambda[i])<- beta0 + w1*beta1*X[i,1] + w2*beta2*X[i,2] + w1*w2*w3*beta3*X[i,2]*X[i,1]
    }
    }
    ",fill=TRUE,file="model.txt")


# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2",
                "beta3", "w1", "w2", "w3")

# Initial values
set.seed(2015)
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0 = runif(1), N = ncap+2,
        w1=1, w2=1, w3=1)
}
# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "model.txt", n.thin = nt,
            n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE, bugs.dir=bd)
print(out, digits = 3)


# Extract the indicator variables for model selection
w1 <- out$sims.list$w1
w2 <- out$sims.list$w2

# Create a new w3 variable which takes on the value 1 if the
#     interaction is in the model, also requires w1 = 1 AND w2=1
w3 <- out$sims.list$w3 * w1 * w2

# Combine into a model indicator string and tabulate posterior frequencies
mod <- paste(w1, w2, w3, sep = "")
table(mod)


# 7.6.3 Poisson formulation of the multinomial mixture model
# ------------------------------------------------------------------------
# Specify model in BUGS language
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    for(i in 1:M){
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    # log-linear model for abundance: UFC + TRBA + UFC:TRBA
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1]
    
    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]
    
    for(j in 1:4){
    y[i,j] ~ dpois(pi[i,j])
    }
    # Generate predictions of N[i]
    N[i] ~ dpois(lambda[i])
    }
    }
    ",fill=TRUE,file="modelP.txt")

# Bundle up the data and inits
data <- list(y = y, M = nrow(y), X = as.matrix(siteCovs(ovenFrame)))
inits <- function(){
  list (p0 = runif(1), alpha1=runif(1), beta0=runif(1), beta1=runif(1), beta2=runif(1), beta3=runif(1))
}

# Define parameters to save and MCMC settings
parameters <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "N")
nc <- 3   ;   ni <- 6000   ;   nb <- 1000   ;   nt <- 1


# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "modelP.txt", n.thin=nt,
            n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE, bugs.dir=bd)
print(out, 3)


# Specify model in BUGS language
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    tau ~ dgamma(0.1,0.1)  # Excess-Poisson variation (precision)
    sigma <- sqrt(1 / tau)
    
    for(i in 1:M){
    # logit-linear model for detection: understory cover effect
    logit(p[i]) <- alpha0 + alpha1 * X[i,1]
    # Normal random effects
    eta[i] ~ dnorm(0, tau)  # 'residuals' for extra-Poisson noise
    # log-linear model for abundance: UFC + TRBA + UFC:TRBA + eta
    log(lambda[i])<- beta0 + beta1*X[i,1] + beta2*X[i,2] + beta3*X[i,2]*X[i,1] + eta[i]     # note 'extra-residual' for overdispersion
    
    # Poisson parameter = multinomial cellprobs x expected abundance
    pi[i,1] <- p[i] * lambda[i]
    pi[i,2] <- p[i] * (1-p[i]) * lambda[i]
    pi[i,3] <- p[i] * (1-p[i]) * (1-p[i]) * lambda[i]
    pi[i,4] <- p[i] * (1-p[i]) * (1-p[i]) * (1-p[i]) * lambda[i]
    
    for(j in 1:4){
    y[i,j] ~ dpois(pi[i,j])
    }
    # Generate predictions of N[i]
    N[i] ~ dpois(lambda[i])
    }
    }
    ",fill=TRUE,file="modelP.txt")

# Inits
inits <- function(){
  list (p0 = runif(1), alpha1=runif(1), beta0=runif(1), beta1=runif(1),
        beta2=runif(1), beta3=runif(1), tau = 1)  }

# Define parameters to save and MCMC settings
parameters <- c("p0", "alpha1", "beta0", "beta1", "beta2", "beta3", "sigma")
nc <- 3   ;   ni <- 32000   ;   nb <- 2000   ;   nt <- 1

# Call WinBUGS from R and summarize marginal posteriors
out <- bugs(data, inits, parameters, "modelP.txt", n.thin=nt,
            n.chains = nc, n.burnin = nb, n.iter = ni, debug = TRUE, bugs.dir = bd)
print(out, digits=3)

plot(density(out$sims.list$sigma), frame = F, main = "") # Fig. 7-4



# 7.7 Building custom multinomial models in unmarked
# ------------------------------------------------------------------------


# Removal model: capture probs for 5 sites, with 3 removal periods
(pRem <- matrix(0.5, nrow=5, ncol=3))

removalPiFun(pRem)   # Multinomial cell probabilities for each site

# Double observer model: capture probs for 5 sites, with 2 observers
(pDouble <- matrix(0.5, 5, 2))

doublePiFun(pDouble)  # Multinomial cell probabilities for each site


instRemPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  p[,1] <- pi[,1] <- 1 - (1 - p[,1])^2
  p[,2] <- 1 - (1 - p[,2])^3
  p[,3] <- 1 - (1 - p[,3])^5
  for(i in 2:J) {
    pi[,i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
  }
  return(pi)
}


instRemPiFun(pRem)

o2y <- matrix(1, 2, 3)
o2y



# 7.8 Spatially Stratified Capture-Recapture Models
# ------------------------------------------------------------------------


crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
  cbind("001" = (1 - p1) * (1 - p2) *      p3,
        "010" = (1 - p1) *      p2  * (1 - p3),
        "011" = (1 - p1) *      p2  *      p3,
        "100" =      p1  * (1 - p2) * (1 - p3),
        "101" =      p1  * (1 - p2) *      p3,
        "110" =      p1  *      p2  * (1 - p3),
        "111" =      p1  *      p2  *      p3)
}


p <- matrix(0.4, 2, 3)
crPiFun(p)


# To compute pi0 we do this:
(pi0 <- 1 - rowSums(crPiFun(p)))


# 7.8.1 Example 2: Fitting Models M0, Mt, and Mx to Chandler???s Flycatcher data
# ------------------------------------------------------------------------
alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl, 5)

alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",package="unmarked"),
                      row.names=1)
head(alfl.covs)

alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3,
                             sep="")
alfl$captureHistory <- factor(alfl$captureHistory,
                              levels=c("001", "010", "011", "100", "101", "110", "111"))
alfl$id <- factor(alfl$id, levels=rownames(alfl.covs))

alfl.v1 <- alfl[alfl$survey==1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
head(alfl.H1, 5)

intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
class(alfl.H1) <- "matrix"
o2y <- matrix(1, 3, 7)
umf.cr1 <- unmarkedFrameMPois(y=alfl.H1,
                              siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
                              obsCovs=list(interval=intervalMat), obsToY=o2y, piFun="crPiFun")
summary(umf.cr1)

M0 <- multinomPois(~ 1 ~ 1, umf.cr1)
Mt <- multinomPois(~ interval - 1 ~ 1, umf.cr1)
Mx <- multinomPois(~ time.1 ~ 1, umf.cr1)

(M0.woody <- multinomPois(~ 1 ~ woody, umf.cr1))

fl <- modSel(fitList(M0, Mt, Mx, M0.woody))

nd <- data.frame(woody=seq(0, 0.8, length=50))
E.abundance <- predict(M0.woody, type="state", newdata=nd, appendData=TRUE)
plot(Predicted ~ woody, E.abundance, type="l", ylim=c(0, 6),
     ylab="Alder flycatchers / plot", xlab="Woody vegetation cover",
     frame = F)
lines(lower ~ woody, E.abundance, col=gray(0.7))
lines(upper ~ woody, E.abundance, col=gray(0.7))


backTransform(M0.woody, type="det")

round(getP(M0.woody), 2)[1,]


# 7.8.2 Models with behavioral response
# ------------------------------------------------------------------------
crPiFun.Mb <- function(p) {
  pNaive <- p[,1]
  pWise <- p[,3]
  cbind("001" = (1 - pNaive) * (1 - pNaive) *      pNaive,
        "010" = (1 - pNaive) *      pNaive  * (1 - pWise),
        "011" = (1 - pNaive) *      pNaive  *      pWise,
        "100" =      pNaive  * (1 - pWise)  * (1 - pWise),
        "101" =      pNaive  * (1 - pWise)  *      pWise,
        "110" =      pNaive  *      pWise   * (1 - pWise),
        "111" =      pNaive  *      pWise   *      pWise)
}

behavior <- matrix(c('Naive', 'Naive', 'Wise'), 50, 3, byrow=TRUE)
umf.cr1Mb <- unmarkedFrameMPois(y=alfl.H1,
                                siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
                                obsCovs=list(behavior=behavior),  obsToY=o2y, piFun="crPiFun.Mb")
M0 <- multinomPois(~1 ~1, umf.cr1Mb)
(Mb <- multinomPois(~behavior-1 ~1, umf.cr1Mb))


# 7.8.3 Models with individual heterogeneity model
# ------------------------------------------------------------------------
parID <- matrix(c('p','sig','sig'), 50, 3, byrow=TRUE)
umf.cr2 <- unmarkedFrameMPois(y=alfl.H1,
                              siteCovs=alfl.covs[,c("woody", "struct", "time.1")],
                              obsCovs=list(parID=parID), obsToY=o2y, piFun="MhPiFun")



# 7.8.4 Bayesian analysis using data augmentation (DA): heterogeneity models in BUGS
# ----------------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400
y <- as.matrix(alfl[,c("interval1","interval2","interval3")] )
nind <- nrow(y)
M <- 400
y <- rbind(y, matrix(0, nrow=(M-nind), ncol=3))

# Site ID
# Make site ID into an integer: This only works ok here because sites are in
#  alphabetical order in the data set!
site <- as.numeric(alfl$id)
site <- c(site, rep(NA, M-nind))

# Next we extract the covariates and standardize them
sitecovs <- scale(as.matrix(alfl.covs[,c("woody", "struct", "time.1")]))

# Bundle data for BUGS
data <- list(y = y, J = 3, M = M , nsites = 50, X = sitecovs, group=site)
str(data)

# Specify model in BUGS language
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- log(p0 / (1-p0))    # same as logit(p0)
    alpha1 ~ dnorm(0, 0.01)
    alpha2 ~ dnorm(0,0.01)
    beta0 ~ dnorm(0,0.01)
    beta1 ~ dnorm(0,0.01)
    psi <- sum(lambda[]) / M   # psi is a derived parameter
    
    # log-linear model for abundance: lambda depends on WOODY
    for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X[s,1]
    probs[s] <- lambda[s] / sum(lambda[])
    }
    
    # Model for individual encounter histories
    for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends on 2 covariates: STRUCT + TIME
    for(j in 1:J){
    logit(p[i,j]) <- alpha0 + alpha1 * X[group[i],2] + alpha2*X[group[i],3]
    pz[i,j] <- p[i,j] * z[i]
    y[i,j] ~ dbern(pz[i,j])
    }
    }
    }
    ",fill=TRUE,file="model.txt")

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi")

# Initial values
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), alpha2 = rnorm(1),beta0=runif(1),
        beta1=rnorm(1), z= c( rep(1,100), rep(0, 300))) }

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

# Call JAGS from R and summarize marginal posteriors
library("jagsUI")
out <- jags(data, inits, parameters, "model.txt", n.thin = nt,
            n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, digits = 3)


# Specify model in BUGS language
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- log(p0/(1-p0))
    alpha1 ~ dnorm(0, 0.01)
    alpha2 ~ dnorm(0, 0.01)
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    psi <- sum(lambda[])/M
    tau ~ dgamma(0.1,0.1) # New parameter, precision of ind. random effects
    sigma <- 1/sqrt(tau)
    
    # log-linear model for abundance: lambda depends on WOODY
    for(s in 1:nsites){
    log(lambda[s])<- beta0 + beta1*X[s,1]
    probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Model for individual encounter histories
    for(i in 1:M){
    eta[i] ~ dnorm(alpha0, tau)  # Individual random effect
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    # Observation model: p depends on STRUCT + TIME + ind. heterogeneity
    for(j in 1:J){
    logit(p[i,j]) <- alpha1 * X[group[i],2] + alpha2*X[group[i],3] + eta[i]
    pz[i,j] <- p[i,j] * z[i]
    y[i,j] ~ dbern(pz[i,j])
    }
    }
    }
    ",fill=TRUE,file="model.txt")

# Parameters monitored: add sigma
parameters <- c("p0", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "psi", "sigma")

# Initial values: add tau
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), alpha2 = rnorm(1),beta0=runif(1),
        beta1=rnorm(1),z= c( rep(1,100), rep(0, 300)), tau = 1)  }

# MCMC settings (others as before)
ni <- 50000   ;   nb <- 10000   ;   nt <- 2   ;   nt <- 10

# Call JAGS from R and summarize marginal posteriors
out <- jagsUI(data, inits, parameters, "model.txt", n.thin = nt,
              n.chains = nc, n.burnin = nb, n.iter = ni)
print(out, digits = 3)


7.8.5 Example 2: analysis of data simulated with the simNmix data function
# ------------------------------------------------------------------------
# Fit model in unmarked
library(unmarked)
time <- matrix(as.character(1:3), data$nsite, 3, byrow = T)

# Define pifun for J=3 occasion capture-recapture protocol
crPiFun <- function(p) {
  p1 <- p[,1] # Extract the columns of the p matrix, one for
  p2 <- p[,2] # each of J = 3 sample occasions
  p3 <- p[,3]
  cbind( # define multinomial cell probabilities:
    "100" = p1 * (1-p2) * (1-p3),
    "010" = (1-p1) * p2 * (1-p3),
    "001" = (1-p1) * (1-p2) * p3,
    "110" = p1 * p2 * (1-p3),
    "101" = p1 * (1-p2) * p3,
    "011" = (1-p1) * p2 * p3,
    "111" = p1 * p2 * p3)
}

# Define mapping function for missing values
o2y <- matrix(1, 3, 7)

# Create unmarked frame and fit couple of models
umf <- unmarkedFrameMPois(y = dhfreq, siteCovs = data.frame(cov3 = data$site.cov[,3]), obsCovs = list(time = time), obsToY = o2y, piFun = "crPiFun")
fm1 <- multinomPois(~1 ~1, umf)    # detection model before abundance model
fm2 <- multinomPois(~time-1 ~1, umf)
fm3 <- multinomPois(~1 ~cov3, umf)
fm4 <- multinomPois(~time-1 ~cov3, umf)

# Assemble the models into a fitList and rank using AIC
ms <- fitList(
  "lam(.)p(.)" = fm1,
  "lam(.)p(time)" = fm2,
  "lam(cov3)p(.)" = fm3,
  "lam(cov3)p(time)" = fm4)

(AICtable <- modSel(ms))

summary(fm4)
p.true <- qlogis(data$p[min(which(data$N>0)),,1])
tmp <- cbind(rbind(lam0 = log(data$mean.lam), beta3 = data$beta3.lam, logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3]), coef(fm4))
colnames(tmp) <- c("Truth", "MLEs")
tmp


# No temporal variation in p and effect of site-covariate on lambda
set.seed(24)
data <- simNmix(mean.lam = exp(1), beta2.lam = 1, beta3.lam = 1, mean.p = plogis(0.2), beta3.p = -1, beta.p.survey = -1)
str(data$DH)

# Get occasions with first detection of each individual
f <- apply(data$DH, c(1,3), function(x) min(which(x != 0)))
head(f)   ;   str(f)

# Produce removal counts (for any number of occasions)
y <- array(NA, dim = c(data$nsite, data$nvisit), dimnames = list(NULL, as.factor(1:data$nvisit)))
for(i in 1:data$nsite){
  y[i,] <- table(factor(f[i,], levels = as.character(1:data$nvisit)))
}
y                # Look at removal data set

# Fit models in unmarked
summary(umf <- unmarkedFrameMPois(y = y, siteCovs = data.frame(cov2 = data$site.cov[,2], cov3 = data$site.cov[,3]), obsCovs = list(obs.cov = data$survey.cov), type = "removal"))  # Create and look at um data frame

fm1 <- multinomPois(~1 ~1, umf)    # Detection model before abundance model
(fm4 <- multinomPois(~cov3 + obs.cov ~cov2 + cov3, umf))

print(c(log(data$mean.lam), data$beta2.lam, data$beta3.lam,
        logit(data$mean.p), data$beta3.p, data$beta.p.survey))


# Simulate detection frequency data  from back in section 7.3.1
set.seed(24)
data <- simNmix(mean.lam = exp(0), beta2.lam = 1, beta3.lam =-1,
                beta4.lam = 0.2, dispersion = 1, mean.p = plogis(0),
                beta3.p = 1, sigma.p.visit = 1, Neg.Bin = TRUE)
dhfreq <- array(NA, dim = c(data$nsite, 7),
                dimnames = list(NULL, c("100", "010", "001", "110", "101", "011", "111")))

for(i in 1:data$nsite){
  dhfreq[i,] <- table(factor(paste(data$DH[i,1,], data$DH[i,2,],
                                   data$DH[i,3,], sep = ""),
                             levels = c("100", "010", "001", "110", "101", "011", "111")))
}
dhfreq                     # Look at resulting data set

# Bundle data in unmarked frame
time <- matrix(as.character(1:3), data$nsite, 3, byrow = T)
summary(umf <- unmarkedFrameGMM(y = dhfreq, numPrimary = 1,
                                siteCovs = data.frame(cov2 = data$site.cov[,2],
                                                      cov3 = data$site.cov[,3], cov4 = data$site.cov[,4]),
                                obsCovs = list(time = time), obsToY = o2y, piFun = "crPiFun"))

# Fit a couple of models, first for detection
fm1 <- gmultmix(lambdaformula = ~1, phiformula = ~1, pformula = ~1, mix = "NB", data = umf)
fm2 <- gmultmix(~1, ~1, ~time-1, mix = "NB", data = umf)
fm3 <- gmultmix(~1, ~1, ~time-1+cov3, mix = "NB", data = umf)

# ... then for abundance,
fm4 <- gmultmix(~cov2, ~1, ~1, mix = "NB", data = umf)
fm5 <- gmultmix(~cov2+cov3, ~1, ~1, mix = "NB", data = umf)
fm6 <- gmultmix(~cov2+cov3+cov4, ~1, ~1, mix = "NB", data = umf)

# ... and the data-generating model
fm7 <- gmultmix(~cov2+cov3+cov4, ~1, ~time-1+cov3, mix = "NB", data = umf)

# Compare models with AIC
ms <- fitList(
  "lam(.)p(.)" = fm1,
  "lam(.)p(time)" = fm2,
  "lam(.)p(time+cov3)" = fm3,
  "lam(cov2)p(.)" = fm4,
  "lam(cov2+cov3)p(.)" = fm5,
  "lam(cov2+cov3+cov4)p(.)" = fm6,
  "lam(cov2+cov3+cov4)p(time+cov3)" = fm7)

(AICtable <- modSel(ms))

# Compare data-generation truth with estimates
p.true <- qlogis(data$mean.p) + data$eta.p.visit
Truth <- rbind(lam0 = log(data$mean.lam), beta2.lam = data$beta2.lam, beta3.lam = data$beta3.lam, beta4.lam = data$beta4.lam, logit.p1 = p.true[1], logit.p2 = p.true[2], logit.p3 = p.true[3], beta3.p = data$beta3.p, log.dispersion = log(data$dispersion))
MLEs <- coef(fm7)
print(cbind(Truth, MLEs), 3)



# 7.9 Example 3: Jays in the Swiss MHB
# ------------------------------------------------------------------------


data(jay)        # Load data
str(jay)         # Inspect data list
dim(jay$caphist) # Look at detection history data

jay$caphist[1:4,]

# 7.9.1 Setting up the data and preparing for analysis
# ------------------------------------------------------------------------
crPiFun <- function(p) {
  p1 <- p[,1] # Extract the columns of the p matrix, one for
  p2 <- p[,2] #   each of J = 3 sample occasions
  p3 <- p[,3]
  cbind(      # define multinomial cell probabilities:
    "100" = p1 * (1-p2) * (1-p3),
    "010" = (1-p1) * p2 * (1-p3),
    "001" = (1-p1) * (1-p2) * p3,
    "110" = p1 * p2 * (1-p3),
    "101" = p1 * (1-p2) * p3,
    "011" = (1-p1) * p2 * p3,
    "111" = p1 * p2 * p3,
    "10x" = p1*(1-p2),
    "01x" = (1-p1)*p2,
    "11x" = p1*p2)
}


o2y <- matrix(1, 3, 10)

# Grab the data objects
covinfo <- jay$covinfo
gridinfo <- jay$gridinfo
sitecovs <- jay$sitecovs
caphist <- jay$caphist

# Get observation covariates to use in model
# Day of year, sample intensity and survey duration.
# Standardize them.
day <- scale(covinfo[,c("date1", "date2", "date3")])
dur <- as.matrix(covinfo[,c("dur1", "dur2", "dur3")])
dur[is.na(dur)] <- mean(dur, na.rm=TRUE) # Pad the 6 missing values
intensity <- dur / sitecovs[,"length"] # Sample rate = duration/length
dur <- scale(dur)
intensity <- scale(intensity)


reps <- apply(!is.na(day), 1, sum)
day[reps==2,3] <- 0
dur[reps==2,3] <- 0
intensity[reps==2, 3] <- 0

# Store the observation covariates in a list
obscovs <- list(intensity = intensity, dur = dur, day = day)

# Standardize site covariates
sitecovs[,"elev"] <- scale(sitecovs[,"elev"])
sitecovs[,"forest"] <- scale(sitecovs[,"forest"])
# NOTE: length is NOT standardized, but over-written with its reciprocal
sitecovs[,"iLength"] <- 1 / sitecovs[,"length"]

# Create unmarkedFrame (need crPiFun above and unmarked loaded)
caphist <- as.matrix(caphist)
mhb.umf <- unmarkedFrameMPois(y=caphist, siteCovs=as.data.frame(sitecovs),
                              obsCovs=obscovs, obsToY=o2y, piFun="crPiFun")



# 7.9.2. Fitting some models
# ------------------------------------------------------------------------
# Fit a series of models
fm1 <- multinomPois(~1 ~1, mhb.umf)
fm2 <- multinomPois(~day ~1, mhb.umf)
fm3 <- multinomPois(~day + I(day^2) ~1, mhb.umf)
fm4 <- multinomPois(~intensity ~1, mhb.umf)
fm5 <- multinomPois(~intensity + I(intensity^2) ~1, mhb.umf)
fm6 <- multinomPois(~day + intensity ~1, mhb.umf)
fm7 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~1, mhb.umf)

# Assemble the models into a fitList and rank them by AIC
mspart1 <- fitList(
  "lam(.)p(.)" = fm1,
  "lam(.)p(day)" = fm2,
  "lam(.)p(day+day^2)" = fm3,
  "lam(.)p(intensity)" = fm4,
  "lam(.)p(intensity + intensity^2)" = fm5,
  "lam(.)p(day + rate)" = fm6,
  "lam(.)p(data + day^2 + intensity + intensity^2)" = fm7)

(mspart1 <- modSel(mspart1))

# Fit series of models
fm7 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~1, mhb.umf)
fm8 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev, mhb.umf)
fm9 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest, mhb.umf)
fm10 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~iLength, mhb.umf)
fm11 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev, mhb.umf)
fm12 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + iLength, mhb.umf)
fm13 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev + iLength, mhb.umf)
fm14 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + iLength, mhb.umf)
fm15 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~elev + I(elev^2), mhb.umf)
fm16 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + I(elev^2), mhb.umf)
fm17 <- multinomPois(~day + I(day^2) + intensity + I(intensity^2) ~forest + elev + I(elev^2) + iLength, mhb.umf)


# Assemble the models into a fitList
mspart2 <- fitList(
  "lam(.)p(best)"                               = fm7,
  "lam(elev)p(best)"                            = fm8,
  "lam(forest)p(best)"                          = fm9,
  "lam(length)p(best)"                          = fm10,
  "lam(forest + elev)p(best)"                   = fm11,
  "lam(forest + iLength)p(best)"                = fm12,
  "lam(elev + iLength)p(best)"                  = fm13,
  "lam(forest + elev + length)p(best)"          = fm14,
  "lam(elev + elev^2)p(best)"                   = fm15,
  "lam(forest + elev + elev^2)p(best)"          = fm16,
  "lam(forest + elev + elev^2 + iLength)p(best)"= fm17)

# Rank them by AIC
(mspart2 <- modSel(mspart2))

fm17

mhb.umf2 <- unmarkedFrameGMM(y=caphist, numPrimary = 1,
                             siteCovs=as.data.frame(sitecovs), obsCovs=obscovs, obsToY=o2y, piFun="crPiFun")

fm1NB <- gmultmix(~1, ~1, ~1, mix = "NB", data = mhb.umf2)

fm17P <- gmultmix(~forest + elev + I(elev^2) + iLength, ~1, ~day + I(day^2) + intensity + I(intensity^2), mix = "P",   data = mhb.umf2)

fm17NB <- gmultmix(~forest + elev + I(elev^2) + iLength, ~1, ~day + I(day^2) + intensity + I(intensity^2), mix = "NB",   data = mhb.umf2)


fm17P    # AIC-best Poisson mixture model


fm17NB  # NegBin version of AIC-best Poisson mixture model



# 7.9.3 Analysis of model fit
# ------------------------------------------------------------------------
# Define new fitstats function
fitstats2 <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  n.obs<- apply(observed,1,sum,na.rm=TRUE)
  n.pred<- apply(expected,1,sum,na.rm=TRUE)
  sse <- sum(resids^2,na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
  freeTuke.n<- sum((sqrt(n.obs)-sqrt(n.pred))^2,na.rm=TRUE)
  sse.n <- sum( (n.obs -n.pred)^2,na.rm=TRUE)
  chisq.n <- sum((n.obs - n.pred)^2 / expected,na.rm=TRUE)
  
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke,
           SSE.n = sse.n, Chisq.n = chisq.n, freemanTukey.n=freeTuke.n)
  return(out)
}


(pb.mhb <- parboot(fm17, fitstats2, nsim=1000, report=1))


(pb.mhbNB <- parboot(fm17NB, fitstats2, nsim=1000, report=1))

# Compute c-hat
pb.mhbNB@t0[2] / mean(pb.mhbNB@t.star[,2])

n.obs <- apply(caphist, 1, sum, na.rm=TRUE)
n.predP <- apply(fitted(fm17P), 1, sum, na.rm=TRUE)
n.predNB <- apply(fitted(fm17NB), 1, sum, na.rm=TRUE)
plot(n.obs, n.predP, frame = F)   # Fig. 7-7
abline(0,1)
points(smooth.spline(n.predP ~ n.obs, df = 4), type = "l", lwd = 2, col = "blue")
points(smooth.spline(n.predNB ~ n.obs, df=4), type= "l", lwd = 2, col = "red")



# 7.9.4 Summary analyses of jay models
# ------------------------------------------------------------------------
range(siteCovs(mhb.umf)$elev)

elev.mean <- attr(siteCovs(mhb.umf)$elev, "scaled:center")
elev.sd <- attr(siteCovs(mhb.umf)$elev, "scaled:scale")
elev.orig <- elev.sd*seq(-1.5, 2.42,,500)  + elev.mean


# Remember length = 0 is saturation sampling because length = 1/L
newL <- data.frame(elev = seq(-1.5,2.42,,500), elev.orig, forest = 0, iLength = 1/5.1)           # 'Low' prediction
newH <- data.frame(elev = seq(-1.5,2.42,,500), elev.orig, forest = 0, iLength = 0)           # 'High' prediction
predL <- predict(fm17NB, type="lambda", newdata=newL, appendData=TRUE)
predH <- predict(fm17NB, type="lambda", newdata=newH, appendData=TRUE)
head(cbind(low = predL[,1:2], high = predH[,1:2]))

plot(Predicted ~ elev.orig, predL, type="l", lwd = 3, xlab="Elevation",
     ylab="Expected # territories", ylim=c(0, 13), frame=F, col = "red")
points(Predicted ~ elev.orig, predH, type="l", lwd = 3, col = "blue")
matlines(elev.orig, predL[,3:4], lty = 1, lwd = 1, col = "red")
matlines(elev.orig, predH[,3:4], lty = 1, lwd = 1, col = "blue")

b <- coef(fm17NB)[3]
c <- coef(fm17NB)[4]
elev.opt <- -b / (2*c)

(elev.opt <- elev.opt*elev.sd + elev.mean)

require(AICcmodavg)
model.list <- list(fm17NB) # candidate model list with single model
model.names <- c("AIC-best model")

# Compute model-averaged predictions of abundance for values of elevation, with uncertainty (SE, CIs) adjusted for overdispersion (c.hat), with latter estimated from bootstrapped Chisquare
pred.c.hatL <- modavgPred(cand.set = model.list, modnames = model.names, newdata = newL, parm.type = "lambda", type = "response", c.hat = 1.11)

# Compare predictions and SE without and with c.hat adjustment
head(cbind(predL[1:2], pred.c.hatL[2:3]), 10) ## see errata


rlength <- seq(1, 30, 0.01)         # Vary route length from 1 to 30 kms
newData <- data.frame(elev=0, forest=0, iLength=1/rlength)
pred <- predict(fm17NB, type="lambda", newdata=newData, c.hat = 1.11)
par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.3)
plot(rlength, pred[[1]], type = "l", lwd = 3, col = "blue", frame = F, xlab = "Transect length (km)", ylab = "Exposed population (lambda)", ylim = c(0, 16), axes = F)
axis(1, at = seq(2,30,2))       ;      axis(2)
abline(v = c(1.2, 5.135, 9.4), lwd = 2)
matlines(rlength, cbind(pred[[1]]-pred[[2]], pred[[1]]+pred[[2]]), type = "l", lty = 1, lwd = 2, col = "gray")

sat.pred <- predict(fm17NB, type="lambda", newdata= data.frame(elev=0, forest=0, iLength=0), c.hat = 1.11)
abline(as.numeric(sat.pred[1]),0, lwd = 2, lty = 2)


pred[round(rlength,2)==1.2,]/as.numeric(sat.pred[1])

pred[round(rlength,2)==5.14,]/as.numeric(sat.pred[1])

pred[round(rlength,2)==9.4,]/as.numeric(sat.pred[1])


# 7.9.5 Spatial prediction
# ------------------------------------------------------------------------
library(raster)
library(rgdal)

# Swiss landscape data and shape files
data(Switzerland)         # Load Swiss landscape data from unmarked
CH <- Switzerland         # this is for 'confoederatio helvetica'
head(CH)
gelev <- CH[,"elevation"] # Median elevation of quadrat
gforest <- CH[,"forest"]
grid <- CH[,c("x", "y")]
lakes <- readOGR(".", "lakes")
rivers <- readOGR(".", "rivers")
border <- readOGR(".", "border")

# Draw two maps of Swiss elevation and forest cover (Fig. 7-10)
par(mfrow = c(1,2), mar = c(1,2,3,5))
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))
mapPalette2 <- colorRampPalette(c("grey", "lightgreen", "darkgreen"))
r1 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = CH$elevation))
r2 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = CH$forest))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Elevation (m a.s.l.)", zlim = c(0, 4000))
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)

plot(r2, col = mapPalette2(100), axes = FALSE, box = FALSE, main = "Forest cover (%)", zlim = c(0, 100))
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)


# Standardize elevation for all grid cells using the mean at sample plots
elev.mean <- attr(siteCovs(mhb.umf)$elev, "scaled:center")
elev.sd <- attr(siteCovs(mhb.umf)$elev, "scaled:scale")
gelev <- (gelev - elev.mean) / elev.sd

# Standardize forest cover also using the mean at sample plots
forest.mean <- attr(siteCovs(mhb.umf)$forest, "scaled:center")
forest.sd <- attr(siteCovs(mhb.umf)$forest, "scaled:scale")
gforest <- (gforest - forest.mean) / forest.sd


# Form predictions for Swiss landscape
newL <- data.frame(elev=gelev, forest=gforest, iLength=1/5.1)
newH <- data.frame(elev=gelev, forest=gforest, iLength=0)
pred.mhb.NB.Low <- predict(fm17NB, type="lambda", newdata=newL, appendData=T)
pred.mhb.NB.High <- predict(fm17NB, type="lambda", newdata=newH, appendData=T)

# Create rasters and mask for elevation (mask areas > 2250 m)
r1 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = pred.mhb.NB.Low[,1]))
r2 <- rasterFromXYZ(cbind(x = CH$x, y = CH$y, z = pred.mhb.NB.Low[,2]))
elev <- rasterFromXYZ(cbind(CH$x, CH$y, gelev))
elev[elev > 2250] <- NA
r1 <- mask(r1, elev)
r2 <- mask(r2, elev)

# Draw maps of jay density and standard error of density (Fig. 7???11)
par(mfrow = c(1,2), mar = c(1,2,3,5))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Density of European Jay", zlim = c(0, 10))
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Standard errors of density", zlim = c(0, 1.5))
plot(rivers, col = "dodgerblue", add = TRUE)
plot(border, col = "transparent", lwd = 1.5, add = TRUE)
plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)



# 7.9.6 Population size of jays in Switzerland
# ------------------------------------------------------------------------
out <- which(CH$water > 50 | CH$elevation > 2250)
(N.jay <- sum(pred.mhb.NB.Low[-out,1]))      # 'low prediction'

(N.jay <- sum(pred.mhb.NB.High[-out,1]))     # 'high prediction'


Nhat <- function(fm) {
  betavec <- coef(fm)[1:5]
  Xg <- cbind(rep(1, length(gforest)), gforest, gelev, gelev*gelev, 1/5.1)
  predLow <- as.numeric(exp(Xg%*%(betavec)))
  predHigh <- as.numeric(exp(Xg[,-5]%*%(betavec[-5])))
  out <- which(CH$water > 50 | CH$elevation > 2250)
  Nlow <- sum(predLow[-out])
  Nhigh <- sum(predHigh[-out])
  return(c(Nlow = Nlow, Nhigh = Nhigh))
}

set.seed(2015)
(pb.N <- parboot(fm17NB, Nhat, nsim=100, report=1))


# 7.9.7 Bayesian Analysis of the MHB Data
# ------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400
y <- as.matrix(getY(mhb.umf))

# Now  we have to stretch out the encounter frequencies into individuals...
# There were 439 unique individuals observed during the survey
eh<- unlist(dimnames(y)[2])
ehid<- col(y)[y>0 & !is.na(y)]  # Column ids, index to encounter history
eh<- eh[ehid]
siteid<- row(y)[y>0 & !is.na(y)] # Site ids
y<- y[y>0 & !is.na(y)]   # Positive counts
eh<- rep(eh, y)
siteid<- rep(siteid, y)

eh.mat<- matrix(NA,nrow=length(eh),ncol=3)
for(i in 1:length(eh)){
  eh.mat[i,]<- as.numeric(unlist(strsplit(eh[i],split="")))
}


# Define some things and do the data augmentation
nsites = nrow(sitecovs)
nind <- nrow(eh.mat)
M <- 800
y <- rbind(eh.mat, matrix(0, nrow=(M-nind), ncol=3))

# Augment site ID
site <- c(siteid, rep(NA, M-nind))
sitecovs <- siteCovs(mhb.umf)

obscovs<- obsCovs(mhb.umf)
Intensity<- matrix(obscovs[,"intensity"], nrow=nsites, ncol=3,byrow=TRUE)

# Bundle data for BUGS
data <- list(y = y, J = 3, M = M , nsites=nsites, X = as.matrix(sitecovs), Intensity= Intensity, group=site)
str(data)

# Specify model in BUGS language
cat("
    model {
    
    # Prior distributions
    p0 ~ dunif(0,1)
    alpha0 <- log(p0 / (1-p0))
    alpha1 ~ dnorm(0, 0.01)
    
    beta0 ~ dnorm(0,0.01)
    beta1 ~ dnorm(0,0.01)
    beta2 ~ dnorm(0,0.01)
    beta3 ~ dnorm(0,0.01)
    psi <- sum(lambda[]) / M   # psi is a derived parameter
    
    # Model for abundance: lambda depends on Elev, Length, Forest
    for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X[s,1] + (beta2/X[s,2]) + beta3*X[s,3]
    probs[s] <- lambda[s] / sum(lambda[])
    }
    
    # Model for individual encounter histories
    for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    
    # Observation model: p depends on Intensity
    for(j in 1:J){
    logit(p[i,j]) <- alpha0 + alpha1 * Intensity[group[i],j]
    
    pz[i,j] <- p[i,j] * z[i]
    y[i,j] ~ dbern(pz[i,j])
    }
    }
    }
    ",fill=TRUE,file="model.txt")

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3", "psi")

# Initial values
zst <- c(rep(1,M-100), rep(0,100)) ## see errata
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0=runif(1),
        beta1=rnorm(1), z= zst ) }

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

# Call JAGS from R and summarize marginal posteriors

out <- jags(data, inits, parameters, "model.txt", n.thin = nt,
            n.chains = nc, n.burnin = nb, n.iter = ni)

print(out, digits = 2)

fm17P    # AIC-best Poisson mixture model



# 7.10 Summary and Outlook
# ------------------------------------------------------------------------











# =========================================================================
#
# 8. Modeling abundance using hierarchical distance sampling (HDS)
#
# =========================================================================





# 8.1 Introduction
# ------------------------------------------------------------------------


# 8.2 Conventional Distance Sampling (CDS)
# ------------------------------------------------------------------------


# 8.2.1 The full likelihood
# ------------------------------------------------------------------------


# 8.2.2 Models of detection probability
# ------------------------------------------------------------------------


# 8.2.3 Simulating distance sampling data
# ------------------------------------------------------------------------
strip.width <- 100  # one side of the transect, really half-width
sigma <- 30         # Scale parameter of half-normal detection function

# Define half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2)) # Function definition
g(30, sig=sigma)    # Detection probability at a distance of 30m

# Plot the detection function
par(mfrow=c(1,2))
curve(g(x, sig=30), 0, 100, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, frame = F)
curve(g(x, sig=60), 0, 100, add=TRUE, lty = 2, lwd = 2)

# Define function to simulate non-hierarchical line transect data
sim.ldata <- function(N = 200, sigma = 30){
  # Function to simulate line transect data under CDS.
  # Function arguments:
  #    N: number of individuals along transect with distance u(-100, 100)
  #    sigma: scale parameter of half-normal detection function
  # Function subjects N individuals to sampling, and then retains the value
  # of x=distance only for individuals that are captured
  par(mfrow = c(1,2))
  # Plot the detection function
  curve(exp(-x^2/(2*sigma^2)), 0, 100, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, main = "Detection function", ylim = c(0,1))
  text(80, 0.9, paste("sigma:", sigma))
  xall <- runif(N, -100,100) # Distances of all N individuals
  hist(abs(xall), nclass=10, xlab = "Distance (x)", col = "grey", main = "True (grey) \nand observed distances (blue)")
  g <- function(x, sig) exp(-x^2/(2*sig^2))
  p <- g(xall, sig=sigma) # detection probability
  y <- rbinom(N, 1, p) # some inds. are detected and their distance measured
  x <- xall[y==1]      # this has direction (right or left transect side)
  x <- abs(x)          # now it doesn't have direction
  hist(x, col = "blue", add = TRUE)
  return(list(N = N, sigma = sigma, xall = xall, x = x))
}

# Obtain a data set for analysis
set.seed(2015)               # If you want to get same results
tmp <- sim.ldata(sigma = 30) # Execute function and assign results to 'tmp'
attach(tmp)


# Conditional likelihood
Lcond <- function(lsigma){  # Define conditional nll
  sigma<- exp(lsigma)
  -1*sum(log(g(x,sig=sigma)/integrate(g, 0, 100, sig=sigma)$value/100))
}

# Full likelihood
Lfull <- function(parm){    # Define full nll
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  N <- length(x)+ n0
  pbar <- integrate(g, 0, 100, sig=sigma)$value/100
  -1*( lgamma(N+1) - lgamma(n0+1) + sum(log(g(x,sig=sigma)/100)) + n0*log(1-pbar) )
}

# Call optim to maximize full likelihood
optim(c(log(30), log(4)), Lfull, hessian=TRUE)


pbar<- integrate(g, 0, 100, sig=exp(3.26))$value/100
n<- length(tmp$x)

(Nhat.condl<- n/pbar)
(Dhat.condl<- Nhat.condl/(10*.2))

n0hat<- exp(5.01)
(Nhat.full<- n + n0hat)
(Dhat.full<- Nhat.full/(10*.2))



# 8.2.4 Binned data
# ------------------------------------------------------------------------

# 8.2.4.1 Conditional and other likelihoods for binned data
# ------------------------------------------------------------------------

# 8.2.4.2 Simulating binned distance sampling data
# ------------------------------------------------------------------------
set.seed(2015)
# Design settings and truth (population size N and detection function g)
interval.width <- 10
strip.width <- 100    # half-width really (one side of transect)
nbins<-strip.width%/%interval.width
sigma <- 30           # Scale parameter of half-normal detection function
g <- function(x, sig) exp(-x^2/(2*sig^2)) # Half-normal detection function
N <- 200              # Population size

# Method 1: simulate continuous distances and put into intervals
x <- runif(N, -strip.width, strip.width) # Distance all animals
p <- g(x, sig=sigma)  # Detection probability
y <- rbinom(N, 1, p)  # only individuals with y=1 are detected
x <- x[y==1]          # this has direction (right or left side of transect)
x <- abs(x)           # now it doesn't have direction

# Compute the distance category of each observation
xbin <- x %/% interval.width + 1   # note integer division function %/%

# Multinomial frequencies, may have missing levels
y.obs <- table(xbin)

# Pad the frequencies to include those with 0 detections
y.padded <- rep(0,nbins)
names(y.padded) <- 1:nbins
y.padded[names(y.obs)] <- y.obs
y.obs <- y.padded
y.true <- c(y.obs, N-length(xbin)) # Last category is "Not detected"

# Relative frequencies by binning continuous data (pi). These should compare
#  with the cell probabilities computed below when N is very large
(y.rel <- y.true/N)      # Last category is pi(0) from above
(pi0.v1 <- y.rel[nbins+1])

# Compute detection probability in each distance interval
dist.breaks <- seq(0, strip.width, by=interval.width)
p <- rep(NA, length(dist.breaks)-1)
for(j in 1:length(p)){
  p[j] <- integrate(g, dist.breaks[j], dist.breaks[j+1],
                    sig=sigma)$value / (dist.breaks[j+1]-dist.breaks[j])
}
round(p, 2)

# Compute the multinomial cell probabilities analytically. These are exact.
# psi = probability of occurring in each interval
interval.width <- diff(dist.breaks)
psi <- interval.width/strip.width
pi <- p * psi
sum(pi)                 # This is 1 ??? pi(0) from above
(pi0.exact <- 1-sum(pi))

# Method 2: Use rmultinom to simulate binned observations directly
# This includes 0 cells AND n0
pi[length(p)+1] <- 1 - sum(pi)
(y.obs2 <- as.vector(rmultinom(1, N, prob=pi)))
(y.obs2 <- y.obs2[1:nbins]) # Discard last cell for n0 (because not observed)

Lik.binned <- function(parm, data, dist.breaks){
  # Note that the parameters are parm[1] = log(sigma), parm[2] = log(n0)
  
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  p <- rep(NA, length(dist.breaks)-1)
  for(j in 1:length(p)) {
    p[j] <- integrate(g, dist.breaks[j], dist.breaks[j+1],
                      sig=sigma)$value / (dist.breaks[j+1]-dist.breaks[j])
  }
  psi <- interval.width/strip.width
  pi <- p * psi
  pi0 <- 1-sum(pi)
  
  N <- sum(data) + n0
  -1*(lgamma(N+1)-lgamma(n0+1) + sum(c(data,n0)*log(c(pi,pi0))))
}

# Evaluate likelihood for some particular value of the parameters
Lik.binned(c(2,0), data=y.obs, dist.breaks=dist.breaks)

# Obtain the MLEs for the simulated data
optim(c(2,0), Lik.binned, data=y.obs, dist.breaks=dist.breaks)


# 8.2.5 Point transect data
# ------------------------------------------------------------------------
# Define function to compute cell probs for binned distance sampling
cp.ri <-function(radius1, radius2, sigma){
  Pi <- 3.141593
  a <- Pi*radius2^2 - Pi*radius1^2
  integrate(function(r, s=sigma) exp(-r^2 / (2 * s^2)) * r, radius1,radius2)$value *(2*Pi/a)
}

# Define distance intervals and compute multinomial probabilities
delta <- 0.5                   # Width of distance bins
B <- 3                         # Max count distance
dist.breaks <-seq(0, B, delta) # Make the interval cut points
nD <-length(dist.breaks)-1
sigma <- 1
p.x <-rep(NA,nD)               # Conditional detection probabilities
for(i in 1:nD){
  p.x[i] <- cp.ri(dist.breaks[i], dist.breaks[i+1], sigma =1)
}
area <- 3.141593 * dist.breaks[-1]^2
ring.area <- diff(c(0, area))
# Pr(detection| in ring)*Pr(in ring)
cp <- p.x* ring.area/sum(ring.area)


# 8.2.5.1 Simulating point transect data
# ------------------------------------------------------------------------
sim.pdata <- function(N=1000, sigma=1, B=3, keep.all=FALSE) {
  # Function simulates coordinates of individuals on a square
  # Square is [0,2*B] x[0,2*B], with a count location on the center
  # point (B,B)
  # Function arguments:
  #    N: total population size in the square
  #    sigma: scale of half-normal detection function
  #    B: circle radias
  #    keep.all: return the data for y = 0 individuals or not
  
  # Plot the detection function
  par(mfrow = c(1,2))
  curve(exp(-x^2/(2*sigma^2)), 0, B, xlab="Distance (x)", ylab="Detection prob.", lwd = 2, main = "Detection function", ylim = c(0,1))
  text(0.8*B, 0.9, paste("sigma:", sigma))
  
  # Simulate and plot simulated data
  library(plotrix)
  u1 <-runif(N, 0, 2*B)           # (u1,u2) coordinates of N individuals
  u2 <- runif(N, 0, 2*B)
  d <- sqrt((u1 - B)^2 + (u2 - B)^2) # distance to center point of square
  plot(u1, u2, asp = 1, pch = 1, main = "Point transect")
  N.real <- sum(d<= B)           # Population size inside of count circle
  
  # Can only count indidividuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them):
  p <- ifelse(d < B, 1, 0) * exp(-d*d/(2*(sigma^2)))
  # Now we decide whether each individual is detected or not
  y <- rbinom(N, 1, p)
  points(u1[d <= B], u2[d <= B], pch = 16, col = "black")
  points(u1[y==1], u2[y==1], pch = 16, col = "blue")
  points(B, B, pch = "+", cex = 3, col = "red")
  draw.circle(B, B, B)
  
  # Put all of the data in a matrix:
  #      (note we don't care about y, u, or v normally)
  
  if(!keep.all){
    u1 <- u1[y==1]
    u2 <- u2[y==1]
    d <- d[y==1]
  }
  return(list(N=N, sigma=sigma, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real))
}

# obtain a data set by distance sampling a population of N=1000
set.seed(1234)
tmp <-sim.pdata(N=1000, sigma=1, keep.all=FALSE, B=3)
attach(tmp)


# Bin the data and tabulate the bin frequencies. Be sure to pad the 0s!
delta <- 0.5                   # width of distance bins
dist.breaks <-seq(0, B, delta) # make the interval cut points
dclass <- tmp$d %/% delta +1   # Convert distances to categorical distances
nD<-length(dist.breaks)-1     # How many intervals do we have ?
y.obs <- table(dclass)         # next pad the frequency vector
y.padded <- rep(0, nD)
names(y.padded) <- 1:nD
y.padded[names(y.obs)] <- y.obs
y.obs <- y.padded


cp <- c(cp, 1-sum(cp))  # Compute the last cell and add it to the vector

as.vector(rmultinom(n=1, size=1000, prob=cp))


# 8.2.5.2 Likelihood analysis of point transect data
# ----------------------------------------------------------------------------------------
# (1) Define multinomial likelihood for binned data
Lik.binned.point <- function(parm, data, dist.breaks){
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  p.x <-rep(NA, nD)
  for(i in 1:nD){
    p.x[i] <- cp.ri(dist.breaks[i], dist.breaks[i+1], sigma =sigma)
  }
  area <- 3.141593 * dist.breaks[-1]^2
  ring.area <- diff(c(0, area))
  cp <- p.x* ring.area/sum(ring.area) # Pr(detection| in ring)*Pr(in ring)
  pi0 <- 1-sum(cp)
  N <- sum(data) + n0
  negLL <- -1*(lgamma(N+1)-lgamma(n0+1) + sum(c(data,n0)*log(c(cp,pi0))))
  return(negLL)
}

# Fit model
mle1 <- optim(c(2,0), Lik.binned.point, data=y.obs, dist.breaks=dist.breaks)

# (2) Define full likelihood for continuous data
Lik.cont.point <- function(parm, data, B){
  sigma <- exp(parm[1])
  n0 <- exp(parm[2])
  n <- length(data)
  N <- n + n0
  p <- exp(-data*data/(2*sigma*sigma))
  f <- 2*data/(B^2)
  pbar <- integrate(function(r, s=sigma) exp(-r^2 / (2 * s^2)) * r, 0, B)$value*2/(B^2)
  negLL <- -1*sum( log(p*f/pbar )) -1*(lgamma(N+1) - lgamma(n0+1) +
                                         n * log(pbar) + n0*log(1-pbar))
  return(negLL)
}

# Fit model
mle2 <- optim(c(0, 5), Lik.cont.point, data=tmp$d, B=B, hessian=TRUE)

# Compare two solutions and with realized true value of N
(Nhat.binned <- length(tmp$d) + exp(mle1$par[2]))
(Nhat.cont <- length(tmp$d) + exp(mle2$par[2]))
tmp$N.real


# (3) Define conditional likelihood for continuous data
Lik.cond.point <- function(parm, data, B){
  sigma <- exp(parm)
  p <- exp(-data*data/(2*sigma*sigma))
  f <- 2*data/(B^2)
  pbar <- integrate(function(r, s=sigma) exp(-r^2 / (2 * s^2)) * r, 0, B)$value*2/(B^2)
  negLL <- -1*sum( log(p*f/pbar ))
  return(negLL)
}

# Fit the model
mle3 <- optim(c(0), Lik.cond.point, data=tmp$d, B=B, method="Brent", hessian=TRUE, lower=-10, upper=10)

# Inspect the output
mle3

# Estimated sigma
(sigma.hat <- exp(mle3$par))


# Estimated average detection probability and conditional estimator of N
pbar <- integrate(
  function(r, s=sigma.hat) exp(-r^2 / (2 * s^2)) * r, 0, B)$value*2/(B^2)

(Nhat.condl <- length(d) / pbar)


# 8.2.6 Sensitivity to bin width
# ------------------------------------------------------------------------
set.seed(1234)
simrep <- 1000              # Number of sim reps
simout <- matrix(NA, nrow=simrep, ncol=3)
colnames(simout) <- c("N.real", "N.binned", "N.continuous")
delta <- 0.5                # Set width of bins

# Begin simulation loop
for(sim in 1:simrep){
  tmp <- sim.pdata(N=1000, sigma=1, keep.all=FALSE, B=3)
  B <- tmp$B
  d <- tmp$d
  N.real <- tmp$N.real
  
  # Bin data, tabulate frequencies and pad 0s if necessary
  dist.breaks <- seq(0, B, delta)
  dclass <- d%/%delta + 1      # Convert distances to categorical distances
  nD <- length(dist.breaks) -1 # How many intervals do we have ?
  y.obs <- table(dclass) # Next pad the frequency vector
  y.padded <- rep(0, nD)
  names(y.padded) <- 1:nD
  y.padded[names(y.obs)] <- y.obs
  y.obs <- y.padded
  
  # Obtain the MLEs using both models
  binned.est <- optim(c(2,0), Lik.binned.point, data=y.obs,
                      dist.breaks=dist.breaks)
  cont.est <- optim(c(1, 6), Lik.cont.point, data=d, B=B, hessian=TRUE)
  Nhat.binned <- length(d) + exp(binned.est$par[2])
  Nhat.cont <- length(d) + exp(cont.est$par[2])
  
  # Store things in a matrix
  simout[sim,] <- c(N.real, Nhat.binned, Nhat.cont)
}

# Now summarize the output
apply(simout, 2, mean)

sqrt(apply(simout, 2, var))



# 8.2.7 Spatial sampling
# ------------------------------------------------------------------------



# 8.3 Bayesian conventional distance sampling
# ------------------------------------------------------------------------


# 8.3.1 Bayesian analysis of line transect data
# ------------------------------------------------------------------------
# Get data and do data-augmentation
# Observed distances (meters)
x <- c(71.93, 26.05, 58.47, 92.35, 163.83, 84.52, 163.83, 157.33,
       22.27, 72.11, 86.99, 50.8, 0, 73.14, 0, 128.56, 163.83, 71.85,
       30.47, 71.07, 150.96, 68.83, 90, 64.98, 165.69, 38.01, 378.21,
       78.15, 42.13, 0, 400, 175.39, 30.47, 35.07, 86.04, 31.69, 200,
       271.89, 26.05, 76.6, 41.04, 200, 86.04, 0, 93.97, 55.13, 10.46,
       84.52, 0, 77.65, 0, 96.42, 0, 64.28, 187.94, 0, 160.7, 150.45,
       63.6, 193.19, 106.07, 114.91, 143.39, 128.56, 245.75, 123.13,
       123.13, 153.21, 143.39, 34.2, 96.42, 259.81, 8.72)

B <- 500 # Strip half-width. Larger than max distance
nind <- length(x)

# Analysis of continuous data using data augmentation (DA)
nz <- 200 # Augment observed data with nz = 200 zeroes
y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
x <- c(x, rep(NA, nz)) # Value of distance are missing for the augmented

# Bundle and summarize data set
str( win.data <- list(nind=nind, nz=nz, x=x, y=y, B=B) )

# Save text file with BUGS model
cat("
    model {
    
    # Priors
    sigma ~ dunif(0,1000)  # Half-normal scale
    psi ~ dunif(0,1)       # DA parameter
    
    # Likelihood
    for(i in 1:(nind+nz)){
    # Process model
    z[i] ~ dbern(psi)   # DA variables
    x[i] ~ dunif(0, B)  # Distribution of distances
    # Observation model
    logp[i] <- -((x[i]*x[i])/(2*sigma*sigma)) # Half-normal detection fct.
    p[i] <- exp(logp[i])
    mu[i] <- z[i] * p[i]
    y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
    }
    # Derived quantities
    N <- sum(z[1:(nind + nz)]) # Population size
    D <- N / 60                # Density, with A = 60 km^2 when B = 500
    }
    ",fill=TRUE,file="model1.txt")

# Inits
zst <- y
inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }

# Params to save
params <- c("N", "sigma", "D")

# Experience the raw power of BUGS and summarize marginal posteriors
library(R2WinBUGS)
bd <- "c:/Program Files/WinBUGS14/"    # May have to adapt
out1 <- bugs(win.data, inits, params, "model1.txt", n.thin=2,n.chains=3,
             n.burnin=1000, n.iter=11000, debug=TRUE, DIC=FALSE, bugs.dir=bd)
print(out1, 3)


# Analysis of binned data using data augmentation
delta <- 50                # Width of distance bins
xg <- seq(0, B, delta)     # Make the interval cut points
dclass <- x %/% delta + 1  # Convert distances to distance category
nD <- length(xg) -1        # N intervals = length(xg) if max(x) = B

# Bundle data
# Note data changed to include dclass, nG, bin-width delta and midpt
midpt <- xg[-1] - delta/2  # Interval mid-points
str( win.data <- list (nind=nind, nz=nz, dclass=dclass, y=y, B=B,
                       delta=delta, nD=nD, midpt=midpt) )   # Bundle and summarize

# BUGS model specification
cat("
    model{
    # Priors
    psi ~ dunif(0, 1)
    sigma ~ dunif(0, 1000)
    
    # Likelihood
    # construct conditional detection probability and Pr(x) for each bin
    for(g in 1:nD){        # midpt = mid point of each cell
    log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
    pi[g] <- delta / B  # probability of x in each interval
    }
    
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)             # model for individual covariates
    dclass[i] ~ dcat(pi[])        # population distribution of distance class
    mu[i] <- z[i] * p[dclass[i]]  # p depends on distance class
    y[i] ~ dbern(mu[i])
    }
    # Derived quantities: Population size and density
    N <- sum(z[])
    D <- N / 60
    }
    ",fill=TRUE, file = "model2.txt")

# Inits function
zst <- y # DA variables start at observed value of y
inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }

# Parameters to save
params <- c("N", "sigma", "D")

# Unleash WinBUGS and summarize posteriors
bd <- "c:/Program Files/WinBUGS14/"
out2 <- bugs(win.data, inits, params, "model2.txt", n.thin=2, n.chains=3,
             n.burnin=1000, n.iter=11000, debug=TRUE, DIC=FALSE, bugs.dir = bd)
print(out2, 2)


# 8.3.2 Other formulations of the distance sampling model
# ------------------------------------------------------------------------

# 8.3.3 A treatise on the integration of mathematical functions in one dimension
# ------------------------------------------------------------------------
sigma <- 2                        # normal scale (standard deviation)
curve(exp(-x^2 / (2*sigma^2)), 0, 10, frame = F)

delta <- 1                        # bin width
mid <- seq(0.5, 9.5, delta)       # 10 rectangles
f.mid <- exp(-mid^2 / (2*sigma^2))
barplot(f.mid, add=T, space=0, col="grey", width=delta)
curve(exp(-x^2 / (2*sigma^2)), 0, 10, add = TRUE, col = "blue", lwd = 3)

# Integral done using the integrate function
integrate( function(x){ exp(-x^2/(2*sigma^2)) }, lower=0, upper=100)

# Summing up the 10 rectangular areas:
areas <- f.mid * delta
sum(areas)


# 8.3.4 Bayesian analysis of point transect data
# ------------------------------------------------------------------------
### Version 1: Point count data in BUGS (conditional likelihood)
# Simulate a data set and harvest the output
set.seed(1234)
tmp <- sim.pdata(N=200, sigma=1, keep.all=FALSE, B=3)
attach(tmp)

# Chop the data into bins
delta <- 0.1           # width of distance bins for approximation
xg <- seq(0, B, delta) # Make the mid points and chop up the data
midpt <- xg[-1] - delta/2

# Convert distances to categorical distances (which bin?)
dclass <- d %/% delta + 1
nD <- length(midpt) # how many intervals
nind <- length(dclass)

# Bundle and summarize data set
str( win.data <- list(midpt=midpt, delta=delta, B=B, nind=nind, nD=nD, dclass=dclass) )


# BUGS model specification, conditional version
cat("
    model{
    
    # Prior for single parameter
    sigma ~ dunif(0, 10)
    
    # Construct cell probabilities for nG cells (rectangle approximation)
    for(g in 1:nD){    # midpt[g] = midpoint of each distance band
    log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
    pi[g] <- (( 2 * midpt[g] ) / (B*B)) * delta
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap
    }
    pcap <- sum(f[]) # capture prob. is the sum of all rectangular areas
    
    # Categorical observation model
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[])
    }
    # Derived quantity: population size
    N <- nind / pcap
    D<- N/(3.141*B*B)
    }
    ",fill=TRUE, file="model3.txt")

# Inits function
inits <- function(){list (sigma=runif(1, 1, 10)) }

# Params to save
params <- c("sigma", "N","D")

# MCMC settings
ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3

# Run BUGS and summarize posteriors
bd <- "c:/Program Files/WinBUGS14/"
out3 <- bugs(win.data, inits, params, "model3.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, debug=FALSE, bugs.dir = bd)


## Version 2: point count data (full likelihood with data augmentation)
# Do data augmentation (for same simulated data set)
M <- 400
nz <- M - nind
y <- c(rep(1, nind), rep(0, nz))
dclass <- c(dclass, rep(NA, nz))

# Bundle and summarize data set
str( win.data <- list(midpt=midpt, delta=delta, B=B, nind=nind, nD=nD, dclass=dclass, y=y, nz=nz) )

# BUGS model
cat("
    model{
    
    # Priors
    sigma ~ dunif(0, 10)
    psi ~ dunif(0, 1)
    
    # Construct cell probabilities for nG cells (rectangle approximation)
    for(g in 1:nD){           # midpt[g] = midpoint of each distance band
    log(p[g]) <- -midpt[g] * midpt[g] / (2*sigma*sigma)
    pi[g] <- ((2 * midpt[g]) / (B * B)) * delta
    pi.probs[g] <- pi[g] / norm
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap   # conditional probabilities
    }
    pcap <- sum(f[])# capture prob. is the sum of all rectangular areas
    norm <- sum(pi[])
    
    # Categorical observation model
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    dclass[i] ~ dcat(pi.probs[])
    mu[i] <- p[dclass[i]] * z[i]
    y[i] ~ dbern(mu[i])
    }
    
    # Derived quantity: population size
    N <- sum(z[])
    D<- N/(3.141*B*B)
    
    }
    ",fill=TRUE,file="model4.txt")

# Inits
inits <- function(){list (sigma=runif(1,1,10), psi=runif(1) ) }

# Params to save
params <- c("sigma", "N","D","psi")

# MCMC settings
ni <- 62000   ;   nb <- 2000   ;   nt   <-   2   ;   nc <- 3

# Run BUGS and summarize posteriors
out4 <- bugs(win.data, inits, params, "model4.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, debug=FALSE, bugs.dir = bd)


# Compare posterior summaries
print(out3,2)

print(out4,2)



# 8.4 Hierarchical Distance Sampling
# ------------------------------------------------------------------------



# 8.4.1 HDS data structure and model
# ------------------------------------------------------------------------


# 8.4.2 HDS in unmarked
# ------------------------------------------------------------------------


# 8.4.3 Example: Estimating the global population size of the Island Scrub Jay
# ------------------------------------------------------------------------
# Load, view and format the ISSJ data
library(unmarked)
data(issj)
round(head(issj), 2)

# Package things up into an unmarkedFrame
covs <- issj[,c("elevation", "forest", "chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameDS(y=as.matrix(issj[,1:3]),
                          siteCovs=data.frame(covs, area),
                          dist.breaks=c(0, 100, 200, 300),
                          unitsIn="m", survey="point")


# Fit model 1
(fm1 <- distsamp(~chaparral ~chaparral + elevation + offset(log(area)),
                 jayumf, keyfun="halfnorm", output="abund"))

# Fit model 2
(fm2 <- distsamp(~1 ~chaparral + elevation + offset(log(area)),
                 jayumf, keyfun="halfnorm", output="abund"))


(pb <- parboot(fm1, fitstats, nsim=1000, report=5))
(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed
# and mean of expected value of Chi2 (under H0)
# (see, e.g., Johnson et al., Biometrics, 2010)

residuals(fm1)             # Can inspect residuals
plot(pb)                   # Not shown
print(pb)


# Standardize the covariates
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
siteCovs(jayumf) <- sc.s
summary(jayumf)


# Fit a bunch of models and produce a model selection table.
fall <- list()   # make a list to store the models

# With the offset output=abund is the same as output = density
fall$Null <- distsamp(~1 ~offset(log(area)), jayumf, output="abund")
fall$Chap. <- distsamp(~1 ~chaparral + offset(log(area)), jayumf,
                       output="abund")
fall$Chap2. <- distsamp(~1 ~chaparral+I(chaparral^2)+offset(log(area)),
                        jayumf, output="abund")
fall$Elev. <- distsamp(~1 ~ elevation+offset(log(area)), jayumf,
                       output="abund")
fall$Elev2. <- distsamp(~1 ~ elevation+I(elevation^2)+offset(log(area)),
                        jayumf, output="abund")
fall$Forest. <- distsamp(~1 ~forest+offset(log(area)), jayumf,
                         output="abund")
fall$Forest2. <- distsamp(~1 ~forest+I(forest^2)+offset(log(area)),
                          jayumf, output="abund")
fall$.Forest <- distsamp(~forest ~offset(log(area)), jayumf,
                         output="abund")
fall$.Chap <- distsamp(~chaparral ~offset(log(area)), jayumf,
                       output="abund")
fall$C2E. <- distsamp(~1 ~ chaparral + I(chaparral^2) + elevation +
                        offset(log(area)),jayumf, output="abund")
fall$C2F2. <- distsamp(~1 ~chaparral + I(chaparral^2) + forest +
                         I(forest^2)+offset(log(area)), jayumf,  output="abund")
fall$C2E.F <- distsamp(~forest ~chaparral+I(chaparral^2)+elevation+
                         offset(log(area)), jayumf, output="abund")
fall$C2E.C <- distsamp(~chaparral ~chaparral + I(chaparral^2) + elevation +
                         offset(log(area)), jayumf, output="abund")

# Create a fitList and a model selection table
(msFall <- modSel(fitList(fits=fall)))


# Check out the best model
fall$C2E.C


# Check out the goodness-of-fit of this model
(pb.try2 <- parboot(fall$C2E.C, fitstats, nsim=1000, report=5))
Call: parboot(object = fall$C2E.C, statistic = fitstats, nsim = 1000, report = 5)

# Express the magnitude of lack of fit by an overdispersion factor
(c.hat <- pb.try2@t0[2] / mean(pb.try2@t.star[,2]))  #    Chisq

covs <- issj[,c("elevation", "forest", "chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                           siteCovs=data.frame(covs, area), numPrimary=1,
                           dist.breaks=c(0, 100, 200, 300),
                           unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
siteCovs(jayumf) <- sc.s
summary(jayumf)


# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation +
                          offset(log(area)), ~1, ~chaparral, data =jayumf, output="abund",
                        mixture="NB", K = 150))

gdistsamp(lambdaformula = ~chaparral + I(chaparral^2) + elevation +
            offset(log(area)), phiformula = ~1, pformula = ~chaparral,
          data = jayumf, output = "abund", mixture = "NB", K = 150)

(pb.try3 <- parboot(nb.C2E.C, fitstats, nsim=1000, report=5))

(c.hat <- pb.try3@t0[2] / mean(pb.try3@t.star[,2]))  #


# *Expected* population size for the sample points
getN <- function(fm, newdata=NULL)
  sum(predict(fm, type="lambda", newdata=newdata)[,1])
getN(nb.C2E.C)

# This does the same thing as the following three commands
X <- model.matrix(~chaparral+I(chaparral^2)+elevation+log(offset(area)),
                  siteCovs(jayumf))
head(X) # The design matrix

# Prediction of total expected population size at the sample points
sum(exp(X %*% c(coef(nb.C2E.C, type="lambda"), 1)))

# Empirical Bayes estimates of posterior distribution:
# Pr(N=x | y, lambda, sigma) for x=0,1,...,K
re.jay <- ranef(nb.C2E.C, K = 150)

# *Realized* population size
sum(bup(re.jay, "mean"))


summary(jayumf) # Note the range of chaparral which we need to know

# Create a new data frame with area 28.27 ha, the area of a 300 m circle
chap.orig <- seq(0, 1, 0.01)    # Values from 0 to 1 prop. chaparral
chap.pred <- (chap.orig - mean(issj$chaparral)) / sd(issj$chaparral)
newdat <- data.frame(chaparral = chap.pred, elevation = 0, area=28.27)

# Expected values of N for covariate values in "newdat"
E.N <- predict(fall$C2E.C, type="state", newdata=newdat, appendData=TRUE)
head(E.N)

# Make a plot of the response curve for the grid of chaparral values
plot(chap.orig, E.N[,"Predicted"], xlab="Proportion chaparral", ylab="Predicted jay abundance", type="l", ylim = c(0, 20), frame = F, lwd = 2)
matlines(chap.orig, E.N[,3:4], lty = 1, col = "grey", lwd = 1)

attributes(sc.s) # means are "scaled:center". SDs are "scaled:scale"

cruz.s <- cruz   # Created a new data set for the scaled variables
cruz.s$elevation <- (cruz$elevation-202)/125
cruz.s$chaparral <- (cruz$chaparral-0.270)/0.234
cruz.s$area <- (300*300)/10000 # The grid cells are 300x300m=9ha
EN <- predict(nb.C2E.C, type="lambda", newdata=cruz.s)

# Total population size (by summing predictions for all pixels)
getN(nb.C2E.C, newdata=cruz.s)

# Parametric bootstrap for CI
# A much faster function could be written to doing the sum
set.seed(2015)
(EN.B <- parboot(nb.C2E.C, stat=getN, nsim=1000, report=5))


library(raster)
cruz.raster <- stack(rasterFromXYZ(cruz.s[,c("x","y","elevation")]),
                     rasterFromXYZ(cruz.s[,c("x","y","chaparral")]),
                     rasterFromXYZ(cruz.s[,c("x","y","area")]))
names(cruz.raster) # These should match the names in the formula

plot(cruz.raster)                      #  not shown
# Elevation map on the original scale (not shown)
plot(cruz.raster[["elevation"]]*125 + 202, col=topo.colors(20), main="Elevation (in feet) and Survey Locations", asp = 1)
points(issj[,c("x","y")], cex=0.8, pch = 16)


EN.raster <- predict(nb.C2E.C, type="lambda", newdata=cruz.raster)
plot(EN.raster, col = topo.colors(20), asp = 1)   # See Fig. 8-8



# 8.5 Bayesian HDS
# ------------------------------------------------------------------------


# 8.5.1 Simulating some HDS data
# ------------------------------------------------------------------------
set.seed(1234)
tmp1 <- simHDS("point")   # Point transect
tmp2 <- simHDS()          # Line transect (this is the default)
str(tmp1)                 # Look at function output


# 8.5.2 Bayesian HDS using data augmentation
# ------------------------------------------------------------------------
# Recreate line transect data set
set.seed(1234)
tmp <- simHDS()                  # Line transect (default)
attach(tmp)

# Data augmentation: add a bunch of "pseudo-individuals"
nz <- 500                        # Augment by 500
nind <- nrow(data)
y <- c(data[,2], rep(0, nz))     # Augmented detection indicator y
site <- c(data[,1], rep(NA, nz)) # Augmented site indicator,
# unknown (i.e., NA) for augmented inds.
d <- c(data[,5], rep(NA,nz))     # Augmented distance data (with NAs)

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, habitat=habitat, wind=wind, B=B, nind=nind, nz=nz, y=y, d=d, site=site) )
win.data$site                    # unknown site cov. for augmented inds.


# BUGS model for line transect HDS (NOT point transects!)
cat("
    model{
    # Prior distributions
    beta0 ~ dunif(-10,10)   # Intercept of lambda-habitat regression
    beta1 ~ dunif(-10,10)   # Slope of log(lambda) on habitat
    alpha0 ~ dunif(-10,10)  # Intercept of log(sigma) (half-normal scale)
    alpha1 ~ dunif(-10,10)  # Slope of log(sigma) on wind
    
    # psi is a derived parameter under DA for stratified populations
    psi <- sum(lambda[]) / (nind+nz)
    
    # 'Likelihood' (sort of...)
    for(i in 1:(nind+nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                    # Data augmentation variables
    d[i] ~ dunif(0, B)                   # distance uniformly distributed
    p[i] <- exp(-d[i]*d[i]/(2*sigma[site[i]]*sigma[site[i]])) # Det. function
    mu[i] <- z[i]* p[i]                  # 'straw man' for WinBUGS
    y[i] ~ dbern(mu[i])                  # basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    }
    
    # Linear models for abundance and for detection
    for(s in 1:nsites){                    # s is index for sites
    # Model for abundance
    # next line not necessary, but allows to make predictions
    N[s] ~ dpois(lambda[s])              # Realized abundance at site s
    log(lambda[s]) <- beta0 + beta1*habitat[s] # Linear model abundance
    site.probs[s] <- lambda[s] / sum(lambda[])
    
    # Linear model for detection
    log(sigma[s]) <- alpha0 + alpha1*wind[s]
    }
    # Derived parameter: total population size across all sites
    Ntotal <- sum(z[])
    area<- nsites*1*2*B   # Unit length == 1, half-width = B
    D<- Ntotal/area
    }
    ",fill=TRUE , file = "model1.txt")


# Inits
zst <- c(rep(1, sum(y)), rep(0, nz)) # ... and for DA variables
inits <- function(){list(beta0=0, beta1=0, alpha0=0, alpha1=0, z=zst)}

# Parameters to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Call BUGS (ART 33 min) ...
bd <- "c:/Program Files/WinBUGS14/" # Never forget this for WinBUGS
out1 <- bugs(win.data, inits, params, "model1.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)

# ... or try JAGS for a change (ART 6 min)
library(jagsUI)       # never forget to load jagsUI
out1 <- jags(win.data, inits, params, "model1.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni)

# Summarize posterior output
print(out1, 2)
sum(tmp$N.true)

# Prepare data
delta <- 0.1                    # width of distance bins for approx.
midpt <- seq(delta/2, B, delta) # make mid-points and chop up data
dclass <- d %/% delta + 1       # convert distances to cat. distances
nD <- length(midpt)             # Number of distance intervals

# Bundle and summarize data set
str( win.data <- list (y=y, dclass=dclass, site=site, midpt=midpt, delta=delta, B=B, nind=nind, nz=nz, nsites=nsites, nD=nD, habitat=habitat, wind=wind) )


# BUGS model specification for line-transect HDS (NOT point transects!)
cat("
    model{
    # Prior distributions
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    
    psi <- sum(lambda[])/(nind+nz)     # psi is a derived parameter
    
    for(i in 1:(nind+nz)){             # Loop over individuals
    z[i] ~ dbern(psi)               # DA variables
    dclass[i] ~ dcat(pi[site[i],])  # Population distribution of dist class
    mu[i] <- z[i] * p[site[i],dclass[i]] # p depends on site AND dist class
    y[i] ~ dbern(mu[i])             # Basic Bernoulli response in DS model
    site[i] ~ dcat(site.probs[1:nsites]) # Site membership of inds
    }
    
    for(s in 1:nsites){                # Loop over sites
    # Construct cell probabilities for nG cells
    for(g in 1:nD){                    # midpt = mid point of each cell
    log(p[s,g]) <- -midpt[g]*midpt[g]/(2*sigma[s]*sigma[s])
    pi[s,g] <- delta/B              # probability of x per interval
    f[s,g] <- p[s,g]*pi[s,g]        # pdf of observed distances
    }
    
    # not necessary   N[s]~dpois(lambda[s]) except for prediction
    N[s] ~ dpois(lambda[s])        # predict abundance at each site
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model for N
    site.probs[s] <- lambda[s]/sum(lambda[])
    log(sigma[s]) <- alpha0 + alpha1*wind[s] # linear model for sigma
    }
    
    # Derived parameter
    Ntotal <- sum(z[])   # Also sum(N[]) which is size of a new population
    area<- nsites*1*2*B  # Unit length == 1, half-width = B
    D<- Ntotal/area
    }
    ",fill=TRUE, file = "model2.txt")

# Inits
zst <- c(rep(1, sum(y)), rep(0, nz))
inits <- function(){list (alpha0=0, alpha1=0, beta0=0, beta1=0, z=zst) }

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run JAGS with parallel processing (ART 1 min)
library(jagsUI)
out2 <- jags(win.data, inits, params, "model2.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
print(out2,2)


# 8.5.3 Bayesian HDS using the 3-part conditional multinomial model
# ------------------------------------------------------------------------
set.seed(1234)
tmp <- simHDS(type="line", discard0=FALSE)
attach(tmp)

# Get number of individuals detected per site
# ncap = 1 plus number of detected individuals per site
ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing detected
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)

# Prepare other data
site <- data[!is.na(data[,2]),1]   # site ID of each observation
delta <- 0.1                       # distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # make mid-points and chop up data
dclass <- data[,5] %/% delta + 1   # convert distances to cat. distances
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt, delta=delta, ncap=ncap, habitat=habitat, wind=wind, dclass=dclass, site=site) )

# BUGS model specification for line-transect HDS (NOT point transects!)
cat("
    model{
    # Priors
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
    }
    
    for(s in 1:nsites){
    # Construct cell probabilities for nD multinomial cells
    for(g in 1:nD){                 # midpt = mid-point of each cell
    log(p[s,g]) <- -midpt[g] * midpt[g] / (2*sigma[s]*sigma[s])
    pi[s,g] <- delta / B          # probability per interval
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas
    
    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model abundance
    log(sigma[s])<- alpha0 + alpha1*wind[s]      # linear model detection
    }
    # Derived parameters
    Ntotal <- sum(N[])
    area<- nsites*1*2*B  # Unit length == 1, half-width = B
    D<- Ntotal/area
    }
    ",fill=TRUE, file = "model3.txt")


# Inits
Nst <- ncap + 1
inits <- function(){list(alpha0=0, alpha1=0, beta0=0, beta1=0, N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Run JAGS (ART 1 min) and summarize posteriors
library(jagsUI)
out3 <- jags(win.data, inits, params, "model3.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni)
print(out3, 2)


# 8.5.4 Point transect HDS using the conditional multinomial formulation
# ------------------------------------------------------------------------
# Simulate a data set using our simHDS function
set.seed(1234)
tmp <- simHDS(type="point", discard0=FALSE)
attach(tmp)

# Prepare data
# Number of individuals detected per site
ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing was seen
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)            # Number of individuals detected per site

# Other data
site <- data[!is.na(data[,2]),1]   # Site ID of each observation
delta <- 0.1                       # Distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # Make mid-points and chop up data
dclass <- data[,5] %/% delta + 1   # Convert distance to distance category
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
                      delta=delta, ncap=ncap, habitat=habitat, wind=wind, dclass=dclass,
                      site=site) )

# BUGS model specification for point transect data
cat("
    model{
    # Priors
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
    }
    for(s in 1:nsites){
    # Construct cell probabilities for nD distance bands
    for(g in 1:nD){                # midpt = mid-point of each band
    log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
    pi[s,g] <- ((2 * midpt[g] ) / (B * B)) * delta # prob. per interval
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas
    
    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1 * habitat[s] # linear model abundance
    log(sigma[s]) <- alpha0 + alpha1*wind[s]     # linear model detection
    }
    
    # Derived parameters
    Ntotal <- sum(N[])
    area <- nsites*3.141*B*B
    D <- Ntotal/area
    }
    ",fill=TRUE, file="model4.txt")


# Inits
Nst <- ncap + 1
inits <- function(){list(alpha0=0, alpha1=0, beta0=0, beta1=0, N=Nst)}

# Params to save
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntotal","D")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Run BUGS (not STAN !) (ART 2.3 min) and summarize posteriors
out4 <- bugs(win.data, inits, params, "model4.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)
print(out4, 2)

sum(tmp$N.true)                  # True realized population size
sum(!is.na(tmp$data[,"y"]))      # Observed population size


# 8.5.5 Analysis of the ISSJ data
# ------------------------------------------------------------------------
# Load the ISSJ data
library(unmarked)
data(issj)

# Prepare some data
nD <- 3                          # Number of intervals
delta <- 100                     # Interval width
B <- 300                         # Upper bound (max. distance)
midpt <- c(50, 150, 250)         # mid points

# Convert vector frequencies to individual distance class
H <- as.matrix(issj[,1:3])
nsites <- nrow(H)
ncap <- apply(H, 1, sum)         # Number of individuals detected per site
dclass <- rep(col(H), H)         # Distance class of each individual
nind <- length(dclass)           # Number of individuals detected
elevation <- as.vector(scale(issj[,c("elevation")])) # Prepare covariates
forest <- as.vector(scale(issj[,"forest"]))
chaparral <- as.vector(scale(issj[,"chaparral"]))


# Bundle and summarize data set
str( win.data <- list(nsites=nsites, nind=nind, B=B, nD=nD, midpt=midpt,
                      delta=delta, ncap=ncap, chaparral=chaparral, elevation=elevation, dclass=dclass) )

# BUGS model specification
cat("
    model{
    # Priors
    sigma ~ dunif(0,1000)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    beta2 ~ dunif(-10,10)
    beta3 ~ dunif(-10,10)
    sigma.site ~ dunif(0,10)
    tau <- 1/(sigma.site*sigma.site)
    
    # Specify hierarchical model
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[]) # Part 1 of HM
    }
    
    # construct cell probabilities for nG cells
    for(g in 1:nD){                # midpt = mid-point of each cell
    log(p[g]) <- -midpt[g] * midpt[g] / (2 * sigma * sigma)
    pi[g] <- ((2 * midpt[g]) / (B * B)) * delta # prob. per interval
    f[g] <- p[g] * pi[g]
    fc[g] <- f[g] / pcap
    }
    pcap <- sum(f[])               # Pr(capture): sum of rectangular areas
    
    for(s in 1:nsites){
    ncap[s] ~ dbin(pcap, N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])      # Part 3 of HM
    log(lambda[s]) <- beta0 + beta1*elevation[s] + beta2*chaparral[s] + beta3*chaparral[s]*chaparral[s] + site.eff[s] # linear model for abundance
    site.eff[s] ~ dnorm(0, tau)   # Site log-normal 'residuals'
    }
    # Derived params
    Ntotal <- sum(N[])
    area<- nsites*3.141*300*300/10000   # Total area sampled, ha
    D<- Ntotal/area
    }
    ",fill=TRUE, file="model5.txt")


# Inits
Nst <- ncap + 1
inits <- function(){list (sigma = runif(1, 30, 100), beta0 = 0, beta1 = 0, beta2 = 0, beta3 = 0, N = Nst, sigma.site = 0.2)}

# Params to save
params <- c("sigma", "beta0", "beta1", "beta2", "beta3", "sigma.site", "Ntotal","D")

# MCMC settings
ni <- 52000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run BUGS (ART 0.9 min) and summarize posteriors
out5 <- bugs(win.data, inits, params, "model5.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni, debug=TRUE, bugs.dir = bd)
out5 <- jags(win.data, inits, params, "model5.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb, n.iter=ni)

# Run JAGS (ART 0.5 min) and summarize posteriors
print(out5, 3)



# 8.6 Summary
# ------------------------------------------------------------------------










# =========================================================================
#
# 9. Advanced Hierarchical Distance Sampling
#
# =========================================================================





# 9.1 Introduction
# ------------------------------------------------------------------------



# 9.2 Distance sampling with clusters, groups, or other individual covariates
# ------------------------------------------------------------------------


# 9.2.1 Simulating HDS data with group size
# ------------------------------------------------------------------------
# Function to simulate data under HDS protocol with groups
simHDSg <- function(type = "line", nsites = 100, lambda.group = 0.75, alpha0 = 0, alpha1 = 0.5, beta0 = 1, beta1 = 0.5, B = 4, discard0 = TRUE){
  #
  # Function simulates hierarchical distance sampling (HDS) data for groups under
  #   either a line (type = "line") or a point (type = "point") transect protocol
  #   and using a half-normal detection function (Buckland et al. 2001).
  #   Other function arguments:
  #     nsites: Number of sites (spatial replication)
  #     lambda.group: Poisson mean of group size
  #     alpha0, alpha1: intercept and slope of log-linear model relating sigma of
  #        half-normal detection function to group size
  #     beta0, beta1: intercept and slope of log-linear model relating the Poisson
  #        mean of the number of groups per unit area to habitat
  #     B: strip half width
  #
  # Get covariates
  habitat <- rnorm(nsites)              # Simulated covariate
  
  # Simulate abundance model for groups (Poisson GLM for N)
  lambda <- exp(beta0 + beta1*habitat)  # Density of groups per "square"
  N <- rpois(nsites, lambda)            # site-specific number of groups
  N.true <- N                           # for point: inside of B
  
  # Simulate observation model
  data <- groupsize <- NULL
  
  for(i in 1:nsites){
    if(N[i]==0){
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # save site, y=1, u, v, d
      next
    }
    
    if(type=="line"){
      # Simulation of distances, uniformly, for each individual in the population
      d <- runif(N[i], 0, B)
      gs <- rpois(N[i],lambda.group) +1  # Observable group sizes >= 1
      groupsize<-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))  # Subtract 1 for interpretation
      # Detection probability for each group
      p <- exp(-d*d/(2*(sigma.vec^2)))
      # Determine if individuals are captured or not
      y <- rbinom(N[i], 1, p)
      u1 <- u2 <- rep(NA,N[i])
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    
    if(type=="point"){
      # Simulation of data on a circle of radius B (algorithm of Wallin)
      angle <- runif(N[i], 0, 360)
      r2 <- runif(N[i], 0, 1)
      r <-  B*sqrt(r2)
      u1 <-  r*cos(angle)  + B
      u2 <-  r*sin(angle)  + B
      
      d <- sqrt((u1 - B)^2 + (u2-B)^2)
      N.true[i] <- sum(d<= B)    # Population size inside of count circle, should be N[i] here.
      gs <- rpois(N[i], lambda.group) + 1
      groupsize <-c(groupsize,gs)
      sigma.vec <- exp(alpha0 + alpha1*(gs-1))
      # For counting individuals on a circle so we truncate p here
      p <- ifelse(d<(B), 1, 0)*exp(-d*d/(2*(sigma.vec^2)))
      y <- rbinom(N[i], 1, p)
      # Subset to "captured" individuals only
      d <- d[y==1] ; u1 <- u1[y==1] ; u2 <- u2[y==1] ; gs <- gs[y==1] ; y <- y[y==1]
    }
    # Now compile things into a matrix and insert NA if no individuals were
    # captured at site i. Coordinates (u,v) are preserved.
    if(sum(y) > 0)
      data <- rbind(data,cbind(rep(i, sum(y)), y, u1, u2, d, gs))
    else
      data <- rbind(data,c(i,NA,NA,NA,NA,NA)) # make a row of missing data
  }
  # Subset to sites at which individuals were captured. You may or may not
  #  do this depending on how the model is formulated so be careful.
  if(discard0)
    data <- data[!is.na(data[,2]),]
  
  # Visualisation
  if(type=="line"){       # For line transect
    par(mfrow = c(1, 3))
    hist(data[,"d"], col = "lightblue", breaks = 20, main =
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  if(type=="point"){       # For point transect
    par(mfrow = c(2,2))
    plot(data[,"u1"], data[,"u2"], pch = 16, main =
           "Located groups in point transects", xlim = c(0, 2*B),
         ylim = c(0, 2*B), col = data[,1], asp = 1)
    points(B, B, pch = "+", cex = 3)
    library(plotrix)
    draw.circle(B, B, B)
    hist(data[,"d"], col = "lightblue", breaks = 20, main =
           "Frequency of distances to groups", xlab = "Distance")
    ttt <- table(data[,1])
    n <- rep(0, nsites)
    n[as.numeric(rownames(ttt))] <- ttt
    plot(habitat, n, main = "Observed group counts (n) vs. habitat", frame = F)
    plot(table(data[,"gs"]), main = "Observed group sizes", ylab = "Frequency", frame = F)
  }
  
  # Output
  list(type = type, nsites = nsites, lambda.group = lambda.group, alpha0 = alpha0, alpha1 = alpha1, beta0 = beta0, beta1 = beta1, B = B, data=data, habitat=habitat, N = N, N.true = N.true, groupsize=groupsize)
}


data <- simHDSg(type = "line")     # Defaults for line transect data
data <- simHDSg(type = "point")    # Default for point transect data
data <- simHDSg(lambda.group = 5)  # Much larger groups
data <- simHDSg(lambda.group = 5, alpha1 = 0) # No effect of groups size on p


# 9.2.2 Analysis in BUGS
# ------------------------------------------------------------------------
set.seed(1234)                 # we all create same data set
temp <- simHDSg(type="line")   # Execute function
data <- temp$data              # harvest data
B <- temp$B                    # Get strip half width
habitat <- temp$habitat        # habitat covariate
nsites <- temp$nsites          # Number of spatial replicates
groupsize <- data[,"gs"] -1    # Input groupsize-1 as data

M <- 400                        # Size of augmented data set is M
nz <- M-nrow(data)              # Number of "pseudo-groups" added
y <- c(data[,2],rep(0,nz))      # Indicator of capture (== 1 for all obs. groups)
nind <- nrow(data)              # Number of observed groups
site <- c(data[,1], rep(NA,nz)) # Site they belong to is unknown
d <- c(data[,5], rep(NA,nz))    # Their distance data are missing ...
groupsize <- c(groupsize, rep(NA,nz)) # .... as is their size
zst <- y                        # Starting values for data augmentation variable

# Bundle data and produce summary
str(bugs.data <- list (y=y, B=B, nind=nind, nsites=nsites, d=d, habitat=habitat,
                       site=site, nz=nz, groupsize=groupsize))


# Define model in BUGS langauge
cat("
    model{
    
    # Prior distributions for model parameters
    alpha0 ~ dunif(-10,10)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-10,10)
    beta1 ~ dunif(-10,10)
    lambda.group ~ dgamma(0.1, 0.1)
    # psi is a derived parameter
    psi <- sum(lambda[])/(nind+nz)
    
    # Individual level model: observations and process
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)                   # Data augmentation variables
    d[i] ~ dunif(0, B)                  # Distance is uniformly distributed
    groupsize[i] ~ dpois(lambda.group)  # Group size is Poisson
    
    log(sigma[i]) <- alpha0 +  alpha1*groupsize[i]
    mu[i] <- z[i]*exp(-d[i]*d[i]/(2*sigma[i]*sigma[i])) #p dep on dist class
    # here using the half normal detection function (Buckland et al. 2001)
    y[i] ~ dbern(mu[i])
    
    site[i] ~ dcat(site.probs[1:nsites]) # Population distribution among sites
    zg[i]<- z[i]*(groupsize[i] + 1)      # Number of individuals in that group
    }
    
    for(s in 1:nsites){
    # Model for population size of groups
    N[s] ~ dpois(lambda[s])
    log(lambda[s])<- beta0 + beta1*habitat[s]
    site.probs[s]<- lambda[s]/sum(lambda[])
    }
    
    # Derived quantities
    G <- sum(z[])        # Total number of groups
    Ntotal <- sum(zg[])  # Total population size (all groups combined)
    }
    ",fill=TRUE, file="model1.txt")

# Load some libraries, define MCMC settings, inits function and parameters to save
library("R2WinBUGS")
library("jagsUI")  #
ni <- 6000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3
inits <- function(){list(alpha0=0, alpha1=0.5, beta0=0, beta1=0, z=zst)}
params <- c("alpha0", "alpha1", "beta0", "beta1", "psi", "Ntotal", "G",
            "lambda.group")

# Call JAGS (ART 1.4 min), check convergence and summarize posterior distributions
out1 <- jags(bugs.data, inits, params, "model1.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)
traceplot(out1)   ;   print(out1, 3)


# 9.2.3 Imperfect observation of cluster size
# ------------------------------------------------------------------------



# 9.3 Time-removal and distance sampling combined
# ------------------------------------------------------------------------


# 9.3.1 The four-part hierarchical model
# ------------------------------------------------------------------------


# 9.3.2 Simulating some time-removal/DS data
# ------------------------------------------------------------------------
# Obtain a data set and harvest the results
set.seed(1235)                 # so we all create the same data set
temp <- simHDStr(type="point") # Simulate point count-removal data set
data <- temp$data              # harvest data
B <- temp$B                    # upper limit of counting (maximum count distance)
nsites <- temp$nsites          # Number of sites
habitat <- temp$habitat        # habitat covariate
K <- temp$K                    # Number of removal periods


# Create the observed encounter frequencies per site (include the zeros! )
data <- data[!is.na(data[,2]),]   # Sites where detections did occur
n <- rep(0,nsites)                # The full site vector
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])  # Put in the counts
site <- data[,1]
nobs <- nrow(data)

# Create the distance class data
nD <- 10             # Number of distance classes
delta <- B/nD        # bin size or width
mdpts <- seq(delta/2,B,delta) # midpoint distance of bins up to max distance
dclass <- data[,"d"] # distance class for each observation
dclass <- dclass%/%delta  +1
tint <- data[,"aux"]

# Bundle data and summarize
str( win.data<-list(n=n, site=site, dclass=as.numeric(dclass),nsites=nsites, nobs=nobs, delta=delta, nD=nD,mdpts=mdpts,B=B, K=K, tint=tint, habitat=habitat) )


cat("
    model {
    # Prior distributions for basic parameters
    # Intercepts
    beta.a0 ~ dnorm(0,0.01)    # intercept for availability
    alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
    alpha1 ~ dnorm(0,0.01)     # slope on sigma covariate
    # Coefficients
    # beta.a1 ~ dnorm(0,0.01)  # slope for availability covariate
    beta0 ~ dnorm(0,0.01)      # intercept for lambda
    beta1 ~dnorm(0,0.01)       # slope for lambda covariate
    
    for(s in 1:nsites){
    # Add covariates to scale parameter DISTANCE (perceptibility)
    log(sigma[s]) <- alpha0 +  alpha1*habitat[s]
    # Add covariates for availability here TIME-REMOVAL (availability)
    p.a[s] <- exp(beta.a0) / (1+exp(beta.a0))
    # Optional covariates on availability
    # exp(beta.a0 + beta.a1*date[s])/(1+exp(beta.a0+beta.a1*date[s]))
    # Distance sampling detection probability model
    for(b in 1:nD){
    log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # Half-normal
    f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Radial density function
    pi.pd[b,s] <- g[b,s]*f[b,s]  #  Product Pr(detect)*Pr(distribution)
    pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]  # Conditional probabilities
    }
    pdet[s] <- sum(pi.pd[,s])  # Probability of detection at all
    
    # Time-removal probabilities
    for (k in 1:K){
    pi.pa[k,s] <- p.a[s] * pow(1-p.a[s], (k-1))
    pi.pa.c[k,s] <- pi.pa[k,s]/phi[s] # Conditional probabilities of availability
    }
    phi[s] <- sum(pi.pa[,s]) # Probability of ever available
    }
    # Conditional observation model for categorical covariates
    for(i in 1:nobs){
    dclass[i] ~ dcat(pi.pd.c[,site[i]])
    tint[i] ~ dcat(pi.pa.c[,site[i]])
    }
    # Abundance model
    for(s in 1:nsites){
    # Binomial model for # of captured individuals
    # n[s] ~ dbin(pmarg[s], M[s]) # Formulation b, see text
    # pmarg[s] <- pdet[s]*phi[s]
    n[s] ~ dbin(pdet[s], N[s])    # Formulation a, see text
    N[s] ~ dbin(phi[s],M[s])      # Number of available individuals
    M[s] ~ dpois(lambda[s])       # Abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s]) <- beta0 + beta1*habitat[s]
    }
    # Derived quantities
    Mtot <- sum(M[])  # Total population size
    Ntot <- sum(N[])  # Total available population size
    PDETmean <- mean(pdet[]) # Mean perceptibility across sites
    PHImean <- mean(phi[]) # Mean availability across sites
    }
    ", fill=TRUE, file="tr-ds.txt")

# Create initial values (including for M and N) and list parameters to save
Mst <- Nst <- n + 1
inits <- function(){list(M=Mst, N=Nst, alpha0=1, beta0=runif(1,-1,1),
                         beta.a1=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
                         beta.a0=runif(1,-1,1))}
params <- c("beta.a0", "beta.a1", "alpha0", "alpha1", "beta0", "beta1", "PDETmean", "PHImean", "Mtot", "Ntot")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 7.3 min), check convergence and summarize posteriors
out2a <- jags(data=win.data, inits=inits, parameters=params,
              model.file ="tr-ds.txt",n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
              parallel = TRUE)
traceplot(out2a)   ;   print(out2a, 3)

sum(temp$M)

print(out2b,3)



# 9.4 Mark-Recapture/Double observer Distance Sampling
# ------------------------------------------------------------------------


# 9.4.1 Simulating MRDS data
# ------------------------------------------------------------------------
# Simulate a double-observer sampling data set
set.seed(1235)
temp <- simHDStr(type="point", method="double") # simulate double observer point count data set
data <- temp$data         # harvest data
B <- temp$B               # upper limit of counting (maximum count distance)
nsites <-temp$nsites      # number of sites
habitat <-temp$habitat    # habitat covariate


# Processing of the data: pad the count vector with 0s etc.
data <- data[!is.na(data[,2]),]
n <- rep(0,nsites)
names(n) <- 1:nsites
n[names(table(data[,1]))] <- table(data[,1])
site <- data[,1]
dclass <- data[,"d"]      # categorical distance class for each observation
aux <- data[,"aux"]       # the auxiliary variable is capture history

# Create the categorical distance variable, use 10 classes here.
nD <- 10
delta <- B/nD # bin width
mdpts <-seq(delta/2,B,delta) # midpoint of bins up to max distance
nobs <- nrow(data)
dclass <- dclass%/%delta  +1

# Bundle data and look at overview of data
str( win.data <-list(n=n,site=site, dclass=as.numeric(dclass), nsites=nsites,
                     nobs=nobs, delta=delta, nD=nD, mdpts=mdpts, B=B, aux=aux, habitat=habitat) )


# 9.4.2 Analysis in BUGS
# ------------------------------------------------------------------------
# Define model in BUGS langauge
cat("
    model {
    
    #Priors for fixed detection parameters
    # 2 observer detection probability parameters
    logitp1 ~ dnorm(0, 0.01)
    logitp2 ~ dnorm(0, 0.01)
    # Intercepts
    alpha0 ~ dnorm(0, 0.01)    # intercept for sigma
    alpha1 ~ dnorm(0, 0.01)    # slope on sigma covariate
    # Coefficients
    beta0 ~ dnorm(0,0.01)      # intercept for lambda
    beta1 ~ dnorm(0,0.01)      # slope for lambda covariate
    
    # Detection scale parameter model
    for(s in 1:nsites){
    # Covariates on scale parameter (perceptibility)
    log(sigma[s]) <- alpha0 + alpha1*habitat[s]
    # Double observer cell probabilities here if there are covariates
    logit(pobs[1,s]) <- logitp1 # + covariates
    logit(pobs[2,s]) <- logitp2 # + covariates
    
    # Distance sampling model and cell probabilities
    for(b in 1:nD){
    log(g[b,s]) <- -mdpts[b]*mdpts[b]/(2*sigma[s]*sigma[s])  # half-normal
    f[b,s] <- ( 2*mdpts[b]*delta )/(B*B) # Scaled radial density function
    pi.pd[b,s] <- g[b,s]*f[b,s]          #  Product Pr(detect)*Pr(distribution)
    pi.pd.c[b,s] <- pi.pd[b,s]/pdet[s]   # Conditional cell probabilities
    }
    pdet[s] <- sum(pi.pd[,s])              # Marginal probability of detection
    
    # Double observer cell probabilities and conditional probabilities
    doprobs[1,s] <- pobs[1,s]*(1-pobs[2,s])
    doprobs.condl[1,s] <- doprobs[1,s]/sum(doprobs[,s])
    doprobs[2,s] <- (1-pobs[1,s])*pobs[2,s]
    doprobs.condl[2,s] <- doprobs[2,s]/sum(doprobs[,s])
    doprobs[3,s] <- pobs[1,s]*pobs[2,s]
    doprobs.condl[3,s] <- doprobs[3,s]/sum(doprobs[,s])
    pavail[s] <- sum(doprobs[,s])  # probability of availability AT ALL
    }
    
    # Observation model for two categorical covariates
    for(i in 1:nobs){
    dclass[i] ~ dcat(pi.pd.c[,site[i]])
    aux[i] ~ dcat(doprobs.condl[,site[i]])
    }
    
    # Abundance model
    for(s in 1:nsites){
    # Binomial model for # of captured individuals
    n[s] ~ dbin(pdet[s], N[s])
    N[s] ~ dbin(pavail[s], M[s])   # binomial availability model
    # Abundance model
    M[s] ~ dpois(lambda[s])        # predicted abundance per survey/site/point
    # Add site-level covariates to lambda
    log(lambda[s])<- beta0 + beta1*habitat[s]
    }
    # Derived parameters
    Mtot <- sum(M[])
    Ntot <- sum(N[])
    logit(p1) <- logitp1
    logit(p2) <- logitp2
    sigma0 <- exp(alpha0)            # Baseline sigma
    }
    ", fill=TRUE,file="do_model.txt")

# Inits function
Nst <- n + 1         # inits for N
inits <- function(){list(M=Nst+1, N=Nst, alpha0=runif(1,1,2),
                         beta0=runif(1,-1,1), beta1=runif(1,-1,1), alpha1=runif(1,-1,1),
                         logitp1=0, logitp2=0)}

# Parameters to monitor
params <- c("alpha0", "alpha1", "beta0", "beta1", "Ntot", "Mtot", "logitp1",
            "logitp2", "p1", "p2", "sigma0")

# MCMC settings
ni <- 50000   ;   nb <- 10000   ;   nt <- 4   ;   nc <- 3

# Run JAGS in parallel (ART 6.8 min), check convergence and summarize the results
out3 <- jags(data=win.data, inits=inits, parameters.to.save=params,
             model.file="do_model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni,
             parallel = TRUE)
traceplot(out3)   ;   print(out3, 3)


# Put true values into a vector
truth <- temp$parms
psi <- 1-(1-truth["p.double1"])*(1-truth["p.double2"])# Compute availability
truth <- c(truth[c("p.double1", "p.double2")], exp(truth["alpha0"]),
           truth["beta0"], "Mtot" = sum(temp$M),
           "Ntot" = sum(temp$M)*as.numeric(psi),
           truth[c("alpha0","alpha1","beta1")])

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- out3$summary[c("p1", "p2", "sigma0", "beta0", "Mtot", "Ntot",
                       "alpha0", "alpha1", "beta1"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.4.3 Remarks
# ------------------------------------------------------------------------



# 9.5 Open HDS models: temporary emigration
# ------------------------------------------------------------------------


# 9.5.1 Data and model structure
# ------------------------------------------------------------------------


# 9.5.2 Cautionary note on temporary emigration processes
# ------------------------------------------------------------------------


# 9.5.3 Modeling temporary emigration with distance sampling in unmarked using the function gdistsamp
# ------------------------------------------------------------------------
# Load the wagtail data, investigate NA patterns in detection data
data("wagtail")
str(wagtail)
Y <- wagtail$Y
table(n.missing <- rowSums(is.na(Y))) # Frequency distribution of number of NAs per site
n.missing
keep <- which(n.missing == 0)     # Sites with complete distance data
Y <- Y[keep,]                     # restrict analysis to those

# Harvest other data for sites with complete distance data
potato <- wagtail$potato[keep]   ;   grass <- wagtail$grass[keep]
lscale <- wagtail$lscale[keep]   ;   hour <- wagtail$hour[keep,]
date <- wagtail$date[keep,]   ;   rep <- wagtail$rep[keep,]
breaks <- wagtail$breaks

# Look at the distance data
str(Y)
tmp <- apply(Y, 2, sum, na.rm = T)
matplot(1:6, t(matrix(tmp, nrow = 4, byrow= T)), type = "b", ylim = c(0, 90), xlab = "Distance class", ylab = "Number of wagtails", frame = F, lwd = 3, lty = 1)

# Standardize all continuous covariates
mn.potato <- mean(potato)   ;   sd.potato <- sd(potato)
mn.grass <- mean(grass)   ;   sd.grass <- sd(grass)
mn.lscale <- mean(lscale)   ;   sd.lscale <- sd(lscale)
mn.date <- mean(date)    ;   sd.date <- sd(c(date))
mn.hour <- mean(hour)   ;   sd.hour <- sd(c(hour))
POTATO <- (potato - mn.potato) / sd.potato
GRASS <- (grass - mn.grass) / sd.grass
LSCALE <- (lscale - mn.lscale) / sd.lscale
DATE <- (date - mn.date) / sd.date
HOUR <- (hour - mn.hour) / sd.hour

# Package into unmarked GDS data frame and inspect the data
umf <- unmarkedFrameGDS(y = Y[,1:24], survey="point", unitsIn="m",
                        dist.breaks=breaks, numPrimary = 4,
                        siteCovs = data.frame(POTATO, GRASS, LSCALE),
                        yearlySiteCovs=list(rep = rep, DATE = DATE, HOUR = HOUR))
str(umf)
summary(umf)


# Model fitting: Null models fm0
# exponential detection function
summary(fm0.exp <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                             keyfun = "exp", output = "density", unitsOut = "ha",
                             mixture = "P", K = 100, se = TRUE, data = umf) )

# hazard detection function
summary(fm0.haz <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                             keyfun = "haz", output = "density", unitsOut = "ha",
                             mixture = "P", K = 100, se = TRUE, data = umf ) )


# half-normal detection function
summary(fm0.hn <- gdistsamp(lambdaformula = ~1, phiformula = ~1, pformula = ~1,
                            keyfun = "halfnorm", output = "density", unitsOut = "ha",
                            mixture = "P", K = 100, se = TRUE, data = umf,control=list(trace=TRUE, REPORT=1)) )

# Compare AIC scores for 3 detection functions
rbind('AIC exp' = fm0.exp@AIC, 'AIC haz' = fm0.haz@AIC, 'AIC hn' = fm0.hn@AIC)

backTransform(fm0.haz, type="lambda")

backTransform(fm0.haz, type="det")

plot(1:300, gxhaz(1:300, shape = exp(5.13), scale=1.59), frame = F, type = "l", xlab = "Distance Wagtail???Observer (metres)", ylab = "Detection probability", lwd=3)


# Model with time-dependent phi
fm1 <- gdistsamp(lambdaformula = ~1, phiformula = ~rep-1,
                 pformula = ~1, keyfun = "haz", output = "density", unitsOut = "ha",
                 mixture = "P", K = 100, se = TRUE, data = umf)

# Compare AIC for models with phi constant and phi time-dependent
rbind('AIC phi constant' = fm0.haz@AIC, 'AIC phi time-dep' = fm1@AIC)


# Add covariates on lambda:  2-phase fitting to assist convergence using
#   first K = 20
summary( fm2.init <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                               phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                               unitsOut = "ha",   control=list(trace=TRUE, REPORT=1),
                               mixture = "P", K = 20, se = TRUE, data = umf))

starts <- coef(fm2.init)
summary( fm2  <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                           phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                           unitsOut = "ha",  starts=starts, control=list(trace=TRUE,
                                                                         REPORT=1), mixture = "P", K = 100, se = TRUE, data = umf))


# Add covariates on lambda: optimisation with default, no inits
summary( fm2a <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                           phiformula = ~rep-1, pformula = ~1, keyfun = "haz", output = "density",
                           unitsOut = "ha", control=list(trace=TRUE, REPORT=1),
                           mixture = "P", K = 100, se = TRUE, data = umf))


# Estimates of availability probabilities
plogis(fm2@estimates@estimates$phi@estimates)

# Models with time-dependent phi, AIC-best key functions and
#    three covariates on phi. Use previous estimates as starting values.
(tmp <- coef(fm2))
starts <- c(tmp[1:4], tmp[5:8], 0,0,0, tmp[9], tmp[10])

summary(fm3 <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                         phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~1,
                         keyfun = "haz", output = "density", unitsOut = "ha",
                         mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                         se = TRUE, data = umf, starts = starts))

# Models with time-dependent phi, AIC-best key function, 3 covariates on phi
#     and, in addition, date and hour on detection
# linear effects on detection
tmp <- fm3@estimates@estimates
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, tmp$det@estimates, 0, 0, tmp$scale@estimates)
summary(fm4A <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                          phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE, pformula = ~ DATE + HOUR,
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

# quadratic effects on detection
tmp <- fm4A@estimates@estimates
p.start <- tmp$det@estimates
p.start <- c(p.start[1:2], 0, p.start[3], 0)
starts <- c(tmp$lambda@estimates, tmp$phi@estimates, p.start, tmp$scale@estimates)

summary(fm4B <- gdistsamp(lambdaformula = ~ POTATO+GRASS+LSCALE,
                          phiformula = ~(rep-1)+ POTATO+GRASS+LSCALE,
                          pformula = ~ DATE + I(DATE^2) + HOUR + I(HOUR^2),
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

starts <- coef(fm4B)[-16]   # Drop coef for HOUR^2
summary(fm4C <- gdistsamp(~ POTATO+GRASS+LSCALE,
                          ~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
                          keyfun = "haz", output = "density", unitsOut = "ha",
                          mixture = "P", K = 100, control=list(trace=TRUE, REPORT=1),
                          se = TRUE, data = umf, starts = starts) )

starts <- c(coef(fm4C), 0)
summary(fm5 <- gdistsamp(~ POTATO+GRASS+LSCALE,
                         ~(rep-1)+ POTATO+GRASS+LSCALE, ~ DATE + I(DATE^2) + HOUR,
                         keyfun = "haz", output = "density", unitsOut = "ha",
                         mixture = "NB", K = 100, control=list(trace=TRUE, REPORT=1),
                         se = TRUE, data = umf , starts = starts) )

# Now we create a model selection table of these various models
modSel(fitList(fm0.haz, fm1, fm2, fm3, fm4A, fm4B, fm4C, fm5) )

summary(fm5)

# Bootstrap Goodness-of-fit assessment: ART ~ 20 hours
set.seed(1234)
(pb <- parboot(fm5, fitstats, nsim=100, report=5))

# Compute magnitude of "overdispersion" c.hat as ratio of observed to expected
#    chisquare test statistic
(c.hat <- pb@t0[2] / mean(pb@t.star[,2]))  # c-hat as ratio of observed/expected


# Predictions of lambda for POTATO, GRASS and LSCALE
newdat1 <- data.frame(POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat2 <- data.frame(POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat3 <- data.frame(POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred1 <- predict(fm5, type="lambda", newdata=newdat1, append = T)
pred2 <- predict(fm5, type="lambda", newdata=newdat2, append = T)
pred3 <- predict(fm5, type="lambda", newdata=newdat3, append = T)

# Predictions of phi for POTATO, GRASS and LSCALE and for rep = 1
newdat4 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=0, GRASS=0, LSCALE = seq(-1.8,4.33,,100))
newdat5 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=seq(-0.75,3,,100), GRASS=0, LSCALE = 0)
newdat6 <- data.frame(rep = factor('1', levels = c('1','2','3','4')), POTATO=0, GRASS=seq(-0.4, 3.6,,100), LSCALE = 0)
pred4 <- predict(fm5, type="phi", newdata=newdat4, append = T)
pred5 <- predict(fm5, type="phi", newdata=newdat5, append = T)
pred6 <- predict(fm5, type="phi", newdata=newdat6, append = T)

# Predictions of detection function sigma for DATE and HOUR
newdat7 <- data.frame(DATE = seq(-1.51,1.69,,100), HOUR = 0)
newdat8 <- data.frame(DATE=0, HOUR = seq(-1.92,3.1,,100))
pred7 <- predict(fm5, type="det", newdata=newdat7, append = T)
pred8 <- predict(fm5, type="det", newdata=newdat8, append = T)

par(mfrow = c(1,3), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
plot(newdat1$LSCALE, pred1[,1], xlab="Standardized covariate", ylab="Density (birds/ha)", lwd=3,type="l", frame = F)
lines(newdat2$POTATO, pred2[,1], lwd=3, col="red")
lines(newdat3$GRASS, pred3[,1], lwd=3, col="blue")
legend(-1.6, 1.65, c("LSCALE", "POTATO", "GRASS"), col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2)

plot(newdat4$LSCALE, pred4[,1], xlab="Standardized covariate", ylab="Availability (phi)", lwd=3,type="l", frame = F)
lines(newdat5$POTATO, pred5[,1], lwd=3, col="red")
lines(newdat6$GRASS, pred6[,1], lwd=3, col="blue")
legend(2, 0.65, c("LSCALE", "POTATO", "GRASS"), col=c("black", "red", "blue"), lty=1, lwd=3, cex=1.2)

plot(newdat7$DATE, pred7[,1], xlab="Standardized covariate", ylab="Detection function (sigma)", lwd=3,type="l", frame = F, ylim = c(100, 200))
lines(newdat8$HOUR, pred8[,1], lwd=3, col="red")
legend(0.5, 140, c("DATE", "HOUR"), col=c("black", "red"), lty=1, lwd=3, cex=1.2)



# 9.5.4 Fitting temporary emigration HDS models in BUGS
# ------------------------------------------------------------------------


# 9.5.4.1 Simulating a temporary emigration system
# ------------------------------------------------------------------------
simHDSopen(type="line", nsites = 100, mean.lam = 2, beta.lam = 0, mean.sig = 1, beta.sig = 0, B = 3, discard0=TRUE, nreps=2, phi=0.7, nyears=5, beta.trend = 0)

# Obtain a temporary emigration data set
set.seed(1234)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100)  )
attach(tmp)


apply(tmp$M.true,2,sum)

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0,dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels= 1:ndclass)
    y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
  }
}

y3d <- y4d[,,,1]


# Bundle and summarize the data set
nobs <- apply(y3d, c(1,3), sum)  # Total detections per site and occasion
str( data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs) )


# Define model in BUGS
cat("
    model {
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient of lambda on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0.01,5)   # Distance function parameter
    
    # Detection probs for each distance interval and related things
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal
    f[b] <- (2*midpt[b]*delta)/(B*B)    # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1]<- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    pdet[s,k] <- sum(cellprobs[1:nD])   # Distance class probabilities
    pmarg[s,k] <- pdet[s,k]*phi         # Marginal probability
    
    # Model part 4: distance class frequencies
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k])
    # Model part 3: total number of detections:
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])
    # nobs[s,k] ~ dbin(pdet[s,k], Navail[s,k]) # Alternative formulation
    # Model part 2: Availability. Not used in this model but simulated.
    Navail[s,k] ~ dbin(phi, M[s])
    }  # end k loop
    # Model part 1: Abundance model
    M[s] ~ dpois(lambda[s])
    log(lambda[s]) <- beta0 + beta1*habitat[s]
    }  # End s loop
    
    # Derived quantities
    Mtot <- sum(M[])
    for(k in 1:K){
    Ntot[k]<- sum(Navail[,k])
    }
    } # End model
    ",file="model.txt")


# Assemble the initial values and parameters to save for JAGS
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max)  +2
inits <- function(){
  list(M=Mst, sigma = 1.0, phi=0.9, beta0=log(2), beta1=0.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 60000   ;   nb <- 10000   ;   nt <- 5   ;   nc <- 3

# Run WinBUGS or JAGS
library("R2WinBUGS")
library("jagsUI")  # JAGS works but WinBUGS does not!
# bd <- "c:/WinBUGS14/"
# out1 <- bugs(data, inits, parameters, "model.txt", n.thin=nthin,
#       n.chains=nc, n.burnin=nb,n.iter=ni,debug=TRUE, bugs.dir = bd)
# We get this error: vector valued relation y3d must involve consecutive
# elements of variable

# Run JAGS: This fails quite often ('invalid parent node'), just keep trying
outTE1 <- jags(data, inits, params, "model.txt", n.thin=nt,n.chains=nc,
               n.burnin=nb,n.iter=ni, parallel = TRUE)
traceplot(outTE1)   ;    print(outTE1, 3)            # ART 4 min


# Put true values into a vector
truth <- c(tmp$parms[c(1:3,5)], Mtot = sum(tmp$M[,1]),
           Ntot = (apply(tmp$Na.real[,,1],2,sum)))

# Get posterior means and 2.5% and 97.5% percentiles (95% CRI)
post <- outTE1$summary[c("mean.lam", "beta1", "sigma", "phi", "Mtot", "Ntot[1]", "Ntot[2]", "Ntot[3]" ,"Ntot[4]", "Ntot[5]", "Ntot[6]", "Ntot[7]"), c(1,3,7)]

# Table compares truth with posterior mean and 95% CRI from JAGS
cbind(truth, posterior = round(post, 3))


# 9.5.4.2 Bayesian analysis of the Wagtail data
# ------------------------------------------------------------------------
y3d <- array(NA,dim=c(nrow(Y), 6, 4) )          # Create 3d array
y3d[,,1] <- Y[,1:6]  ;  y3d[,,2] <- Y[,7:12]    # Fill the array
y3d[,,3] <- Y[,13:18]  ;  y3d[,,4] <- Y[,19:24]

K <- 4                          # Number of primary occasions
nsites <- nrow(Y)               # Number of sites
nD <- 6                         # Number of distance classes
midpt <- seq(25,275,50)         # Class midpoint distance
delta <- 50                     # Class width
B <- 300                        # Maximum distance
nobs <- apply(y3d, c(1,3), sum) # Total detections per site and occasion

# Bundle and summarize data set
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta,
                 B=B, nobs=nobs, area=area))

# Write out the BUGS model file
cat("
    model {
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    
    # Availability parameter
    phi ~ dunif(0,1)
    
    # Detection parameter
    sigma ~ dunif(0,500)
    
    # Multinomial cell probabilities
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma)  # Half-normal model
    f[b] <- (2*midpt[b]*delta)/(B*B) # Scaled radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    # Conditional 4-part version of the model
    pdet[s,k] <- sum(cellprobs[1:nD])
    pmarg[s,k] <- pdet[s,k]*phi
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k]) # Part 4: distance
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: number of detected individuals
    Navail[s,k] ~ dbin(phi,M[s])        # Part 2: Number of available individuals
    }  # end k loop
    
    M[s] ~ dpois(lambda[s])    #  Part 1: Abundance model
    log(lambda[s]) <- beta0    #  Habitat variables would go here
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
    Davail[k] <- phi*exp(beta0)/area
    }
    Mtotal <- sum(M[])
    Dtotal<- exp(beta0)/area
    } # end model
    ",fill=TRUE,file="wagtail.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma = 100.0)

# Parameters to save
params <- c("sigma", "phi", "beta0", "beta1", "Mtotal", "Davail", "Dtotal")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 3

# Run JAGS (ART 3 min)
library("jagsUI")
wag1 <- jags(data, inits, params, "wagtail.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
par(mfrow = c(3,3))    ;   traceplot(wag1)
summary(wag1)

exp(-2.06)


# Bundle and summmarize data set for BUGS
rep <- matrix(as.numeric(rep), ncol=4)
area <- pi*(300^2)/10000
str(data <- list(y3d=y3d, nsites=nsites, K=K, nD=nD, midpt = midpt, delta=delta, B=B, nobs=nobs, POTATO=POTATO, GRASS=GRASS, LSCALE=LSCALE, rep=rep, DATE=DATE,
                 HOUR=HOUR, area=area))

# Define model in BUGS
cat("
    model {
    
    # Priors
    # Abundance parameters
    beta0 ~ dnorm(0, 0.01)
    beta1 ~ dnorm(0, 0.01)
    beta2 ~ dnorm(0, 0.01)
    beta3 ~ dnorm(0, 0.01)
    
    # Availability parameters
    phi0 ~ dunif(0,1)
    logit.phi0 <- log(phi0/(1-phi0))
    for(k in 1:4){
    gamma1[k] ~ dunif(0, 1) # Availability effects of surveys 1 - 4
    logit.gamma1[k]<- log(gamma1[k]/(1-gamma1[k]))
    }
    gamma2 ~ dnorm(0, 0.01)
    gamma3 ~ dnorm(0, 0.01)
    gamma4 ~ dnorm(0, 0.01)
    
    # Detection parameters
    sigma0 ~ dunif(0.1,500)   # Intercept
    alpha2 ~ dnorm(0, 0.01)   # effect of DATE (linear)
    alpha3 ~ dnorm(0, 0.01)   # effect of DATE (squared)
    alpha4 ~ dnorm(0, 0.01)   # effect of HOUR
    theta ~ dgamma(0.1, 0.1)
    r ~ dunif(0, 10)
    
    for (s in 1:nsites) {
    for (k in 1:K) {
    # Availability parameter
    logit.phi[s,k] <- logit.gamma1[k] + gamma2*POTATO[s] + gamma3*GRASS[s] + gamma4*LSCALE[s]
    phi[s,k] <- exp(logit.phi[s,k])/(1+ exp(logit.phi[s,k]))
    # Distance sampling parameter
    log(sigma[s,k]) <- log(sigma0) + alpha2*DATE[s,k] + alpha3*pow(HOUR[s,k],2) + alpha4*HOUR[s,k]
    # Multinomial cell probability construction
    for(b in 1:nD){
    #log(g[s,b,k]) <- -midpt[b]*midpt[b]/(2*sigma[s,k]*sigma[s,k]) # half-normal
    cloglog(g[s,b,k]) <- theta*log(sigma[s,k])  - theta*log(midpt[b])  # hazard
    f[s,b,k] <- (2*midpt[b]*delta)/(B*B)
    cellprobs[s,b,k] <- g[s,b,k]*f[s,b,k]
    cellprobs.cond[s,b,k] <- cellprobs[s,b,k]/sum(cellprobs[s,1:nD,k])
    }
    cellprobs[s,nD+1,k]<- 1-sum(cellprobs[s,1:nD,k])
    
    #  Conditional 4-part hierarchical model
    pdet[s,k] <- sum(cellprobs[s,1:nD,k])
    pmarg[s,k] <- pdet[s,k]*phi[s,k]
    y3d[s,1:nD,k] ~ dmulti(cellprobs.cond[s,1:nD,k], nobs[s,k]) # Part 4
    nobs[s,k] ~ dbin(pmarg[s,k], M[s])  # Part 3: Number of detected individuals
    Navail[s,k] ~ dbin(phi[s,k],M[s])   # Part 2: Number of available individuals
    } # end k loop
    
    M[s] ~ dnegbin(prob[s], r)
    prob[s] <- r/(r+lambda[s])
    # M[s] ~ dpois(lambda[s])         # Part 1: Abundance model
    log(lambda[s]) <- beta0 + beta1*POTATO[s] + beta2*GRASS[s] + beta3*LSCALE[s]
    }  # end s loop
    
    # Derived quantities
    for(k in 1:K){
    Davail[k] <- mean(phi[,k])*exp(beta0)/area
    }
    Mtotal <- sum(M[])
    Dtotal <- exp(beta0)/area
    } # End model
    ",fill=TRUE,file="wagtail2.txt")

# Inits
Navail.st <- apply(y3d, c(1,3),sum)
Mst <- apply(Navail.st, c( 1), max,na.rm=TRUE) + 2
inits <- function() list(M=Mst, sigma0 = 100, alpha2=0, alpha3=0, alpha4=0, gamma2=0, gamma3=0, gamma4=0, beta1=0,beta2=0,beta3=0, r = 1)

# Parameters to save
params <- c("r","sigma0", "beta0", "beta1", "beta2", "beta3", "Mtotal", "alpha2", "alpha3",  "alpha4", "theta", "Dtotal", "Davail", "phi0", "gamma1", "gamma2", "gamma3", "gamma4" ,"logit.gamma1")

# MCMC settings
ni <- 32000   ;   nb <- 2000   ;   nt <- 2   ;   nc <- 5

# Run JAGS (ART 79 min), check convergence and summarize posteriors
wag2 <- jags(data, inits, params, "wagtail2.txt", n.thin=nt,n.chains=nc,
             n.burnin=nb,n.iter=ni, parallel = TRUE)

# Compare posterior means to MLEs obtained from unmarked
mle <- coef(fm5)
mle[12] <- exp(mle[12]) # convert to distance units
mle[16] <- exp(mle[16]) # back-transform the hazard parameter
mle[17] <- exp(mle[17]) # back-transform the NB dispersion parameter
bayes <- wag2$summary[,1]
bayes <- c(bayes[3:6],bayes[c(25:28,22:24)],bayes[c(2,8:10)],bayes[11],bayes[1])
bayes[1] <- log(exp(bayes[1])/area)# Convert from N per site to log(density) per ha

round( cbind(mle,bayes), 3)


# 9.5.4.3 Robust Design: Replicates within and among years
# ------------------------------------------------------------------------



# 9.6 Open HDS models: Implicit Dynamics
# ------------------------------------------------------------------------


# Obtain a data set
set.seed(1236)
str(tmp <- simHDSopen("point", nreps=7, nyears=5, nsites=100, beta.trend=0.2) )
attach(tmp)

apply(tmp$M.true,2,sum)  # True population size per year

# Define distance class information
delta <- 0.5
nD <- B%/%delta                 # Number of distance classes
midpt <- seq(delta/2, B, delta) # Mid-point of distance intervals

# Create the 4-d array
y4d <- array(0, dim=c(nsites, nD, K, nyears))
for(yr in 1:nyears){
  for(rep in 1:K){
    data <- tmp$data[[yr]][[rep]]
    site <- data[,1]
    dclass <- data[,"d"]%/%delta + 1
    ndclass <- B%/%delta
    dclass <- factor(dclass, levels=  1:ndclass)
    y4d[1:nsites,1:nD,rep,yr] <- table(site, dclass)
  }
}


# Bundle and summarize the data set
nobs <- apply(y4d, c(1,3,4), sum)  # Total detections per site and occasion
str( data <- list(y4d=y4d, nsites=nsites, K=K, nD=nD, midpt=midpt, delta=delta, habitat=habitat, B=B, nobs = nobs, T=tmp$nyears) )

# Define model in BUGS
cat("
    model {
    
    # Prior distributions
    beta0 ~ dnorm(0, 0.01)  # Intercept for log(lambda)
    mean.lam <- exp(beta0)
    beta1 ~ dnorm(0, 0.01)  # Coefficient on habitat
    phi ~ dunif(0,1)        # Probability of availability
    sigma ~ dunif(0,5)      # Detection function parameter
    beta.trend ~ dnorm(0, 0.01)
    
    # Construct the multinomial cell probabilities
    for(b in 1:nD){
    log(g[b]) <- -midpt[b]*midpt[b]/(2*sigma*sigma) # half-normal
    f[b] <- (2*midpt[b]*delta)/(B*B)                # radial density function
    cellprobs[b] <- g[b]*f[b]
    cellprobs.cond[b] <- cellprobs[b]/sum(cellprobs[1:nD])
    }
    cellprobs[nD+1] <- 1-sum(cellprobs[1:nD])
    for (s in 1:nsites) {
    for (k in 1:K) {
    pdet[s,k] <- sum(cellprobs[1:nD]) # Distance class probabilities
    pmarg[s,k] <- pdet[s,k]*phi       # Marginal probability
    }
    }
    
    for(t in 1:T){                        # Years
    for (s in 1:nsites) {               # Sites
    for (k in 1:K) {                  # Replicates
    # Model part 4: distance class frequencies
    y4d[s,1:nD,k,t] ~ dmulti(cellprobs.cond[1:nD], nobs[s,k,t])
    # Model part 3: total number of detections:
    nobs[s,k,t] ~ dbin(pmarg[s,k], M[s,t])
    # Model part 2: Availability. Not used in this model but simulated.
    Navail[s,k,t] ~ dbin(phi, M[s,t])
    }  # end k loop
    # Model part 1: Abundance model
    M[s,t] ~ dpois(lambda[s,t])
    log(lambda[s,t]) <- beta0 + beta1*habitat[s] + beta.trend*(t-2.5)
    }  # end s loop
    } # end t loop
    
    # Derived quantities
    for(t in 1:T){
    Mtot[t] <- sum(M[,t])
    for(k in 1:K){
    Ntot[k,t] <- sum(Navail[,k,t])
    }
    }
    } # End model
    ",file="tempemig4d.txt")

# Inits and parameters to save
Navail.st <- apply(y4d, c(1,3,4),sum)
Mst <- apply(Navail.st, c( 1,3), max) +2
inits <- function(){
  list(M=Mst, Navail = Navail.st, sigma = 1.0, phi=.9,beta0=log(2),beta1=.5)
}
params <- c("sigma", "phi", "beta0", "mean.lam", "beta.trend",
            "beta1", "Mtot", "Ntot")

# MCMC settings
ni <- 12000   ;   nb <- 2000   ;   nt <- 5   ;   nc <- 3

# Run JAGS (ART 9 min), look at trace plots and summarize
outRD <- jags(data, inits, params, "tempemig4d.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
par(mfrow = c(3,3))   ;   traceplot(outRD)
summary(outRD)



# 9.7 Open HDS models: modelling population dynamics
# ------------------------------------------------------------------------


# 9.7.1 Simulating the ISSJ data over multiple years
# ------------------------------------------------------------------------
# We load the ISSJ data analyzed in chapter 8, package into an unmarked frame
library(unmarked)
library(rjags)
data(issj)
covs <- issj[,c("elevation","forest","chaparral")]
area <- pi*300^2 / 100^2             # Area in ha
jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                           siteCovs=data.frame(covs, area), numPrimary=1,
                           dist.breaks=c(0, 100, 200, 300),
                           unitsIn="m", survey="point")
sc <- siteCovs(jayumf)
sc.s <- scale(sc)
sc.s[,"area"] <- pi*300^2 / 10000  # Don't standardize area
covs<- siteCovs(jayumf) <- sc.s
summary(jayumf)

# Fit the model using gdistsamp and look at the fit summary
(nb.C2E.C <- gdistsamp( ~chaparral + I(chaparral^2) + elevation , ~1, ~chaparral,
                        data =jayumf, output="abund", mixture="NB", K = 150))

# Get coefficient estimates to be used in data simulation
beta <- coef(nb.C2E.C)
betaFall <- beta[c("lambda(Int)", "lambda(chaparral)",
                   "lambda(elevation)", "lambda(I(chaparral^2))")]

# Predict expected abundance per point count on log-scale for simulation
Xmat <- cbind(rep(1,307),covs[,3],covs[,3]^2,covs[,1]) # Order: chap, chap^2, elev
loglam <- Xmat%*%(betaFall)
lamnew <- exp(loglam)

# Parameters of the detection function
dparm <- beta[c("p(Int)", "p(chaparral)")]
sigma <- exp(Xmat[,c(1, 2)]%*%dparm)
J <- nsites <- 307 # number of sampling points

# Number of years
nyrs <- 6

# Set dynamics parameters to achieve a target growth rate of 0.95
phi <- 0.6       # Survival probability
gamma <- 0.35    # Recruitment rate

# Distance category info
db <- c(0,50, 100, 150, 200, 250, 300)
midpt <- c(25, 75, 125, 175, 225, 275)
nD <- length(midpt)
delta <- 50       # Distance interval width
B <- 300

# Simulate an ISSJ data set and harvest the data objects
set.seed(2015)
dat <- issj.sim(B=300, db = db, lam=lamnew, sigma=sigma, phi=phi, gamma=gamma, npoints=nsites, nyrs=nyrs)

y <- dat$y
dclass <- dat$dclass
site <- dat$site

# Bundle and summarize the data set
str(data1<-list(nsites=nsites, chap=as.vector(covs[,"chaparral"])[dat$cell],
                chap2=as.vector(covs[,"chaparral"]^2)[dat$cell],
                elev=as.vector(covs[,"elevation"])[dat$cell], T=nyrs, nD=nD, midpt=midpt,
                B=B, delta=delta, y=y, dclass=dclass, site=site, nind=sum(y)) )



# 9.7.2 Fitting a slurry of open population models
# ------------------------------------------------------------------------


# 9.7.2.1 The independence model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    beta4 ~ dunif(-20,20) # Population trend parameter
    r ~ dunif(0,5)        # NegBin dispersion parameter
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
    pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    for (t in 1:T){
    log(lambda[s,t]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s] + beta4*(t - t/2)  # Note trend parameter here
    y[s,t] ~ dbin(pcap[s], N[s,t])
    N[s,t] ~ dnegbin(prob[s,t], r)
    prob[s,t] <- r/(r+lambda[s,t])
    } # End loop over years
    } # End loop over sites
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
    }
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)   # 300 m point = 28.27 ha
    }
    }
    ", file="Sollmann1.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1  # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1),
                         beta3=runif(1), beta4=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), r = 1)}
params <-c('beta0', 'beta1', 'beta2', 'beta3', 'beta4', 'alpha0', 'alpha1', 'Ntot', 'D', 'r')
ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# JAGS setting b/c otherwise JAGS cannot build a sampler, rec. by M. Plummer
set.factory("bugs::Conjugate", FALSE, type="sampler")

# Execute JAGS, look at convergence and summarize the results
library(jagsUI)
open1 <- jags (data1, inits, params, "Sollmann1.txt", n.thin=nt, n.chains=nc,
               n.burnin=nb, n.iter=ni)
par(mfrow = c(3,3))   ;   traceplot(open1)   ;   print(open1, 2)




# 9.7.2.2 The reduced-dynamics model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    theta ~ dunif(0,5)
    # NegBin dispersion parameter
    r ~ dunif(0,5)
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    # Compute detection probability
    for(k in 1:nD){
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    # Abundance model for Yr1 as in Sillett et al 2012
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])
    
    # Population dynamics model for subsequent years
    for (t in 2:T){
    N[s,t] ~ dpois(N[s, t-1] * theta)
    y[s,t] ~ dbin(pcap[s], N[s,t])
    }
    }
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
    }
    
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)  # 300 m point = 28.27 ha
    }
    }
    ", file="Sollmann2.txt")

# Set up initial values, parameters vector and MCMC settings
Nst <- y+1 # this is for trend model
inits <- function(){list(N=Nst, beta0=runif(1), beta1=runif(1), beta2=runif(1),
                         beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), theta=runif(1,0.6,0.99))}
params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'theta',
            'rout', 'Ntot', 'D', 'r')
ni <- 22000   ;   nb <- 2000   ;   nt <- 1   ;   nc <- 3

# Execute JAGS, look at convergence and summarize the results
open2 <- jags (data1, inits, params, "Sollmann2.txt", n.thin=nt, n.chains=nc,
               n.burnin=nb, n.iter=ni)
par(mfrow = c(3,3))   ;   traceplot(open2)   ;   print(open2, 2)



# 9.7.2.3 The Glorious Integrated HDS/Dail-Madsen Model
# ------------------------------------------------------------------------
# Write out the BUGS model file
cat("
    model{
    # Prior distributions
    # Regression parameters
    alpha0 ~ dunif(0,20)
    alpha1 ~ dunif(-10,10)
    beta0 ~ dunif(-20,20)
    beta1 ~ dunif(-20,20)
    beta2 ~ dunif(-20,20)
    beta3 ~ dunif(-20,20)
    
    # Priors for dynamics parameters: here they are constant across years
    # We could add covariate models for logit(phi) and log(gamma)
    phi ~ dunif(0,1)
    gamma ~ dunif(0,5)
    
    # NegBin dispersion parameter
    r ~ dunif(0,5)
    rout <- log(r)
    
    # 'Likelihood'
    for (s in 1:nsites){
    # Linear model for detection function scale
    log(sigma[s]) <- alpha0+alpha1*chap[s]
    
    # Compute detection probability
    for(k in 1:nD){
    log(p[k,s]) <- -midpt[k]*midpt[k]/(2*sigma[s]*sigma[s])
    f[k,s] <- p[k,s]*pi[k,s]
    fc[k,s] <- f[k,s]/pcap[s]
    fct[k,s] <- fc[k,s]/sum(fc[1:nD,s])
    pi[k,s] <- (2*midpt[k]*delta )/(B*B)
    }
    pcap[s]<-sum(f[1:nD,s])  # Overall detection probability
    
    # Process model
    # Abundance model for year 1
    log(lambda[s,1]) <- beta0 + beta1*chap[s] + beta2*chap2[s] + beta3*elev[s]
    y[s,1] ~ dbin(pcap[s], N[s,1])
    N[s,1] ~ dnegbin(prob[s,1], r)
    prob[s,1] <- r/(r+lambda[s,1])
    
    # Population dynamics model for subsequent years
    for (t in 2:T){                      # Loop over years
    S[s,t] ~ dbinom(phi, N[s, t-1])   # Survivors
    R[s,t] ~ dpois(gamma * N[s, t-1]) # Recruits
    N[s,t] <- S[s,t] + R[s,t]         # N = Survivors + Recruits
    y[s,t]~ dbin(pcap[s],N[s,t])       # Measurement error
    }
    }
    
    # Distance sampling observation model for observed (binned) distance data
    for(i in 1:nind){
    dclass[i] ~ dcat(fct[1:nD,site[i]])
    }
    
    # Derived parameters
    for(t in 1:6){
    Ntot[t] <- sum(N[,t])
    D[t] <- Ntot[t] / (28.27*nsites)     # 300 m point = 28.27 ha
    }
    }
    ", file="Sollmann3.txt")


# Set up some sensible starting values for S and R
yin <- y+1
yin[,2:6] <- NA
Sin <- Rin <- matrix(NA, nrow=nsites, ncol=nyrs)
y1 <- y + 1
for(s in 1:nsites){
  for (t in 2:6){
    Sin[s,t] <- rbinom(1,y1[s,t-1], phi )
    Rin[s,t] <- ifelse((y1[s,t]-Sin[s,t])>0, y1[s,t]-Sin[s,t], 0)
  }
}

# Set up initial values, parameters vector and MCMC settings
inits <-function(){list(N=yin, beta0=runif(1), beta1=runif(1), beta2=runif(1),
                        beta3=runif(1), alpha0=runif(1,3,5), alpha1=runif(1), phi=0.6, gamma=0.3,
                        R=Rin, S=Sin) }
params <- c('beta0', 'beta1', 'beta2', 'beta3', 'alpha0', 'alpha1', 'phi',
            'gamma', 'Ntot', 'D', 'r')
ni <- 152000   ;   nb <- 2000   ;   nt <- 10   ;   nc <- 5

# Run JAGS, look at convergence and summarize the results
library(jagsUI)
open3  <- jags (data1, inits, params, "Sollmann3.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, parallel=TRUE)
par(mfrow = c(3,3))   ;   traceplot(open3)   ;   print(open3, 2)


# Compare inferences in graph .... (Fig. 9-6)
plot(apply(dat$N,2,sum),ylim=c(600,1300),xlab="Year",ylab="Population size (307 sample units)")
lines(apply(open1$sims.list$Ntot,2,mean), lty=1, col="blue", lwd=2)
lines(apply(open2$sims.list$Ntot,2,mean), lty=1, col="red", lwd=2)
lines(apply(open3$sims.list$Ntot,2,mean), lty=1, col="green", lwd=2)
open1.95cri<- apply(open1$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open2.95cri<- apply(open2$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
open3.95cri<- apply(open3$sims.list$N,2,function(x) quantile(x, c(0.025,0.975)))
legend(1,750,legend=c("Independence","Reduced dynamics","Full dynamics"), lty=1, col=c("blue","red","green"))

matlines(1:6, t(open1.95cri), type="l", lty=2, lwd=2, col="blue")
matlines(1:6, t(open2.95cri), type="l", lty=2, lwd=2, col="red")
matlines(1:6, t(open3.95cri), type="l", lty=2, lwd=2, col="green")

# .... and table
parms <- c(betaFall, dparm)
round(post <- cbind(parms, Independent=open1$summary[c(1:4,6,7),1],
                    Partial=open2$summary[1:6,1], Full=open3$summary[1:6,1]), 3)

parms Independent Partial   Full
lambda(Int)             0.827       0.990   0.971  0.949
lambda(chaparral)       1.432       1.755   1.783  1.806
lambda(elevation)      -0.227      -0.509  -0.519 -0.523
lambda(I(chaparral^2)) -0.376      -0.343  -0.317 -0.323
p(Int)                  4.679       4.682   4.672  4.680
p(chaparral)           -0.199      -0.194  -0.187 -0.192



# 9.7.2.4 Summary remarks on modelling populations over time
# ------------------------------------------------------------------------



# 9.8 Spatial Distance Sampling: Modelling within-unit variation in density
# ------------------------------------------------------------------------


# 9.8.1 Distance sampling with location of encounter
# ------------------------------------------------------------------------
# Simulate a data set and harvest the output
set.seed(1234)
str(tmp <- sim.pdata(N=200,sigma=1,keep.all=FALSE,B=3))

List of 8
$ N     : num 200
$ sigma : num 1
$ B     : num 3
$ u1    : num [1:37] 3.73 3.66 3.27 1.72 3.32 ...
$ u2    : num [1:37] 3.17 1.9 1.68 3.83 2.54 ...
$ d     : num [1:37] 0.753 1.276 1.345 1.529 0.557 ...
$ y     : int [1:200] 0 1 1 0 0 0 0 0 0 0 ...
$ N.real: int 152

# Harvest some data objects
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
nind <- length(d)

# Data augmentation
M <- 400                       # max of 400 individuals
nz <- M-nind                   # augment by nz individuals
y <- c(rep(1,nind), rep(0,nz)) # augmented data augmentation variable
u <- cbind(u1,u2)               # Locations of individuals detected
u <- rbind(u, matrix(NA, nrow=nz, ncol=2))


# Bundle and summarize the data set
str(data <- list (B=B, nind=nind, u=u, y=y, nz=nz))


# Write out the BUGS model file
cat("
    model{
    
    # Priors
    sigma ~ dunif(0,10)
    psi ~ dunif(0,1)
    
    # Categorical observation model
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    u[i,1] ~ dunif(0, 2*B)  # Here is the uniformity assumption made explicit
    u[i,2] ~ dunif(0, 2*B)
    # Compute distance as a derived quantity
    d[i] <- pow( pow( u[i,1]-B,2) + pow(u[i,2]-B,2), 0.5) # Pythagoras
    p[i] <- exp(-d[i]*d[i] / (2*sigma*sigma))
    mu[i] <- p[i] * z[i]
    y[i] ~ dbern(mu[i])
    }
    # Other derived quantity
    N <- sum(z[])
    D <- N / (B*B)
    }
    ",fill=TRUE,file="model1.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")   ;   library(jagsUI)
ni <- 22000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Inits and parameters
inits <- function(){
  list(sigma=runif(1,1,10), psi=runif(1),z = c(rep(1,nind),rep(0,nz)) ) }
params <- c("sigma", "N", "psi")

# Execute jags and summarize the posterior distributions
out1 <- jags (data, inits, parameters, "model1.txt", n.thin=nthin,
              n.chains=nc, n.burnin=nb,n.iter=ni, parallel = FALSE)
par(mfrow = c(2,2))   ;   traceplot(out1)
print(out1, 2)

# 9.8.2 The line transect case
# ------------------------------------------------------------------------


# 9.8.3 Modelling spatial covariates
# ------------------------------------------------------------------------
# Simulator function for spatial distance sampling data
sim.spatialDS <-
  function(N=1000, beta = 1, sigma=1, keep.all=FALSE, B=B, model="halfnorm"){
    # Function simulates coordinates of individuals on a square
    # Square is [0,2B] x [0,2B], with a count location on the point (B, B)
    #   N: total population size in the square
    #   beta: coefficient of SOEMTHING on spatial covariate x
    #   sigma: scale of half-normal detection function
    #   B: circle radius
    #   keep.all: return the data for y=0 individuals or not
    library(raster)      # Load required packages
    library(plotrix)
    
    # Create coordinates for 30 x 30 grid
    delta <- (2*B-0)/30                # '2D bin width'
    grx <- seq(delta/2, 2*B - delta/2, delta) # mid-point coordinates
    gr <- expand.grid(grx,grx)         # Create grid coordinates
    
    # Create spatially correlated covariate x and plot it
    V <- exp(-e2dist(gr,gr)/1)
    x <- t(chol(V))%*%rnorm(900)
    par(mar=c(3,3,3,6))
    image(rasterFromXYZ(cbind(gr,x)), col=topo.colors(10))
    draw.circle(3, 3, B)
    points(3, 3, pch="+", cex=3)
    image_scale(x, col=topo.colors(10))  ## function renamed
    
    # Simulate point locations as function of habitat covariate x
    probs <- exp(beta*x)/sum(exp(beta*x)) # probability of point in pixel (sum = 1)
    pixel.id <- sample(1:900, N, replace=TRUE, prob=probs)
    # could simulate randomly within the pixel but it won't matter so place centrally
    u1 <- gr[pixel.id,1]
    u2 <- gr[pixel.id,2]
    points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
    title("This is so cool !")         # express your appreciation of all this
    
    d <- sqrt((u1 - B)^2 + (u2-B)^2)   # distance to center point of square
    #plot(u1, u2, pch = 1, main = "Point transect")
    N.real <- sum(d<= B)               # Population size inside of count circle
    
    # Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
    # p <- ifelse(d< B, 1, 0) * exp(-d*d/(2*(sigma^2)))
    # We do away with the circle constraint here.
    if(model=="hazard")
      p <- 1-exp(-exp(-d*d/(2*sigma*sigma)))
    if(model=="halfnorm")
      p <- exp(-d*d/(2*sigma*sigma))
    # Now we decide whether each individual is detected or not
    y <- rbinom(N, 1, p)                                           # detected or not
    points(u1[d<= B], u2[d<= B], pch = 16, col = "black", cex = 1) # not detected
    points(u1[y==1], u2[y==1], pch = 16, col = "red", cex = 1)     # detected
    
    # Put all of the data in a matrix
    if(!keep.all){
      u1 <- u1[y==1]
      u2 <- u2[y==1]
      d <- d[y==1]   }
    # Output
    return(list(model=model, N=N, beta=beta, B=B, u1=u1, u2=u2, d=d, y=y, N.real=N.real, Habitat=x, grid=gr))
  }

# Generate one data set and harvest the output
set.seed(1234)
str(tmp <- sim.spatialDS(N=200, beta=1, sigma=1.5, keep.all=FALSE, B=3)) # Fig. 9-7


# Harvest data
B <- tmp$B
d <- tmp$d
u1 <- tmp$u1
u2 <- tmp$u2
Habitat <- as.vector(tmp$Habitat)
Habitat <- Habitat - mean(Habitat)
Habgrid <- tmp$grid
nind <- length(d)
G <- nrow(Habgrid)

# Do data augmentation, including for pixel ID
M <- 400
nz <- M-nind
pixel <- rep(NA, M)   # We use discrete "pixel ID" here instead of "s"
y <- c(rep(1,nind), rep(0,nz))

# Pick some starting values and figure out the pixel of each observation
s <- cbind(u1,u2)
s <- rbind(s, matrix(NA,nrow=nz,ncol=2))
D <- e2dist(s[1:nind,], Habgrid)
for(i in 1:nind){
  pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
}

# Bundle and summarize the data for BUGS
str(data <- list (B=B, nind=nind, y=y, nz=nz, Habitat=Habitat, Habgrid=Habgrid, G=G, pixel=pixel))


# Write BUGS model
cat("
    model{
    
    # Prior distributions
    sigma ~ dunif(0,10)
    psi ~ dunif(0,1)
    beta ~ dnorm(0,0.01)
    
    for(g in 1:G){   # g is the pixel index, there are G total pixels
    probs.num[g] <- exp(beta*Habitat[g])
    probs[g] <- probs.num[g]/sum(probs.num[])
    }
    
    # Models for DA variables and location (pixel)
    for(i in 1:(nind+nz)){
    z[i] ~ dbern(psi)
    pixel[i] ~ dcat(probs[])
    s[i,1:2] <- Habgrid[pixel[i],]   # location = derived quantity
    # compute distance = derived quantity
    d[i] <- pow(   pow( s[i,1]-B,2) + pow(s[i,2]-B,2), 0.5)
    p[i] <- exp(-d[i]*d[i]/(2*sigma*sigma))  # Half-normal detetion function
    mu[i]<- p[i]*z[i]
    y[i] ~ dbern(mu[i])                      # Observation model
    }
    # Derived parameters
    N <- sum(z[])                      # N is a derived parameter
    D <- N/9                           # area = 9 ha
    }
    ",fill=TRUE, file="spatialDS.txt")


# Load libraries and specify MCMC settings
library("R2WinBUGS")
library(jagsUI)
ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3

# Create inits and define parameters to monitor
inits <- function(){  list (sigma=runif(1,1,10),psi=runif(1),
                            z = c(rep(1,nind), rep(0,nz)) ) }
params <- c("sigma", "N", "psi", "beta", "D")

# Run JAGS, check convergence and summarize posteriors
out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
              n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
par(mfrow = c(2,3)   ;   traceplot(out2)
    print(out2, 2)
    
    
    # Add pixel in order to make a density map
    params <- c("sigma", "N", "psi", "beta", "D", "pixel")
    
    # Run JAGS, check convergence and summarize posteriors
    out2 <- jags (data, inits, params, "spatialDS.txt", n.thin=nthin,
                  n.chains=nc, n.burnin=nb, n.iter=ni, parallel = FALSE)
    
    # Plot density maps
    library(raster)
    par(mfrow=c(1,2))
    pixel <- out2$sims.list$pixel
    post <- table(pixel)/nrow(pixel)   # Average number of locations in each pixel
    prior.mean <- mean(out2$sims.list$beta)*as.vector(Habitat)
    prior.mean <- mean(out2$sims.list$psi)*M*exp(prior.mean)/sum(exp(prior.mean))
    plot(rast.data <- rasterFromXYZ(cbind(Habgrid,prior.mean)), axes=FALSE,
         col=topo.colors(10) )
    title("Prior mean density (estimated)")
    plot(rast.post <- rasterFromXYZ(cbind(Habgrid,as.vector(post))),axes=FALSE,
         col=topo.colors(10) )
    title("Posterior mean density")
    
    
    # 9.8.4 Spatial HDS models in unmarked using the pcount function
    # ------------------------------------------------------------------------
    ## see errata; here we use the "logit" detection function; the results will be different
    ##  but not 'nonsense'
    # Simulate a data set, N = 600 for the population size
    set.seed(1234)
    tmp <-sim.spatialDS(N=600, sigma=1.5, keep.all=FALSE, B=3, model= "logit")
    
    # Harvest stuff
    B <- tmp$B
    d <- tmp$d
    u1 <- tmp$u1
    u2 <- tmp$u2
    Habitat <- as.vector(tmp$Habitat)
    Habitat <- Habitat - mean(Habitat)
    Habgrid <- tmp$grid
    nind <- length(d)
    G <- nrow(Habgrid)
    
    # Find which pixel each observation belongs to
    s <- cbind(u1,u2)
    D <- e2dist(s[1:nind,], Habgrid)
    pixel <-rep(NA,nind)
    for(i in 1:nind){
      pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
    }
    
    # Create a vector of counts in each pixel and pad it with zeros
    pixel.count <- rep(0, G)
    names(pixel.count) <- 1:G
    pixel.count[names(table(pixel))] <- table(pixel)
    # Create a covariate: distance between observer and pixel center
    dist <- sqrt( (Habgrid[,1]-3)^2 + (Habgrid[,2]-3)^2  )
    # Construct an unmarkedFrame
    umf <- unmarkedFramePCount(y=matrix(pixel.count,ncol=1),
                               siteCovs=data.frame(dist=dist,Habitat=Habitat))
    summary(umf)
    
    
    # Fit an N-mixture model with no intercept and distance squared using
    #   the function pcount.spHDS now in package unmarked
    
    (fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20))
    
    lam <- exp( coef(fm1)[1] + coef(fm1)[2]*Habitat )
    pred <- predict(fm1, type='state')
    sum(lam)
    sum(pred[,1])  # Same
    
    
    # 9.8.5 Hierarchical spatial distance sampling
    # ------------------------------------------------------------------------
    sim.spatialHDS(lam0 = 4, sigma = 1.5, B = 3, nsites = 100)
    # lam0 = expected population size per site
    # nsites = number of point count locations
    # B = count radius. Function simulates coordinates of individuals on a square
    #       [0,2*B] x[0,2*B], with a count location on the point (B,B)
    # sigma = scale of half-normal detection function
    
    library(raster)
    
    # Simulate a data set and harvest the output
    set.seed(1234)
    str(tmp <-sim.spatialHDS(lam0 = 3, sigma = 1.5, B = 3, nsites = 100))
    
    # Process the simulated data set
    data <- tmp$data
    # To make it a ???real??? data set:
    data <- data[!is.na(data[,2]),] # Get rid of the 0 sites
    data <- data[data[,"y"]==1,]    # Only keep detected individuals
    
    # Now zero-pad the observed counts
    nsites <- tmp$nsites
    nobs <- rep(0, nsites)
    names(nobs) <- 1:nsites
    nobs[names(table(data[,1]))] <- table(data[,1])
    
    # Extract data elements that we need
    site <- data[,"site"]
    s <- data[,c("u1","u2")]
    B <- tmp$B
    Habitat <- (tmp$Habitat) # Raster values
    Habgrid <- tmp$grid   # Raster coordinates
    nind <- nrow(data)
    G <- nrow(Habgrid)
    
    # We have to convert observed locations to pixels
    pixel <- rep(NA,nrow(data))
    D <- e2dist(s[1:nind,], Habgrid)
    for(i in 1:nind){
      pixel[i] <- (1:ncol(D))[D[i,]==min(D[i,])]
    }
    
    # Do data augmentation of data from each site ???S-fold??? DA. Three objects need
    # to have DA applied to them: Ymat, umat, pixmat
    Msite <- 2*max(nobs)  # Perhaps use a larger value
    Ymat <- matrix(0,nrow=Msite,ncol=nsites)
    umat <- array(NA, dim=c(Msite,nsites,2))
    pixmat <- matrix(NA,nrow=Msite,ncol=nsites)
    for(i in 1:nsites){
      if(nobs[i]==0) next
      Ymat[1:nobs[i],i]<- data[data[,1]==i,"y"]
      umat[1:nobs[i],i,1:2]<- data[data[,1]==i,c("u1","u2")]
      pixmat[1:nobs[i],i]<- pixel[data[,1]==i]
    }
    
    # Bundle the data for BUGS
    str(data <- list (y=Ymat, pixel=pixmat, Habitat=Habitat, Habgrid=Habgrid, G = G,
                      nsites=nsites, M = Msite, B = B))
    
    # Write out the BUGS model file
    cat("
        model{
        
        # Prior distributions
        sigma ~ dunif(0,10)
        beta1 ~ dnorm(0,0.01)
        beta0 ~ dnorm(0,0.01)
        lam0 <- exp(beta0)*G           # Baseline lambda in terms of E(N) per sample unit
        
        # For each site, construct the DA parameter as a function of lambda
        for(s in 1:nsites){
        lamT[s] <- sum(lambda[,s])   # total abundance at a site
        psi[s] <- lamT[s]/M
        for(g in 1:G){             # g is the pixel index, there are G total pixels
        lambda[g,s] <- exp(beta0 + beta1*Habitat[g,s])
        probs[g,s] <- lambda[g,s]/sum(lambda[,s])
        }
        }
        
        # DA variables and spatial location variables:
        for(s in 1:nsites){
        for(i in 1:M){
        z[i,s] ~ dbern(psi[s])
        pixel[i,s] ~ dcat(probs[,s])
        u[i,s,1:2] <- Habgrid[pixel[i,s],]   # location = derived quantity
        # distance = derived quantity
        d[i,s] <- pow(   pow( u[i,s,1]-B,2) + pow(u[i,s,2]-B,2), 0.5)
        p[i,s] <- exp(-d[i,s]*d[i,s]/(2*sigma*sigma))  # Half-normal model
        mu[i,s] <- p[i,s]*z[i,s]
        y[i,s] ~ dbern(mu[i,s])    # Observation model
        }
        # Derived parameters
        N[s]<- sum(z[,s])            # Site specific abundance
        }
        Ntotal <- sum(N[])             # Total across all sites
        D <- Ntotal/(9*nsites)         # Density: point area = 9 ha
        }
        ",fill=TRUE, file="spatialHDS.txt")
    
    
    # Inits and parameters saved
    zst <- Ymat
    inits <- function(){ list (sigma=1.5, z = zst, beta0 = -5, beta1=1 ) }
    params <- c("sigma", "Ntotal", "beta1", "beta0", "D", "lam0")
    
    # MCMC settings
    ni <- 12000   ;   nb <- 2000   ;   nthin <- 2   ;   nc <- 3
    
    # Call JAGS, check convergence and summarize the results
    out3 <- jags(data, inits, params, "spatialHDS.txt", n.thin=nthin, n.chains=nc, n.burnin=nb, n.iter=ni, parallel = TRUE)
    par(mfrow = c(2,3))   ;   traceplot(out3)
    print(out3, 3)
    
    Ntrue <- tmp$N
    Nhat <- out3$summary[7:106,1]
    plot(Ntrue, Nhat, xlab="Local population size, N, for each site", ylab="Posterior mean N for each site", pch=20)
    abline(0, 1, lwd=2)
    
    
    
    # 9.9 Summary
    # ------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    # ====================================================================================
    #
    # 10. Modeling static occurrence and species distributions using site-occupancy models
    #
    # ====================================================================================
    
    
    
    
    
    
    # 10.1. Introduction to the modeling of occurrence, including species distributions
    # ---------------------------------------------------------------------------------
    
    
    alpha <- seq(-10, 10, by = 1)
    curve(plogis(-10 +  1 * x), -2, 2, lwd = 3, ylim = c(0,1), xlab = "Covariate X", ylab = "Occupancy prob.", frame = F, main = "", cex.lab = 1.5, cex.axis = 1.5)
    for(i in 2:21){  curve(plogis(alpha[i] +  1 * x), -2, 2, lwd = 3, add = T)  }
    (min <- (plogis(-10 +  1 * 2) - plogis(-10 +  1 * -2)) )
    (max <- (plogis(0 +  1 * 2) - plogis(0 +  1 * -2)) )
    
    
    
    # 10.2. Another exercise in hierarchical modeling: Derivation of the site-occupancy model
    # ----------------------------------------------------------------------------------------
    
    
    
    # 10.3 Simulation and analysis of the simplest possible site-occupancy model
    # --------------------------------------------------------------------------
    
    
    # Choose sample sizes and prepare observed data array y
    set.seed(24)                  # So we all get same data set
    M <- 100                      # Number of sites
    J <- 2                        # Number of presence/absence measurements
    y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data
    
    # Parameter values
    psi <- 0.8                    # Probability of occupancy or presence
    p <- 0.5                      # Probability of detection
    
    # Generate presence/absence data (the truth)
    z <- rbinom(n = M, size = 1, prob = psi)  # R has no Bernoulli
    
    # Generate detection/nondetection data (i.e. presence/absence measurements)
    for(j in 1:J){
      y[,j] <- rbinom(n = M, size = 1, prob = z*p)
    }
    
    # Look at data
    sum(z)                        # True number of occupied sites
    sum(apply(y, 1, max))         # Observed number of occupied sites
    head(cbind(z=z, y))           # Truth and measurements for first 6 sites
    
    library(unmarked)
    umf <- unmarkedFrameOccu(y = y)  # Create unmarked data frame
    summary(umf)                     # Summarize data frame
    (fm1 <- occu(~1 ~1, data = umf)) # Fit model
    
    backTransform(fm1, "state")      # Get estimates on probability scale
    backTransform(fm1, "det")
    
    
    # Bundle data and summarize data bundle
    str( win.data <- list(y = y, M = nrow(y), J = ncol(y)) )
    
    # Specify model in BUGS language
    sink("model.txt")
    cat("
        model {
        # Priors
        psi ~ dunif(0, 1)
        p ~ dunif(0, 1)
        # Likelihood
        for (i in 1:M) {    # Loop over sites
        z[i] ~ dbern(psi)         # State model
        for (j in 1:J) { # Loop over replicate surveys
        y[i,j] ~ dbern(z[i]*p)  # Observation model (only JAGS !)
        #        y[i,j] ~ dbern(mu[i])  # For WinBUGS define 'straw man'
        }
        #   mu[i] <- z[i]*p          # Only WinBUGS
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    zst <- apply(y, 1, max)       # Avoid data/model/inits conflict
    inits <- function(){list(z = zst)}
    
    # Parameters monitored
    params <- c("psi", "p")
    
    # MCMC settings
    ni <- 5000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3
    
    # Call JAGS and summarize posteriors
    library(jagsUI)
    fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
                n.thin = nt, n.iter = ni, n.burnin = nb)
    print(fm2, dig = 3)
    
    
    
    # 10.4 A slightly more complex site-occupancy model with covariates
    # ------------------------------------------------------------------------
    
    
    # Choose sample sizes and prepare obs. data array y
    set.seed(1)                   # So we all get same data set
    M <- 100                      # Number of sites
    J <- 3                        # Number of presence/absence measurements
    y <- matrix(NA, nrow = M, ncol = J) # to contain the obs. data
    
    # Create a covariate called vegHt
    vegHt <- sort(runif(M, -1, 1)) # sort for graphical convenience
    
    # Choose parameter values for occupancy model and compute occupancy
    beta0 <- 0                    # Logit-scale intercept
    beta1 <- 3                    # Logit-scale slope for vegHt
    psi <- plogis(beta0 + beta1 * vegHt) # Occupancy probability
    # plot(vegHt, psi, ylim = c(0,1), type = "l", lwd = 3) # Plot psi relationship
    
    # Now visit each site and observe presence/absence perfectly
    z <- rbinom(M, 1, psi)        # True presence/absence
    
    # Look at data so far
    table(z)
    
    # Plot the true system state
    par(mfrow = c(1, 3), mar = c(5,5,2,2), cex.axis = 1.5, cex.lab = 1.5)
    plot(vegHt, z, xlab="Vegetation height", ylab="True presence/absence (z)", frame = F, cex = 1.5)
    plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")
    
    # Create a covariate called wind
    wind <- array(runif(M * J, -1, 1), dim = c(M, J))
    
    # Choose parameter values for measurement error model and compute detectability
    alpha0 <- -2                        # Logit-scale intercept
    alpha1 <- -3                        # Logit-scale slope for wind
    p <- plogis(alpha0 + alpha1 * wind) # Detection probability
    # plot(p ~ wind, ylim = c(0,1))     # Look at relationship
    
    # Take J = 3 presence/absence measurements at each site
    for(j in 1:J) {
      y[,j] <- rbinom(M, z, p[,j])
    }
    sum(apply(y, 1, max))               # Number of sites with observed presences
    
    # Plot observed data and true effect of wind on detection probability
    plot(wind, y, xlab="Wind", ylab="Observed det./nondetection data (y)", frame = F, cex = 1.5)
    plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
    
    # Look at the data: occupancy, true presence/absence (z), and measurements (y)
    cbind(psi=round(psi,2), z=z, y1=y[,1], y2=y[,2], y3=y[,3])
    
    # Create factors
    time <- matrix(rep(as.character(1:J), M), ncol = J, byrow = TRUE)
    hab <- c(rep("A", 33), rep("B", 33), rep("C", 34))  # Must have M = 100
    
    
    # Load unmarked, format data and summarize
    library(unmarked)
    umf <- unmarkedFrameOccu(
      y = y,                                            # Pres/Abs measurements
      siteCovs = data.frame(vegHt = vegHt, hab = hab),  # site-specific covs.
      obsCovs = list(wind = wind, time = time))         # obs-specific covs.
    summary(umf)
    
    
    # Fit model and extract estimates
    # Detection covariates follow first tilde, then occupancy covariates
    summary(fm.occ1 <- occu(~wind ~vegHt, data=umf))
    
    # Predict occupancy and detection as function of covs (with 95% CIs)
    # Add truth from data simulation (below for full code to produce fig. 10-2)
    newdat <- data.frame(vegHt=seq(-1, 1, 0.01))
    pred.occ <- predict(fm.occ1, type="state", newdata=newdat)
    newdat <- data.frame(wind=seq(-1, 1, 0.1))
    pred.det <- predict(fm.occ1, type="det", newdata=newdat)
    
    # Predictions for specified values of vegHt, say 0.2 and 2.1
    newdat <- data.frame(vegHt=c(0.2, 2.1))
    predict(fm.occ1, type="state", newdata=newdat, append = T)
    
    # ... for values of wind of -1 to 1
    newdat <- data.frame(wind=seq(-1, 1, , 5))
    predict(fm.occ1, type="det", newdata=newdat, append = T)
    
    
    # Fit detection-naive GLM to observed occurrence and plot comparison
    summary(fm.glm <- glm(apply(y, 1, max) ~ vegHt, family=binomial))
    plot(vegHt, apply(y, 1, max), xlab="Vegetation height", ylab="Observed occurrence ('ever observed ?')", frame = F, cex = 1.5)
    plot(function(x) plogis(beta0 + beta1*x), -1, 1, add=T, lwd=3, col = "red")
    lines(vegHt, predict(fm.glm,,"response"), type = "l", lwd = 3)
    lines(vegHt, predict(fm.occ1, type="state")[,1], col = "blue", lwd = 3)
    legend(-1, 0.9, c("Truth", "'LR' with p", "LR without p"), col=c("red", "blue", "black"), lty = 1, lwd=3, cex = 1.2)
    
    
    ranef(fm.occ1)
    
    
    (psi1 <- predict(fm.occ1, type="state")[1,1])
    (p1 <- predict(fm.occ1, type="det")[c(1:3),1])
    
    (z1 <- (psi1 * prod(1-p1)) / ((1 - psi1) + psi1 * prod(1-p1)))
    
    # Define function for finite-sample number and proportion of occupied sites
    fs.fn <- function(fm){
      Nocc <- sum(ranef(fm)@post[,2,])
      psi.fs <- Nocc / nrow(fm@data@y)
      out <- c(Nocc = Nocc, psi.fs = psi.fs)
      return(out)
    }
    
    # Bootstrap the function
    fs.hat <- fs.fn(fm.occ1)           # Point estimate
    pb.fs <- parboot(fm.occ1, fs.fn, nsim=10000, report=2) # Takes a while (33 min)
    # system.time(pb.fs <- parboot(fm.occ1, fs.fn, nsim=100, report=10)) # quicker
    
    # Summarize bootstrap distributions
    summary(pb.fs@t.star)
    
    # Get 95% bootstrapped confidence intervals
    (tmp1 <- quantile(pb.fs@t.star[,1], prob = c(0.025, 0.975)))
    
    (tmp2 <- quantile(pb.fs@t.star[,2], prob = c(0.025, 0.975)))
    
    # Plot bootstrap distribution of number of occupied sites (Fig. 10-3 left)
    par(mfrow = c(1,2), mar = c(5,4,3,2))
    hist(pb.fs@t.star[,1], col = "grey", breaks = 80, xlim = c(20, 100), main = "", freq = F)
    abline(v = fs.hat[1], col = "blue", lwd = 3)    # add point estimate
    abline(v = tmp1, col = "grey", lwd = 3)         # add 95% CI
    abline(v = sum(apply(y, 1, max)), lty = 2, lwd = 3) # observed #occ sites
    abline(v = sum(z), col = "red", lwd = 3)        # true #occ sites
    
    
    # Fit model p(time+wind), psi(hab+vegHt)
    summary(fm2.occ <- occu(~time+wind-1 ~hab+vegHt-1, data=umf))
    
    # Predict occupancy for habitat factor levels at average covariate values
    newdat <- data.frame(vegHt=0, hab = c("A", "B", "C"))
    predict(fm2.occ, type="state", newdata = newdat, appendData = TRUE)
    
    # Predict detection for time factor levels at average covariate values
    newdat <- data.frame(wind=0, time = c("1", "2", "3"))
    predict(fm2.occ, type="det", newdata=newdat, appendData = TRUE)
    
    
    # Fit model p(time*wind), psi(hab*vegHt)
    summary(fm3.occ <- occu(~time*wind-1-wind ~hab*vegHt-1-vegHt, data=umf))
    
    # Do likelihood ratio test
    LRT(fm2.occ, fm3.occ)
    
    # Bundle and summarize data set
    str( win.data <- list(y = y, vegHt = vegHt, wind = wind, M = nrow(y), J = ncol(y), XvegHt = seq(-1, 1, length.out=100), Xwind = seq(-1, 1, length.out=100)) )
    
    
    # Specify model in BUGS language
    sink("model.txt")
    cat("
        model {
        
        # Priors
        mean.p ~ dunif(0, 1)         # Detection intercept on prob. scale
        alpha0 <- logit(mean.p)      # Detection intercept
        alpha1 ~ dunif(-20, 20)      # Detection slope on wind
        mean.psi ~ dunif(0, 1)       # Occupancy intercept on prob. scale
        beta0 <- logit(mean.psi)     # Occupancy intercept
        beta1 ~ dunif(-20, 20)       # Occupancy slope on vegHt
        
        # Likelihood
        for (i in 1:M) {
        # True state model for the partially observed true state
        z[i] ~ dbern(psi[i])      # True occupancy z at site i
        logit(psi[i]) <- beta0 + beta1 * vegHt[i]
        for (j in 1:J) {
        # Observation model for the actual observations
        y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j
        p.eff[i,j] <- z[i] * p[i,j]   # 'straw man' for WinBUGS
        logit(p[i,j]) <- alpha0 + alpha1 * wind[i,j]
        }
        }
        
        # Derived quantities
        N.occ <- sum(z[])       # Number of occupied sites among sample of M
        psi.fs <- N.occ/M       # Proportion of occupied sites among sample of M
        for(k in 1:100){
        logit(psi.pred[k]) <- beta0 + beta1 * XvegHt[k] # psi predictions
        logit(p.pred[k]) <- alpha0 + alpha1 * Xwind[k]  # p predictions
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values: must give for same quantities as priors given !
    zst <- apply(y, 1, max)        # Avoid data/model/inits conflict
    inits <- function(){list(z = zst, mean.p = runif(1), alpha1 = runif(1), mean.psi = runif(1), beta1 = runif(1))}
    
    # Parameters monitored
    params <- c("alpha0", "alpha1", "beta0", "beta1", "N.occ", "psi.fs", "psi.pred", "p.pred", "z") # Also estimate z = "conditional occ. prob."
    
    # MCMC settings
    ni <- 25000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART 2 min) and summarize posteriors
    out1B <- bugs(win.data, inits, params, "model.txt", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
    print(out1B, dig = 3)
    
    
    # Compare truth with MLEs and bayesian posterior inference in table ...
    truth <- c(alpha0, alpha1, beta0, beta1, sum(z), sum(z)/M)
    tmp <- summary(fm.occ1)
    MLEs <- rbind(tmp[[2]][1:2,1:2], tmp[[1]][1:2,1:2], sumZ = c(mean(pb.fs@t.star[,1]), sd(pb.fs@t.star[,1])), psi.fs = c(mean(pb.fs@t.star[,2]), sd(pb.fs@t.star[,2])))
    print(cbind(truth, MLEs, out1B$summary[1:6, 1:2]))
    
    # .... and in a graph (Fig. 10-4)
    par(mfrow = c(2, 2), mar = c(4, 5, 2, 2), las = 1, cex.lab = 1, cex = 1.2)
    plot(vegHt, z, xlab="Vegetation height", ylab="Occupancy prob. (MLE)", ylim = c(0, 1), frame = F)   # True presence/absence
    lines(seq(-1,1, 0.01), pred.occ[,1], col = "blue", lwd = 2)
    matlines(seq(-1,1, 0.01), pred.occ[,3:4], col = "grey", lty = 1)
    lines(vegHt, psi, lwd=3, col="red")   # True psi
    plot(wind, y, xlab="Wind", ylab="Detection prob. (MLE)", ylim = c(0,1), frame=F)
    lines(seq(-1, 1, 0.1), pred.det[,1], col = "blue", lwd = 2)
    matlines(seq(-1, 1, 0.1), pred.det[,3:4], col = "grey", lty = 1)
    plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
    plot(vegHt, z, xlab="Vegetation height", ylab="Occupancy prob. (P. mean)", las = 1, frame = F)   # True presence/absence
    lines(vegHt, psi, lwd=3, col="red")   # True psi
    lines(win.data$XvegHt, out1B$summary[7:106,1], col="blue", lwd = 2)
    matlines(win.data$XvegHt, out1B$summary[7:106,c(3,7)], col="grey", lty = 1)
    plot(wind, y, xlab="Wind", ylab="Detection prob. (P. mean)", frame = F)
    plot(function(x) plogis(alpha0 + alpha1*x), -1, 1, add=T, lwd=3, col = "red")
    lines(win.data$Xwind, out1B$summary[107:206,1], col="blue", lwd = 2)
    matlines(win.data$Xwind, out1B$summary[107:206,c(3,7)], col="grey", lty = 1)
    
    
    # Plot posterior distribution of number of occupied sites (see Fig. 10-3, right)
    hist(out1B$sims.list$N.occ, col = "grey", breaks = 60, xlim = c(20, 100),
         main = "", freq = F)
    abline(v = out1B$mean$N.occ, col = "blue", lwd = 3)   # add point estimate
    abline(v = out1B$summary[5,c(3,7)], col = "grey", lwd = 3)   # add 95% CRI
    abline(v = sum(apply(y, 1, max)), lty = 2, lwd = 3)   # observed #occ sites
    abline(v = sum(z), col = "red", lwd = 2)              # true #occ sites
    
    
    
    # 10.5 A general data-simulation function for static occupancy models: simOcc()
    # -----------------------------------------------------------------------------
    
    
    simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, time.effects = c(-1, 1), alpha1 = -1, alpha2 = -3, alpha3 = 0, sd.lp = 0.5, b = 2, show.plot = TRUE)
    
    simOcc()                  # Execute function with default arguments
    simOcc(show.plot = FALSE) #    same, without plots
    simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = -2, beta2 = 2, beta3 = 1, mean.detection = 0.3, time.effects = c(-1, 1), alpha1 = -1, alpha2 = -3, alpha3 = 0, sd.lp = 0.5, b = 2, show.plot = TRUE) # Explicit defaults
    
    # Create a 'fix' data set and look at what we created
    set.seed(24)
    data <- simOcc()          # Assign results to an object called 'data'
    str(data)
    
    
    # Simplest possible occupancy model, with constant occupancy and detection
    tmp <- simOcc(mean.occ=0.6, beta1=0, beta2=0, beta3=0, mean.det=0.3, time.effects=c(0, 0), alpha1=0, alpha2=0, alpha3=0, sd.lp=0, b=0)
    str(tmp)                 # give overview of results
    
    # psi = 1 (i.e., species occurs at every site)
    tmp <- simOcc(mean.occ=1)   ;   str(tmp)
    
    # p = 1 (i.e., species is always detected when it occurs)
    tmp <- simOcc(mean.det=1)   ;   str(tmp)
    
    Other potentially interesting settings include these:
      simOcc(J = 2)                 # Only 2 surveys
    simOcc(M = 1, J = 100)        # No spatial replicates, but 100 measurements
    simOcc(beta3 = 1)             # Including interaction elev-wind on p
    simOcc(mean.occ = 0.96)       # A really common species
    simOcc(mean.occ = 0.05)       # A really rare species
    simOcc(mean.det = 0.96)       # A really easy species
    simOcc(mean.det = 0.05)       # A really hard species
    simOcc(mean.det = 0)          # The dreaded invisible species
    simOcc(alpha1=-2, beta1=2)    # Opposing effects of elev on psi and p
    simOcc(J = 10, time.effects = c(-5, 5)) # Huge time effects on p
    simOcc(sd.lp = 10)            # Huge (random) site effects on p
    simOcc(J = 10, b = 0)         # No behavioural response in p
    simOcc(J = 10, b = 2)         # Trap happiness
    simOcc(J = 10, b = -2)        # Trap shyness
    
    
    
    # 10.6 A model with lots of covariates: use of R function model.matrix with BUGS
    # ------------------------------------------------------------------------------
    
    
    set.seed(148)
    data <- simOcc(time.effects = c(0,0), sd.lp = 0, b = 0)
    str(data)                         # Look at data object
    
    # Create habitat factor
    hab <- c(rep("A", 90), rep("B", 90), rep("C", 87)) # must have M = 267 sites
    
    # Load library, format data and summarize unmarked data frame
    library(unmarked)
    umf <- unmarkedFrameOccu(
      y = data$y,
      siteCovs = data.frame(elev = data$elev, forest = data$forest, hab = hab),
      obsCovs = list(wind = data$wind))
    summary(umf)
    
    summary(fm <- occu(~elev+wind ~elev*forest*hab, data=umf))
    
    # Bundle and summarize data set
    HAB <- as.numeric(as.factor(hab))  # Get numeric habitat factor
    str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y), elev = data$elev, forest = data$forest, wind = data$wind, HAB = HAB) )
    
    # Specify model in BUGS language
    sink("modelA.txt")
    cat("
        model {
        
        # Priors
        mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
        alpha0 <- logit(mean.p)       #   same on logit scale
        mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
        beta0 <- logit(mean.psi)      #   same on logit scale
        for(k in 1:2){                # 2 terms in detection model
        alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
        }
        for(k in 1:11){               # 11 terms in occupancy model
        beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
        }
        
        # Likelihood
        for (i in 1:M) {              # Loop over sites
        z[i] ~ dbern(psi[i])
        logit(psi[i]) <- beta0 +                 # occupancy (psi) intercept
        beta[1] * elev[i] +                    # effect of elev
        beta[2] * forest[i] +                  # effect of forest
        beta[3] * equals(HAB[i],2) +           # effect of habitat 2 (= B)
        beta[4] * equals(HAB[i],3) +           # effect of habitat 3 (= C)
        beta[5] * elev[i] * forest[i] +                     # elev:forest
        beta[6] * elev[i] * equals(HAB[i],2) +              # elev:habB
        beta[7] * elev[i] * equals(HAB[i],3) +              # elev:habC
        beta[8] * forest[i] * equals(HAB[i],2) +            # forest:habB
        beta[9] * forest[i] * equals(HAB[i],3) +            # forest:habC
        beta[10] * elev[i] * forest[i] * equals(HAB[i],2) + # elev:forest:habB
        beta[11] * elev[i] * forest[i] * equals(HAB[i],3)   # elev:forest:habC
        for (j in 1:J) {           # Loop over replicates
        y[i,j] ~ dbern(z[i] * p[i,j])        # WinBUGS would need 'straw man' !
        logit(p[i,j]) <- alpha0 +            # detection (p) intercept
        alpha[1] * elev[i] +              # effect of elevation on p
        alpha[2] * wind[i,j]              # effect of wind on p
        }
        }
        }
        ",fill = TRUE)
    sink()
    
    # Inits
    inits <- function(){list(z = apply(data$y, 1, max), mean.psi = runif(1), mean.p = runif(1), alpha = rnorm(2), beta = rnorm(11))}
    
    # Parameters monitored
    params <- c("alpha0", "alpha", "beta0", "beta")
    
    # MCMC settings
    ni <- 50000   ;   nt <- 10   ;   nb <- 10000   ;   nc <- 3
    
    # Run JAGS (ART 4 min), look at convergence and summarize posteriors
    outA <- jags(win.data, inits, params, "modelA.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    traceplot(outA)    ;    print(outA, 3)
    
    # Compare MLEs and SEs with posterior means and sd's
    tmp <- summary(fm)
    cbind(rbind(tmp$state[1:2], tmp$det[1:2]), Post.mean = outA$summary[c(4:15, 1:3), 1], Post.sd = outA$summary[c(4:15, 1:3), 2])
    
    
    # Create design matrix for occupancy covariates and look at it
    occDM <- model.matrix(~ data$elev * data$forest * hab)[,-1] # Drop first col.
    head(occDM)             # Look at design matrix
    str(occDM)
    
    # Bundle and summarize data set
    str( win.data <- list(y = data$y, M = nrow(data$y), J = ncol(data$y), elev = data$elev, wind = data$wind, occDM = occDM) )
    
    # Specify model in BUGS language
    sink("modelB.txt")
    cat("
        model {
        
        # Priors
        mean.p ~ dunif(0, 1)          # Detection intercept on prob. scale
        alpha0 <- logit(mean.p)       #   same on logit scale
        mean.psi ~ dunif(0, 1)        # Occupancy intercept on prob. scale
        beta0 <- logit(mean.psi)      #   same on logit scale
        for(k in 1:2){                # 2 terms in detection model
        alpha[k] ~ dnorm(0, 0.1)   # Covariates on logit(detection)
        }
        for(k in 1:11){               # 11 terms in occupancy model
        beta[k] ~ dnorm(0, 0.1)    # Covariates on logit(occupancy)
        }
        
        # Likelihood
        for (i in 1:M) {
        z[i] ~ dbern(psi[i])
        logit(psi[i]) <- beta0 + inprod(beta[], occDM[i,])  # slick !
        for (j in 1:J) {
        y[i,j] ~ dbern(z[i] * p[i,j]) # In WinBUGS need 'straw man'
        logit(p[i,j]) <- alpha0 +     # detection (p) intercept
        alpha[1] * elev[i] +       # effect of elevation on p
        alpha[2] * wind[i,j]       # effect of wind on p
        }
        }
        }
        ",fill = TRUE)
    sink()
    
    # Call JAGS from R (ART 3.3 min) and summarize posteriors
    outB <- jags(win.data, inits, params, "modelB.txt", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    traceplot(outB)    ;    print(outB, 3)
    
    
    
    # 10.7 Study design, and bias and precision of site-occupancy estimators
    # ------------------------------------------------------------------------
    
    
    # Do simulation with 1000 reps
    simreps <- 1000
    library(unmarked)
    
    # Define arrays to hold the results
    p <- array(dim = c(simreps,3,3))
    estimates <- array(dim = c(2,simreps,3,3))
    
    # Choose number and levels of simulation factors
    nsites <- c(20, 120, 250)     # number of sites
    nsurveys <- c(2, 5, 10)       # number of repeat surveys
    
    # Start simulation
    system.time(                  # time whole thing
      for(j in 1:3) {               # Loop j over site factor
        for(k in 1:3) {             # Loop k over survey factor
          for(i in 1:simreps){      # Loop i over simreps
            # Counter
            cat("** nsites", j, "nsurveys", k, "simrep", i, "***\n")
            # Generate a data set: pick p and use p in simOcc()
            det.prob <- runif(1, 0.01, 0.99)
            data <- simOcc(M = nsites[j], J = nsurveys[k], mean.occupancy = 0.5,
                           beta1 = 0, beta2 = 0, beta3 = 0, mean.detection = det.prob,
                           time.effects = c(0, 0), alpha1 = 0, alpha2 = 0, alpha3 = 0,
                           sd.lp = 0, b = 0, show.plot = F)
            # Fit model
            umf <- unmarkedFrameOccu(y = data$y)
            tmp <- occu(~1 ~1, umf, se = FALSE)   # Only get MLEs, not SEs
            # Save results (p and MLEs)
            p[i,j,k] <- data$mean.det
            estimates[,i,j,k] <- coef(tmp)
          }
        }
      }
    )
    
    # Plot results
    par(mfrow = c(3,3), mar = c(4,5,3,1), cex.main = 1.2)
    for(j in 1:3){
      for(k in 1:3){
        lab <- paste(nsites[j],"sites,", nsurveys[k],"surveys")
        plot(p[,j,k], plogis(estimates[1,,j,k]), xlab = "Detection prob.",
             ylab = "Occupancy prob.", main = lab, ylim = c(0,1))
        abline(h = 0.5, col = "red", lwd = 2)
        lines(smooth.spline(plogis(estimates[1,,j,k])~ p[,j,k], df = 5),
              col = "blue", lwd = 2)
      }
    }
    
    
    
    # 10.8 Goodness of fit
    # ------------------------------------------------------------------------
    
    
    simOcc(M = 267, J = 3, mean.occupancy = 0.6, beta1 = 0, beta2 = 1, beta3 = 0, mean.detection = 0.3, time.effects = c(0, 0), alpha1 = 1, alpha2 = -1, alpha3 = 0, sd.lp = 0, b = 0)
    
    
    
    # 10.9 Distribution modeling and mapping of Swiss red squirrels
    # ------------------------------------------------------------------------
    ## Code modified to use the file "SwissSquirrels.txt" included in the AHMbook package
    
    # Read in data set, select squirrels and harvest data
    fn <- file.path(find.package("AHMbook"), "extdata", "SwissSquirrels.txt")
    data <- read.table(fn, header = TRUE)
    str(data)
    y <- as.matrix(data[,7:9])     # Grab 2007 squirrel det/nondet data
    elev.orig <- data[,"ele"]      # Unstandardised, original values of covariates
    forest.orig <- data[,"forest"]
    time <- matrix(as.character(1:3), nrow=265, ncol = 3, byrow = T)
    date.orig <- as.matrix(data[,10:12])
    dur.orig <- as.matrix(data[,13:15])
    
    # Overview of covariates
    covs <- cbind(elev.orig, forest.orig, date.orig, dur.orig)
    par(mfrow = c(3,3))
    for(i in 1:8){
      hist(covs[,i], breaks = 50, col = "grey", main = colnames(covs)[i])
    }
    pairs(cbind(elev.orig, forest.orig, date.orig, dur.orig))
    
    # Standardize covariates and mean-impute date and duration
    # Compute means and standard deviations
    (means <- c(apply(cbind(elev.orig, forest.orig), 2, mean), date.orig = mean(c(date.orig), na.rm = TRUE), dur.orig=mean(c(dur.orig), na.rm = TRUE)))
    (sds <- c(apply(cbind(elev.orig, forest.orig), 2, sd), date.orig = sd(c(date.orig), na.rm = TRUE), dur.orig=sd(c(dur.orig), na.rm = TRUE)))
    
    # Scale covariates
    elev <- (elev.orig - means[1]) / sds[1]
    forest <- (forest.orig - means[2]) / sds[2]
    date <- (date.orig - means[3]) / sds[3]
    date[is.na(date)] <- 0
    dur <- (dur.orig - means[4]) / sds[4]
    dur[is.na(dur)] <- 0
    
    # Load unmarked, format data and summarize
    library(unmarked)
    umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(elev = elev, forest = forest), obsCovs = list(time = time, date = date, dur = dur))
    summary(umf)
    
    We want to identify a model that is useful for inference, specifically for prediction of squirrel distribution to the whole of Switzerland. We do some step-wise model selection first on the detection part, then on the occupancy part, while keeping the detection part as identified in the first step.
    
    # Fit a series of models for detection first and do model selection
    summary(fm1 <- occu(~1 ~1, data=umf))
    summary(fm2 <- occu(~date ~1, data=umf))
    summary(fm3 <- occu(~date+I(date^2) ~1, data=umf))
    summary(fm4 <- occu(~date+I(date^2)+I(date^3) ~1, data=umf))
    summary(fm5 <- occu(~dur ~1, data=umf))
    summary(fm6 <- occu(~date+dur ~1, data=umf))
    summary(fm7 <- occu(~date+I(date^2)+dur ~1, data=umf))
    summary(fm8 <- occu(~date+I(date^2)+I(date^3)+dur ~1, data=umf))
    summary(fm9 <- occu(~dur+I(dur^2) ~1, data=umf))
    summary(fm10 <- occu(~date+dur+I(dur^2) ~1, data=umf))
    summary(fm11 <- occu(~date+I(date^2)+dur+I(dur^2) ~1, data=umf))
    summary(fm12 <- occu(~date+I(date^2)+I(date^3)+dur+I(dur^2) ~1, data=umf))
    
    # Put the fitted models in a "fitList" and rank them by AIC
    fms <- fitList("p(.)psi(.)"                     = fm1,
                   "p(date)psi(.)"                      = fm2,
                   "p(date+date2)psi(.)"                = fm3,
                   "p(date+date2+date3)psi(.)"          = fm4,
                   "p(dur)psi(.)"                       = fm5,
                   "p(date+dur)psi(.)"                  = fm6,
                   "p(date+date2+dur)psi(.)"            = fm7,
                   "p(date+date2+date3+dur)psi(.)"      = fm8,
                   "p(dur+dur2)psi(.)"                  = fm9,
                   "p(date+dur+dur2)psi(.)"             = fm10,
                   "p(date+date2+dur+dur2)psi(.)"       = fm11,
                   "p(date+date2+date3+dur+dur2)psi(.)" = fm12)
    
    (ms <- modSel(fms))
    
    # Continue with model fitting for occupancy, guided by AIC as we go
    # Check effects of elevation
    summary(fm13 <- occu(~date+dur+I(dur^2) ~elev, data=umf))
    summary(fm14 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2), data=umf))
    summary(fm15 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+ I(elev^3), data=umf))
    cbind(fm13@AIC, fm14@AIC, fm15@AIC) # model 14 with elev2 best
    
    # Check effects of forest and interactions
    summary(fm16 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest, data=umf))
    summary(fm17 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2), data=umf))
    summary(fm18 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest, data=umf))
    summary(fm19 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2), data=umf))
    summary(fm20 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
    summary(fm21 <- occu(~date+dur+I(dur^2) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest+ I(elev^2):I(forest^2), data=umf))
    cbind(fm16@AIC, fm17@AIC, fm18@AIC, fm19@AIC, fm20@AIC) # fm20 is best
    
    # Check for some additional effects in detection
    summary(fm22 <- occu(~date+dur+I(dur^2)+elev ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
    summary(fm23 <- occu(~dur+I(dur^2)+date*(elev+I(elev^2)) ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
    summary(fm24 <- occu(~dur+I(dur^2)+date*(elev+I(elev^2))+forest ~elev+I(elev^2)+forest+I(forest^2)+elev:forest+elev:I(forest^2)+I(elev^2):forest, data=umf))
    cbind(fm22@AIC, fm23@AIC, fm24@AIC) # None better, hence, stay with model 20
    
    
    library(AICcmodavg)
    system.time(gof.boot <- mb.gof.test(fm20, nsim = 1000))
    gof.boot
    
    
    # Create new covariates for prediction ('prediction covs')
    orig.elev <- seq(200, 2500,,100)    # New covs for prediction
    orig.forest <- seq(0, 100,,100)
    orig.date <- seq(15, 110,,100)
    orig.duration <- seq(100, 550,,100)
    ep <- (orig.elev - means[1]) / sds[1] # Standardise them like actual covs
    fp <- (orig.forest - means[2]) / sds[2]
    dp <- (orig.date - means[3]) / sds[3]
    durp <- (orig.duration - means[4]) / sds[4]
    
    # Obtain predictions
    newData <- data.frame(elev=ep, forest=0)
    pred.occ.elev <- predict(fm20, type="state", newdata=newData, appendData=TRUE)
    newData <- data.frame(elev=0, forest=fp)
    pred.occ.forest <- predict(fm20, type="state", newdata=newData, appendData=TRUE)
    newData <- data.frame(date=dp, dur=0)
    pred.det.date <- predict(fm20, type="det", newdata=newData, appendData=TRUE)
    newData <- data.frame(date=0, dur=durp)
    pred.det.dur <- predict(fm20, type="det", newdata=newData, appendData=TRUE)
    
    # Plot predictions against unstandardized 'prediction covs'
    par(mfrow = c(2,2), mar = c(5,5,2,3), cex.lab = 1.2)
    plot(pred.occ.elev[[1]] ~ orig.elev, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Elevation (m)", frame = F)
    matlines(orig.elev, pred.occ.elev[,3:4], lty = 1, lwd = 1, col = "grey")
    plot(pred.occ.forest[[1]] ~ orig.forest, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. occupancy prob.", xlab = "Forest cover (%)", frame = F)
    matlines(orig.forest, pred.occ.forest[,3:4], lty = 1, lwd = 1, col = "grey")
    plot(pred.det.date[[1]] ~ orig.date, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Date (1 = 1 April)", frame = F)
    matlines(orig.date, pred.det.date[,3:4], lty = 1, lwd = 1, col = "grey")
    plot(pred.det.dur[[1]] ~ orig.duration, type = "l", lwd = 3, col = "blue", ylim = c(0,1), las = 1, ylab = "Pred. detection prob.", xlab = "Survey duration (min)", frame = F)
    matlines(orig.duration, pred.det.dur[,3:4], lty = 1, lwd = 1, col = "grey")
    
    
    # Predict abundance and detection jointly along two separate covariate gradients
    # abundance ~ (forest, elevation) and detection ~ (survey duration, date)
    pred.matrix1 <- pred.matrix2 <- array(NA, dim = c(100, 100)) # Define arrays
    for(i in 1:100){
      for(j in 1:100){
        newData1 <- data.frame(elev=ep[i], forest=fp[j])       # For abundance
        pred <- predict(fm20, type="state", newdata=newData1)
        pred.matrix1[i, j] <- pred$Predicted
        newData2 <- data.frame(dur=durp[i], date=dp[j])        # For detection
        pred <- predict(fm20, type="det", newdata=newData2)
        pred.matrix2[i, j] <- pred$Predicted
      }
    }
    
    par(mfrow = c(1,2), cex.lab = 1.2)
    mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
    image(x=orig.elev, y=orig.forest, z=pred.matrix1, col = mapPalette(100), axes = FALSE, xlab = "Elevation [m]", ylab = "Forest cover [%]")
    contour(x=orig.elev, y=orig.forest, z=pred.matrix1, add = TRUE, lwd = 1.5, col = "blue", labcex = 1.3)
    axis(1, at = seq(min(orig.elev), max(orig.elev), by = 250))
    axis(2, at = seq(0, 100, by = 10))
    box()
    title(main = "Expected squirrel occurrence prob.", font.main = 1)
    points(data$ele, data$forest, pch="+", cex=1)
    
    image(x=orig.duration, y=orig.date, z=pred.matrix2, col = mapPalette(100), axes = FALSE, xlab = "Survey duration [min]", ylab = "Date (1 = April 1)")
    contour(x=orig.duration, y=orig.date, z=pred.matrix2, add = TRUE, lwd = 1.5, col = "blue", labcex = 1.3)
    axis(1, at = seq(min(orig.duration), max(orig.duration), by = 50))
    axis(2, at = seq(0, 100, by = 10))
    box()
    title(main = "Expected squirrel detection prob.", font.main = 1)
    matpoints(as.matrix(data[, 13:15]), as.matrix(data[, 10:12]), pch="+", cex=1)
    
    
    # Load the Swiss landscape data from unmarked
    data(Switzerland)             # Load Swiss landscape data in unmarked
    CH <- Switzerland
    
    # Get predictions of occupancy prob for each 1km2 quadrat of Switzerland
    newData <- data.frame(elev = (CH$elevation - means[1])/sds[1], forest = (CH$forest - means[2])/sds[2])
    predCH <- predict(fm20, type="state", newdata=newData)
    
    # Prepare Swiss coordinates and produce map
    library(raster)
    library(rgdal)
    
    # Define new data frame with coordinates and outcome to be plotted
    PARAM <- data.frame(x = CH$x, y = CH$y, z = predCH$Predicted)
    r1 <- rasterFromXYZ(PARAM)     # convert into raster object
    
    # Mask quadrats with elevation greater than 2250
    elev <- rasterFromXYZ(cbind(CH$x, CH$y, CH$elevation))
    elev[elev > 2250] <- NA
    r1 <- mask(r1, elev)
    
    # Plot species distribution map (Fig. 10-14 left)
    par(mfrow = c(1,2), mar = c(1,2,2,5))
    mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
    plot(r1, col = mapPalette(100), axes = F, box = F, main = "Red squirrel distribution in 2007")
    lakes <- readOGR(".", "lakes")
    rivers <- readOGR(".", "rivers")
    border <- readOGR(".", "border")
    plot(rivers, col = "dodgerblue", add = TRUE)
    plot(border, col = "transparent", lwd = 1.5, add = TRUE)
    plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
    
    # Plot SE of the species distrbution map (Fig. 10-14 right)
    r2 <- rasterFromXYZ(data.frame(x = CH$x, y = CH$y, z = predCH$SE))
    r2 <- mask(r2, elev)
    plot(r2, col = mapPalette(100), axes = F, box = F, main = "Uncertainty map 2007")
    plot(rivers, col = "dodgerblue", add = TRUE)
    plot(border, col = "transparent", lwd = 1.5, add = TRUE)
    plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
    points(data$coordx, data$coordy, pch = "+", cex = 0.8)
    
    
    # Get extent of squirrel occurrence in 2007
    sum(predCH$Predicted)                      # All quadrats
    sum(predCH$Predicted[CH$elevation < 2250]) # Only at elevations < 2250 m
    
    
    # Standardise prediction covariate identical to those in analysis
    pelev <- (CH$elevation - means[1]) / sds[1]
    pforest <- (CH$forest - means[2]) / sds[2]
    
    # Define function that predicts occupancy under model 20
    Eocc <- function(fm) {
      betavec <- coef(fm)[1:8]       # Extract coefficients in psi
      DM <- cbind(rep(1,length(pelev)), pelev, pelev^2, pforest, pforest^2, pelev*pforest, pelev*pforest^2, pelev^2*pforest) # design matrix
      pred <- plogis(DM%*%(betavec)) # Prediction = DM * param. vector
      Eocc <- sum(pred)              # Sum over all Swiss quadrats (no mask)
      Eocc
    }
    
    (estimate.of.occurrence <- Eocc(fm20))    # Same as before, without mask
    system.time(Eocc.boot <- parboot(fm20, Eocc, nsim=1000, report=10)) # 100 sec
    plot(Eocc.boot)         # Plot bootstrap distribution of extent of occurrence
    quantile(Eocc.boot@t.star, c(0.025, 0.975))
    
    
    
    # 10.10 Multi-scale occupancy models
    # ------------------------------------------------------------------------
    
    
    data <- sim3Occ()       # Execute function with default args
    str(data)               # Summary of output with some comments added
    
    # 'Null' model (model 1)
    str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(0, 0), beta.Xtheta = 0, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(0,0), beta.Xp = 0, sd.logit.p = 0))
    
    # No covariate effects, no random variability (model 2)
    str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 0, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 0, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = 0, sd.logit.p = 0))
    
    # All covariate effects, but no random variability (model 3)
    str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 0, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 0, mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 0))
    
    # Most complex model with all effects allowed for by sim function (model 4)
    str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 1, mean.theta = 0.6, theta.time.range = c(-1, 1), beta.Xtheta = 1, sd.logit.theta = 1, mean.p = 0.4, p.time.range = c(-2,2), beta.Xp = -1, sd.logit.p = 1))
    
    
    set.seed(1)
    str(data <- sim3Occ(nunit = 100, nsubunit = 5, nrep = 3, mean.psi = 0.8, beta.Xpsi = 1, mean.theta = 0.3, theta.time.range = c(-1, 1), beta.Xtheta = 1, mean.p = 0.2, p.time.range = c(-1,1), beta.Xp = -1))
    
    # Look at data
    str(data$z)           # True quadrat (pond) occurrence state
    str(data$a)           # True subquadrat (water sample) occurrence state
    str(data$y)           # Observed data
    cbind("pond"=data$z, "sample 1"= data$a[,1], "sample 2"= data$a[,2], "sample 3"= data$a[,3], "sample 4"= data$a[,4], "sample 5"= data$a[,5])
    which(data$z-apply(data$a, 1, max)==1) # Fungus present in pond, but not in examined samples
    
    
    # Bundle and summarize data set
    y <- data$y
    str( win.data <- list(y = y, n.pond = dim(y)[1], n.samples = dim(y)[2], n.pcr = dim(y)[3], covA = data$covA, covB = data$covB, covC = data$covC) )
    
    # Define model in BUGS language
    sink("model.txt")
    cat("
        model {
        
        # Priors and model for params
        int.psi ~ dunif(0,1)         # Intercept of occupancy probability
        for(t in 1:n.samples){
        int.theta[t] ~ dunif(0,1) # Intercepts availability probability
        }
        for(t in 1:n.pcr){
        int.p[t] ~ dunif(0,1)     # Intercepts detection probability (1-PCR error)
        }
        beta.lpsi ~ dnorm(0, 0.1)    # Slopes of three covariates
        beta.ltheta ~ dnorm(0, 0.1)
        beta.lp ~ dnorm(0, 0.1)
        
        # 'Likelihood' (or basic model structure)
        for (i in 1:n.pond){
        # Occurrence in pond i
        z[i] ~ dbern(psi[i])
        logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
        for (j in 1:n.samples){
        # Occurrence in sample j
        a[i,j] ~ dbern(mu.a[i,j])
        mu.a[i,j] <- z[i] * theta[i,j]
        logit(theta[i,j]) <- logit(int.theta[j]) + beta.ltheta * covB[i,j]
        for (k in 1:n.pcr){
        # PCR detection error process in sample k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j,k]
        logit(p[i,j,k]) <- logit(int.p[k]) + beta.lp * covC[i,j,k]
        }
        }
        tmp[i] <- step(sum(a[i,])-0.1)
        }
        
        # Derived quantities
        sum.z <- sum(z[])   # Total # of occupied ponds in sample
        sum.a <- sum(tmp[]) # Total # of ponds with presence in <=1 of the 5 samples
        } # end model
        ",fill=TRUE)
    sink()
    
    # Initial values
    zst <- apply(y, 1, max)        # inits for presence (z)
    ast <- apply(y, c(1,2), max)   # inits for availability (a)
    inits <- function() list(z = zst, a = ast, int.psi = 0.5, beta.lpsi = 0)
    
    # Parameters monitored
    params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta", "beta.lp", "sum.z", "sum.a")
    
    # MCMC setting
    ni <- 5000   ;   nt <- 2   ;   nb <- 1000   ;   nc <- 3
    
    # Call WinBUGS and summarize posterior
    out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.dir = bd) # bd="c:/WinBUGS14/"
    print(out, 3)
    data$sum.z
    
    
    
    # 10.11 Space-for-time substitution
    # ------------------------------------------------------------------------
    
    
    # 10.11.1 A magical covariate
    # ------------------------------------------------------------------------
    set.seed(1)
    data <- sim3Occ(nunit = 500, nsubunit = 5, nrep = 1, mean.psi = 0.8, beta.Xpsi = 1, sd.logit.psi = 0.4, mean.theta = 0.6, theta.time.range = c(0, 0), beta.Xtheta = 1, sd.logit.theta = 0.6, mean.p = 0.4, p.time.range = c(0,0), beta.Xp = -1, sd.logit.p = 0.8)
    
    
    # Bundle and summarize data set
    y <- data$y
    str( win.data <- list(y = y, nunit = dim(y)[1], nsubunit = dim(y)[2], nrep = dim(y)[3], covA = data$covA, covB = data$covB, covC = data$covC) )
    
    # Define model in BUGS langauge
    sink("model.txt")
    cat("
        model {
        
        # Priors
        int.psi ~ dunif(0,1)   # Occupancy probability
        int.theta ~ dunif(0,1) # Availability probability
        int.p ~ dunif(0,1)     # Detection probability
        beta.lpsi ~ dnorm(0, 0.01)
        beta.ltheta ~ dnorm(0, 0.01)
        beta.lp ~ dnorm(0, 0.01)
        
        # Likelihood
        for (i in 1:nunit){
        # Occupancy model for quad i
        z[i] ~ dbern(psi[i])
        logit(psi[i]) <- logit(int.psi) + beta.lpsi * covA[i]
        for (j in 1:nsubunit){
        # Availability in subquad j
        a[i,j] ~ dbern(mu.a[i,j])
        mu.a[i,j] <- z[i] * theta[i,j]
        logit(theta[i,j]) <- logit(int.theta) + beta.ltheta * covB[i,j]
        for (k in 1:nrep){
        # PCR detection error process in replicate k
        y[i,j,k] ~ dbern(mu.y[i,j,k])
        mu.y[i,j,k] <- a[i,j] * p[i,j]
        logit(p[i,j]) <- logit(int.p) + beta.lp * covC[i,j,1]
        }
        }
        tmp[i] <- step(sum(a[i,])-0.1)
        }
        
        # Derived quantities
        sum.z <- sum(z[])      # Total number of occupied quadrats
        sum.a <- sum(tmp[])    # Total number of quads with presence in samples
        p.theta <- int.p * int.theta # What a 2-level model estimates as 'p'
        }
        ",fill=TRUE)
    sink()
    
    # Initial values
    inits <- function() list(z = array(1, dim = data$nunit), a = array(1, dim =c(data$nunit, data$nsubunit)) )    # Set all to 1 to avoid conflict
    
    # Parameters monitored
    params <- c("int.psi", "int.theta", "int.p", "beta.lpsi", "beta.ltheta", "beta.lp","p.theta", "sum.z", "sum.a")
    
    # MCMC settings
    ni <- 25000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call JAGS (ART 15 min) and summarize posterior
    out <- jags(win.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)
    traceplot(out)   ;   print(out, 3)
    
    # Compare truth and estimate in table
    tmp <- cbind(rbind(mean.psi = data$mean.psi, mean.theta = data$mean.theta, mean.p = data$mean.p, beta.lpsi = data$beta.Xpsi, beta.ltheta = data$beta.Xtheta, beta.lp = data$beta.Xp, product.theta.p = data$mean.theta* data$mean.p, sum.z = data$sum.z), out$summary[c(1:8),c(1, 3, 7)])
    colnames(tmp) <- c("Truth", "Post.mean", "LCRL", "UCRL")
    print(tmp, 3)
    
    
    # 10.11.2 No magical covariate known: theta and p are confounded
    # ------------------------------------------------------------------------
    # Load unmarked and define a function to fit the model with unmarked
    library(unmarked)
    occUM.fn <- function(data = data, inits = c(1, 1, -1)){
      umf <- unmarkedFrameOccu(y = data$y[,,1], siteCovs =
                                 data.frame(covA = data$covA))
      tmp1 <- summary(fm <- occu(~1 ~covA, data=umf, starts=inits))
      tmp <- matrix(unlist(c(tmp1$state[,1], tmp1$det[,1], tmp1$state[,2],
                             tmp1$det[,2])), ncol = 2, byrow = F)
      dimnames(tmp) <- list(c("Occ_Int", "Occ_A", "Det_Int"), c("MLE", "SE"))
      return(MLE = tmp)
    }
    
    # Choose number of simulations and create structures to hold results
    simreps <- 10000                # takes about 30 min
    obs.stats <- array(NA, dim = c(simreps, 3))
    dimnames(obs.stats) <- list(NULL, c("sum.z", "obs.sum.z", "sum.z.x"))
    MLEs <- array(NA, dim = c(3,2,simreps))
    dimnames(MLEs) <- list(c("Occ_Int", "Occ_A", "Det_Int"), c("MLE", "SE"), NULL)
    
    # Set timer and launch simulation
    system.time(
      for(i in 1:simreps){
        cat("\n\n*** Simrep Number:", i, "***\n\n")
        # Generate data set
        data <- sim3Occ(nunit = 500, nsubunit = 5, nrep = 1, mean.psi = 0.8,
                        beta.Xpsi = 1, sd.logit.psi = 0.4, mean.theta = 0.6,
                        theta.time.range = c(-1, 1), beta.Xtheta = 0, sd.logit.theta = 0.6,
                        mean.p = 0.4, p.time.range = c(0,0), beta.Xp = 0, sd.logit.p = 0.8)
        # Save stats
        obs.stats[i,] <- unlist(data[23:25])
        # Get MLEs of occupancy model and save them
        UMmle <- try(occUM.fn(data = data, inits = c(1,1,-1)))
        if (class(UMmle) == "try-error") {v<-1} else {
          MLEs[,,i] <- UMmle
        }
        rm(data, UMmle)
      }
    )
    
    # Visualize results
    par(mfrow = c(1,3), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
    # Estimate of occupancy (psi)
    hist(plogis(MLEs[1,1,]), breaks = 40, col = "grey", main = "Quadrat occupancy (psi)")
    abline(v = mean(plogis(MLEs[1,1,]), na.rm = T), col = "blue", lwd = 3)
    abline(v = 0.8, col = "red", lwd = 3)
    # Estimate of occupancy covariate (A)
    hist(MLEs[2,1,], breaks = 40, col = "grey", main = "Quadrat occupancy covariate (covA)")
    abline(v = mean(MLEs[2,1,], na.rm = T), col = "blue", lwd = 3)
    abline(v = data$beta.Xpsi, col = "red", lwd = 3)
    # Estimate of "detection": product of theta and p
    hist(plogis(MLEs[3,1,]), breaks = 40, col = "grey", main = "'Detection probability' = \n theta * p")
    abline(v = mean(plogis(MLEs[3,1,]), na.rm = T), col = "blue", lwd = 3)
    abline(v = data$mean.theta * data$mean.p, col = "red", lwd = 3, lty = 2)
    
    
    
    # 10.12 Models for data along transects: Poisson, exponential,
    #       Weibull and removal observation models
    # ------------------------------------------------------------------------
    
    
    # 10.12.1 Occupancy models with "survival model" observation process:
    #         exponential time-to-detection (TTD) model with simulated data
    # ------------------------------------------------------------------------
    simOccttd(M = 250, mean.psi = 0.4, mean.lambda = 0.3, beta1 = 1, alpha1 = -1, Tmax = 10)
    
    set.seed(1)
    data <- simOccttd()
    str(data)
    
    # Plot response (not shown)
    hist(data$ttd, breaks = 50, col = "grey", main = "Observed distribution of time to detection", xlim = c(0, data$Tmax), xlab = "Measured time to detection")
    abline(v = data$Tmax, col = "grey", lwd = 3)
    
    # Bundle data
    str( win.data <- list(ttd = data$ttd, d = data$d, covA = data$covA,
                          covB = data$covB, nobs = data$M, Tmax = data$Tmax) )
    
    # Define exponential observation model
    cat(file = "model1.txt", "
        model {
        
        # Priors
        int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
        beta1 ~ dnorm(0, 0.001)             # Slope coefficient in logit(occupancy)
        int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter
        alpha1 ~ dnorm(0, 0.001)            # Slope coefficient in log(rate)
        
        # Likelihood
        for (i in 1:nobs){
        # Model for occurrence
        z[i] ~ dbern(psi[i])
        logit(psi[i]) <- logit(int.psi) + beta1 * covB[i]
        
        # Observation model
        # Exponential model for time to detection ignoring censoring
        ttd[i] ~ dexp(lambda[i])
        log(lambda[i]) <- log(int.lambda) + alpha1 * covA[i]
        # Model for censoring due to species absence and ttd>=Tmax
        d[i] ~ dbern(theta[i])
        theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
        }
        # Derived quantities
        n.occ <- sum(z[])                   # Number of occupied sites among M
        }
        ")
    
    # Inits function for some params
    # Initialize with z = 1 throughout and
    #   all missings due to censoring, rather than non-occurrence
    zst <- rep(1, length(win.data$ttd))
    ttdst <-rep(win.data$Tmax+1, data$M)
    ttdst[win.data$d == 0] <- NA
    inits <- function(){list(z =zst, ttd = ttdst, int.psi = runif(1), int.lambda = runif(1))}
    
    # Parameters to estimate
    params <- c("int.psi", "beta1", "int.lambda", "alpha1", "n.occ")
    
    # MCMC settings
    ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART 1.3 min) and summarize posteriors
    out1 <- bugs(win.data, inits, params, "model1.txt", n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt, debug = TRUE, bugs.directory = bd)
    print(out1, dig = 3)
    
    
    # 10.12.2 Time-to-detection analysis with real data:
    #         Weibull occupancy model for the peregrine spring survey
    # ------------------------------------------------------------------------
    # Code modified slightly to use the ttdPeregrine data set in the AHMbook package
    ?ttdPeregrine  # check the description of the data
    data(ttdPeregrine)
    data <- ttdPeregrine
    # Manage data and standardize time of day
    nobs <- length(data$SiteNumber)     # Number of observations
    d <- as.numeric(is.na(data$ttd))    # Censoring indicator
    mean.tod <- mean(data$MinOfDay)
    sd.tod <- sd(data$MinOfDay)
    tod <- (data$MinOfDay -mean.tod) / sd.tod
    
    # Bundle and summarize data set
    str( win.data <- list(M = max(data$SiteNumber), site = data$SiteNumber,
                          tod = tod, male = as.numeric(data$sex)-1, ttd = data$ttd, d = d, nobs = nobs,
                          Tmax = data$Tmax) )
    
    
    # Define model
    cat(file = "model2.txt", "
        model {
        
        # Priors
        psi ~ dunif(0, 1)              # Occupancy intercept
        lambda.int[1] ~ dgamma(0.001, 0.001) # Poisson rate parameter for females
        lambda.int[2] ~ dgamma(0.001, 0.001) # Poisson rate parameter for males
        alpha1 ~ dnorm(0, 0.001)       # Coefficient of time of day (linear)
        alpha2 ~ dnorm(0, 0.001)       # Coefficient of time of day (squared)
        shape ~ dgamma(0.001,0.001)    # Weibull shape
        sexratio ~ dunif(0,1)          # Sex ratio (proportion males)
        
        # Likelihood
        for (i in 1:M){                # Model for occurrence at site level
        z[i] ~ dbern(psi)
        }
        
        for (i in 1:nobs){             # Observation model at observation level
        # Weibull model for time to detection ignoring censoring
        ttd[i] ~ dweib(shape, lambda[i])
        log(lambda[i]) <- (1-male[i])*log(lambda.int[1]) + male[i]*log(lambda.int[2])+ alpha1 * tod[i] + alpha2 * pow(tod[i],2)
        # Model for censoring due to species absence and ttd>=Tmax
        d[i] ~ dbern(theta[i])
        theta[i] <- z[site[i]] * step(ttd[i] - Tmax[i]) + (1 - z[site[i]])
        # Model for sex of unobserved individuals
        male[i] ~ dbern(sexratio)   # Will impute sex for unobserved individuals
        }
        # Derived quantities
        n.occ <- sum(z[])              # Number of occupied sites among M
        }
        ")
    
    # Inits function
    zst <- rep(1, win.data$M)
    ttdst <-rep(win.data$Tmax+1)
    ttdst[win.data$d == 0] <- NA
    inits <- function(){list(z =zst, ttd = ttdst, psi = runif(1), lambda.int = runif(2), alpha1 = rnorm(1), alpha2 = rnorm(1), shape = runif(1))}
    
    # Parameters to estimate
    params <- c("psi", "lambda.int", "alpha1", "alpha2", "n.occ", "z", "sexratio", "shape")
    
    # MCMC settings
    ni <- 15000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART 0.6 min) and summarize posteriors
    out2 <- bugs(win.data, inits, params, "model2.txt", n.chains=nc, n.iter=ni, n.burn = nb, n.thin=nt, debug = T, bugs.directory = bd)
    print(out2, dig = 3)
    
    
    # Predict detection over time of day: prediction cov. runs from 7h to 19h
    minutes <- 60:780
    pred.tod <- (minutes - mean.tod) / sd.tod   # Standardize as real data
    
    # Predict p over time of day, averaging over sex, and for duration of 10 min
    sex.mean <- apply(out2$sims.list$lambda.int, 1, mean)
    p.pred1 <- 1 - exp(-exp(log(mean(sex.mean)) + out2$mean$alpha1 * pred.tod + out2$mean$alpha2 * pred.tod^2) * 10)
    
    # Predict p for durations of 1-60 min, averaging over time of day and sex
    duration <- 1:60
    p.pred2 <- 1 - exp(-exp(log(mean(sex.mean))) * duration)
    
    # Visualize analysis
    par(mfrow = c(2,2), mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
    hist(data$ttd, breaks = 40, col = "grey", xlab = "Time to first detection (min)", main = "")
    plot(table(out2$sims.list$n.occ)/length(out2$sims.list$n.occ), xlab = "Number of occupied sites", ylab = "Density", frame = F)
    plot(minutes, p.pred1, xlab = "Minutes after 6.00 hours (i.e., 7.00 ??? 19.00h)", ylab = "Detection prob.", ylim = c(0.6, 1), type = "l", col = "blue", lwd = 3, frame = F)
    plot(duration, p.pred2, xlab = "Survey duration (min)", ylab = "Detection prob.", ylim = c(0, 1), type = "l", col = "blue", lwd = 3, frame = F)
    
    
    
    # 10.12.3 Occupancy models with removal design observation process
    # ------------------------------------------------------------------------
    
    
    
    # 10.13 Occupancy modeling of a community of species
    # ------------------------------------------------------------------------
    
    
    
    # 10.14 Modeling wiggly covariate relationships: penalized splines in hierarchical models
    # ------------------------------------------------------------------------
    
    
    # Execute the function and inspect file produced
    data <- wigglyOcc(seed = 1)
    str(data)
    
    # Convert matrix data into vectors and prepare stuff
    y <- c(data$y)                      # Detection/nondetection data (response)
    Xsite <- data$Xsite                 # Fine as is
    Xsurvey <- c(data$Xsurvey)          # Survey covariate
    site <- rep(1:data$M, data$J)       # Site index
    
    # tmp1 <-spline.prep(Xsite, 20) # This would choose 20 knots for Xsite
    tmp1 <-spline.prep(Xsite, NA)   # Choose variable default number of knots
    tmp2 <-spline.prep(Xsurvey, NA)
    Xocc <- tmp1$X         # Fixed-effects part of covariate Xsite in occ
    Zocc <- tmp1$Z         # Random-effects part of covariate Xsite in occ
    Xdet <- tmp2$X         # Fixed-effects part of covariate Xsite in det
    Zdet <- tmp2$Z         # Random-effects part of covariate Xsite in det
    nk.occ <- length(tmp1$knots)    # Number of knots in occupancy spline
    nk.det <- length(tmp2$knots)    # Number of knots in detection spline
    
    
    # Bundle and summarize data set
    win.data <- list(y1 = y, site = site, M = data$M, Xocc = Xocc, Zocc = Zocc, nk.occ = nk.occ, Xdet = Xdet, Zdet = Zdet, nk.det = nk.det, nobs = data$M*data$J, y2 = y, onesSite = rep(1, 240), onesSurvey = rep(1, 720), Xsite = Xsite, Xsite2 = Xsite^2,  Xsurvey = Xsurvey, Xsurvey2 = Xsurvey^2)
    str(win.data)    # onesSite and onesSurvey are for occ and det intercepts
    
    # Specify two models in one in BUGS language
    cat(file = "hypermodel.txt",
        "model {
        
        # *** Spline model for the data***
        # --------------------------------
        # Priors
        for(k in 1:3){                 # Regression coefficients
        alpha1[k] ~ dnorm(0, 0.1)   # Detection model
        beta1[k] ~ dnorm(0, 0.1)    # Occupancy model
        }
        for(k in 1:nk.occ){ # Random effects at specified knots (occupancy)
        b.occ[k] ~ dnorm(0, tau.b.occ)
        }
        for(k in 1:nk.det){ # Random effects at specified knots (detection)
        b.det[k] ~ dnorm(0, tau.b.det)
        }
        tau.b.occ ~ dgamma(0.01, 0.01)
        tau.b.det ~ dgamma(0.01, 0.01)
        
        # Likelihood
        # Model for latent occupancy state
        for (i in 1:M) {
        z1[i] ~ dbern(psi1[i])
        logit(psi1[i]) <- fix.terms.occ[i] + smooth.terms.occ[i]
        fix.terms.occ[i] <- beta1[1]*Xocc[i,1] + beta1[2]*Xocc[i,2] + beta1[3]*Xocc[i,2]
        smooth.terms.occ[i] <- inprod(b.occ[], Zocc[i,])
        }
        # Model for observations
        for(i in 1:nobs){
        y1[i] ~ dbern(mu.y1[i])
        mu.y1[i] <- z1[site[i]] * p1[i]
        logit(p1[i]) <- fix.terms.det[i] + smooth.terms.det[i]
        fix.terms.det[i] <- alpha1[1]*Xdet[i,1] + alpha1[2]*Xdet[i,2] +
        alpha1[3]*Xdet[i,2]
        smooth.terms.det[i] <- inprod(b.det[], Zdet[i,])
        }
        
        # Derived quantities
        sum.z1 <- sum(z1[])            # Number of occupied sites in sample
        sd.b.occ <- sqrt(1/tau.b.occ)  # SD of spline random effects variance Occ.
        sd.b.det <- sqrt(1/tau.b.det)  # SD of spline random effects variance Det.
        
        
        # *** Polynomial model for same data ***
        # --------------------------------------
        # Priors
        for(k in 1:3){                 # Regression coefficients
        alpha2[k] ~ dnorm(0, 0.1)  # Detection model
        beta2[k] ~ dnorm(0, 0.1)   # Occupancy model
        }
        
        # Likelihood
        # Model for latent occupancy state
        for (i in 1:M) {
        z2[i] ~ dbern(psi2[i])
        logit(psi2[i]) <- beta2[1]*onesSite[i] + beta2[2]*Xsite[i] + beta2[3]*Xsite2[i]
        }
        # Model for observations
        for(i in 1:nobs){
        y2[i] ~ dbern(mu.y2[i])
        mu.y2[i] <- z2[site[i]] * p2[i]
        logit(p2[i]) <- alpha2[1]*onesSurvey[i] + alpha2[2]*Xsurvey[i] +
        alpha2[3] * Xsurvey2[i]
        }
        
        # Derived quantities
        sum.z2 <- sum(z2[])            # Number of occupied sites in sample
        }
        ")
    
    # Initial values
    zst <- apply(data$y, 1, max)
    inits <- function(){list(z1=zst, alpha1=rnorm(3), beta1=rnorm(3), b.occ = rnorm(nk.occ), b.det = rnorm(nk.det), tau.b.occ = runif(1), tau.b.det = runif(1), z2=zst, alpha2=rnorm(3), beta2=rnorm(3))}
    
    # Parameters monitored
    params <- c("alpha1", "beta1", "psi1", "p1", "fix.terms.occ", "smooth.terms.occ", "b.occ", "fix.terms.det", "smooth.terms.det", "b.det", "sum.z1", "sd.b.occ", "sd.b.det", "alpha2", "beta2", "psi2", "p2", "sum.z2")
    
    # MCMC settings
    ni <- 100000   ;   nb <- 10000   ;   nt <- 90   ;   nc <- 3
    
    # Call JAGS from R (ART 95 min) and summarize posteriors
    system.time(fhm <- jags(win.data, inits, params, "hypermodel.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE))
    traceplot(fhm)   ;   print(fhm, 3)
    print(fhm$summary[c(1:6, 2957:2965, 3926),], 3)  # Compare some key estimands
    
    # Plot prediction of psi and p
    par(mfrow = c(1,2), mar = c(5,4,3,2), cex.lab = 1.5, cex.axis = 1.5)
    plot(Xsite, data$psi, main = "Occupancy probability",type = "l", ylim = c(-0.1, 1.1), col = "red", xlab = "Site covariate (Xsite)", ylab = "", lwd = 3)
    #abline(v = tmp1$knots, col = "grey")
    points(Xsite, jitter(data$z, amount = 0.02))
    lines(Xsite, fhm$mean$psi1, col = "blue", lty = 1, lwd = 3)
    lines(Xsite, fhm$mean$psi2, col = "brown", lty = 1, lwd = 3)
    
    plot(Xsurvey[order(data$x.index)], data$p.ordered, main = "Detection probability ",type = "l", ylim = c(-0.1, 1.1), col = "red", xlab = "Survey covariate (Xsurvey)", ylab = "", lwd = 3)
    points(Xsurvey, jitter(y, amount = 0.02))
    #abline(v = tmp2$knots, col = "grey")
    lines(Xsurvey[order(data$x.index)], fhm$mean$p1[order(data$x.index)], col = "blue", lwd = 3)
    lines(Xsurvey[order(data$x.index)], fhm$mean$p2[order(data$x.index)], col = "brown", lwd = 3)
    
    
    
    # 10.15 Summary and outlook
    # ------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    # =========================================================================
    #
    # 11. Hierarchical models for communities
    #
    # =========================================================================
    
    
    
    
    
    # 11.1 Introduction
    # ------------------------------------------------------------------------
    
    
    
    # 11.2 Simulation of a metacommunity
    # ------------------------------------------------------------------------
    
    
    simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
            mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
            mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
            mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)
    
    
    # Execute function with default arguments
    set.seed(1234)
    data <- simComm(type="det/nondet", nsite=30, nrep=3, nspec=100,
                    mean.psi=0.25, sig.lpsi=1, mu.beta.lpsi=0, sig.beta.lpsi=0,
                    mean.lambda=2, sig.loglam=1, mu.beta.loglam=1, sig.beta.loglam=1,
                    mean.p=0.25, sig.lp=1, mu.beta.lp=0, sig.beta.lp=0, show.plot = TRUE)
    # data <- simComm() # same
    
    str(data)
    
    # Some possibly interesting settings of the function
    data <- simComm(nsite = 267, nspec = 190, mean.psi = 0.25, sig.lpsi = 2, mean.p = 0.12, sig.lp = 2) # similar to Swiss MHB
    data <- simComm(mean.psi = 1)         # all species occur at every site
    data <- simComm(mean.p = 1)           # no measurement error (perfect detection)
    
    # Effect of spatial sample size (nsite) on species richness in sample (Ntotal.fs)
    data <- simComm(nsite=50, nspec = 200) # 1-3 are usually missed in sample
    data <- simComm(nsite=30, nspec = 200) # 4-6 usually missed
    data <- simComm(nsite=10, nspec = 200) # around 30 typically missed
    
    # Check for frequentist characteristics of such statistics
    temp <- rep(NA, 100)
    for(i in 1:100){
      cat("\nSimrep", i)
      temp[i] <- simComm(nsite=10, nspec = 200, show.plot = F)$Ntotal.fs
    }
    hist(200-temp, breaks = 30, main = "Number of species in the metacommunity \nthat do not occur in the 10 sampled sites", col = "gold")
    
    # Simulation 1: effects of psi and sd(logit(psi)) on number of species actually occurring in the 50 sampled sites
    simrep <- 50                   # Run 50 simulation reps
    mpsi <- seq(0.01, 0.25,,10)
    slpsi <- seq(0.1, 5,,10)
    results1 <- array(NA, dim = c(10, 10, simrep))
    for(i in 1:10){      # Loop over levels of factor mean.psi (mpsi)
      for(j in 1:10){    # Loop over levels of factor sig.lpsi (slpsi)
        for(k in 1:simrep){
          cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
          tmp <-  simComm(nsite=50, nspec = 200, show.plot = F, mean.psi = mpsi[i],
                          sig.lpsi = slpsi[j])
          results1[i,j,k] <- tmp$Ntotal.fs
        }
      }
    }
    
    
    # Simulation 2: effects of p and sd(logit(p)) on the proportion of the species occurring in the 50 sampled sites that are detected at least once
    simrep <- 50         # Run 50 simulation reps again
    mp <- seq(0.01, 0.25,,10)
    slp <- seq(0.1, 5,,10)
    results2 <- array(NA, dim = c(10, 10, simrep, 2))
    for(i in 1:10){      # Loop over levels of factor mean.p (mp)
      for(j in 1:10){    # Loop over levels of factor sig.lp (slp)
        for(k in 1:simrep){
          cat("\nDim 1:",i, ", Dim 2:", j, ", Simrep", k)
          tmp <-  simComm(nsite=50, nspec = 200, show.plot = F, mean.p = mp[i],
                          sig.lp = slp[j])
          results2[i,j,k,] <- c(tmp$Ntotal.fs, tmp$Ntotal.obs)
        }
      }
    }
    
    
    # Plot these two prediction matrices
    par(mfrow = c(1, 2), mar = c(5,5,2,2), cex.lab = 1.5, cex.axis = 1.5)
    mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
    
    # Plot proportion of species occurring in sampled sites (Fig. 11-4 left)
    z1 <- apply(results1/200, c(1,2), mean) # Prop species occurring
    image(x=mpsi, y=slpsi, z=z1, col = mapPalette(100), axes = T, xlab = "Occupancy probability (psi)", ylab = "Among-species variability in psi")
    contour(x=mpsi, y=slpsi, z=z1, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
    
    # Plot proportion of species detected in sampled sites (Fig. 11-4 right)
    z2 <- apply(results2[,,,2] / results2[,,,1], c(1,2), mean)
    image(x=mp, y=slp, z=z2, col = mapPalette(100), axes = T, xlab = "Detection probability (p)", ylab = "Among-species variability in p")
    contour(x=mp, y=slp, z=z2, add = T, col = "blue", labcex = 1.5, lwd = 1.5)
    
    
    
    # 11.3 Metacommunity data from the Swiss breeding bird survey MHB
    # ------------------------------------------------------------------------
    
    ## Code modified to use the built-in data set MHB2014 instead of the file "MHB_2014.csv"
    data(MHB2014)
    ?MHB2014
    str(MHB2014)
    # NB some of the data preprocessing on p.644 has already been done.
    
    # Check the detection data in 3D array MHB2014$count: site x rep x species
    ( nsite <- nrow(MHB2014$sites) )    # number of sites in Swiss MHB
    nrep <- 3                           # maximum number of replicate surveys per season
    ( nspec <- nrow(MHB2014$species) )  # 158 species occur in the 2014 data
    dim(MHB2014$count) == c(nsite, nrep, nspec) # check
    
    # Create the detection/nondetection (1/0) array
    y <- MHB2014$count ; y[y > 1] <- 1  ## 'Y' replaced with 'y'
    str(y)
    
    # Check data for one species, here chaffinch, and pull them out from 3D array
    (tmp <- y[, , "Common Chaffinch"])
    
    # Frequency distribution of number of surveys actually carried out per site in 2014
    # NB MHB2014$sites$nsurvey gives the number of surveys *planned*.
    table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))
    
    # Which site has all NA data in 2014 ?
    (NAsites <- which(nsurveys == 0) )
    
    # Observed number of occupied sites
    tmp <- apply(y, c(1,3), max, na.rm = TRUE)
    # For the 'all NA' site, max returns -Inf with a warning
    tmp[tmp == -Inf] <- NA         # Change -Inf to NA
    sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))
    
    # Plot species 'occurrence frequency' distribution (not shown)
    plot(sort(obs.occ), xlab = "Species number", ylab = "Number of quads with detections")
    
    # Drop data from species that were not observed in 2014
    toss.out <- which(obs.occ == 0)
    y <- y[,,-toss.out]
    obs.occ <- obs.occ[-toss.out]
    
    # Redefine nspec as the number of species observed in 2014: 145
    ( nspec <- dim(y)[3] )
    
    str(y)
    
    # Get observed number of species per site
    tmp <- apply(y, c(1,3), max, na.rm = TRUE)
    tmp[tmp == "-Inf"] <- NA
    sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts
    
    plot(table(C), xlim = c(0, 60), xlab = "Observed number of species", ylab = "Number of quadrats", frame = FALSE)
    abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)
    
    
    
    # 11.4 Overview of some models for metacommunities
    # ------------------------------------------------------------------------
    
    
    
    # 11.5 Community models that ignore species identity
    # ------------------------------------------------------------------------
    
    
    
    # 11.5.1 Simple Poisson regression for the observed community size
    # ------------------------------------------------------------------------
    # Get covariates (from AHMbook::MHB2014) and standardise them
    
    # Quadrat elevation and forest cover
    orig.ele <- MHB2014$sites$elev
    (mean.ele <- mean(orig.ele, na.rm = TRUE))
    (sd.ele <- sd(orig.ele, na.rm = TRUE))
    ele <- (orig.ele - mean.ele) / sd.ele
    orig.forest <- MHB2014$sites$forest
    (mean.forest <- mean(orig.forest, na.rm = TRUE))
    (sd.forest <- sd(orig.forest, na.rm = TRUE))
    forest <- (orig.forest - mean.forest) / sd.forest
    
    # Average date and duration of survey
    orig.mdate <- apply(MHB2014$date, 1, mean, na.rm = TRUE)
    (mean.mdate <- mean(orig.mdate[-NAsites]))   # drop unsurved site
    (sd.mdate <- sd(orig.mdate[-NAsites]))
    mdate <- (orig.mdate - mean.mdate) / sd.mdate
    mdate[NAsites] <- 0                 # impute mean for missing
    
    orig.mdur <- apply(MHB2014$dur, 1, mean, na.rm = TRUE)
    (mean.mdur <- mean(orig.mdur[-NAsites]))
    (sd.mdur <- sd(orig.mdur[-NAsites]))
    mdur <- (orig.mdur - mean.mdur) / sd.mdur
    mdur[NAsites] <- 0                  # impute mean for missing
    
    # Bundle data and summarize input data for BUGS
    str( win.data <- list(C = C, nsite = length(C), ele = ele, forest = forest,
                          mdate = mdate, mdur = mdur) )
    
    # Specify model in BUGS language
    sink("model1.txt")
    cat("
        model {
        
        # Priors
        gamma0 ~ dnorm(0, 0.001)           # Regression intercept
        for(v in 1:6){                     # Loop over regression coef's
        gamma[v] ~ dnorm(0, 0.001)
        }
        
        # Likelihood for Poisson GLM
        for (i in 1:nsite){
        C[i] ~ dpois(lambda[i])
        log(lambda[i]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) + gamma[3] * forest[i] + gamma[4] * mdate[i] + gamma[5] * pow(mdate[i],2) + gamma[6] * mdur[i]
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))
    
    # Parameters monitored
    params <- c("gamma0", "gamma")
    
    # MCMC settings
    ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART <1 min)
    out1 <- bugs(win.data, inits, params, "model1.txt", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
    
    # Call JAGS from R (ART <1 min), check convergence and summarize posteriors
    library(jagsUI)
    out1J <- jags(win.data, inits, params, "model1.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    traceplot(out1J)    ;    print(out1J, dig = 3)
    
    
    # Plot posterior distributions for potentially 'ecological' parameters
    par(mfrow = c(1,3))
    hist(out1J$sims.list$gamma[,1], breaks = 50, col = "grey", main = "", xlab = "Slope of elevation (linear)")
    hist(out1J$sims.list$gamma[,2], breaks = 50, col = "grey", main = "", xlab = "Slope of elevation (squared)")
    hist(out1J$sims.list$gamma[,3], breaks = 50, col = "grey", main = "", xlab = "Slope of forest")
    
    
    # Get covariate values for prediction
    orig.pred.ele <- seq(250, 2750,, 500)   # 500 vals spread between 250 and 2750
    p.ele <- (orig.pred.ele - mean.ele) / sd.ele
    orig.pred.forest <- seq(1, 100,, 500)
    p.forest <- (orig.pred.forest - mean.forest) / sd.forest
    
    # Compute predictions
    nsamp <- out1J$mcmc.info$n.samples
    pred.ele <- pred.forest <- array(NA, dim = c(500, nsamp))
    for(i in 1:nsamp){
      pred.ele[,i] <- exp(out1J$sims.list$gamma0[i] + out1J$sims.list$gamma[i,1] * p.ele + out1J$sims.list$gamma[i,2]* p.ele^2)
      pred.forest[,i] <- exp(out1J$sims.list$gamma0[i] + out1J$sims.list$gamma[i,3] * p.forest)
    }
    
    # Plot posterior mean and a random sample of 100 from posterior of regression
    selection <- sample(1:nsamp, 100)
    par(mfrow = c(1,3))
    matplot(orig.pred.ele, pred.ele[,selection], ylab = "Predicted species count", xlab = "Elevation (m a.s.l.)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 50), frame = F)
    lines(orig.pred.ele, apply(pred.ele, 1, mean), lwd = 3, col = "blue")
    matplot(orig.pred.forest, pred.forest[,selection], ylab = "Predicted species count", xlab = "Forest cover (%)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 50), frame = F)
    lines(orig.pred.forest, apply(pred.forest, 1, mean), lwd = 3, col = "blue")
    
    
    # Get observed species richness per site and rep and plot
    CC <- apply(y, c(1,2), sum, na.rm = TRUE)
    CC[CC == 0] <- NA            # 0 means not surveyed
    matplot(t(CC), type = 'l', lty = 1, lwd = 2, xlab = "First to third survey", ylab = "Number of species detected", frame = F)  # Fig. 11???6 right
    
    # Get survey date and survey duration and standardise both
    # Survey date (this is Julian date, with day 1 being April 1)
    orig.DAT <- MHB2014$date
    (mean.date <- mean(orig.DAT, na.rm = TRUE))
    (sd.date <- sd(c(orig.DAT), na.rm = TRUE))
    DAT <- (orig.DAT - mean.date) / sd.date      # scale
    DAT[is.na(DAT)] <- 0                         # impute missings
    # Survey duration (in minutes)
    orig.DUR <- MHB2014$dur
    (mean.dur <- mean(orig.DUR, na.rm = TRUE))
    (sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
    DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
    DUR[is.na(DUR)] <- 0                         # mean impute missings
    
    # Bundle data and summarize
    str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )
    
    
    # Specify model in BUGS language
    sink("model2.txt")
    cat("
        model {
        
        # Priors
        gamma0 ~ dnorm(0, 0.001)
        for(v in 1:6){
        gamma[v] ~ dnorm(0, 0.001)
        }
        
        # Likelihood for Poisson GLM
        for (i in 1:M){             # Loop over sites
        for(j in 1:J){            # Loop over occasions
        CC[i,j] ~ dpois(lambda[i,j])
        log(lambda[i,j]) <- gamma0 + gamma[1] * ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
        gamma[6] * DUR[i,j]
        }
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    inits <- function() list(gamma0 = rnorm(1), gamma = rnorm(6))
    
    # Parameters monitored
    params <- c("gamma0", "gamma")
    
    # MCMC settings
    ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART 1.7 min)
    out2 <- bugs(win.data, inits, params, "model2.txt", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
    
    # Call JAGS from R (ART 0.6 min)
    out2J <- jags(win.data, inits, params, "model2.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    traceplot(out2J)   ;   print(out2J, dig = 2)
    
    
    # 11.5.2 Poisson random effects model for the observed community size
    # ------------------------------------------------------------------------
    # Bundle and summarize data set
    str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )
    
    # Specify model in BUGS language
    sink("model3.txt")
    cat("
        model {
        
        # Priors
        mugamma0 ~ dnorm(0, 0.001)     # Hyperparameters
        taugamma0 <- pow(sd.gamma0,-2)
        sd.gamma0 ~ dunif(0, 10)
        for(v in 1:6){                 # Parameters
        gamma[v] ~ dnorm(0, 0.001)
        }
        
        # Likelihood for Poisson GLMM
        for (i in 1:M){                # Loop over sites
        gamma0[i] ~ dnorm(mugamma0, taugamma0)     # site intercepts random now
        for(j in 1:J){               # Loop over repeated measurements
        CC[i,j] ~ dpois(lambda[i,j])
        log(lambda[i,j]) <- gamma0[i] + gamma[1]*ele[i] + gamma[2] * pow(ele[i],2) +
        gamma[3] * forest[i] + gamma[4] * DAT[i,j] + gamma[5] * pow(DAT[i,j],2) +
        gamma[6] * DUR[i,j]
        }
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    inits <- function() list(gamma0 = rnorm(nrow(CC)), gamma = rnorm(6))
    
    # Parameters monitored
    params <- c("mugamma0", "sd.gamma0", "gamma0", "gamma")
    
    # MCMC settings
    ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3
    
    # Call WinBUGS from R (ART 2.9 min)
    out3 <- bugs(win.data, inits, params, "model3.txt", n.chains = nc,
                 n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
    
    # Call JAGS from R (ART 0.7 min)
    out3J <- jags(win.data, inits, params, "model3.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    traceplot(out3J, c('mugamma0', 'sd.gamma0', 'gamma'))   ;    print(out3J, dig = 2)
    
    
    # 11.5.3 N-mixture model for the observed community size
    # -----------------------------------------------------
    # Bundle and summarize data set
    str( win.data <- list(CC = CC, M = nrow(CC), J = ncol(CC), ele = ele, forest = forest, DAT = DAT, DUR = DUR) )
    
    # Specify model in BUGS language
    sink("model4.txt")
    cat("
        model {
        
        # Priors
        alpha0 ~ dnorm(0, 0.01)      # Base-line community detection probability
        beta0 ~ dnorm(0, 0.01)       # Base-line community size (number of species)
        for(v in 1:3){
        alpha[v] ~ dnorm(0, 0.01) # Covariate effects on detection
        beta[v] ~ dnorm(0, 0.01)  # Covariate effects on community size
        }
        
        # Likelihood
        # Ecological model for true community size
        for (i in 1:M){              # Loop over sites
        N[i] ~ dpois(lambda[i])   # Community size
        lambda[i] <- exp(beta0 + beta[1] * ele[i] + beta[2] * pow(ele[i],2) +
        beta[3] * forest[i])
        
        # Observation model for repeated measurements
        for (j in 1:J){          # Loop over occasions
        CC[i,j] ~ dbin(p[i,j], N[i])
        p[i,j] <- 1 / (1 + exp(-lp[i,j]))
        lp[i,j] <- alpha0 + alpha[1] * DAT[i,j] + alpha[2] * pow(DAT[i,j],2) +
        alpha[3] * DUR[i,j]
        # logit(p) = ... causes undefined real result in WinBUGS (but not JAGS)
        }
        }
        }
        ",fill = TRUE)
    sink()
    
    # Define function to generate random initial values
    Nst <- apply(CC, 1, max, na.rm = TRUE) + 1
    Nst[Nst == -Inf] <- max(Nst, na.rm = T)  # Some nonzero val. for unsurv. sites
    inits <- function() list(N = Nst, alpha0 = rnorm(1), alpha = rnorm(3), beta0 = rnorm(1), beta = rnorm(3))
    
    # Parameters monitored
    params <- c("alpha0", "alpha", "beta0", "beta","N")
    
    # MCMC settings
    ni <- 6000   ;   nt <- 4   ;   nb <- 2000   ;   nc <- 3
    
    # Run JAGS from R (ART 1.5 min) in parallel, look at traceplots
    #   and summarize posteriors
    out4 <- jags(win.data, inits, params, "model4.txt",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    traceplot(out4, c('alpha0', 'alpha', 'beta0', 'beta'))  ;  print(out4, 3)
    
    print(tmp <- cbind(out3$summary[c(1, 270:272),c(1:3,7)], out4$summary[5:8, c(1:3, 7)]), 3)
    
    plot(orig.ele, out4$summary[9:275, 1], pch = 16, xlab = "Elevation (m)", ylab = "Species richness", ylim = c(0, 60), frame = F)
    segments(orig.ele, out4$summary[9:275, 3], orig.ele, out4$summary[9:275, 7])
    points(orig.ele+20, C)      # elevation jittered
    
    
    
    
    # 11.6 Community models that retain species identity
    # ------------------------------------------------------------------------
    
    
    # 11.6.1 Simplest community occupancy model: n-fold single species
    #        occupancy model with species treated as fixed effects
    # ------------------------------------------------------------------------
    # Collapse 3D detection/nondetection data to 2D detection frequencies
    ysum <- apply(y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
    ysum[NAsites,] <- NA                     # Have to NA out sites with NA data
    
    # Bundle and summarize data set
    str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )
    
    # Specify model in BUGS language
    sink("model5.txt")
    cat("
        model {
        
        # Priors
        for(k in 1:nspec){          # Loop over species
        psi[k] ~ dunif(0, 1)
        p[k] ~ dunif(0, 1)
        }
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:nspec){          # Loop over species
        for (i in 1:M) {         # Loop over sites
        z[i,k] ~ dbern(psi[k])
        }
        }
        
        # Observation model for observed data y
        for(k in 1:nspec){          # Loop over species
        for (i in 1:M) {
        mup[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mup[i,k], J[i])
        }
        }
        
        # Derived quantities
        for(k in 1:nspec){          # Loop over species
        Nocc.fs[k] <- sum(z[,k]) # Add up number of occupied sites among the 267
        }
        for (i in 1:M) {            # Loop over sites
        Nsite[i] <- sum(z[i,])   # Add up number of occurring species at each site
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
    zst[is.na(zst)] <- 1
    inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))
    
    # Parameters monitored
    params <- c("psi", "p", "Nsite", "Nocc.fs")
    
    # MCMC settings
    ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3
    
    # Call JAGS from R (ART 2.1 min)
    out5 <- jags(win.data, inits, params, "model5.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(4,4))   ;   traceplot(out5)   ;   print(out5, dig = 3)
    
    # Compare observed and estimated site species richness
    par(cex = 1.3)
    twice <- which(MHB2014$sites$nsurvey == 2)
    plot(C[twice], out5$summary[291:557,1][twice], xlab = "Observed number of species", ylab = "Estimated number of species", frame = F, xlim = c(0, 60), ylim = c(0, 70), col = "red", pch = 16)
    segments(C[twice], out5$summary[291:557,3][twice], C[twice], out5$summary[291:557,7][twice], col = "red")
    points(C[-twice], out5$summary[291:557,1][-twice], col = "blue", pch = 16)
    segments(C[-twice], out5$summary[291:557,3][-twice], C[-twice], out5$summary[291:557,7][-twice], col = "blue")
    abline(0,1)
    
    # Observed and estimated number of occupied sites for each species
    # in a table
    cbind(obs.occu = obs.occ, out5$summary[558:702, c(1,3,7)])
    
    # and in a plot
    plot(obs.occ, out5$summary[558:702, 1], xlab = "Observed number of occupied sites", ylab = "Estimated version of quantity", ylim = c(0, 267), frame = F, pch = 16)
    abline(0,1)
    segments(obs.occ, out5$summary[558:702,3], obs.occ, out5$summary[558:702,7], col = "grey", lwd = 2)
    
    
    # Estimated occupancy and detection probability for each species
    plot(out5$summary[1:145,1], out5$summary[146:290,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
    segments(out5$summary[1:145,3], out5$summary[146:290,1], out5$summary[1:145,7], out5$summary[146:290,1], col = "grey", lwd = 2)
    segments(out5$summary[1:145,1], out5$summary[146:290,3], out5$summary[1:145,1], out5$summary[146:290,7], col = "grey", lwd = 2)
    
    
    # 11.6.2 Community occupancy model with bivariate species-specific random effects
    # -----------------------------------------------------
    # Bundle and summarize data set
    str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2], R = matrix(c(5,0,0,1), ncol = 2), df = 3) )
    
    
    # Specify model in BUGS language
    sink("model6.txt")
    cat("
        model {
        
        # Priors
        for(k in 1:nspec){  # Group lpsi and lp together in array eta
        lpsi[k] <- eta[k,1]
        lp[k] <- eta[k,2]
        eta[k, 1:2] ~ dmnorm(mu.eta[], Omega[,])
        }
        # Hyperpriors
        # Priors for mu.lpsi=mu.eta[1] and mu.lp=mu.eta[2]
        # probs = community means of occupancy and detection probability
        for(v in 1:2){
        mu.eta[v] <- log(probs[v] / (1-probs[v]))
        probs[v] ~ dunif(0,1)
        }
        # Prior for variance-covariance matrix
        Omega[1:2, 1:2] ~ dwish(R[,], df)
        Sigma[1:2, 1:2] <- inverse(Omega[,])
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:nspec){
        logit(psi[k]) <- lpsi[k]   # Must take outside of i loop
        for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
        }
        }
        
        # Observation model for observed data y
        for(k in 1:nspec){
        logit(p[k]) <- lp[k]       # Needs to be outside of i loop
        for (i in 1:M) {
        mu.p[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mu.p[i,k], J[i])
        }
        }
        
        # Derived quantities
        rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2])  # Correlation coefficient
        for(k in 1:nspec){
        Nocc.fs[k] <- sum(z[,k])   # Number of occupied sites among the 267
        }
        for (i in 1:M) {
        Nsite[i] <- sum(z[i,])     # Number of occurring species
        }
        }
        ",fill = TRUE)
    sink()
    
    
    # Initial values
    zst <- apply(y, c(1,3), max) # Observed occurrence as starting values for z
    zst[is.na(zst)] <- 1
    inits <- function() list(z = zst, Omega = matrix(c(1,0,0,1), ncol = 2), eta = matrix(0, nrow = nspec, ncol = 2))
    
    # Parameters monitored
    params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")
    
    # MCMC settings
    ni <- 20000   ;   nt <- 15   ;   nb <- 5000   ;   nc <- 3
    
    # Call JAGS from R (ART 12 min), check traceplots and summarize posteriors
    out6 <- jags(win.data, inits, params, "model6.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(3,3))   ;   traceplot(out6, c('mu.eta', 'probs', 'Sigma', 'rho'))
    print(out6, 3)
    
    # Graphically compare some estimates between fixed- and random-effects model
    par(mfrow = c(2,2))      # not shown
    # Species-specific occupancy (probability scale)
    plot(out5$summary[1:145,1], out6$summary[5:149,1], main = "Species-specific occupancy probability")   ;   abline(0,1)
    # Species-specific detection (probability scale)
    plot(out5$summary[146:290,1], out6$summary[150:294,1], main = "Species-specific detection probability")   ;   abline(0,1)
    # Site-specific species richness
    plot(out5$summary[291:557,1], out6$summary[295:561,1], main = "Site-specific species richness (conditional on list of 145 detected)")   ;   abline(0,1)
    # Species-specific number of presences
    plot(out5$summary[558:702,1], out6$summary[562:706,1], main = "Species-specific number of presences (in 267 sites)")   ;   abline(0,1)
    
    # Estimated occupancy and detection probability for each species
    plot(out6$summary[5:149,1], out6$summary[150:294,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
    segments(out6$summary[5:149,3], out6$summary[150:294,1], out6$summary[5:149,7], out6$summary[150:294,1], col = "grey", lwd = 2)
    segments(out6$summary[5:149,1], out6$summary[150:294,3], out6$summary[5:149,1], out6$summary[150:294,7], col = "grey", lwd = 2)
    
    
    
    # 11.6.3 Modeling species-specific effects in community occupancy models
    # ------------------------------------------------------------------------
    # Look at distribution of body mass among 145 observed species (see errata)
    mass <- MHB2014$species$body.mass[-toss.out] # Get  species mass of observed species
    hist(log10(mass), breaks = 40, col = "grey")      # Look at log10
    gmass <- as.numeric(log10(mass) %/% 1.3 + 1)      # size groups 1, 2 and 3
    gmass[gmass == 4] <- 3                            # Mute swan is group 3, too
    
    # Bundle and summarize data set
    str( win.data <- list(ysum = ysum, g = gmass, M = nrow(ysum), J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )
    
    # Specify model in BUGS language
    sink("model7.txt")
    cat("
        model {
        
        # Priors
        for(k in 1:nspec){      # loop over species
        lpsi[k] ~ dnorm(mu.lpsi[g[k]], tau.lpsi[g[k]]) # note g-dependence now
        lp[k] ~ dnorm(mu.lp[g[k]], tau.lp[g[k]])
        }
        
        # Hyperpriors
        for(g in 1:3){          # loop over groups (g)
        mu.lpsi[g] <- logit(mu.psi[g])      # everything is indexed g now
        mu.lp[g] <- logit(mu.p[g])
        mu.psi[g] ~ dunif(0,1)
        mu.p[g] ~ dunif(0,1)
        tau.lpsi[g] <- pow(sd.lpsi[g], -2)
        sd.lpsi[g] ~ dunif(0,5)
        tau.lp[g] <- pow(sd.lp[g], -2)
        sd.lp[g] ~ dunif(0,5)
        }
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:nspec){      # no change at all down here in model
        logit(psi[k]) <- lpsi[k]
        for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
        }
        }
        
        # Observation model for observed data ysum
        for(k in 1:nspec){      # Loop over species
        logit(p[k]) <- lp[k]
        for (i in 1:M) {
        mu.px[i,k] <- z[i,k] * p[k]  # call mu.px to avoid conflict with above
        ysum[i,k] ~ dbin(mu.px[i,k], J[i])
        }
        }
        
        # Derived quantities
        for(k in 1:nspec){          # Loop over species
        Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
        }
        for (i in 1:M) {            # Loop over sites
        Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    zst <- apply(y, c(1,3), max)
    zst[is.na(zst)] <- 1
    inits <- function() list(z = zst)
    
    # Parameters monitored
    params <- c("mu.psi", "mu.lpsi", "sd.lpsi", "mu.p", "mu.lp", "sd.lp")
    
    # MCMC settings
    ni <- 6000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call JAGS from R (ART 6 min), look at convergence and summarize posteriors
    out7 <- jags(win.data, inits, params, "model7.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(3,3))   ;   traceplot(out7)
    print(out7, dig = 3)
    
    
    # Bundle and summarize data set
    logmass <- as.numeric(log10(mass))         # Take log10 of body mass
    str( win.data <- list(ysum = ysum, logmass = logmass, M = nrow(ysum), nsite = nrow(ysum),
                          J = MHB2014$sites$nsurvey, nspec = dim(ysum)[2]) )
    
    
    # Specify model in BUGS language
    sink("model8.txt")
    cat("
        model {
        
        # Priors
        for(k in 1:nspec){              # loop over species
        lpsi[k] ~ dnorm(mu.lpsi[k], tau.lpsi[k]) # now all indexed by k, not g
        tau.lpsi[k] <- 1/var.lpsi[k]
        lp[k] ~ dnorm(mu.lp[k], tau.lp[k])
        tau.lp[k] <- 1/var.lp[k]
        mu.lpsi[k] <- delta0.lpsi + delta1.lpsi * logmass[k]
        mu.lp[k] <- delta0.lp + delta1.lp * logmass[k]
        log(var.lpsi[k]) <- phi0.lpsi + phi1.lpsi * logmass[k]
        log(var.lp[k]) <- phi0.lp + phi1.lp * logmass[k]
        }
        # Priors for regression params for means
        delta0.lpsi ~ dnorm(0, 0.01)
        delta1.lpsi ~ dnorm(0, 0.01)
        delta0.lp ~ dnorm(0, 0.01)
        delta1.lp ~ dnorm(0, 0.01)
        # Priors for regression params for variances
        phi0.lpsi ~ dnorm(0, 0.01)
        phi1.lpsi ~ dnorm(0, 0.01)
        phi0.lp ~ dnorm(0, 0.01)
        phi1.lp ~ dnorm(0, 0.01)
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:nspec){
        logit(psi[k]) <- lpsi[k]
        for (i in 1:M) {
        z[i,k] ~ dbern(psi[k])
        }
        }
        
        # Observation model for observed data ysum
        for(k in 1:nspec){              # Loop over species
        logit(p[k]) <- lp[k]
        for (i in 1:M) {
        mu.p[i,k] <- z[i,k] * p[k]
        ysum[i,k] ~ dbin(mu.p[i,k], J[i])
        }
        }
        
        # Derived quantities
        for(k in 1:nspec){          # Loop over species
        Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
        }
        for (i in 1:M) {            # Loop over sites ## see errata
        Nsite[i] <- sum(z[i,])   # Number of occurring species at each site
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    zst <- apply(y, c(1,3), max)
    zst[is.na(zst)] <- 1
    inits <- function() list(z = zst, delta0.lpsi = rnorm(1), delta1.lpsi = rnorm(1),
                             delta0.lp = rnorm(1), delta1.lp = rnorm(1), phi0.lpsi = rnorm(1),
                             phi1.lpsi = rnorm(1), phi0.lp = rnorm(1), phi1.lp = rnorm(1))
    
    # Parameters monitored
    params <- c("delta0.lpsi", "delta1.lpsi", "delta0.lp", "delta1.lp", "phi0.lpsi",
                "phi1.lpsi", "phi0.lp", "phi1.lp", "psi", "p", "Nocc.fs", "Nsite")
    
    # MCMC settings
    ni <- 12000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call JAGS from R (ART 12 min), look at convergence and summarize posteriors
    out8 <- jags(win.data, inits, params, "model8.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(3,3))
    traceplot(out8, c('delta0.lpsi', 'delta1.lpsi', 'delta0.lp', 'delta1.lp', 'phi0.lpsi', 'phi1.lpsi', 'phi0.lp', 'phi1.lp'))
    print(out8, dig = 3)
    
    
    # Get covariate values for prediction
    predm <- seq(10, 10000,,500)        # Predict for mass of 10g to 10 kg
    pred.logm <- log10(predm)
    
    # Compute predictions (all in one array)
    tmp <- out8$sims.list               # Grab simulation list
    nsamp <- out8$mcmc.info$n.samples   # Number of MCMC samples
    pred <- array(NA, dim = c(500, nsamp, 4)) # Array for predictions
    # [,,1] mu.psi, [,,2] mu.p, [,,3] var.lpsi, [,,4] var.lp
    for(i in 1:nsamp){                  # Fill array
      pred[,i,1] <- plogis(tmp$delta0.lpsi[i] + tmp$delta1.lpsi[i] * pred.logm)
      pred[,i,2] <- plogis(tmp$delta0.lp[i] + tmp$delta1.lp[i] * pred.logm)
      pred[,i,3] <- exp(tmp$phi0.lpsi[i] + tmp$phi1.lpsi[i] * pred.logm)
      pred[,i,4] <- exp(tmp$phi0.lp[i] + tmp$phi1.lp[i] * pred.logm)
    }
    
    # Plot posterior mean and a random sample of 100 from posterior of regression
    selection <- sample(1:nsamp, 100)   # Choose random sample of MCMC output
    par(mfrow = c(2,2), mar = c(5,5,2,2))
    matplot(predm, pred[,selection,1], ylab = "Occupancy mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.4), frame = F)
    lines(predm, apply(pred[,,1], 1, mean), lwd = 3, col = "blue")
    matplot(predm, pred[,selection,2], ylab = "Detection mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.8), frame = F)
    lines(predm, apply(pred[,,2], 1, mean), lwd = 3, col = "blue")
    matplot(predm, pred[,selection,3], ylab = "Occupancy variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
    lines(predm, apply(pred[,,3], 1, mean), lwd = 3, col = "blue")
    matplot(predm, pred[,selection,4], ylab = "Detection variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
    lines(predm, apply(pred[,,4], 1, mean), lwd = 3, col = "blue")
    
    
    
    # 11.6.4 Modeling species richness in a two-step analysis
    # ------------------------------------------------------------------------
    # Extract estimates of N from model 5
    N.pm <- out5$summary[291:557, 1]       # Posterior means of Nsite
    N.psd <- out5$summary[291:557, 2]      # ... posterior sd's of Nsite
    N.cri <- out5$summary[291:557, c(3,7)] # ... CRL's of Nsite
    
    # Plot estimates as a function of elevation
    elev <- MHB2014$sites$elev
    plot(elev, N.pm, xlab = "Altitude (m a.s.l.)", ylab = "Estimated avian species richness", ylim = c(0, 70), frame = F)
    segments(elev, N.cri[,1], elev, N.cri[,2], col = "grey")
    lines(smooth.spline(N.pm ~ elev, w = 1 / N.psd), col = "grey", lwd = 3)
    
    
    # Bundle and summarize data set
    pred.ele <- (seq(200, 2750,5) - mean.ele) / sd.ele # elevation standardised
    str(win.data <- list(ele = ele, N = N.pm, psd = N.psd, n = length(N.pm), pred.ele = pred.ele, npred = length(pred.ele)))
    
    # Define model in BUGS language
    sink("meta.analysis.txt")
    cat("
        model{
        
        # Priors
        for(v in 1:4){         # Priors for intercept and polynomial coefficients
        beta[v] ~ dnorm(0, 0.0001)
        }
        tau.site <- pow(sd.site, -2)
        sd.site ~ dunif(0,10)
        
        # Likelihood
        for (i in 1:n){
        N[i] ~ dnorm(muN[i], tau.psd[i]) # Measurement error model for estimated N
        tau.psd[i] <- pow(psd[i], -2)    # 'Known' part of residual: meas. error
        muN[i] <- beta[1] + beta[2] * ele[i] + beta[3] * pow(ele[i],2) +
        beta[4] * pow(ele[i],3) + eps.site[i] # add another source of uncertainty
        eps.site[i] ~ dnorm(0, tau.site) # this is the usual 'residual'
        }
        # Get predictions for plot
        for(i in 1:npred){
        Npred[i] <- beta[1] + beta[2] * pred.ele[i] + beta[3] * pow(pred.ele[i],2) + beta[4] * pow(pred.ele[i],3)
        }
        } # end model
        ",fill=TRUE)
    sink()
    
    # Initial values, params monitored, and MCMC settings
    inits <- function() list(beta = rnorm(4))
    params <- c("beta", "sd.site", "Npred")
    ni <- 12000   ;   nt <- 10   ;   nb <- 2000   ;   nc <- 3
    
    # Call JAGS and summarize posterior
    out <- jags(win.data, inits, params, "meta.analysis.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
    print(out, 3)
    
    lines(seq(200, 2750,5), out$mean$Npred, col = "blue", lwd = 3)
    matlines(seq(200,2750,5), out$summary[6:516,c(3, 7)], col = "blue", lwd = 2, lty= "dashed")
    
    
    par(mfrow = c(3, 3), mar = c(5,4,3,2))
    for(i in 1:267){
      for(j in 1:267){
        plot(jitter(out5$sims.list$Nsite[,i]), jitter(out5$sims.list$Nsite[,j]),
             main = paste("Joint posterior sites", i, "and", j))
        #   browser()
      }
    }
    
    
    
    # 11.7 The Dorazio-Royle (DR) community occupancy model with data augmentation
    # ----------------------------------------------------------------------------
    
    
    # 11.7.1 The simplest DR community model with data augmentation
    # ------------------------------------------------------------------------
    # Augment data set (DA part)
    nz <- 150                # Number of potential species in superpopulation
    M <- nspec + nz          # Size of augmented data set ('superpopulation')
    yaug <- cbind(ysum, array(0, dim=c(nsite, nz))) # Add all zero histories
    
    # Bundle and summarize data set
    str( win.data <- list(yaug = yaug, nsite = nrow(ysum), nrep = MHB2014$sites$nsurvey, M = M, nspec = nspec, nz = nz) )
    
    # Specify model in BUGS language
    sink("model9.txt")
    cat("
        model {
        
        # Priors to describe heterogeneity among species in community
        for(k in 1:M){                  # Loop over all species in augmented list
        lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
        lp[k] ~ dnorm(mu.lp, tau.lp)
        }
        
        # Hyperpriors to describe full community
        omega ~ dunif(0,1)              # Data augmentation or 'occupancy' parameter
        mu.lpsi ~ dnorm(0,0.001)        # Community mean of occupancy (logit)
        mu.lp ~ dnorm(0,0.001)          # Community mean of detection (logit)
        tau.lpsi <- pow(sd.lpsi, -2)
        sd.lpsi ~ dunif(0,5)            # Species heterogeneity in logit(psi)
        tau.lp <- pow(sd.lp, -2)
        sd.lp ~ dunif(0,5)              # Species heterogeneity in logit(p)
        
        # Superpopulation process:this is the 'paramater expansion' part of PX-DA
        for(k in 1:M){
        w[k] ~ dbern(omega)           # Metacommunity membership indicator
        }                               # (or data augmentation variable)
        
        # Ecological model for latent occurrence z (process model)
        for(k in 1:M){
        mu.psi[k] <- w[k] * psi[k]    # species not part of community zeroed out for z
        logit(psi[k]) <- lpsi[k]
        for (i in 1:nsite) {
        z[i,k] ~ dbern(mu.psi[k])
        }
        }
        
        # Observation model for observed detection frequencies
        for(k in 1:M){
        logit(p[k]) <- lp[k]
        for (i in 1:nsite) {
        mu.p[i,k] <- z[i,k] * p[k]  # non-occurring species are zeroed out for p
        yaug[i,k] ~ dbin(mu.p[i,k], nrep[i])
        }
        }
        
        # Derived quantities
        for(k in 1:M){
        Nocc.fs[k] <- sum(z[,k])     # Number of occupied sites among the 267
        }
        for (i in 1:nsite) {
        Nsite[i] <- sum(z[i,])       # Number of occurring species at each site
        }
        n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species in metacommunity
        Ntotal <- sum(w[])              # Total metacommunity size (= nspec + n0)
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    wst <- rep(1, nspec+nz)                   # Simply set everybody at 'occurring'
    zst <- array(1, dim = c(nsite, nspec+nz)) # ditto for z
    inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz))
    
    # Parameters monitored
    params <- c("mu.lpsi", "sd.lpsi", "mu.lp", "sd.lp", "psi", "p", "Nsite", "Ntotal", "omega", "n0")
    
    # MCMC settings
    ni <- 22000   ;   nt <- 2   ;   nb <- 2000   ;   nc <- 3
    
    # Call JAGS from R (ART 62 min), check convergence and summarize posteriors
    out9 <- jags(win.data, inits, params, "model9.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(2,2)) ; traceplot(out9, c('mu.lpsi', 'sd.lpsi', 'mu.lp', 'sd.lp'))
    print(out9, dig = 3)
    
    
    # Plot posterior distribution of site-specific species richness (Nsite)
    par(mfrow = c(3,3), mar = c(5,4,3,2))
    for(i in 1:267){
      plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
           xlab = "Local species richness", ylab = "", frame = F,
           xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                    max(out9$sims.list$Nsite[,i]) ))
      abline(v = C[i], col = "grey", lwd = 4)
      browser()
    }
    
    # Plot it only for a selection of sites
    par(mfrow = c(3,3), mar = c(5,4,3,2))
    for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
      plot(table(out9$sims.list$Nsite[,i]), main = paste("Quadrat", i),
           xlab = "Local species richness", ylab = "", frame = F,
           xlim = c((min(C[i], out9$sims.list$Nsite[,i], na.rm = T)-2),
                    max(out9$sims.list$Nsite[,i]) ))
      abline(v = C[i], col = "grey", lwd = 4)
    }
    
    
    # Plot posterior distribution of total species richness (Ntotal)
    plot(table(out9$sims.list$Ntotal), main = "", ylab = "", xlab = "Avian metacommunity size in Swiss MHB survey (267 1km2 quadrats)", frame = F, xlim = c(144, 245))
    abline(v = nspec, col = "grey", lwd = 4)
    
    
    
    # 11.7.2 Dorazio-Royle community model with covariates
    # ------------------------------------------------------------------------
    # Augment data set: choose one of two different priors on Ntotal
    nz <- 250                 # Use for vague prior on Ntotal: M = 395
    nz <- 215 - nspec         # Use for informative prior on Ntotal: M = 215
    yaug <- array(0, dim=c(nsite, nrep, nspec+nz)) # array with only zeroes
    yaug[,,1:nspec] <- y      # copy into it the observed data
    
    # Create same NA pattern in augmented species as in the observed species
    missings <- is.na(yaug[,,1]) # e.g., third survey in high-elevation quads
    for(k in (nspec+1):(nspec+nz)){
      yaug[,,k][missings] <- NA
    }
    
    # Bundle and summarize data
    str(win.data <- list(y = yaug, nsite = dim(y)[1], nrep = dim(y)[2], nspec = dim(y)[3], nz = nz, M = nspec + nz, ele = ele, forest = forest, DAT = DAT, DUR = DUR) )
    
    
    # Specify model in BUGS language
    sink("model10.txt")
    cat("
        model {
        
        # Priors
        omega ~ dunif(0,1)
        # Priors for species-specific effects in occupancy and detection
        for(k in 1:M){
        lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)    # Hyperparams describe community
        betalpsi1[k] ~ dnorm(mu.betalpsi1, tau.betalpsi1)
        betalpsi2[k] ~ dnorm(mu.betalpsi2, tau.betalpsi2)
        betalpsi3[k] ~ dnorm(mu.betalpsi3, tau.betalpsi3)
        lp[k] ~ dnorm(mu.lp, tau.lp)
        betalp1[k] ~ dnorm(mu.betalp1, tau.betalp1)
        betalp2[k] ~ dnorm(mu.betalp2, tau.betalp2)
        betalp3[k] ~ dnorm(mu.betalp3, tau.betalp3)
        }
        
        # Hyperpriors
        # For the model of occupancy
        mu.lpsi ~ dnorm(0,0.01)
        tau.lpsi <- pow(sd.lpsi, -2)
        sd.lpsi ~ dunif(0,8)   # as always, bounds of uniform chosen by trial and error
        mu.betalpsi1 ~ dnorm(0,0.1)
        tau.betalpsi1 <- pow(sd.betalpsi1, -2)
        sd.betalpsi1 ~ dunif(0, 4)
        mu.betalpsi2 ~ dnorm(0,0.1)
        tau.betalpsi2 <- pow(sd.betalpsi2, -2)
        sd.betalpsi2 ~ dunif(0,2)
        mu.betalpsi3 ~ dnorm(0,0.1)
        tau.betalpsi3 <- pow(sd.betalpsi3, -2)
        sd.betalpsi3 ~ dunif(0,2)
        
        # For the model of detection
        mu.lp ~ dnorm(0,0.1)
        tau.lp <- pow(sd.lp, -2)
        sd.lp ~ dunif(0, 2)
        mu.betalp1 ~ dnorm(0,0.1)
        tau.betalp1 <- pow(sd.betalp1, -2)
        sd.betalp1 ~ dunif(0,1)
        mu.betalp2 ~ dnorm(0,0.1)
        tau.betalp2 <- pow(sd.betalp2, -2)
        sd.betalp2 ~ dunif(0,1)
        mu.betalp3 ~ dnorm(0,0.1)
        tau.betalp3 <- pow(sd.betalp3, -2)
        sd.betalp3 ~ dunif(0,1)
        
        # Superpopulation process: Ntotal species sampled out of M available
        for(k in 1:M){
        w[k] ~ dbern(omega)
        }
        
        # Ecological model for true occurrence (process model)
        for(k in 1:M){
        for (i in 1:nsite) {
        logit(psi[i,k]) <- lpsi[k] + betalpsi1[k] * ele[i] +
        betalpsi2[k] * pow(ele[i],2) + betalpsi3[k] * forest[i]
        mu.psi[i,k] <- w[k] * psi[i,k]
        z[i,k] ~ dbern(mu.psi[i,k])
        }
        }
        
        # Observation model for replicated detection/nondetection observations
        for(k in 1:M){
        for (i in 1:nsite){
        for(j in 1:nrep){
        logit(p[i,j,k]) <- lp[k] + betalp1[k] * DAT[i,j] +
        betalp2[k] * pow(DAT[i,j],2) + betalp3[k] * DUR[i,j]
        mu.p[i,j,k] <- z[i,k] * p[i,j,k]
        y[i,j,k] ~ dbern(mu.p[i,j,k])
        }
        }
        }
        
        # Derived quantities
        #for(k in 1:M){
        #   Nocc.fs[k] <- sum(z[,k])       # Number of occupied sites among the 267
        #}
        for (i in 1:nsite){
        Nsite[i] <- sum(z[i,])          # Number of occurring species at each site
        }
        n0 <- sum(w[(nspec+1):(nspec+nz)]) # Number of unseen species
        Ntotal <- sum(w[])                 # Total metacommunity size
        
        # Vectors to save (S for ???save???; discard posterior samples for
        # all minus 1 of the potential species to save disk space)
        # we do this for nz = 250 (i.e., M = 395)
        lpsiS[1:(nspec+1)] <- lpsi[1:(nspec+1)]
        betalpsi1S[1:(nspec+1)] <- betalpsi1[1:(nspec+1)]
        betalpsi2S[1:(nspec+1)] <- betalpsi2[1:(nspec+1)]
        betalpsi3S[1:(nspec+1)] <- betalpsi3[1:(nspec+1)]
        lpS[1:(nspec+1)] <- lp[1:(nspec+1)]
        betalp1S[1:(nspec+1)] <- betalp1[1:(nspec+1)]
        betalp2S[1:(nspec+1)] <- betalp2[1:(nspec+1)]
        betalp3S[1:(nspec+1)] <- betalp3[1:(nspec+1)]
        }
        ",fill = TRUE)
    sink()
    
    
    # Initial values
    wst <- rep(1, nspec+nz)                   # Simply set everybody at occurring
    zst <- array(1, dim = c(nsite, nspec+nz)) # ditto
    inits <- function() list(z = zst, w = wst, lpsi = rnorm(n = nspec+nz), betalpsi1 = rnorm(n = nspec+nz), betalpsi2 = rnorm(n = nspec+nz), betalpsi3 = rnorm(n = nspec+nz), lp = rnorm(n = nspec+nz), betalp1 = rnorm(n = nspec+nz), betalp2 = rnorm(n = nspec+nz), betalp3 = rnorm(n = nspec+nz))
    
    # Set 1
    params1 <- c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal", "Nsite")
    
    # MCMC settings
    ni <- 15000   ;   nt <- 10   ;   nb <- 5000   ;   nc <- 3
    
    # Run JAGS, check convergence and summarize posteriors
    out101 <- jags(win.data, inits, params1, "model10.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    par(mfrow = c(2, 2))
    traceplot(out101, c(c("omega", "mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "mu.lp", "sd.lp", "mu.betalp1", "sd.betalp1", "mu.betalp2", "sd.betalp2", "mu.betalp3", "sd.betalp3", "Ntotal")) )
    
    # Set 2
    params2 <- c("mu.lpsi", "sd.lpsi", "mu.betalpsi1", "sd.betalpsi1", "mu.betalpsi2", "sd.betalpsi2", "mu.betalpsi3", "sd.betalpsi3", "lpsi", "betalpsi1", "betalpsi2", "betalpsi3", "lp", "betalp1", "betalp2", "betalp3", "z", "w")
    ni <- 12000   ;   nt <- 20   ;   nb <- 2000   ;   nc <- 3
    out102 <- jags.basic(win.data, inits, params2, "model10.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    library(coda)
    all10 <- as.matrix(out102) # Put output from 3 chains into a matrix
    summary(out102)            # May take a loooong time
    gelman.diag(out102)        # ditto
    
    
    # Comparison of main hyperparameters when M = 215 and with M = 395
    # (not all code to produce this output is shown)
    print(cbind(out10.215$summary[1:17,c(1:3, 7)], out10.395$summary[1:17, c(1:3, 7)]), 2)
    
    
    out10 <- out101
    
    
    par(mfrow = c(1,2))       # Fig. 11-16
    psi.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lpsi, sd = out10$mean$sd.lpsi))
    p.sample <- plogis(rnorm(10^6, mean = out10$mean$mu.lp, sd = out10$mean$sd.lp))
    hist(psi.sample, freq = F, breaks = 50, col = "grey", xlab = "Species occupancy probability", ylab = "Density", main = "")
    hist(p.sample, freq = F, breaks = 50, col = "grey", xlab = "Species detection probability", ylab = "Density", main = "")
    summary(psi.sample)   ;   summary(p.sample)
    
    par(mfrow = c(2,4))  # Among-species variability in parameters (not shown)
    hist(out10$sims.list$sd.lpsi, breaks = 100, col = "grey", xlim = c(0,6), main = "Occupancy: intercept")
    abline(v = mean(out10$sims.list$sd.lpsi), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalpsi1, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: linear effect of elevation")
    abline(v = mean(out10$sims.list$sd.betalpsi1), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalpsi2, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: quadratic effect of elevation")
    abline(v = mean(out10$sims.list$sd.betalpsi2), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalpsi3, breaks = 100, col = "grey", xlim = c(0,3), main = "Occupancy: linear effect of forest cover")
    abline(v = mean(out10$sims.list$sd.betalpsi3), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.lp, breaks = 100, col = "grey", xlim = c(0,2), main = "Detection: intercept")
    abline(v = mean(out10$sims.list$sd.lp), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalp1, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: linear effect of survey date")
    abline(v = mean(out10$sims.list$sd.betalp1), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalp2, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: quadratic linear effect of survey date")
    abline(v = mean(out10$sims.list$sd.betalp2), col = "blue", lwd = 3)
    hist(out10$sims.list$sd.betalp3, breaks = 100, col = "grey", xlim = c(0,1), main = "Detection: linear effect of survey duration")
    abline(v = mean(out10$sims.list$sd.betalp3), col = "blue", lwd = 3)
    
    
    # Visualize covariate mean relationships for the average species
    o.ele <- seq(200, 2500,,500)               # Get covariate values for prediction
    o.for <- seq(0, 100,,500)
    o.dat <- seq(15, 120,,500)
    o.dur <- seq(100, 420,,500)
    ele.pred <- (o.ele - mean.ele) / sd.ele
    for.pred <- (o.for - mean.forest) / sd.forest
    dat.pred <- (o.dat - mean.date) / sd.date
    dur.pred <- (o.dur - mean.dur) / sd.dur
    
    # Predict occupancy for elevation and forest and detection for date and duration
    # Put all fourpredictions into a single
    str( tmp <- out10$sims.list )              # grab MCMC samples
    nsamp <- length(tmp[[1]])    # number of mcmc samples
    predC <- array(NA, dim = c(500, nsamp, 4)) # "C" for 'community mean'
    for(i in 1:nsamp){
      predC[,i,1] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi1[i] * ele.pred +
                              tmp$mu.betalpsi2[i] * ele.pred^2 )
      predC[,i,2] <- plogis(tmp$mu.lpsi[i] + tmp$mu.betalpsi3[i] * for.pred)
      predC[,i,3] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp1[i] * dat.pred +
                              tmp$mu.betalp2[i] * dat.pred^2 )
      predC[,i,4] <- plogis(tmp$mu.lp[i] + tmp$mu.betalp3[i] * dur.pred)
    }
    
    # Get posterior means and 95% CRIs and plot (Fig. 11???17)
    pmC <- apply(predC, c(1,3), mean)
    criC <- apply(predC, c(1,3), function(x) quantile(x, prob = c(0.025, 0.975)))
    
    par(mfrow = c(2, 2))
    plot(o.ele, pmC[,1], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Elevation (m a.s.l)", ylab = "Community mean occupancy")
    matlines(o.ele, t(criC[,,1]), col = "grey", lty = 1)
    plot(o.for, pmC[,2], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0, 0.05), xlab = "Forest cover", ylab = "Community mean occupancy")
    matlines(o.for, t(criC[,,2]), col = "grey", lty = 1)
    plot(o.dat, pmC[,3], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey date", ylab = "Community mean detection")
    matlines(o.dat, t(criC[,,3]), col = "grey", lty = 1)
    plot(o.dur, pmC[,4], col = "blue", lwd = 3, type = 'l', lty = 1, frame = F, ylim = c(0.2, 0.8), xlab = "Survey duration", ylab = "Community mean detection")
    matlines(o.dur, t(criC[,,4]), col = "grey", lty = 1)
    
    
    # Plot posterior distribution of site-specific species richness (Nsite)
    par(mfrow = c(3,3), mar = c(5,4,3,2))
    for(i in 1:267){
      plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
           xlab = "Local species richness", ylab = "", frame = F,
           xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
                    max(out10$sims.list$Nsite[,i]) ))
      abline(v = C[i], col = "grey", lwd = 4)
      browser()
    }
    
    # Plot it only for a selection of sites (Fig. 11-18)
    par(mfrow = c(3,3), mar = c(5,4,3,2))
    for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
      plot(table(out10$sims.list$Nsite[,i]), main = paste("Quadrat", i),
           xlab = "Local species richness", ylab = "", frame = F,
           xlim = c((min(C[i], out10$sims.list$Nsite[,i], na.rm = T)-2),
                    max(out10$sims.list$Nsite[,i]) ))
      abline(v = C[i], col = "grey", lwd = 4)
    }
    
    # Plot Nsite estimates under models 9 & 10 vs. elevation (Fig. 11-19)
    offset <- 30    # Set off elevation for better visibility
    plot(elev, out9$mean$Nsite, xlab = "Elevation (metres)", ylab = "Community size estimate (Nsite)", frame = F, ylim = c(0,60), pch = 16) # black: model 9
    lines(smooth.spline(out9$mean$Nsite ~ elev), lwd = 3)
    points(elev+offset, out10$mean$Nsite, pch = 16, col = "blue") # red: model 10
    lines(smooth.spline(out10$mean$Nsite ~ elev), lwd = 3, col = "blue")
    
    
    str(all10)                    # look at the MCMC output
    pm <- apply(all10, 2, mean)    # Get posterior means and 95% CRIs
    cri <- apply(all10, 2, function(x) quantile(x, prob = c(0.025, 0.975))) # CRIs
    
    
    # Effects of date (linear and quadratic) and of duration on detection
    #par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # Can put all three in one
    par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
    # Date linear (Fig. 11 ??? 20 left)
    plot(pm[1:145], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (linear) on detection", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 1:145], 1:145, cri[2, 1:145], 1:145, col = "grey", lwd = 1)
    sig1 <- (cri[1, 1:145] * cri[2, 1:145]) > 0
    segments(cri[1, 1:145][sig1 == 1], (1:145)[sig1 == 1], cri[2, 1:145][sig1 == 1], (1:145)[sig1 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[11,1], lwd = 3, col = "red")
    abline(v = out101$summary[11,c(3,7)], lwd = 2, col = "red", lty = 2)
    
    
    
    # Date quadratic (not shown)
    plot(pm[216:360], 1:145, xlim = c(-1.5, 1.5), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of date (quadratic) on detection", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 216:360], 1:145, cri[2, 216:360], 1:145, col = "grey", lwd = 1)
    sig2 <- (cri[1, 216:360] * cri[2, 216:360]) > 0
    segments(cri[1, 216:360][sig2 == 1], (1:145)[sig2 == 1], cri[2, 216:360][sig2 == 1], (1:145)[sig2 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[13,1], lwd = 3, col = "red")
    abline(v = out101$summary[13, c(3,7)], lwd = 3, col = "red", lty = 2)
    
    
    # Survey duration (Fig. 11-20 right)
    plot(pm[431:575], 1:145, xlim = c(-0.5, 1), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of survey duration on detection", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 431:575], 1:145, cri[2, 431:575], 1:145, col = "grey", lwd = 1)
    sig3 <- (cri[1, 431:575] * cri[2, 431:575]) > 0
    segments(cri[1, 431:575][sig3 == 1], (1:145)[sig3 == 1], cri[2, 431:575][sig3 == 1], (1:145)[sig3 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[15,1], lwd = 3, col = "red")
    abline(v = out101$summary[15, c(3,7)], lwd = 3, col = "red", lty = 2)
    
    
    # Effects of elevation (linear and quadratic) and of forest on occupancy
    # par(mfrow = c(1,3), cex.lab = 1.3, cex.axis = 1.3) # can do all in one
    # Effect of elevation (linear) on occupancy probability (Fig. 11-21)
    plot(pm[646:790], 1:145, xlim = c(-8, 8), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (linear) on occupancy", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 646:790], 1:145, cri[2, 646:790], 1:145, col = "grey", lwd = 1)
    sig4 <- (cri[1, 646:790] * cri[2, 646:790]) > 0
    segments(cri[1, 646:790][sig4 == 1], (1:145)[sig4 == 1], cri[2, 646:790][sig4 == 1], (1:145)[sig4 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[3,1], lwd = 3, col = "red")
    abline(v = out101$summary[3,c(3,7)], lwd = 3, col = "red", lty = 2)
    
    
    # Effect of elevation (quadratic) on occupancy probability (Fig. 11-22)
    plot(pm[861:1005], 1:145, xlim = c(-4, 2), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of elevation (quadratic) on occupancy", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 861:1005], 1:145, cri[2, 861:1005], 1:145, col = "grey", lwd=1)
    sig5 <- (cri[1, 861:1005] * cri[2, 861:1005]) > 0
    segments(cri[1, 861:1005][sig5 == 1], (1:145)[sig5 == 1], cri[2, 861:1005][sig5 == 1], (1:145)[sig5 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[5,1], lwd = 3, col = "red")
    abline(v = out101$summary[5,c(3,7)], lwd = 3, col = "red", lty = 2)
    
    
    # Effect of forest (linear) on occupancy probability (Fig. 11-23)
    plot(pm[1076:1220], 1:145, xlim = c(-3, 4), xlab = "Parameter estimate", ylab = "Species number", main = "Effect of forest cover on occupancy", pch = 16)
    abline(v = 0, lwd = 2, col = "black")
    segments(cri[1, 1076:1220], 1:145, cri[2, 1076:1220],1:145, col = "grey", lwd=1)
    sig6 <- (cri[1, 1076:1220] * cri[2, 1076:1220]) > 0
    segments(cri[1, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], cri[2, 1076:1220][sig6 == 1], (1:145)[sig6 == 1], col = "blue", lwd = 2)
    abline(v = out101$summary[7,1], lwd = 3, col = "red")
    abline(v = out101$summary[7,c(3,7)], lwd = 3, col = "red", lty = 2)
    negsig6 <- (cri[1, 1076:1220] < 0 & cri[2, 1076:1220] < 0) == 1 # sig negative
    possig6 <- (cri[1, 1076:1220] > 0 & cri[2, 1076:1220] > 0) == 1 # sig positive
    
    
    # Predict detection for date and duration and occupancy for elevation and forest
    # for each of the 145 observed species
    predS <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "S" for 'species'
    p.coef <- cbind(lp=pm[1292:1436], betalp1 = pm[1:145], betalp2 = pm[216:360], betalp3 = pm[431:575])
    psi.coef <- cbind(lpsi=pm[1507:1651], betalpsi1 = pm[646:790], betalpsi2 = pm[861:1005], betalpsi3 = pm[1076:1220])
    
    for(i in 1:nspec){          # Loop over 145 observed species
      predS[,i,1] <- plogis(p.coef[i,1] + p.coef[i,2] * dat.pred +
                              p.coef[i,3] * dat.pred^2 )     # p ~ date
      predS[,i,2] <- plogis(p.coef[i,1] + p.coef[i,4] * dur.pred) # p ~ duration
      predS[,i,3] <- plogis(psi.coef[i,1] + psi.coef[i,2] * ele.pred +
                              psi.coef[i,3] * ele.pred^2 )     # psi ~ elevation
      predS[,i,4] <- plogis(psi.coef[i,1] + psi.coef[i,4] * for.pred) # psi ~ forest
    }
    
    # Plots for detection probability and survey date and duration (Fig. 11-24)
    par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
    plot(o.dat, predS[,1,1], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
         ylab = "Detection probability")
    for(i in 2:145){
      lines(o.dat, predS[,i,1], col = i, lwd = 3)
    }
    
    plot(o.dur, predS[,1,2], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Survey duration (min)",
         ylab = "Detection probability")
    for(i in 2:145){
      lines(o.dur, predS[,i,2], col = i, lwd = 3)
    }
    
    
    # Plots for occupancy probability and elevation and forest cover (Fig. 11-25)
    par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
    plot(o.ele, predS[,1,3], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Elevation (m a.s.l.)",
         ylab = "Occupancy probability")
    for(i in 2:145){
      lines(o.ele, predS[,i,3], col = i, lwd = 3)
    }
    
    plot(o.for, predS[,1,4], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Forest cover (%)", ylab = "Occupancy probability")
    for(i in 2:145){
      lines(o.for, predS[,i,4], col = i, lwd = 3)
    }
    
    
    
    # 11.8 Inferences based on the estimated Z matrix: similarity among sites and species
    # -----------------------------------------------------------------------------------
    
    
    # Plug MCMC samples for full z matrix into 3D array
    str(all10)
    nsite <- 267
    nspec <- 215
    nsamp <- dim(all10)[1]        # 1200 MCMC samples
    z <- array(NA, dim = c(nsite, nspec, nsamp))
    Jacc <- array(NA, dim = c(nsite, nspec, nsamp))
    for(j in 1:nsamp){    # Fill z matrix by column (default)
      cat(paste("\nMCMC sample", j, "\n"))
      z[,,j] <- all10[j, 1937:59341]
    }
    
    # Restrict computations to observed species
    zobs <- z[,1:145,]      # Species 1 to 145
    
    # Compute Jaccard index for sites and for species
    Jsite <- array(NA, dim = c(nsite, nsamp))
    Jspec <- array(NA, dim = c(145, nsamp))
    
    
    # Choose reference site and species for Jaccard indices
    ref.site <- 1         # Just choose first site
    ref.species <- 13     # European Sparrowhawk (check object 'obs.occ')
    
    # Get posterior distributions for Jsite and Jspec (for references)
    for(k in 1:nsamp){
      for(i in 1:nsite){ # Jaccard index for sites (in terms of shared species)
        Jsite[i,k] <- sum(zobs[ref.site,,k] * zobs[i,,k]) /
          (sum(zobs[ref.site,,k]) + sum(zobs[i,,k]) -
             sum(zobs[ref.site,,k] * zobs[i,,k]))
      }
      for(i in 1:(nspec-nz)){ # Jacc. index for species (in terms of shared sites)
        Jspec[i,k] <- sum(zobs[,ref.species,k] * zobs[,i,k]) /
          (sum(zobs[,ref.species,k]) + sum(zobs[,i,k]) -
             sum(zobs[,ref.species,k] * zobs[,i,k]))
      }
    }
    # NA's arise when a site has no species or a species no sites
    
    # Get posterior means, standard deviations and 95% CRI
    # Jaccard index for sites, compared to reference site 1
    pm <- apply(Jsite, 1, mean, na.rm = TRUE)  # Post. mean of Jsite wrt. site 1
    psd <- apply(Jsite, 1, sd, na.rm = TRUE)   # Post. mean of Jsite wrt. site 1
    cri <- apply(Jsite, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
    cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
    
    
    # Make a map of Jaccard site indices (Fig. 11-26)
    x <- 3        # poportional size of plotting symbol
    plot(MHB2014$sites$coordx, MHB2014$sites$coordy, xlab = "x coordinate", ylab = "y coordinate", cex = x*pm, asp = 1, pch = 16)
    points(MHB2014$sites$coordx[which(pm == 1)], MHB2014$sites$coordy[which(pm == 1)], cex = x*pm, col = "red", pch = 16)
    
    
    # Jaccard index for species, compared to reference species
    # (species 13, European Sparrowhawk)
    pm <- apply(Jspec, 1, mean, na.rm = TRUE)  # Post. mean of Jspec wrt. species 1
    psd <- apply(Jspec, 1, sd, na.rm = TRUE)   # Post. mean of Jspec wrt. species 1
    cri <- apply(Jspec, 1, function(x) quantile(x, prob = c(0.025, 0.975), na.rm = TRUE)) # CRI
    tmp <- cbind('post. mean' = pm, 'post. sd' = psd, '2.5%' = cri[1,], '97.5%' = cri[2,])
    rownames(tmp) <- names(obs.occ)
    print(tmp])          # print in systematic order
    print(tmp[rev(order(tmp[,1])),]) # print in order of decreasing Jacc. values
    plot(1:145, tmp[rev(order(tmp[,1])),1])   # can also plot
    
    
    
    
    # 11.9 Species richness maps and species accumulation curves
    # ------------------------------------------------------------------------
    
    
    # Get Swiss landscape data and standardise covariates as for model 10
    library(unmarked)
    data(Switzerland)
    ch <- Switzerland
    ELE <- (ch$elevation - mean.ele) / sd.ele
    FOREST <- (ch$forest - mean.forest) / sd.forest
    
    nsamp <- nrow(all10)            # 1200   ..... far too many
    nkm2 <- length(ch[[1]])         # 42275, that's a LOT!
    select.samp <- sort(sample(1:nsamp, 50)) # Chose random sample of 50
    nsamp <- length(select.samp)    # new sample size 50
    
    # Create posterior predictive distribution for Z for Swiss landscape
    str( zCH <- array(NA, dim = c(nkm2, 215, nsamp)) ) # BIG array !
    W <- all10[,1722:1936]          # Grab MCMC samples from w
    LPSI <- all10[,1507:1721]       # Grab MCMC samples from logit(psi)
    BETALPSI1 <- all10[,646:860]    # Grab MCMC samples from betalpsi1
    BETALPSI2 <- all10[,861:1075]   # Grab MCMC samples from betalpsi2
    BETALPSI3 <- all10[,1076:1290]  # Grab MCMC samples from betalpsi3
    for(i in 1:nkm2){               # takes about 5 mins !
      cat(paste("\nQuadrat", i, "\n"))
      for(u in 1:length(select.samp)){
        psi <- W[select.samp[u],] * plogis(LPSI[select.samp[u],] +
                                             BETALPSI1[select.samp[u],] * ELE[i] +
                                             BETALPSI2[select.samp[u],] * ELE[i]^2 +
                                             BETALPSI3[select.samp[u],] * FOREST[i] )
        zCH[i,,u] <- rbinom(215, 1, psi)
      }
    }
    
    # Compute posterior distribution of species richness by collapsing z array
    SR <- apply(zCH, c(1,3), sum)   # posterior distribution
    pmSR <- apply(SR, 1, mean)      # posterior mean
    sdSR <- apply(SR, 1, sd)        # posterior standard deviation
    
    
    library(raster)
    library(rgdal)
    par(mfrow = c(1,2), mar = c(2,2,3,5))
    # Posterior mean map
    r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = pmSR))
    elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
    elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
    r1 <- mask(r1, elev)
    mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
    plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
    lakes <- readOGR(".", "lakes")
    rivers <- readOGR(".", "rivers")
    border <- readOGR(".", "border")
    plot(rivers, col = "dodgerblue", add = TRUE)
    plot(border, col = "transparent", lwd = 1.5, add = TRUE)
    plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
    
    # Posterior standard deviation map
    r1 <- rasterFromXYZ(data.frame(x = ch$x, y = ch$y, z = sdSR))
    elev <- rasterFromXYZ(cbind(ch$x, ch$y,ch$elevation))
    elev[elev > 2250] <- NA         # Mask areas > 2250 m a.s.l.
    r1 <- mask(r1, elev)
    mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
    plot(r1, col = mapPalette(100), axes = F, box = FALSE, main ="")
    lakes <- readOGR(".", "lakes")
    rivers <- readOGR(".", "rivers")
    border <- readOGR(".", "border")
    plot(rivers, col = "dodgerblue", add = TRUE)
    plot(border, col = "transparent", lwd = 1.5, add = TRUE)
    plot(lakes, col = "skyblue", border = "royalblue", add = TRUE)
    
    
    # Get 3,000 posterior samples of omega, and the mean and sd hyperparameters
    omega <- out101$sims.list$omega
    mu.lpsi <- out101$sims.list$mu.lpsi
    str( sd.lpsi <- out101$sims.list$sd.lpsi )    # Confirms we have 3,000 draws
    
    # compute posterior predictions of species occurrence probabilities
    nsites <- 100
    ndraws <- length(omega)
    Nmax <- 215
    psi <- matrix(NA, nrow=ndraws, ncol=Nmax)
    for (i in 1:ndraws) {
      w <- rbinom(215, 1, omega[i])
      psi[i,] <- w * plogis(rnorm(Nmax, mean = mu.lpsi[i], sd=sd.lpsi[i]))
    }
    
    # compute posterior predictions of species presence at each site
    z <- array(NA, dim=c(ndraws, Nmax, nsites))
    for (i in 1:ndraws) {
      for (j in 1:Nmax) {
        z[i,j, ] <- rbinom(nsites, size=1, prob=psi[i,j])
      }
    }
    
    # compute posterior predictions of cumulative number of species present
    Ntot <- matrix(NA, nrow=ndraws, ncol=nsites)
    for (i in 1:ndraws) {
      for (j in 1:nsites) {
        zsum <- rep(NA, Nmax)
        if (j>1) {
          zsum <- apply(z[i, , 1:j], 1, sum)
        }
        else {
          zsum <- z[i, , 1]
        }
        Ntot[i,j] <- sum(zsum>0)
      }
    }                        # takes about 4 min
    
    # compute summary stats of species accumulation curve
    nSpeciesPresent <- matrix(NA, nrow=3, ncol=nsites)
    for (j in 1:nsites) {
      x <- Ntot[,j]
      nSpeciesPresent[1, j] <- mean(x)
      nSpeciesPresent[2:3, j] <- quantile(x, probs=c(0.05, 0.95))
    }
    
    # Plot species accumulation curve
    ylim = c(min(nSpeciesPresent[2,]), max(nSpeciesPresent[3,]))
    plot(1:nsites, nSpeciesPresent[1,], pch=16, ylim=ylim, type="b",
         xlab="Number of sample locations", ylab="Number of occurring species",
         las=1, cex.axis=1.2, cex.lab=1.5, cex=1.2, frame = F)
    segments(1:nsites, nSpeciesPresent[2,], 1:nsites, nSpeciesPresent[3,])
    
    
    
    # 11.10 Community N-mixture models
    # ------------------------------------------------------------------------
    
    # Organize counts in 3D array: site x rep x species
    Yc <- MHB2014$counts
    str(Yc)
    
    # Observed maximum and mean maximum count per species
    tmp <- apply(Yc, c(1,3), max, na.rm = TRUE)
    tmp[tmp == -Inf] <- NA         # 1 quadrat with NA data in 2014
    sort(round(meanmax <- apply(tmp, 2, mean, na.rm = TRUE), 3)) # mean of max
    sort(obs.max.C <- apply(tmp, 2, max, na.rm = TRUE))          # max
    
    # Plot observed species abundance distribution
    plot(sort(meanmax), xlab = "Species number", ylab = "Mean maximum count")
    
    # Spatio-temporal patterns in counts (mean over sites)
    tmp <- apply(Yc, c(2,3), mean, na.rm = TRUE)
    matplot(log10(tmp+0.1), type = "l", lty = 1, lwd = 3, xlab = "MHB survey 1 - 3", ylab = "log10 of mean count over sites", frame = F, cex.lab = 1.3, cex.axis = 1.3)
    
    # Drop data from 13 species not observed in 2014
    toss.out <- which(obs.max.C == 0)   # list of species not seen
    Yc <- Yc[,,-toss.out]               # toss them out
    obs.max.C <- obs.max.C[-toss.out]
    ( nspec <- dim(Yc)[3] )             # Redefine nspec as 145
    
    # So here are our data
    str(Yc)
    plot(table(Yc))   # Extremely skewed distribution of observed counts
    
    
    # Bundle and summarize data set
    str(win.data <- list(Yc = Yc, nsite = dim(Yc)[1], nrep = dim(Yc)[2],
                         nspec = dim(Yc)[3], ele = ele, forest = forest, DAT = DAT, DUR = DUR))
    
    
    # Specify model in BUGS language
    sink("model11.txt")
    cat("
        model {
        
        # Community priors (with hyperparameters) for species-specific parameters
        for(k in 1:nspec){
        phi[k] ~ dunif(0,1)                              # Zero-inflation
        alpha0[k] ~ dnorm(mu.alpha0, tau.alpha0)         # Detection intercepts
        beta0[k] ~ dnorm(mu.beta0, tau.beta0)            # Abundance intercepts
        for(v in 1:3){
        alpha[k, v] ~ dnorm(mu.alpha[v], tau.alpha[v]) # Slopes detection
        beta[k, v] ~ dnorm(mu.beta[v], tau.beta[v])    # Slopes abundance
        }
        }
        
        # Hyperpriors for community hyperparameters
        # abundance model
        mu.beta0 ~ dunif(-1, 2)
        tau.beta0 <- pow(sd.beta0, -2)
        sd.beta0 ~ dunif(0, 3)
        for(v in 1:3){
        mu.beta[v] ~ dunif(-1.5, 1)
        tau.beta[v] <- pow(sd.beta[v], -2)
        }
        sd.beta[1] ~ dunif(0, 3)
        sd.beta[2] ~ dunif(0, 1.5)
        sd.beta[3] ~ dunif(0, 1)
        
        # detection model
        mu.alpha0 ~ dunif(-2, 0)
        tau.alpha0 <- pow(sd.alpha0, -2)
        sd.alpha0 ~ dunif(0, 2)
        for(v in 1:3){
        mu.alpha[v] ~ dunif(-0.5, 0.5)
        tau.alpha[v] <- pow(sd.alpha[v], -2)
        }
        sd.alpha[1] ~ dunif(0, 0.8)
        sd.alpha[2] ~ dunif(0, 0.5)
        sd.alpha[3] ~ dunif(0, 0.3)
        
        # Ecological model for true abundance (process model)
        for(k in 1:nspec){
        for (i in 1:nsite){
        a[i,k] ~ dbern(phi[k])   # zero-inflation
        N[i,k] ~ dpois(a[i,k] * lambda[i,k])
        log(lambda[i,k]) <- beta0[k] + beta[k,1] * ele[i] +
        beta[k,2] * pow(ele[i],2) + beta[k,3] * forest[i]
        # Compute presence/absence matrix z (for N > 0) from latent abundance
        z[i,k] <- step(N[i,k]-1)  # returns TRUE if N >= 0
        }
        }
        
        # Observation model for replicated counts
        for(k in 1:nspec){
        for (i in 1:nsite){
        for (j in 1:nrep){
        Yc[i,j,k] ~ dbin(p[i,j,k], N[i,k])
        logit(p[i,j,k]) <- alpha0[k] + alpha[k,1] * DAT[i,j] +
        alpha[k,2] * pow(DAT[i,j],2) + alpha[k,3] * DUR[i,j]
        }
        }
        }
        
        # Other derived quantities
        for(k in 1:nspec){
        mlambda[k] <- phi[k] * exp(beta0[k]) # Expected abundance on natural scale
        logit(mp[k]) <- alpha0[k]     # Mean detection on natural scale
        Nocc.fs[k] <- sum(z[,k])      # Number of occupied sites among the 267
        }
        for (i in 1:nsite) {
        Nsite[i] <- sum(z[i,])        # Number of occurring species at each site
        }
        }
        ",fill = TRUE)
    sink()
    
    # Initial values
    ast <- matrix(rep(1, nspec*nsite), nrow = nsite)
    some.more <- 5          # May have to play with this until JAGS is happy
    Nst <- apply(Yc, c(1,3), max, na.rm = T) + some.more
    Nst[Nst == '-Inf'] <- 20          # May have to play with this, too
    Nst <- Nst
    inits <- function()list(a = ast, N = Nst)
    
    # OR: use inits at earlier solutions (greatly speeds up convergence)
    pm <- out11$mean     # Pull out posterior means from earlier run
    inits <- function() list(a = ast, N = Nst, alpha0 = rnorm(nspec), beta0 = rnorm(nspec), alpha = matrix(rnorm(n = nspec*3), ncol = 3), beta = matrix(rnorm(n = nspec*3), ncol = 3), mu.beta0 = pm$mu.beta0, sd.beta0 = pm$sd.beta0, mu.beta = pm$mu.beta, sd.beta = pm$sd.beta, mu.alpha0 = pm$mu.alpha0, sd.alpha0 = pm$sd.alpha0, mu.alpha = pm$mu.alpha, sd.alpha = pm$sd.alpha )
    
    # Parameters monitored
    params <- c("phi", "mp", "mlambda", "alpha0", "beta0", "alpha", "beta", "mu.beta0", "sd.beta0", "mu.beta", "sd.beta", "mu.alpha0", "sd.alpha0", "mu.alpha", "sd.alpha", "Nsite")
    
    # MCMC settings
    ni <- 60000   ;   nt <- 30   ;   nb <- 30000   ;   nc <- 3
    
    # Call JAGS from R (BRT XXX min), check convergence and summarize posteriors
    out11 <- jags(win.data, inits, params, "model11.txt", n.chains = nc,
                  n.thin = nt, n.iter = ni, n.burnin = nb, parallel = FALSE)
    par(mfrow = c(3,3))   ;    traceplot(out11, c("mu.beta0", "sd.beta0", "mu.beta", "sd.beta", "mu.alpha0", "sd.alpha0", "mu.alpha", "sd.alpha") )
    print(out11, 2)
    
    
    summary(p.sample <- plogis(rnorm(10^6, mean = -1.170, sd = 0.980)) )
    hist(p.sample, breaks = 50, col = "grey", xlab = "Per-individual detection probability", freq = FALSE)
    
    
    # Predict detection for date and duration and occupancy for elevation and forest
    # for each of the 145 observed species
    predI <- array(NA, dim = c(500, nspec, 4))   # covariate value x species x response, "I" for 'individual' (as opposed to 'species' in model 10)
    pm <- out11$mean            # Grab posterior means from model 11
    for(i in 1:nspec){          # Loop over 145 observed species
      predI[,i,1] <- plogis(pm$alpha0[i] + pm$alpha[i,1] * dat.pred +
                              pm$alpha[i,2] * dat.pred^2 )     # p ~ date
      predI[,i,2] <- plogis(pm$alpha0[i] + pm$alpha[i,3] * dur.pred) # p ~ duration
      predI[,i,3] <- pm$phi[i] * exp(pm$beta0[i] + pm$beta[i,1] * ele.pred +
                                       pm$beta[i,2] * ele.pred^2 )     # psi ~ elevation
      predI[,i,4] <- pm$phi[i] * exp(pm$beta0[i] + pm$beta[i,3] * for.pred) # psi ~ forest
    }
    
    # Plots for detection probability and survey date and duration (Fig. 11-29)
    par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
    plot(o.dat, predI[,1,1], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Survey date (1 = 1 April)",
         ylab = "Per-individual detection probability")
    for(i in 2:145){
      lines(o.dat, predI[,i,1], col = i, lwd = 3)
    }
    
    plot(o.dur, predI[,1,2], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 1), xlab = "Survey duration (min)",
         ylab = "Per-individual detection probability")
    for(i in 2:145){
      lines(o.dur, predI[,i,2], col = i, lwd = 3)
    }
    
    
    # Plots for expected abundance and elevation and forest cover (Fig. 11-30)
    par(mfrow = c(1,2), cex.lab = 1.3, cex.axis = 1.3)
    plot(o.ele, predI[,1,3], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 60), xlab = "Elevation (m a.s.l.)",ylab = "Expected abundance")
    for(i in 2:145){
      lines(o.ele, predI[,i,3], col = i, lwd = 3)
    }
    
    plot(o.for, predI[,1,4], lwd = 3, type = 'l', lty = 1, frame = F,
         ylim = c(0, 60), xlab = "Forest cover (%)", ylab = "Expected abundance")
    for(i in 2:145){
      lines(o.for, predI[,i,4], col = i, lwd = 3)
    }
    
    
    
    # 11.11 Summary and outlook
    # -----------------------------------------------------
    
    
    
    
    
    # =========================================================================
    #
    # Summary and conclusion
    #
    # =========================================================================
    
    # Objectives of Applied hierarchical modeling in ecology (AHM)
    # -----------------------------------------------------
    
    # What's in the two AHM books?
    # -----------------------------------------------------
    
    # Outlook: AHM volume 2
    # -----------------------------------------------------
    