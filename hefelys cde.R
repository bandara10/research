#Appendix S1. Annotated R code to implement marked inhomogeneous Poisson point process 
#species distribution model and two-phase nonparametric bootstrap algorithm.

###############################################################################
#	Ensure you are using R Version 3.0.1 or newer
#	Required Packages: install packages 'VGAM', 'mosaic' and 'spatstat'
###############################################################################

library(VGAM)
library(mosaic)
library(spatstat)

###############################################################################
#	Generate single simulated data example 
#
#     Note: the simulated data example is the same as the large sample size used  
#	in our simulation study. See manuscript for details.
###############################################################################

set.seed(1234)
x=rnorm(10^6,0,1)

eta.pres=8.5+1*x
mat=matrix(exp(eta.pres), nrow=1000, ncol=1000)
lam=im(mat, xrange=c(0,1), yrange=c(0,1))
points=rpoispp(lam) # complete spatial ramdomness.  homogenous or inhomogenous.

x.ipp=(log(lam[points])-8.5)/1   #num [1:8156] 1.312 1.107 1.316 0.617 2.492
x.int=rnorm(1000,0,1)	  #montecarlo points#        #num [1:1000]   -0.0786 -1.2963 -0.2973 0.2261 1.2265

y.IWLR=c(rep(1,length(x.ipp)),rep(0,length(x.int)))#num [1:9156] 0 1
IPP.data=data.frame(y=y.IWLR,x=c(x.ipp,x.int))   #9156 obs. of  2 variables: y , x.

eta.group=1+0.5*x.ipp  #num [1:8156] 1.66 1.55 1.66 1.31 2.25
y.ZTGLM=rpospois(length(eta.group),exp(eta.group)) #num [1:8156] 1 8 3 2 8 2 4 3 6 3
ZTGLM.data=data.frame(group.size=y.ZTGLM,x=x.ipp) #8156 obs. of  2 variables:
#                                                 $ group.size: num  1 8 3 2 8 2 4 3 6 3 ...
#                                                 $ x         : num  1.312 1.107 1.316 0.617 2.492 

eta.det=-2+-1*x.ipp+0.5*scale(y.ZTGLM)  #num [1:8156, 1] -3.92 -2.7 -3.63 -3.08 -4.09 ...
                                        #- attr(*, "scaled:center")= num 5.2
                                        #- attr(*, "scaled:scale")= num 3.45
detected=rbinom(length(eta.det),1,ilogit(eta.det)) #int [1:8156] 0 0 0 0 0 0 0 0 0 0
keep=c(1:round(0.20*length(x.ipp)))  #20% selected           #int [1:1631] 1 2 3 4 5 6 7 8 9 10

Detection.data=data.frame(y=detected[keep],x=x.ipp[keep],group.size=y.ZTGLM[keep]) # 20% data #1631 obs. of  3 variables: y,x, group.size

IPP.data=IPP.data[c(which(detected==1),which(y.IWLR==0)),]  ##1499 obs. of  2 variables y, x 

ZTGLM.data=ZTGLM.data[which(detected==1),]      ##499 obs. of  2 variables: group.size, x.

###############################################################################
#	Examine the example data 
#
#     'IPP.data': Column "y" is the response variable used to estimate the inhomogeneous Poisson point process model
#     parameters using infinitely weighted logistic regression (Fithian & Hastie 2013). A value of 1 corresponds to a  
#     presence-only location and a value of 0 corresponds to a Monte Carlo integration point. Column "x" is the habitat covariate.
#
#   	'ZTGLM.data':Column "group.size" is the size of the group at each presence-only location. 
#	Column "x" is the habitat covariate.
# 
#	'Detection.data': Column "y" is the response variable. A value of 1 corresponds to detection and a value of 0 corresponds
#	to non-detection. Column "x" is the detection covariate and "group.size" is the group size of each detection/non-detection.
###############################################################################

names(IPP.data)
names(ZTGLM.data)
names(Detection.data)


###############################################################################
#	Estimate the coefficient for the habitat covariate "x" for the inhomogeneous Poisson point process model
#	ignoring non-detection sampling bias using infinitely weighted logistic regression (Fithian & Hastie 2013)
#
#     Note: Because of non-detection sampling bias the coefficient estimate (0.29911) is not close to the true value of 1.
###############################################################################

IPP.ignored=glm(y~x,family="binomial",weights=10000^(1-y),data=IPP.data)
summary(IPP.ignored)
library(latex2exp)
xtable(IPP.ignored, auto = TRUE)
###############################################################################
#	Estimate the coefficient for the habitat covariate "x" for the zero truncated Poisson generalized linear model 
#	ignoring non-detection sampling bias.
#
#	Note: Because of non-detection sampling bias the estimated intercept (1.15756) is biased. The true value is 1.
#	Non-detection sampling bias does not affect the estimated covariate "x" (0.51988) by a large amount. The true value is 0.5.
###############################################################################

ZTGLM.ignored=vglm(group.size~x,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.ignored)


################################################################################
#	Correcting for non-detection sampling bias. 
#
#	Note: Steps numbered 1-7 are in reference to the two-phase nonparametric bootstrap algorithm in the manuscript.
#	Bootstrap sampling (i.e., sampling with replacement) and steps 1, 3, and 7 are not needed to get
#	corrected coefficient estimates and are only needed to calculate standard errors and confidence intervals. 
#	We will skip bootstrap sampling and steps 1, 3 and 7 until later on.
#
#	Note: True values of the intercept, "x" and "group.size" used to generate detection data are -2, -1 and 0.5 respectively. 
#	
#	Note: Because we corrected for non-detection sampling bias the estimate of the 
#     coefficient "x" (0.94315) in step 5 is very close to the true value of 1. 
#	The standard error (0.01185), however, should not be trusted as it is estimated without 
#	taking into account the variablity in the estimated probability of detection("p.det").
#	Similarly the estimated intercept (0.98370) and coefficient "x" (0.53195) in step 6 is 
#	close to the true value of 1 and 0.5 respectively. The standard errors for both estimates should
#	not be trusted.
###############################################################################

##Step 2 - Fit an appropriate model to the detection data set
Detection.model=glm(y~x+scale(group.size),family="binomial",data=Detection.data)
summary(Detection.model)
Detecion <- xtable(Detection.model, auto=TRUE)
autoformat(Detecion)
##Step 4 - Estimate the probability of detection for each presence-only location.
p.det=ilogit(predict(Detection.model,new=ZTGLM.data))
IPP.data$p.det=c(p.det,rep(1,length(x.int)))
ZTGLM.data$p.det=p.det


##Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det
IPP.corrected=glm(y~x,family="binomial",weights=(1/p.det)*10000^(1-y),data=IPP.data)
summary(IPP.corrected)

ipp <- xtable(IPP.corrected)
autoformat(ipp)
##Step 6 - Fit an zero-truncated Poisson generalized linear model that weights the log-likelihood by 1/p.det
ZTGLM.corrected=vglm(group.size~x,weights=1/p.det,family="pospoisson",data=ZTGLM.data)
summary(ZTGLM.corrected)

###############################################################################
#	Implementing the two-phase nonparametric bootstrap algorithm (i.e., steps 1-7)
#
#     'set.seed' establishes a random starting seed so simulation results are the same each time the bootstrap code is run.
#     If the set.seed code is ignored, the mean, standard deviation and confidence intervals will vary slightly each time 
#	the code is run.
#
#     The function "tpnbs" is the two-phase nonparametric bootstrap and implements steps 1-7. 
#
#     'n.bootstraps' is the number of bootstrap samples to be included in the analysis. We used 1,000 in our study.
#
#     The 'do()' function in package mosaic implements the bootstrap.
###############################################################################
set.seed(1234)

tpnbs=function()	{
  bss=resample(1:dim(ZTGLM.data)[1])
  IPP.data.bss=rbind(IPP.data[bss,1:2],IPP.data[which(IPP.data$y==0),1:2])
  ZTGLM.data.bss=ZTGLM.data[bss,]
  Detection.data.bss=resample(Detection.data)
  Detection.model=glm(y~x+scale(group.size),family="binomial",data=Detection.data.bss)
  p.det.bss=ilogit(predict(Detection.model,new=ZTGLM.data.bss))
  IPP.data.bss$p.det=c(p.det.bss,rep(1,length(x.int)))
  ZTGLM.data$p.det=p.det.bss
  options(warn=-1)
  IPP.model=glm(y~x,family="binomial",weights=(1/p.det)*10000^(1-y),data=IPP.data.bss)
  ZTGLM.corrected=vglm(group.size~x,weights=1/p.det,family="pospoisson",data=ZTGLM.data.bss)
  c(coef(IPP.model)[2],coef(ZTGLM.corrected))
}
n.bootstraps=1000
bootstrap.sample = mosaic::do(n.bootstraps) *tpnbs()
str(bootstrap.sample)
bootstrap.sample=data.matrix(bootstrap.sample)
###############################################################################
#	Calculating mean, standard deviation and 95%, equal-tailed confidence intervals
#	from the empirical distribution. See "Introduction to the Bootstrap" (Efron & Tibshirani 1994) 
#	for more details.
#
#     Estimates of the coefficient for the covariate "x" (0.9520575) of the inhomogeneous Poisson point 
#	process model and the intercept(0.9832485) and covariate "x" (0.531719) of the zero-truncated 
#	Poisson generalized linear model are close to the true values (1, 1, 0.5).   
#
#	Standard errors when detection bias is accounted for(0.1141398; 0.03696048; 0.0224879), and thus 95%  
#	confidence intervals, and standard errors (0.01185; 0.0099313; 0.0060487) are larger when the variability  
#	in the probability of detection is not accounted for. Note the estimated standard errors(0.01185; 0.0099313; 0.0060487)
#	are obtained from steps 5 & 6 above.	
###############################################################################

colMeans(bootstrap.sample)[1]
sd(bootstrap.sample)[1]
data(c(.025, .975),bootstrap.sample[,1])
colMeans(bootstrap.sample)[2]
sd(bootstrap.sample)[2]
data(c(.025, .975),bootstrap.sample[,2])

colMeans(bootstrap.sample)[3]
sd(bootstrap.sample)[3]
data(c(.025, .975),bootstrap.sample[,3])


######histogrames for each colomns. (sample; means, intercept;means, covariate;means,  )
for (i in 1:ncol(bootstrap.sample)) {
  
  hist(bootstrap.sample[,i], breaks=50)
}


# keep objects for reuse as we run the model with same names.
plot(Detection.model) # this is hefley detectiono model
Hefley.Detection.model.plot <- recordPlot(Detection.model) # keep this and delete Detection.model
#     save as a plot object for resue. 
#     as we want to change detection model in other abalysis with same name.
replayPlot(Hefley.Detection.model.plot) #  replay from global env. 
dev.off()
par(mfrow= c(2,2))
plot(IPP.ignored)
Hefley.IPP.ignored.plot <- recordPlot(IPP.ignored)
plot(IPP.corrected)
Hefley.IPP.corrected.plot <- recordPlot(IPP.corrected)
plot(ZTGLM.ignored)
plot(ZTGLM.corrected)
Hefley.ZTGLM.plot <- recordPlot(ZTGLM.corrected)


# below are plots for models
replayPlot(Hefley.Detection.model.plot)
replayPlot(Hefley.IPP.ignored.plot)
replayPlot(Hefley.IPP.corrected.plot)
replayPlot(Hefley.ZTGLM.plot)

a <- 5
b <- 3
substitute(expression(a * e^(b * x)), list(a=a, b=b))
