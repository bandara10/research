fileNameRoot="Jags-ExampleScript" # For output file names.

# Load the data:
myData = read.csv("z15N50.csv") # Read data file; must be in curr. work. dir.
y = myData$y        # The y values are in the column named y.
Ntotal = length(y)  # Compute the total number of flips.
dataList = list(    # Put the information into a list.list tells JAGS the names and values of variables that define the dat
  y = y ,
  Ntotal = Ntotal 
)

# Define the model:
modelString = "
model {
for ( i in 1:Ntotal ) {
y[i] ~ dbern( theta )
}
theta ~ dbeta( 1 , 1 )
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Initialize the chains based on MLE of data.
# Option: Use single initial value for all chains:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: Use function that generates random values for each chain:
initsList = function() {
  resampledY = sample( y , replace=TRUE )
  thetaInit = sum(resampledY)/length(resampledY)
  thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
  return( list( theta=thetaInit ) )
}

# Run the chains:
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=5 , n.adapt=500 )
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , variable.names=c("theta") ,
                            n.iter=3334 )
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

# Examine the chains:
# Convergence diagnostics:
source("DBDA2E-utilities.R")
diagMCMC( codaObject=codaSamples , parName="theta" )
saveGraph( file=paste0(fileNameRoot,"ThetaDiag") , type="eps" )
# Posterior descriptives:
openGraph(height=3,width=4)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) )
saveGraph( file=paste0(fileNameRoot,"ThetaPost") , type="eps" )
# Re-plot with different annotations:
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )
saveGraph( file=paste0(fileNameRoot,"ThetaPost2") , type="eps" )
