
# Postive values so, can apply any transformation
for( i in c("norm", "rstd", "std", "stretch", "nl", "slog", "sr")) {
  density( raster.transformation(myfullstack[[43]], trans = i) )
}
#("norm", "rstd", "std", "stretch", "nl", "slog", "sr"
density(r1) # aspect normalization
density(r2)#DP_QLD_FPC20141 #log
density(r3) # waves no trasformation is good
density(r5) #rainfallmeanannual 3log
library(spatialEco)

for( i in c("norm", "rstd", "std", "stretch", "nl", "slog", "sr")) {
  density( raster.transformation(myfullstack[[75]], trans = i) )
}
transformed <- abs(r6 - mean(r6)) # for bomodel data
r7 <- myfullstack[[14]]



#https://stats.stackexchange.com/questions/209241/what-transformation-should-i-use-for-a-bimodal-distribution
# transformed <- abs(binomial - mean(binomial))
# shapiro.test(transformed)
# hist(transformed)
# which produces something close to a slightly censored normal distribution and (depending on your seed)
# 
# Shapiro-Wilk normality test
# 
# data:  transformed
# W = 0.98961, p-value = 0.1564

# it is not a standard transformation, but looking at the mixture distribution,
# it was a simple way in this case of superimposing the two modes to make a
# single mode. It would not have worked so well if there there had been 
# differing numbers of observations from the two original normals or if they 
# had different standard deviations