# 2_Creates a random set of negatives
#############################################################
plot(myProvShape, axes = T)
# Sets the number of points to generate
n = 2000
#Generates random coordinates
x = (runif(n)*(myMask@extent@xmax - myMask@extent@xmin))+myMask@extent@xmin
y = (runif(n)*(myMask@extent@ymax - myMask@extent@ymin))+myMask@extent@ymin
# Keep only points according to different conditions
myTDF = as.data.frame(cbind(x,y))
myTDF$IsLand = extract(myMask,myTDF[,1:2])
#####################
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

# Function to estimate the distance to the nearest point of the other vector
Lib_DistEstimates <- function (myX1Vec, myY1Vec, myX2Vec,myY2Vec){
  MDist = myX1Vec * 0
  for (i in 1:length(myX1Vec)) {
    myDistVec = (((myX1Vec[[i]] - myX2Vec)^2) + ((myY1Vec[[i]] - myY2Vec)^2))^0.5
    MDist[i] = min(myDistVec, na.rm = T)
  }
  return(MDist)
}


# Function to estimate the distance to the nearest point of the same vector
Lib_DistEstimatesToItself <- function (myXVec, myYVec){
  myTemp_DF = as.data.frame(cbind(myXVec, myYVec))
  myTemp_DF$MDist1 = 0
  for (i in 1:nrow(myTemp_DF)) {
    myDistVec = (((myTemp_DF$myXVec[[i]] - myTemp_DF$myXVec)^2) + ((myTemp_DF$myYVec[[i]] - myTemp_DF$myYVec)^2))^0.5
    myDistVec = myDistVec[-i]
    myTemp_DF$MDist1[i] = min(myDistVec)
  }
  return(myTemp_DF$MDist1)
}



koalafull <- read.csv("wartondata\\mydatasighting_cleaned.csv", header = TRUE)
koalaxyT <- koalafull[c(117,127,128)]
koalaxyT0 <- subset(koalaxyT, yearnew==2010)
koalaxyT0 <- koalaxyT0[c(2,3)]
koalaxyT1 <- subset(koalaxyT, yearnew==2009)
koalaxyT1 <- koalaxyT1[c(2,3)]

# add here base year sighting  and new year sighting data.
koalaxyT0$DistToPos = Lib_DistEstimates(koalaxyT0$X, koalaxyT0$Y, koalaxyT1$X, koalaxyT1$Y)
koalaxyT2 = subset(koalaxyT0, koalaxyT0$DistToPos > 1000)
koalaxyT2 <- koalaxyT2[c(1,2)]
# Merge the two datasets (base year and new year)
#############################################################
koalaxyTn = koalaxyT1[,c("X","Y")]
koalaxyTnn = koalaxyT2[,c("X","Y")]
koalaxyTb2 = rbind(koalaxyTn,koalaxyTnn)


############## repeat step one for each year. #####
##############
koalaxyT <- koalafull[c(117,127,128)]
koalaxyTr1 <- subset(koalaxyT, yearnew==2011)
koalaxyTr1<- koalaxyTr1[c(2,3)]


# add here base year sighting  and new year sighting data.
koalaxyTb2$DistToPos = Lib_DistEstimates(koalaxyTb2$X, koalaxyTb2$Y, koalaxyTr1$X, koalaxyTr1$Y)
koalaxyT3 = subset(koalaxyTb2, koalaxyTb2$DistToPos > 1000)
koalaxyT3 <- koalaxyT3[c(1,2)]
# Merge the two datasets (base year and new year)
#############################################################
koalaxyT3 = koalaxyT3[,c("X","Y")]
koalaxyTb2 = koalaxyTb2[,c("X","Y")]
koalaxyTb4 = rbind(koalaxyTb3,koalaxyT3)


############## repeat step one for each year. #####
##############
koalaxyT <- koalafull[c(117,127,128)]
koalaxyTr2 <- subset(koalaxyT, yearnew==2012)
koalaxyTr2<- koalaxyTr2[c(2,3)]


# add here base year sighting  and new year sighting data.
koalaxyTb4$DistToPos = Lib_DistEstimates(koalaxyTb4$X, koalaxyTb4$Y, koalaxyTr2$X, koalaxyTr2$Y)
koalaxyT5 = subset(koalaxyTb4, koalaxyTb4$DistToPos > 1000)
koalaxyT5 <- koalaxyT5[c(1,2)]
# Merge the two datasets (base year and new year)
#############################################################
koalaxyT5 = koalaxyT5[,c("X","Y")]
koalaxyTb4 = koalaxyTb4[,c("X","Y")]
koalaxyTb5 = rbind(koalaxyTb4,koalaxyT5)

coordinates(koalaxyT0) <- ~X+Y
plot(koalaxyT0)




