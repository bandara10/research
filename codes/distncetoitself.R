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