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