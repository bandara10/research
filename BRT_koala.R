#3_BRT_Analysis
###################
library(dismo)
library(gbm)

myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack = stack(myfullstack.a)
#### Step 2: import koala .csv data from the full study land region. Full square study are consists of land and sea. ####

hefleydata <- read.csv("hefley_fishnet_rastermatch2011.csv", header = TRUE) # centroids for full study area as some sros are required.
names(hefleydata)
# plot koala presence locations
hefleydata.presence <-subset(hefleydata, presence==1)
coordinates(hefleydata.presence) <- ~x+y
plot(hefleydata.presence, add=TRUE)
# get only coordinates from dataset.
hefleydata.s <- hefleydata[c("x","y")]
hefleydata.s = as.data.frame(hefleydata.s)
#select presence records

#### Step 3:  extract X=vector.boot from raster anad combine wth hefleydata, presence and group varibels.####
myFD = cbind(hefleydata.s,raster::extract(myfullstack,hefleydata.s))
# get presence and group vriables
myFD1 = cbind(myFD,hefleydata [6]) # get presence and group data. 
myFD1 = na.omit(myFD1) # remove all NA valuves

#### Step 4: select presence and abnsences data and combine with. #### 
ZTGLM.myFD1=myFD1[which(myFD1$presence==1),] # select presence data. THis will be reorganised as 0 and 1 later.
ZTGLM.myFD2=myFD1[which(myFD1$presence==0),] # select all absence data 

#####Step 5: select only 1000 absences (monticarlo points as in hefleys method??)####
set.seed(12345)
ZTGLM.myFD3 <- sample(seq_len(nrow(ZTGLM.myFD2)), size = 1000,replace=FALSE)#select only 1000 absences use for ipp.data
ZTGLM.myFD4 <- ZTGLM.myFD2[ZTGLM.myFD3, ] #x.int data frame now
#This is  similar to hefley`s IWLR data set.
ZTGLM.myFD5=rbind(ZTGLM.myFD1,ZTGLM.myFD4) 

# Creates the BRT object
nVar = length(myfullstack)
 # select varibles 
newZTGLM5 <- ZTGLM.myFD5[c(1,2,53,8:23)]
myBRT <- gbm.step(newZTGLM5,gbm.x = 4:19,gbm.y = 3,family = "bernoulli",tree.complexity = 2,learning.rate = 0.01,bag.fraction = 0.75)
plot(myBRT)
par(mgp=c(3,1,0),mar=c(10,12,3,2)+0.1) # set mrgins for the plot
summary(myBRT, las=2, asp = 1)
dev.off()
par(mfrow=c(4,4))
gbm.plot(myBRT, n.plots = 16, write.title = FALSE)
#get predictions and plot full model.
predictions <- predict(myfullstack, myBRT, n.trees=myBRT$gbm.call$best.trees, type="response")
plot(predictions)
#reduce model
set.seed(125)
myBRT.2 <- gbm.step(newZTGLM5, gbm.x = c(7, 8, 16, 18), gbm.y = 3,family = "bernoulli",tree.complexity = 2,learning.rate = 0.01,bag.fraction = 0.75) #Build initial model
gbm.plot(myBRT.2, n.plots = 16, write.title = FALSE)
summary(myBRT.2, las=2, asp = 1)
predictions.2 <- predict(myfullstack, myBRT.2, n.trees=myBRT$gbm.call$best.trees, type="response")
#plot predictions
plot(predictions.2,main="Detection proabilities ")
plot(hefleydata.presence,  cex = 0.3,add=TRUE)
writeRaster(predictions.2,"Detection proabilities_BRT.asc")
gbm.plot.fits(myBRT.2)
find.int <- gbm.interactions(myBRT.2)
find.int$interactions
DD<-find.int$rank.list
#rainfall_paddy
par(mgp=c(10,10,0),mar=c(0,3,3,2)+0.1)
gbm.perspec(myBRT,4, 3)
myBRT$cv.loss.matrix
