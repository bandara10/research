#1_Load Data
# Sets the general path
#############################################################

setwd("C:/Users/compaq/Dropbox/Lubies-FaoNepal/Data")



# Loads data
#############################################################
#myD = read.dbf("case_series.dbf", as.is=TRUE)
#myD = read.dbf("case_seriesHPAI_Oct.dbf", as.is=TRUE)
#myD = read.dbf("alllab.dbf", as.is=TRUE)
#myD = read.dbf("Hninelabdata.dbf", as.is=TRUE)
#myD = read.dbf("Hfive.dbf", as.is=TRUE)
myD = read.dbf("H5N1_novem.dbf", as.is=TRUE)

#myD = subset(myD, TYPE == "OUT")

# Reproject the data
#############################################################

#H5N1_novem
myND = myD[,c("Longitude","Latitude","Year","Start_Date")]
names(myND) = c("x_ll","y_ll","Year","Start_Date")

myInPtDF = SpatialPointsDataFrame(myND[c("x_ll","y_ll")], myND)
# Set the projection of the input spatial dataframe
myInProj = CRS("+proj=longlat +ellps=WGS84")
proj4string(myInPtDF) = myInProj
# Verify by plotting the coordinates
plot(myInPtDF, axes = T)
# Set the projection of the output spatial dataframe
myOutProj = CRS("+init=epsg:32645") # UTM zone 45
# Makes the conversion
myOutPtDF = spTransform(myInPtDF, myOutProj)
# Add new projected coordinates to the table
myND$x_utm45 = myOutPtDF@coords[,1]
myND$y_utm45 = myOutPtDF@coords[,2]

# Loads State map  & map positives
#############################################################
myProvShape = readShapePoly("Shapefile/Nepal_DUTM")
# Plot Nepal maps
plot(myProvShape, axes = T)
# Add outbreak locations
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red")



###Loads the file list new: after removing climati cones, production system,distance to indian border,Elevation and forest:
#All variables avaibake in the folder
#(an_mt)Annual mean temperature	Fahrenheit
#(an_pp)Annual precipitation	mm
#(du_dn)	Duck density	Duck density for each district (feature per km2)
#(el_mt)	Elevation	Elevation  from mean sea level (m)
#(fr_pc)	Forest 	Percentage coverage of forest for each district (%)
#(fw_dn)	Fowl density 	fowl density for each district (feature per km2)
#(hh_dn)	Household density	Household density for each districts (feature per km2 )
#(hp_dn)	Population density	Population density for each district (persons per km2)
#(hw_ds)	Highway distance	Distance from east west highway
#(ib_ds)	Indian border	Distance from Indian border
#(in_ds)	Indian border	Distance from centroid of each district to Indian border (km)
#(ir_pc)	Irrigatable land	Percentage coverage of irrigatable land for each district (%)
#(iw_pc)	Inland water	Percentage coverage of inland water for each district (%)
#(ld_dn)	Layer duck	Layer duck density for each district (feature per km2)
#(lh_dn)	Layer hen density	Layer hen density for each district (feature per km2)
#(mb_ds)	Migratory birds	Distance from migratory  bird habitat areas
#(mu_ds)	Municipalities	Distance from municipalities
#(pd_pc)	Paddy land	Percentage coverage of paddy land for each district (%)
#(po_dn)	Population density	Not aggregated at district level
#(ps_cb)	Production system	livestock production systems(FAO)
#(rd _ds)	Distance to roads	Distance  from roads centre  (m)
#(rd_dn)	Road density	Density roads for each districts(persons per km2)
#(rn_fl)	Rain fall 	Annual rain fall for each climatic zone (mm)
#(rv_ds)	Rivers	Distance to rivers (m)
#(wn_sp)	Wind speed	km/h
#(ndvi)	Normalise diff. veg. index	
#(tm_jn)	Temperature	Fahrenheit
#(tm_fb)	Temperature	Fahrenheit
#(tm_ma)	Temperature	Fahrenheit
#(tm_ap)	Temperature	Fahrenheit
#(tm_my)	Temperature	Fahrenheit
#(tm_jn)	Temperature	Fahrenheit
#(tm_ju)	Temperature	Fahrenheit
#(tm_au)	Temperature	Fahrenheit
#(tm_sp)	Temperature	Fahrenheit
#(tm_oc)	Temperature	Fahrenheit
#(tm_no)	Temperature	Fahrenheit
#(tm_de)	Temperature	Fahrenheit


#myPredictors = c(
  # Poultry, layer duck density, fowl density, layer hen density
  "ld_dn", "fw_dn", "lh_dn",
  
  # Anthropogenic, # Topography
    # rd_dn (road density), mb_ds(disance from migratory bird areas), mu_ds(distance from municipalities), hw_ds (distice from highways) , iw_pc (inland water percentage)
  "rd_dn", "mb_ds", "mu_ds", "hw_ds", "iw_pc")


myPredictors = c(
  # Poultry
  "fw_dn","du_dn","lh_dn",
  # Anthropogenic
  "rd_dn","ir_pc",
  # Topography
  "hp_dn",
  # Others
  "an_mt","an_pp","hh_dn", "ld_dn","iw_pc", "pd_pc", "fr_pc", "el_mt", "rn_fl", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds", "po_dn", "ps_cb", "ib_ds" )




# add extensions
myPredictorsASC = paste(myPredictors,".asc", sep = "")
myPredList = as.list(paste(myPredictors,".asc", sep = ""))                


myStack = stack(myPredList)
# Plot the predictors
plot(myStack)
# Creates a mask layer
myMask = raster("mask.asc") >= 0
plot(myMask)
plot(myProvShape, axes = T, add = T)

# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

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
myTDF$DistToPos = Lib_DistEstimates(myTDF$x, myTDF$y, myND$x_utm45, myND$y_utm45)
myTDF = subset(myTDF, IsLand == 1 & myTDF$DistToPos > 3000)

# Merge the two datasets (positives and negatives
#############################################################
myFD1 = myND[,c("x_utm45","y_utm45")]
names(myFD1) = c("x","y")
myFD1$HasHPAI = 1
myFD2 = myTDF[,1:2]
names(myFD2) = c("x","y")
myFD2$HasHPAI = 0
myFD = rbind(myFD1,myFD2)


# Remove potential duplicates falling in same pixel
#############################################################
myFD$DistToItself = Lib_DistEstimatesToItself(myFD$x,myFD$y)
myFD = subset(myFD, DistToItself > 1000)
myFD = myFD[,-4]
# displays the output
plot(myProvShape, axes = T)
points(myFD$x, myFD$y, pch = 15, col=c("blue","red")[myFD$HasHPAI+1])

# Extract all predictors from raster stack
############################################################# 
myFD = cbind(myFD,extract(myStack,myFD[,1:2]))
str(myFD)
# removes lines with no data

myFD = na.omit(myFD)

#write.csv(myFD, file="myFDb.csv", row.names=FALSE) 

#see corelation of varibales# oir else go to line=199#
#library(corrgram)
#corrgram(myFD[,c("ld_dn",  "fw_dn",  "lh_dn",	"rd_dn",	"mb_ds",	"mu_ds",	"hw_ds",	"iw_pc")], order=TRUE, lower.panel=panel.ellipse,
         upper.panel=panel.pts, text.panel=panel.txt,
         diag.panel=panel.minmax, 
         main="myFDin PC2/PC1 Order")
library(corrgram)
#corrgram(myFD[,c("ld_dn",  "fw_dn",  "lh_dn",  "rd_dn",	"mb_ds",	"mu_ds",	"hw_ds",	"iw_pc")], order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie, text.panel=panel.txt,
         main="Car Milage Data in PC2/PC1 Order")

#corrgram(myFD[,c("fw_dn","du_dn","lh_dn","rd_dn","ir_pc","hp_dn","an_mt","an_pp","hh_dn", "ld_dn","iw_pc", "pd_pc", "fr_pc", "el_mt", "rn_fl", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds", "po_dn", "ps_cb", "ib_ds")], order=NULL, lower.panel=panel.shade,
         upper.panel=NULL, text.panel=panel.txt,
         main="Predictors")

#L<-cor(myFD[,c("fw_dn","du_dn","lh_dn","rd_dn","ir_pc","hp_dn","an_mt","an_pp","hh_dn", "ld_dn","iw_pc", "pd_pc", "fr_pc", "el_mt", "rn_fl", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds", "po_dn", "ps_cb", "ib_ds")], use="complete.obs", method="kendall") 

#capture.output(L,file="test.xls")

fit <- princomp(myFD, cor=TRUE)
summary(fit)
loadings(fit)
plot(fit,type="lines")
fit$scores
biplot(fit)

library(psych)
ct <- corr.test(myFD) 


#############################################################
nBoot = 10
SensitivityTest = 0

# Runs a standard Logistic regression
#############################################################
# Build the GLM object
myStr = "HasHPAI ~"
for (i in 1:length(myPredictors)){myStr = paste(myStr, "+", myPredictors[i])}  
print(myStr)
print("##########################################")
# Runs the GLM  
myGLM = glm(as.formula(myStr), data = myFD, family = "binomial")
#myGLM = glm(HasHPAI ~ brdn+ dsdn+dudn+hpdn+rddn+access+urban_r+dem+ncrop, data = myFD, family = "binomial")
# Show summary data
summary(myGLM)
# Stores the residuals
myFD$res = residuals(myGLM)
# Create a spatial correlogram of the residuals
myResCorr <- correlog(myFD$x, myFD$y, myFD$res,na.rm=T, increment=1000, resamp=0, latlon = F)
#mySigVec = ifelse(myResCorr$p<0.05,1,0)  
#plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", col = mySigVec[1:20]+1, pch=16, lwd=1.5, cex = 1.2,
#xlab="distance", ylab="Moran's I")
plot(myResCorr$mean.of.class[1:20], myResCorr$correlation[1:20] ,type="b", pch=16, lwd=1.5, cex = 1.2,
     xlab="distance", ylab="Moran's I")
abline(h=0)

# Crase approach to account for SA
#############################################################
# Build the standard GLM object
myFD$AR1 = myLib.AutoRegressiveMean(residuals(myGLM), myFD$x, myFD$y, 18000)
# Changes the formular string
myStrAR = paste(myStr,"+ AR1")
# Re-run the model  
myGLMR1 = glm(as.formula(myStrAR), data = myFD, family = "binomial")
summary(myGLMR1)
plot(myGLMR1)  



# Estimate correlogram of new residuals

myFD$R1 = residuals(myGLMR1)
Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=1000, resamp=0, latlon = F)              
#Corr <- correlog(myFD$x, myFD$y, myFD$R1,na.rm=T, increment=10000, resamp=0, latlon = F)
#mySigVec = ifelse(Corr$p<0.01,1,0)  

lines(Corr$mean.of.class[1:20], Corr$correlation[1:20], col = "red")
#points(Corr$mean.of.class[1:20], Corr$correlation[1:20], pch = 16, col = mySigVec[1:20]+1)
# Plot the correlogram
plot(Corr)              

# Estimates DAIC for each predictor
#############################################################
myChi2 = rep(0,length(myPredictors))
myPValue = rep(0,length(myPredictors))

for (i in 1:length(myPredictors)){
  myShortPred = myPredictors[-i]
  myStrN = "HasHPAI ~"
  for (j in 1:length(myShortPred)){myStrN = paste(myStrN, "+", myShortPred[j])}
  myStrARN = paste(myStrN,"+ AR1")
  myGLMR2 = glm(as.formula(myStrARN), data = myFD, family = "binomial")
  mylrtest = lrtest(myGLMR1,myGLMR2)
  myChi2[i] = mylrtest[[4]]
  myPValue[i] = mylrtest[[6]]
}
myLRTable = data.frame(myChi2,myPValue)
myLRTable = cbind(myPredictors,myLRTable)
myLRTable

# ROC & predictions
#############################################################
# Map predictions
AR1 = myMask * 0
AR1@file@name = "AR1"
# Add it to the stack
myARStack = addLayer(myStack, AR1)
names(myARStack)[nlayers(myARStack)]="AR1"

myPred = predict(myARStack, myGLMR1, type = "response")
#  myPred = predict(myStack, myGLM, type = "response")

# plot the prediction
plot(myPred, xlab = "x-coord (UTM45)", ylab= "y-coord (UTM45)",main="Predictive risk map for H5N1")
plot(myProvShape, axes = T, add = T)
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red", cex = 0.6)

plot(myPred)
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red", cex = 0.6)
#Subtitle in two rows
#sub="  Map 2. Predicitve risk map for HPAI in Nepal. Reported cases \nused for the model are from 2009-2013. Data source:WAHID ",col.lab="blue", cex.lab=0.75,col.sub="blue")
# ROC curve & GOF metrics
myPred = prediction(predict(myGLM, type = "response"), myFD$HasHPAI)
perf <- performance(myPred, measure = "tpr", x.measure = "fpr")
plot(perf, colorize = T)
myPredPerfs = myLib.EstimateLogisticGOF(predict(myGLM, type = "response"), myFD$HasHPAI, 0.5)
myPredPerfs

#test = raster("C:\Dropbox\5_Supervision\Lubies-FaoIndia\Outputs\SatscanElliptical\

######################average false positive  rate#############################
pred <- prediction(predict(myGLM, type = "response"), myFD$HasHPAI)
perf <- performance(pred, "tpr", "fpr")

plot(perf, avg="threshold",
     spread.estimate="boxplot")

#3_BRT_Analysis
###################
# Creates the BRT object
nVar = length(myPredictors)

myBRT <- gbm.step(myFD,gbm.x = 4:(nVar+3),gbm.y = 3,family = "bernoulli",tree.complexity = 5,learning.rate = 0.004,bag.fraction = 0.75)
summary(myBRT)
gbm.plot(myBRT)
gbm.plot(myBRT,n.plots=21,write.title = FALSE)
myPDF = as.data.frame(values(myMask))

# Creates a dataframe with the predictors
for (i in 1:length(myPredictors)){
  myVec = values(raster(paste(myPredictors[i],".asc",sep="")))
  myPDF = cbind(myPDF,myVec)
}
names(myPDF) = c("Mask",myPredictors)
# Get the predicted values
myPredValues = predict.gbm(myBRT, n.trees = myBRT$gbm.call$best.trees,tree.complexity = 5, learning.rate = 0.005,newdata = myPDF, type = "response")
# Convert the predicted values to a raster
myPredR = setValues(myMask,myPredValues)*myMask
plot(myPredR,main="Predictive risk map for H5N1")
plot(myProvShape, axes = T, add = T)
points(myND$x_utm45, myND$y_utm45, pch = 15, col = "red", cex = 0.6)
writeRaster(myPredR,"BRTPredictionfao_18_06_14.asc")
gbm.plot.fits(myBRT)
find.int <- gbm.interactions(myBRT)
find.int$interactions
find.int$rank.list
#rainfall_paddy
gbm.perspec(myBRT,17, 15, y.range=c(9,15), z.range=c(0,0.3))

