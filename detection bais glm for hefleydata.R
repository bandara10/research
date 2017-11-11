
#detection data comes from Hefley method. 
# select varibales use for this analysis from hefley data preparation.
#raster data cannot be used for this method.
d <-  Detection.data[c(10,21,34, 51)] # presence, hpop+dis_visitor + distance_tertiaryandlink
s <-dropLayer(myfullstack,c( 1:8, 10:19, 21:32, 34:50))

######standarise variables 

stand.distance_tertiaryandlink=scale.default(d$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
d$distance_tertiaryandlink = stand.distance_tertiaryandlink

stand.dis_visitor=scale.default(d$dis_visitor, center = TRUE, scale = TRUE) #standarise
d$dis_visitor = stand.dis_visitor

Detection.model=glm(presence~ hpop+dis_visitor, family= "binomial", data=d)
unclass(summary(Detection.model))
myPred = predict(s, Detection.model, type = "response")

plot(myPred, xlab = "x", ylab= "y", main="detection model")



################# Correct for distance 
d2 <- d # change s.
d2$dis_visitor = min(d$dis_visitor)# standarise this raster instead of df as it needed for prediction.

# scale visitor raster.
ss <- scale(s[[1]])
#fill this with minimum distance.
dis_visitor <- ss*0+-1.98
sd <- dropLayer(s,1)
s <- addLayer(sd,dis_visitor) # this is similar to new data in warton method.
myPred3 = predict(s, Detection.model2, type = "response")

#rank-deficient fit may be misleading:This warning checks if the rank of your data matrix is at least equal
#to the number of parameters you want to fit
plot(myPred3, xlab = "x", ylab= "y",main="detection model2") # bias corrected


#####Step 4: Estimate the probability of detection for each presence-only location.####
#IPP data from myFD1. Select relavant colomns.



p.det=ilogit(predict(Detection.model2,new=myFD1))# chnaged myD to ZTGLM.data length =461. 3 X=vector.boot 
hist(p.det, breaks=70)
p.det <- as.data.frame(p.det)
myFD1$p.det=p.det  # or cbind(IPP.data,p.det)
#now select group is greater than 0.
myFD12 <- subset(myFD1, group > 0)
#IPP model is corected. here group is the count of koala in each grid cell. use distance method if possible. 
IPP.corrected= glm(group ~ twi + tpo + temp + aspect + elev+habit2pc+hpop+lot_density+sbd,
                   weights=1/p.det, family="poisson", data=myFD1)


