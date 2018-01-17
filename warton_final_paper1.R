##############

load("C://Users//uqrdissa//ownCloud//Covariates_analysis//Mark_S//Data_raw_koala//mydatasighting_cleaned.RData")
mydata <- mydatasighting_cleaned
mydata <- mydata[c("X","Y","yearnew")]
names(mydata) <- tolower(names(mydata))

# Analyse data for 2011:
all.locations= subset(mydata,yearnew >2000:2015, select=x:y)
plot(all.locations)
# set a minimum distance betweek koalas
source("Lib_DistEstimatesToItself.r")
all.locations$disttoitself = Lib_DistEstimatesToItself(all.locations$x, all.locations$y)
select.locations = subset(all.locations, all.locations$disttoitself > 400)
selected.locations = select.locations[,1:2]

# keep a buffer distance of 2000m as required by this method.
# get  raster dimentions:441829, 541829, 6901098, 7001098  (xmin, xmax, ymin, ymax)

selected.locations2<- subset(selected.locations, x > 443829 & x < 539829)
selected.locations <- subset(selected.locations2, y > 6903098 & y < 6999098) # xy only within the study area.
coordinates(selected.locations) <- ~x+y



#check which distance variables are associated with sightings
myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder
myfullstack = stack(myfullstack.a) 
Distance_primary <- myfullstack[[12]] 
plot(Distance_primary, main="primary roads"); plot(selected.locations, add=TRUE)
plot(Distance_motorway <- myfullstack[[15]],main="Distance to motorway"); plot(selected.locations, add=TRUE)

# this is the presence data set. extract?.
#get the full raster data set.
#projection(myfullstack) <- gsub("units=m", "units=km", projection(myfullstack))
myfullstack.a <- list.files(path="rasters_cropped_renner_methods",pattern="\\.tif$") #Mark_s folder
myfullstack = stack(myfullstack.a) 
extent(myfullstack) <- extent(c(xmin(myfullstack), xmax(myfullstack), ymin(myfullstack), ymax(myfullstack))/1000)
#quadrature points
habitat.r<- subset(myfullstack, c(2,3,12,15,22,23,34,45,46,47,48,49)) # habitat covariates
bigquad <- as.data.frame(habitat.r, xy=TRUE, na.rm=T)

# to predict using model based control of obser bias at minimum distance.
bigquad.2 <- bigquad
bigquad.2$distance_primaryandlink = min(bigquad.2$distance_primaryandlink) 

bigquad.2$distance_motorwayandlink = min(bigquad.2$distance_motorwayandlink)
#stt <- na.omit(stt)
colnames(bigquad)[1] <- 'X'; colnames(bigquad)[2] <- 'Y'
# stt[is.na(stt)] <- 0
xydatan <- bigquad[c(1,2)]
# stt requires xy as integers.
xydata <- as.data.frame(lapply(xydatan, as.integer)) # stt[] <- lapply(stt, as.integer)#this line edited on 04/01# make only xy integer in line with dadta shared with Mark S. 
bigquad <- cbind(xydata, bigquad[c(-1,-2)])

# get species selected data as a dataframe
selected.locations=as.data.frame(selected.locations)/1000

sp.xy = data.frame(selected.locations)
colnames(sp.xy)[1] <- 'X'; colnames(sp.xy)[2] <- 'Y'
sp.xy <- as.data.frame(lapply(sp.xy, as.integer))

ppm.form.e = ~ poly(awc,clay,elev,fpcnew, nitro,sbd,temp_max,temp_min,tpo,twi,degree = 2, raw = TRUE)
scales = c( 0.5, 1, 2, 4, 8)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form.e)

ppmFit.e = ppmlasso(ppm.form.e, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
#Predict and plot
pred.fit.e = predict.ppmlasso(ppmFit.e, newdata=bigquad)

predictions.fit.e <- cbind(xydata, pred.fit.e) # xydatan was chnaged to xydata.
pred.final0.e<- rasterFromXYZ(as.data.frame(predictions.fit.e )[, c("X", "Y", "pred.fit.e")])
plot(pred.final0.e, main=" koala density-warton method/ env only")
#### Env and distance both

ppm.form = ~ poly(awc,clay,elev,fpcnew, nitro,sbd,temp_max,temp_min,tpo,twi, degree = 2, raw = TRUE)+ poly(distance_primaryandlink,distance_motorwayandlink, degree = 2, raw = TRUE)

scales = c( 0.5, 1, 2, 4, 8)
findres(scales, coord = c("X", "Y"), sp.xy = sp.xy, env.grid = bigquad, formula = ppm.form)
#4.2 Fitting a regularisation path of point process models
#a LASSO penalty that optimises non-linear GCV
ppmFit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = bigquad, sp.scale = 1, n.fits = 200)
#Predict and plot
pred.fit = predict.ppmlasso(ppmFit, newdata=bigquad)

predictions.fit <- cbind(xydata, pred.fit) # xydatan was chnaged to xydata.
pred.final0<- rasterFromXYZ(as.data.frame(predictions.fit )[, c("X", "Y", "pred.fit")])
plot(pred.final0, main=" koala density-WM env & dist/ bias not corrected")

# now correct for bias.
pred.fit.correct = predict.ppmlasso(ppmFit, newdata=bigquad.2)

predictions.fit.correct <- cbind(xydata, pred.fit.correct) # xydatan was chnaged to xydata.
pred.final0.correct<- rasterFromXYZ(as.data.frame(predictions.fit.correct )[, c("X", "Y", "pred.fit.correct")])
plot(pred.final0.correct, main=" koala density-warton method/ bias corrected")

### residuals:
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas bias corrected model")
### with lurking varibale plots.
diagnose.ppmlasso(ppmFit)
#K-envelop
kenv = envelope(ppmFit, fun = Kinhom) # simulated envelop for summary function
plot(kenv,main= "Inhomogeneous K-function with 95% simulation envelope")

#A regularisation path of Poisson point process models
quad.1k = sample.quad(bigquad, 1)
ppm.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, criterion = "nlgcv")

#4.3 Block cross-validation
#block cross-validation as a method for choosing the LASSO penalty
#area interaction model with radius 2k and lasso penalty by5-fold cross validation.
final.fit = ppmlasso(ppm.form, sp.xy = sp.xy, env.grid = quad.1k,sp.scale = 1, 
                     criterion = "blockCV", n.blocks = 5, block.size = 10)
#Predict and plot
pred.final.fit = predict.ppmlasso(final.fit, newdata=bigquad.2)

predictions.final.fit <- cbind(xydata, pred.final.fit) # xydatan was chnaged to xydata.
pred.final<- rasterFromXYZ(as.data.frame(predictions.final.fit )[, c("X", "Y", "pred.final.fit")])
plot(pred.final, main=" Koala density-warton method bias corrected blockCV")

#diagnostic plots:residuals and lurking variable plots.
diagnose.ppmlasso(final.fit)
diagnose.ppmlasso(final.fit, which = "x", type = "Pearson", compute.sd = TRUE)
diagnose.ppmlasso(final.fit, which = "y", type = "Pearson", compute.sd = TRUE)
resid.plot = diagnose(final.fit, which = "smooth", type = "Pearson", main="smoothed pesrson residulas biascorrected interaction model")
###### runs this code without error.

### another way to fit area interaction model. can add: criterion = "blockCV", n.blocks = 5, block.size = 10 
final.fita = ppmlasso(final.fit,sp.xy, env.grid = bigquad,criterion = "blockCV", n.blocks = 5,
                      block.size = 10,sp.scale = 1, n.fits = 100, family = "area.inter", r = .5)
diagnose(final.fita, which = "smooth", type = "Pearson")
pred.interaction = predict(final.fita, newdata=bigquad.2)
pred.inter.action <- cbind(xydata, pred.interaction)
pred.ct.inter <- rasterFromXYZ(as.data.frame(pred.inter.action)[, c("X", "Y", "pred.interaction")])
plot(pred.ct.inter )
diagnose.ppmlasso(final.fita)

