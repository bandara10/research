library(ppmlasso)
library(raster)
library(sp)
library(rgdal)
library(spatial.tools)
# Pre-standardise observer bias variables
setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")
myenv <- list.files( path="wartondata", pattern="\\.tif$", full.names = TRUE) 
myenv.stack <- stack(myenv)
#crs(myenv.stack) <- NA
extent(myenv.stack) <- extent(c(xmin(myenv.stack), xmax(myenv.stack), ymin(myenv.stack), ymax(myenv.stack))/1000)
plot(myenv.stack[[1]])
stt <- as.data.frame(myenv.stack, xy=TRUE, na.rm=T)
colnames(stt)[1] <- "X"
colnames(stt)[2] <- "Y"
# stt[is.na(stt)] <- 0
xydatan <- stt[c(1,2)]
stt[] <- lapply(stt, as.integer)


# test resolution on back transformatiotn
# sbd <- rasterFromXYZ(as.data.frame(stt)[, c("X", "Y", "temp")])
# plot(sbd)
# scales = c(0.5, 1, 1.2, 2, 4, 8, 16,32) # scales = sp.scale=1
#stt = sample.quad(env.grid =stt , sp.scale= 1, file = "Quad") # this is quadrature points to be use for the analysis.
# 
# 
# 
# #A matrix containing locations of species presences in the first two columns and the interpolated
# #environmental data in the remaining columns
# 
# 
# species.env = env.var(kolaxyT, env.grid = stt, env.scale = .5,
#                       file.name = "Sp Env Data") # simialr to extract. isn`t it?

# A matrix dat.ppm with columns representing the latitude and longitude of presence locations and
# quadrature points along with the associated environmental data, as well as a column Pres indicating
# whether either point corresponds to a presence location or a quadrature point, and a column wt of
# observation weights

# determines observation weights
# and sets up the design matrix required for fitting a regularisation path
# species.ppm = ppmdat(sp.xy = kolaxyT, back.xy = stt,
#                      sp.scale = 1, file.name = "Sp PPM Data")
# summary(species.ppm)

#########
#########Pre-standardise observer bias variables
stand.distance_tertiaryandlink=scale.default(stt$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
stt$distance_tertiaryandlink = stand.distance_tertiaryandlink

stand.dis_visitor=scale.default(stt$dis_visitor, center = TRUE, scale = TRUE) #standarise
stt$dis_visitor = stand.dis_visitor

stand.dis_city=scale.default(stt$dis_city, center = TRUE, scale = TRUE) #standarise
stt$dis_city = stand.dis_city



#koala data
kolaxyT <- read.csv("wartondata\\koalaxy.csv", header = TRUE) # in km.XY| go to ppmFrom

# kolaxy2 <- subset(kolaxy, X > 442 & X < 540)
# kolaxyT <- subset(kolaxy2, Y > 6902 & Y < 7000) # xy within the area only.



#########
ppmForm1 = ~  poly(temp, elev,uhabit3, degree = 2) + poly(dis_visitor, degree = 2)
ppmFit1 = ppmlasso(ppmForm1, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200)

pred.biasCorrectnot = predict(ppmFit1, newdata=stt)
predictions <- cbind(xydatan, pred.biasCorrectnot)

xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrectnot")])
plot(xy.rr)

# bias corrected
## To predict using model-based control of observer bias (at min value for distance):
newstt <- stt
newstt$stand.dis_visitor = min(stand.dis_visitor) # based on minimum distance

ppmForm2 = ~  poly(temp,elev,uhabit3, degree = 2)+ poly(dis_visitor, degree = 2)
ppmFit2 = ppmlasso(ppmForm2, sp.xy = kolaxyT, env.grid = stt, sp.scale = 1, n.fits = 200)
pred.biasCorrect = predict(ppmFit2, newdata=newstt)
predictions.correct <- cbind(xydatan, pred.biasCorrect)
xy.rr2 <- rasterFromXYZ(as.data.frame(predictions.correct)[, c("X", "Y", "pred.biasCorrect")])
plot(xy.rr2)

# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, coord = c("X", "Y"), sp.xy = kolaxyT, env.grid = stt, formula = ppmForm1)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit1, fun = Kinhom)
resid.plot = diagnose(ppmFit1, which = "smooth", type = "Pearson")

kenv = envelope(ppmFit2, fun = Kinhom)
resid.plot = diagnose(ppmFit2, which = "smooth", type = "Pearson")
#newdata' had 8378 rows but variables found have 8655 rowsdata. check ppmfit data.

#####species interaction at 5km. NO avaiability grid is supplied here. 
#species.int = point.interactions(dat.ppm, 5)
ai.fit = ppmlasso(ppmForm1, data = dat.ppm, family = "area.inter", r = 1.5)
diagnose(ai.fit, which = "smooth", type = "Pearson")
pred.inter = predict(ai.fit, newdata=newstt)
pred.inter <- cbind(xydatan, pred.inter)
xy.inter <- rasterFromXYZ(as.data.frame(pred.inter)[, c("X", "Y", "pred.inter")])
plot(xy.inter)
