library(ppmlasso)
data(BlueMountains)
# Pre-standardise observer bias variables
backg.env = BlueMountains$env
#head(backg.env)
xydatan <- backg.env[c(1,2)] # get x y colomns to be cbind later with predicted valuves.
#xydata <- backg.env[c(1,2)]
#coordinates(xydata) <- ~X+Y
#plot(xydata)

#stand.D_MAIN_RDS = backg.env$D_MAIN_RDS
#stand.D_URBAN = backg.env$D_URBAN
#stand.D_MAIN_RDS = standardise.X(backg.env$D_MAIN_RDS)$X
stand.D_MAIN_RDS=scale.default(backg.env$D_MAIN_RDS, center = TRUE, scale = TRUE) #standarise
stand.D_URBAN=scale.default(backg.env$D_URBAN, center = TRUE, scale = TRUE) #standarise



#stand.D_MAIN_RDS # set distance to zero
#stand.D_URBAN
#backg.env$D_MAIN_RDS = stand.D_MAIN_RDS
#backg.env$D_URBAN = stand.D_URBAN
# To fit a Poisson point process model at a spatial resolution of 1km
ppmForm = ~ poly(FC, TMP_MAX,TMP_MIN, RAIN_ANN, degree = 2) 
ppmFit = ppmlasso(ppmForm, sp.xy = BlueMountains$eucalypt,
                  env.grid = backg.env, sp.scale = 1)

# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = BlueMountains$env
#newEnv$D_MAIN_RDS = min(stand.D_MAIN_RDS)
pred.biasCorrect.1 = predict(ppmFit, newdata=newEnv)
predictions <- cbind(xydatan, pred.biasCorrect.1)
xy.rr <- rasterFromXYZ(as.data.frame(predictions)[, c("X", "Y", "pred.biasCorrect.1")])
plot(xy.rr, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = BlueMountains$eucalypt,
        env.grid = BlueMountains$env, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv = envelope(ppmFit, fun = Kinhom)
resid.plot = diagnose(ppmFit, which = "smooth", type = "Pearson")


########
# Pre-standardise observer bias variables
backg.env = BlueMountains$env
head(backg.env)
#stand.D_MAIN_RDS = backg.env$D_MAIN_RDS
#stand.D_URBAN = backg.env$D_URBAN
#stand.D_MAIN_RDS = standardise.X(backg.env$D_MAIN_RDS)$X
stand.D_MAIN_RDS=scale.default(backg.env$D_MAIN_RDS, center = TRUE, scale = TRUE) #standarise
stand.D_URBAN=scale.default(backg.env$D_URBAN, center = TRUE, scale = TRUE) #standarise
#stand.D_MAIN_RDS # set distance to zero
#stand.D_URBAN
backg.env$D_MAIN_RDS = stand.D_MAIN_RDS
backg.env$D_URBAN = stand.D_URBAN
# To fit a Poisson point process model at a spatial resolution of 1km
ppmForm.2 = ~ poly(FC,TMP_MAX,TMP_MIN, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppmFit.2 = ppmlasso(ppmForm.2, sp.xy = BlueMountains$eucalypt,
                  env.grid = backg.env, sp.scale = 1)
summary(ppmForm.2)
summary(ppmFit.2)
# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = BlueMountains$env
#newEnv$D_MAIN_RDS = min(stand.D_MAIN_RDS, D_URBAN)
pred.biasCorrect.2 = predict(ppmFit.2, newdata=newEnv)
head(pred.biasCorrect.2)
predictions.2 <- cbind(xydatan, pred.biasCorrect.2)
xy.rr.2 <- rasterFromXYZ(as.data.frame(predictions.2)[, c("X", "Y", "pred.biasCorrect.2")])
plot(xy.rr.2, las=0)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = BlueMountains$eucalypt,
        env.grid = BlueMountains$env, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv.2 = envelope(ppmFit.2, fun = Kinhom)
resid.plot.2 = diagnose(ppmFit.2, which = "smooth", type = "Pearson")

#######
# To fit a Poisson point process model at a spatial resolution of 1km
library(ppmlasso)
data(BlueMountains)
backg.env = BlueMountains$env
stand.D_MAIN_RDS = backg.env$D_MAIN_RDS
stand.D_URBAN = backg.env$D_URBAN
#stand.D_MAIN_RDS = standardie.X(backg.env$D_MAIN_RDS)$X
stand.D_MAIN_RDS=scale.default(backg.env$D_MAIN_RDS, center = TRUE, scale = TRUE) #standarise
stand.D_URBAN=scale.default(backg.env$D_URBAN, center = TRUE, scale = TRUE) #standarise
backg.env$D_MAIN_RDS = stand.D_MAIN_RDS
backg.env$D_URBAN = stand.D_URBAN

ppmForm.3 = ~ poly(FC,TMP_MAX,TMP_MIN, RAIN_ANN, degree = 2) + poly(D_MAIN_RDS, D_URBAN, degree = 2)
ppmFit.3 = ppmlasso(ppmForm.3, sp.xy = BlueMountains$eucalypt,
                  env.grid = backg.env, sp.scale = 1)

# To predict using model-based control of observer bias (at min value for D_MAIN_RDS):
newEnv = BlueMountains$env
newEnv$D_MAIN_RDS = min(stand.D_MAIN_RDS)
newEnv$D_URBAN = min(stand.D_MAIN_RDS)
pred.biasCorrect.3 = predict(ppmFit.3, newdata=newEnv)
predictions.3 <- cbind(xydatan, pred.biasCorrect.3)
xy.rr.3 <- rasterFromXYZ(as.data.frame(predictions.3)[, c("X", "Y", "pred.biasCorrect.3")])
plot(xy.rr.3,asp=1)
# To find the resolution (in the range from 0.5 to 16 km):
scales = c(0.5, 1, 2, 4, 8, 16)
findres(scales, sp.xy = BlueMountains$eucalypt,
        env.grid = BlueMountains$env, formula = ppmForm)
#which returns the log-likelihood at each scale, difference < 2 at 1km scale
# Diagnostic plots as in Fig 5:
kenv.3 = envelope(ppmFit.3, fun = Kinhom)
resid.plot.3 = diagnose(ppmFit.3, which = "smooth", type = "Pearson")





