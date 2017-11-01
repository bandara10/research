d <-  Detection.data[c(51,15,20)]
s <-dropLayer  (myfullstack, 1:3, 15:18, 20:50)

Detection.model=glm(presence~ distance_pedestrian + distance_tertiaryandlink, family= "binomial", data=d)



stand.distance_tertiaryandlink=scale.default(d$distance_tertiaryandlink, center = TRUE, scale = TRUE) #standarise
d$distance_tertiaryandlink = stand.distance_tertiaryandlink

d2 <- d
d2$distance_tertiaryandlink = min(stand.distance_tertiaryandlink)

Detection.model2=glm(presence~  distance_pedestrian + distance_tertiaryandlink, family= "binomial", d2)
unclass(summary(Detection.model2))

myPred = predict(s, Detection.model, type = "response")
plot(myPred, xlab = "x", ylab= "y",main="detection model")

myPred2 = predict(s, Detection.model2, type = "response")
plot(myPred2, xlab = "x", ylab= "y",main="detection model2")
