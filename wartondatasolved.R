

a <- backg.env[c(1,2,3)]
b <- backg.env[c(1,2,4)]
c <- backg.env[c(1,2,5)]
d <- backg.env[c(1,2,6)]
e <- backg.env[c(1,2,7)]
f <- backg.env[c(1,2,8)]

ar <- rasterFromXYZ(as.data.frame(a)[, c("X", "Y", "FC")])
br <- rasterFromXYZ(as.data.frame(b)[, c("X", "Y", "D_MAIN_RDS")])
cr <- rasterFromXYZ(as.data.frame(c)[, c("X", "Y", "D_URBAN")])
dr <- rasterFromXYZ(as.data.frame(d)[, c("X", "Y", "RAIN_ANN")])
er <- rasterFromXYZ(as.data.frame(e)[, c("X", "Y", "TMP_MAX")])
fr <- rasterFromXYZ(as.data.frame(f)[, c("X", "Y", "TMP_MIN")])

s <- stack(ar,br,cr,dr,er,fr)
ss <- as.data.frame(s, xy = TRUE)

colnames(ss)[1] <- "X"
colnames(ss)[2] <- "Y"
quad.1 = sample.quad(env.grid = ss, sp.scale = 1, file = "Quad")
