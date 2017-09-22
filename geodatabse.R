require(rgdal)

# The input file geodatabase
fgdb <- "C:\\Users\\uqrdissa\\BACKUPS\\KoaladataNew.gdb"

# List all feature classes in a file geodatabase
ogrListLayers(fgdb)

# Read the feature class
LGA <- readOGR(fgdb,layer="LGA_study_areashp")
# Determine the FC extent, projection, and attribute information
summary(LGA)
LGA <- spTransform (LGA, crs("+proj=longlat +ellps=WGS84"))# trasform to lat lng before see in google.
# View the feature class
plot(LGA)
LGA <- as.data.frame(LGA)
#save(LGA, obs, file="LGA.RData")
library(leaflet)
lnd84 <- readRDS('data/LGA')
leaflet() %>%
  addTiles() %>%
  addPolygons(data = LGA)
