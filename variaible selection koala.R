setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack <- stack(myfullstack.a)
plot(myfullstack)
#v <- Variogram(myfullstack[[20]])# check variograme for variables.http://r-gis.net/?q=positional_uncertainty
#plot(v)
#variable selection
#check the coliniarity of variables 
v1 <- vifcor(myfullstack, th=0.9) #  identify collinear variables that should be excluded
myfullstack <- exclude(myfullstack,v1) # exclude the collinear variables that were identified in
vifstep(myfullstack, th=10) # select variables which have VIF less than 10.

#select the best detection model.

Detection.model.s =glm (presence~Dis_habitat_suitable_1+Dis_habitat_suitable_2+Dis_habitat_suitable_3+distance_bridleway+distance_motorwayandlink+
                               distance_path+distance_pedestrian+distance_primaryandlink+distance_residentil+distance_secondaryandlink+
                               distance_tertiaryandlink+distance_trunkandlink+distance_unclassified+s1_residential_dist+
                               s1_unclassified_dist+s2_residential_dist+s2_unclassified_dist+s3_residential_dist+habit1+habit2+
                          habit3+scale(group),family= "binomial", data=Detection.data)
summary(Detection.model.s)
######Step 5 - Fit an inhomogeneous Poisson point process  that weights the log-likelihood by 1/p.det . ####
#use step function here then use significant variable in the next model. or else go to line 79. I asume this is correct way to do it.
IPP.corrected.1= step(glm(presence~twi+tpo+temp +aspect+awc+clay+elev+fpc+habit2pc+hpop+
                            lot_density+nitro+roadk+sbd+suitable_1+suitable_3,
                          family="binomial",weights=(1/p.det)*10000^(1-presence),data=IPP.data)) 
summary(IPP.corrected.1)