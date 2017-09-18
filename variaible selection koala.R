setwd("C:\\Users\\uqrdissa\\ownCloud\\Covariates_analysis\\Mark_S\\raster_stack")

#####  Step 1: read raster data from the folder and create a stack. ####

myfullstack.a <- list.files(pattern="\\.tif$", full.names = TRUE) 
myfullstack <- stack(myfullstack.a)
plot(myfullstack)

#variable selection
#check the coliniarity of variables 
v1 <- vifcor(myfullstack, th=0.9) #  identify collinear variables that should be excluded
myfullstack <- exclude(myfullstack,v1) # exclude the collinear variables that were identified in
#vifstep(myfullstack, th=10) # select variables which have VIF less than 10.

#select bestd etection model.
Detection.model=step(glm(presence~Dis_habitat_suitable_1+Dis_habitat_suitable_2+Dis_habitat_suitable_3+
                           distance_bridleway+distance_motorwayandlink+distance_path+distance_pedestrian+
                           distance_primaryandlink+distance_residentil+distance_secondaryandlink+distance_tertiaryandlink+
                           distance_trunkandlink+ distance_unclassified+s1_residential_dist+s1_unclassified_dist+s2_residential_dist+
                           s2_unclassified_dist+s3_residential_dist +scale(group), family= "binomial", data=Detection.data))