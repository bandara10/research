dat.ppm00 <- ppm(dat.ppp, trend = ~ 1, 
              covariates = list( 
                CLY_005_015= CLY_005_015.im, dem=dem.im,distance_residentil = distance_residentil.im,  
                 foliagePC = foliagePC.im, habit3 = habit3.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,
                PTO_000_005=PTO_000_005.im,rainfallmeanannual = rainfallmeanannual.im, residential = residential.im,SAWC_030_060=SAWC_030_060.im,
                 suitable_0=suitable_0.im,TopoWI=TopoWI.im))
summary(dat.ppm00)$coefs.SE.CI
AIC(dat.ppm00) # 3488.409
windows();plot(dat.ppp)

dat.ppm01 <- ppm(trend =  ~ aspect91+CLY_005_015+dem+distance_bridleway+distance_cycleway+
                   distance_footway+distance_motorwayandlink +distance_pedestrian +distance_primaryandlink+
                   distance_residentil+distance_secondaryandlink+distance_steps +
                   distance_tertiaryandlink+distance_track+distance_trunkandlink+foliagePC+ 
                   habit3+habit1+habit2+hpop+lot_density+NTO_000_005+PTO_000_005+
                   rainfallmeanannual+residential+roads_motor+SAWC_030_060+Slope_percent_61 +suitable_0+
                   suitable_1+suitable_3+TopoWI+ unclassified,
                 covariates = list(aspect91= aspect91.im, CLY_005_015= CLY_005_015.im, dem=dem.im, distance_bridleway=distance_bridleway.im, distance_cycleway=distance_cycleway.im,
                                   distance_footway=distance_footway.im, distance_motorwayandlink = distance_motorwayandlink.im, distance_pedestrian = distance_pedestrian.im,  distance_primaryandlink = distance_primaryandlink.im,
                                   distance_residentil = distance_residentil.im, distance_secondaryandlink = distance_secondaryandlink.im, distance_steps = distance_steps.im, 
                                   distance_tertiaryandlink = distance_tertiaryandlink.im, distance_track = distance_track.im,distance_trunkandlink = distance_trunkandlink.im, foliagePC = foliagePC.im, 
                                   habit3 = habit3.im,habit1=habit1.im,habit2=habit2.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,PTO_000_005=PTO_000_005.im,
                                   rainfallmeanannual = rainfallmeanannual.im, residential = residential.im, roads_motor=roads_motor.im, SAWC_030_060=SAWC_030_060.im, Slope_percent_61 = Slope_percent_61.im, suitable_0=suitable_0.im,
                                   suitable_1=suitable_1.im, suitable_3=suitable_3.im,TopoWI=TopoWI.im, unclassified=unclassified.im))
# Extract the quadrature scheme from dat.ppm01:
Qz <- quad.ppm(dat.ppm01, drop = TRUE)



dat.ppm02 <- ppm(Qz,trend =  ~ aspect91+CLY_005_015+dem+distance_bridleway+distance_cycleway+
                   distance_footway+distance_motorwayandlink +distance_pedestrian +distance_primaryandlink+
                   distance_residentil+distance_secondaryandlink+distance_steps +
                   distance_tertiaryandlink+distance_track+distance_trunkandlink+foliagePC+ 
                   habit3+habit1+habit2+hpop+lot_density+NTO_000_005+PTO_000_005+
                   rainfallmeanannual+residential+roads_motor+SAWC_030_060+Slope_percent_61 +suitable_0+
                   suitable_1+suitable_3+TopoWI+ unclassified,
                 covariates = list(aspect91= aspect91.im, CLY_005_015= CLY_005_015.im, dem=dem.im, distance_bridleway=distance_bridleway.im, distance_cycleway=distance_cycleway.im,
                                   distance_footway=distance_footway.im, distance_motorwayandlink = distance_motorwayandlink.im, distance_pedestrian = distance_pedestrian.im,  distance_primaryandlink = distance_primaryandlink.im,
                                   distance_residentil = distance_residentil.im, distance_secondaryandlink = distance_secondaryandlink.im, distance_steps = distance_steps.im, 
                                   distance_tertiaryandlink = distance_tertiaryandlink.im, distance_track = distance_track.im,distance_trunkandlink = distance_trunkandlink.im, foliagePC = foliagePC.im, 
                                   habit3 = habit3.im,habit1=habit1.im,habit2=habit2.im,hpop=hpop.im, lot_density=lot_density.im, NTO_000_005=NTO_000_005.im,PTO_000_005=PTO_000_005.im,
                                   rainfallmeanannual = rainfallmeanannual.im, residential = residential.im, roads_motor=roads_motor.im, SAWC_030_060=SAWC_030_060.im, Slope_percent_61 = Slope_percent_61.im, suitable_0=suitable_0.im,
                                   suitable_1=suitable_1.im, suitable_3=suitable_3.im,TopoWI=TopoWI.im, unclassified=unclassified.im))
summary(dat.ppm02)$coefs.SE.CI
tdat <- data.frame(summary(dat.ppm02)$coefs.SE.CI); tdat[order(tdat$Zval),]
AIC(dat.ppm02) # 3464.227