# script to complete seismic calcuations

# reading-in seismic earthquake prediction and error
seis_eq_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kmp-')
seis_eq_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kme-')

# reading-in seismic stress angle prediction and error
seis_stress_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2p-')
seis_stress_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2e-')

# setting values to NA where there is no data
seis_eq_pred[(seis_eq_pred %in% -9999)] <- NA
seis_stress_pred[(seis_stress_pred %in% -9999)] <- NA

seis_eq_err[(seis_eq_err %in% -9999)] <- NA

# modiviation of values for two special cases
seis_eq_err[(seis_eq_err %in% 3000)] <- 2550
seis_eq_err[(seis_eq_err %in% 2750)] <- 2515


# thresholds for stress angle
seis_stress_min <- 0.001
seis_stress_max <- 25
seis_stress_thresh3 <- c(seis_stress_min,c(8,16),seis_stress_max)
seis_stress_thresh5<- c(seis_stress_min,c(5,10,15,20),seis_stress_max)

# thresholds for earthquake distance
seis_eq_min <- 0.001
seis_eq_max <- 25
seis_eq_thresh3 <- 10^3*c(seis_eq_min,c(8,16),seis_eq_max)
seis_eq_thresh5<- 10^3*c(seis_eq_min,c(5,10,15,20),seis_eq_max)

# histograms
setwd(wd_image)
makeHist(rast=seis_eq_pred[(seis_eq_pred<10^5)]
         ,thresh3=seis_eq_thresh3
         ,thresh5=seis_eq_thresh5
         ,rev_sc=FALSE
         ,plotnm='seEq_hist.png'
         ,yloc=-1.5*10^-5
         ,yshift=3*10^-6
         ,title='Seismic Risk for Proximity to Earthquake (m)')

setwd(wd_image)
makeHist(rast=seis_stress_pred[(seis_stress_pred<100)]
         ,thresh3=seis_stress_thresh3
         ,thresh5=seis_stress_thresh5
         ,rev_sc=FALSE
         ,plotnm='seSt_hist.png'
         ,yloc=-0.017
         ,yshift=0.003
         ,title='Seismic Risk for Angle in Stress (diff degree)')


# coverting into the play fairway scheme EARTHQUAKES-----
# three color
seEq_3_0_3_NA <- convRastPFRank(rast=seis_eq_pred
                                ,thresholds=seis_eq_thresh3
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seEq_3_0_3_NA
         ,wd=wd_raster
         ,rastnm='seEq_3_0_3_NA.tif')
makeMap(rast=seEq_3_0_3_NA
        ,plotnm='seEq_3_0_3_NA.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

# five color
seEq_5_0_5_NA <- convRastPFRank(rast=seis_eq_pred
                                ,thresholds=seis_eq_thresh5
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seEq_5_0_5_NA
         ,wd=wd_raster
         ,rastnm='seEq_5_0_5_NA.tif')
makeMap(rast=seEq_5_0_5_NA
        ,plotnm='seEq_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

# converting into the play fairway scheme STRESS-----
# three color
seSt_3_0_3_NA <- convRastPFRank(rast=seis_stress_pred
                                ,thresholds=seis_stress_thresh3
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seSt_3_0_3_NA
         ,wd=wd_raster
         ,rastnm='seSt_3_0_3_NA.tif')
makeMap(rast=seSt_3_0_3_NA
        ,plotnm='seSt_3_0_3_NA.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

# five color
seSt_5_0_5_NA <- convRastPFRank(rast=seis_stress_pred
                                ,thresholds=seis_stress_thresh5
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seSt_5_0_5_NA
         ,wd=wd_raster
         ,rastnm='seSt_5_0_5_NA.tif')
makeMap(rast=seSt_5_0_5_NA
        ,plotnm='seSt_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)


# combining stress and earthquakes into play fairway seismic map
se_3_0_3_a <- calc(stack(c(seEq_3_0_3_NA,seSt_3_0_3_NA)),fun=mean)
se_5_0_5_a  <- calc(stack(c(seEq_5_0_5_NA,seSt_5_0_5_NA)),fun=mean)


saveRast(rast=se_3_0_3_a 
         ,wd=wd_raster
         ,rastnm='se_3_0_3_a.tif')
makeMap(rast=se_3_0_3_a
        ,plotnm='se_3_0_3_a.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

saveRast(rast=se_5_0_5_a
         ,wd=wd_raster
         ,rastnm='se_5_0_5_a.tif')
makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)


## calcualting uncertainty maps for the play fairway scheme
# reading-in tables for interpolated values
se_stress_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/se_stress_pfvar3.xlsx',1,header=FALSE))
se_stress_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/se_stress_pfvar5.xlsx',1,header=FALSE))

se_eq_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/se_eq_pfvar3.xlsx',1,header=FALSE))
se_eq_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/se_eq_pfvar5.xlsx',1,header=FALSE))

# values to interpolate variance
se_stress_means <- seis_stress_pred@data@values
se_stress_sds <- seis_stress_err@data@values

se_eq_means <- seis_eq_pred@data@values
se_eq_sds <- seis_eq_err@data@values

# values used in making the interpolation table
mean_seSt <- seq(0,70,by=3) # range of means
std_seSt <- seq(0,260,by=10) # range of standard deviations

mean_seEq <- seq(0,26000,by=200) # range of means
std_seEq <- seq(0,2550,by=100) # range of standard deviations


# substituting in values for NAs so interpolation algorithm will not crash
se_stress_means[which(se_stress_means %in% NA)] <- max(mean_seSt)
se_stress_sds[which(se_stress_sds %in% NA)] <- min(std_seSt)

se_eq_means[which(se_eq_means %in% NA)] <- max(mean_seEq)
se_eq_sds[which(se_eq_sds %in% NA)] <- min(std_seEq)

# for STRESS
# interpolating for the 3 color scheme
seStvecPFvar3 <- interp2(x=std_seSt
                       ,y=mean_seSt
                       ,Z=seSt_interp_tab3
                       ,xp=se_stress_sds
                       ,yp=se_stress_means
                       ,method='linear'
)
# interpolating for the 5 color scheme
seStvecPFvar5 <- interp2(x=std_seSt
                         ,y=mean_seSt
                         ,Z=seSt_interp_tab5
                         ,xp=se_stress_sds
                         ,yp=se_stress_means
                         ,method='linear'
)


# setting values back to NAs
seStvecPFvar3[which(seis_stress_pred@data@values %in% NA)] <- NA
seStvecPFvar5[which(seis_stress_pred@data@values %in% NA)] <- NA

# initializing raster for the stored values
seSt_pfa_var3 <- seis_stress_err
seSt_pfa_var5 <- seis_stress_err

# updating values of the raster
# note: setValues() did not work
values(seSt_pfa_var3) <- seStvecPFvar3
values(seSt_pfa_var5) <- seStvecPFvar5

# saving rasters and making maps
saveRast(rast=seSt_pfa_var3
         ,wd=wd_raster
         ,rastnm='seSt_pfa_var3.tif')
makeMap (rast=seSt_pfa_var3
         ,plotnm='seSt_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=seSt_pfa_var5
         ,wd=wd_raster
         ,rastnm='seSt_pfa_var5.tif')
makeMap (rast=seSt_pfa_var5
         ,plotnm='seSt_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)



# for EARTHQUAKE
# interpolating for the 3 color scheme
seEqvecPFvar3 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=seEq_interp_tab3
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)
# interpolating for the 5 color scheme
seEqvecPFvar5 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=seEq_interp_tab5
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)

# setting values back to NAs
seEqvecPFvar3[which(seis_eq_pred@data@values %in% NA)] <- NA
seEqvecPFvar5[which(seis_eq_pred@data@values %in% NA)] <- NA

# initializing raster for the stored values
seEq_pfa_var3 <- seis_eq_err
seEq_pfa_var5 <- seis_eq_err

# updating values of the raster
# note: setValues() did not work
values(seEq_pfa_var3) <- seEqvecPFvar3
values(seEq_pfa_var5) <- seEqvecPFvar5

# saving rasters and making maps
saveRast(rast=seEq_pfa_var3
         ,wd=wd_raster
         ,rastnm='seEq_pfa_var3.tif')
makeMap (rast=seEq_pfa_var3
         ,plotnm='seEq_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=seEq_pfa_var5
         ,wd=wd_raster
         ,rastnm='seEq_pfa_var5.tif')
makeMap (rast=seEq_pfa_var5
         ,plotnm='seEq_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)


# for COMBINED
sevecPFvar3s <- stack(c(seEqvecPFvar3,seStvecPFvar3))
sevecPFvar5s <- stack(c(seEqvecPFvar5,seStvecPFvar5))

sevecPFvar3 <- calc(sevecPFvar3s,fun=sum)
sevecPFvar5 <- calc(sevecPFvar5s,fun=sum)

# 0.25 because in the averaging of the two each raster is multiplied 
# by 0.5, so the variance will be multiplied by 0.25
values(sevecPFvar3) <- sevecPFvar3s@data@values*0.25
values(sevecPFvar5) <- sevecPFvar5s@data@values*0.25


# saving rasters and making maps
saveRast(rast=sevecPFvar3
         ,wd=wd_raster
         ,rastnm='se_pfa_var3.tif')
makeMap (rast=sevecPFvar3
         ,plotnm='se_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=sevecPFvar5
         ,wd=wd_raster
         ,rastnm='se_pfa_var5.tif')
makeMap (rast=sevecPFvar5
         ,plotnm='se_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

