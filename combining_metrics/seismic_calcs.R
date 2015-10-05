# script to complete seismic calcuations

# reading-in seismic earthquake prediction and error
seis_eq_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kmp-')
seis_eq_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kme-')

# reading-in seismic stress angle prediction and error
seis_stress_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2p-')
seis_stress_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2e-')

# setting values to NA where there is no data for STRESS
preds2 <- values(seis_stress_pred)
errs2 <- values(seis_stress_err)

preds2[which(preds2 %in% -9999)] <- NA
errs2[which(errs2 %in% -9999)] <- NA

values(seis_stress_pred) <- preds2
values(seis_stress_err) <- errs2

# setting values to NA where there is no data for EARTHQUAKE
preds1 <- values(seis_eq_pred)
errs1 <- values(seis_eq_err)

preds1[which(preds1 %in% -9999)] <- NA
errs1[which(errs1 %in% -9999)] <- NA

# special values 
errs1[which(errs1 %in% 3000)] <-  2550
errs1[which(errs1 %in% 2750)] <-  2515

values(seis_eq_pred) <- preds1
values(seis_eq_err) <- errs1

## thresholds
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

## values to interpolate variance
se_stress_means <- values(seis_stress_pred)
se_stress_sds <- values(seis_stress_err)

se_stress_means[se_stress_means > 70] <- 70
se_stress_sds[se_stress_means == 70] <- 0

se_eq_means <- values(seis_eq_pred)
se_eq_sds <- values(seis_eq_err)

# values used in making the interpolation table
mean_seSt <- seq(0,72,by=3) # range of means
std_seSt <- seq(0,310,by=10) # range of standard deviations

mean_seEq <- seq(0,26000,by=200) # range of means
std_seEq <- seq(0,2600,by=100) # range of standard deviations

# substituting in values for NAs so interpolation algorithm will not crash
se_stress_means[which(se_stress_means %in% NA)] <- min(mean_seSt)
se_stress_sds[which(se_stress_sds %in% NA)] <- min(std_seSt)

se_eq_means[which(se_eq_means %in% NA)] <- min(mean_seEq)
se_eq_sds[which(se_eq_sds %in% NA)] <- min(std_seEq)
se_eq_means[which(se_eq_means > 1234566)] <- min(mean_seEq)
se_eq_sds[which(se_eq_means > 1234566)] <- 0

# for STRESS
# interpolating for the 3 color scheme
seStvecPFvar3 <- interp2(x=std_seSt
                       ,y=mean_seSt
                       ,Z=se_stress_interp_tab3
                       ,xp=se_stress_sds
                       ,yp=se_stress_means
                       ,method='linear'
)
# interpolating for the 5 color scheme
seStvecPFvar5 <- interp2(x=std_seSt
                         ,y=mean_seSt
                         ,Z=se_stress_interp_tab5
                         ,xp=se_stress_sds
                         ,yp=se_stress_means
                         ,method='linear'
)


# setting values back to NAs
seStvecPFvar3[which(values(seis_stress_pred)%in% NA)] <- NA
seStvecPFvar5[which(values(seis_stress_pred) %in% NA)] <- NA

seStvecPFvar3[which(values(seis_stress_pred) %in% 70)] <- 0
seStvecPFvar5[which(values(seis_stress_pred) %in% 70)] <- 0

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


# saving rasters and making maps
saveRast(rast=calc(seSt_pfa_var3,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='seSt_pfa_sd3.tif')
makeMap (rast=calc(seSt_pfa_var3,fun=sqrt)
         ,plotnm='seSt_pfa_sd3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(seSt_pfa_var5,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='seSt_pfa_sd5.tif')
makeMap (rast=calc(seSt_pfa_var5,fun=sqrt)
         ,plotnm='seSt_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)


# for EARTHQUAKE
# interpolating for the 3 color scheme
seEqvecPFvar3 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab3
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)
# interpolating for the 5 color scheme
seEqvecPFvar5 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab5
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)

# setting values back to NAs
seEqvecPFvar3[which(values(seis_eq_pred) %in% NA)] <- NA
seEqvecPFvar5[which(values(seis_eq_pred) %in% NA)] <- NA

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


# saving rasters and making maps
saveRast(rast=calc(seEq_pfa_var3,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='seEq_pfa_sd3.tif')
makeMap (rast=calc(seEq_pfa_var3,fun=sqrt)
         ,plotnm='seEq_pfa_sd3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(seEq_pfa_var5,fun=sd)
         ,wd=wd_raster
         ,rastnm='seEq_pfa_sd5.tif')
makeMap (rast=calc(seEq_pfa_var5,fun=sqrt)
         ,plotnm='seEq_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# for COMBINED
se_pfa_var3s <- stack(c(seEq_pfa_var3,seSt_pfa_var3))
se_pfa_var5s <- stack(c(seEq_pfa_var5,seSt_pfa_var5))

se_pfa_var3 <- calc(se_pfa_var3s,fun=sum)
se_pfa_var5 <- calc(se_pfa_var5s,fun=sum)

# 0.25 because in the averaging of the two each raster is multiplied 
# by 0.5, so the variance will be multiplied by 0.25
values(se_pfa_var3) <- values(se_pfa_var3)*0.25
values(se_pfa_var5) <- values(se_pfa_var5)*0.25

# saving rasters and making maps
saveRast(rast=se_pfa_var3
         ,wd=wd_raster
         ,rastnm='se_pfa_var3.tif')
makeMap (rast=se_pfa_var3
         ,plotnm='se_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=se_pfa_var5
         ,wd=wd_raster
         ,rastnm='se_pfa_var5.tif')
makeMap (rast=se_pfa_var5
         ,plotnm='se_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# saving rasters and making maps
saveRast(rast=calc(se_pfa_var3,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='se_pfa_sd3.tif')
makeMap (rast=calc(se_pfa_var3,fun=sqrt)
         ,plotnm='se_pfa_sd3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(se_pfa_var5,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='se_pfa_sd5.tif')
makeMap (rast=calc(se_pfa_var5,fun=sqrt)
         ,plotnm='se_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# deleting files
rm(se_3_0_3_a,se_5_0_5_a,seSt_3_0_3_NA,seSt_5_0_5_NA,seEq_3_0_3_NA,seEq_5_0_5_NA
   ,se_eq_means,se_eq_sds,se_stress_means,se_stress_sds
   ,se_pfa_var3,se_pfa_var3s,se_pfa_var5,se_pfa_var5s
   ,seSt_pfa_var3,seSt_pfa_var5,seEq_pfa_var3,seEq_pfa_var5
   ,seis_eq_pred,seis_eq_err,seis_stress_pred,seis_stress_err
   ,sevecPFvar3s,sevecPFvar5s,seStvecPFvar3,seStvecPFvar5)