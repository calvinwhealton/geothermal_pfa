# code to run the reservoir calcuations
# to convert to play fairway scheme

# reading-in rasters
res_pred_l1000 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_l1000p')
res_pred_1015 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1015p2')
res_pred_1520 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1520p')
res_pred_2025 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2025p')
res_pred_2530 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2530p')
res_pred_3035 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3035p')
res_pred_3540 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3540p')

# reservoir error
res_err_l1000 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_l1000e')
res_err_1015 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1015e')
res_err_1520 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1520e')
res_err_2025 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2025e')
res_err_2530 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2530e')
res_err_3035 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3035e')
res_err_3540 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3540e')

# substituting NA for -9999
preds1 <- values(res_pred_l1000)
preds2 <- values(res_pred_1015)
preds3 <- values(res_pred_1520)
preds4 <- values(res_pred_2025)
preds5 <- values(res_pred_2530)
preds6 <- values(res_pred_3035)
preds7 <- values(res_pred_3540)

errs1 <- values(res_err_l1000)
errs2 <- values(res_err_1015)
errs3 <- values(res_err_1520)
errs4 <- values(res_err_2025)
errs5 <- values(res_err_2530)
errs6 <- values(res_err_3035)
errs7 <- values(res_err_3540)

preds1[which(preds1 %in% -9999)] <- NA
preds2[which(preds2 %in% -9999)] <- NA
preds3[which(preds3 %in% -9999)] <- NA
preds4[which(preds4 %in% -9999)] <- NA
preds5[which(preds5 %in% -9999)] <- NA
preds6[which(preds6 %in% -9999)] <- NA
preds7[which(preds7 %in% -9999)] <- NA

errs1[which(preds1 %in% -9999)] <- NA
errs2[which(preds2 %in% -9999)] <- NA
errs3[which(preds3 %in% -9999)] <- NA
errs4[which(preds4 %in% -9999)] <- NA
errs5[which(preds5 %in% -9999)] <- NA
errs6[which(preds6 %in% -9999)] <- NA
errs7[which(preds7 %in% -9999)] <- NA

values(res_pred_l1000) <- preds1
values(res_pred_1015) <- preds2
values(res_pred_1520) <- preds3
values(res_pred_2025) <- preds4
values(res_pred_2530) <- preds5
values(res_pred_3035) <- preds6
values(res_pred_3540) <- preds7

values(res_err_l1000) <- errs1
values(res_err_1015) <- errs2
values(res_err_1520) <- errs3
values(res_err_2025) <- errs4
values(res_err_2530) <- errs5
values(res_err_3035) <- errs6
values(res_err_3540) <- errs7

# stacking rasters and taking maximum, assuming going for highest quality reservoir
# ignoring reservoirs shallower than 1000 m but all others taken
res_pred <- stack(c(res_pred_1015,res_pred_1520,res_pred_2025,res_pred_2530,res_pred_3035,res_pred_3540))
res_pred_max <- calc(res_pred,fun=max)

# stacking reservoir errors
# shallower than 1000m ignored again
res_err <- stack(c(res_err_1015,res_err_1520,res_err_2025,res_err_2530,res_err_3035,res_err_3540))

# making reservoir error term
# true/false raster of whether the max value is equal to the value in the layer
tf_rast1 <- (res_pred_max == res_pred[[1]])
tf_rast2 <- (res_pred_max == res_pred[[2]])
tf_rast3 <- (res_pred_max == res_pred[[3]])
tf_rast4 <- (res_pred_max == res_pred[[4]])
tf_rast5 <- (res_pred_max == res_pred[[5]])
tf_rast6 <- (res_pred_max == res_pred[[6]])

tf_stack <- stack(tf_rast1,tf_rast2,tf_rast3,tf_rast4,tf_rast5,tf_rast6)
tf_stacksum <- calc(tf_stack,fun=sum,na.rm=TRUE)

# stacking raster with error raster
# values iterate starting with 1000-1500m raster
# product used because 1 (true)
err_tf1 <- calc(stack(c(tf_rast1,res_err[[1]])),fun=prod)
err_tf2 <- calc(stack(c(tf_rast2,res_err[[2]])),fun=prod)
err_tf3 <- calc(stack(c(tf_rast3,res_err[[3]])),fun=prod)
err_tf4 <- calc(stack(c(tf_rast4,res_err[[4]])),fun=prod)
err_tf5 <- calc(stack(c(tf_rast5,res_err[[5]])),fun=prod)
err_tf6 <- calc(stack(c(tf_rast6,res_err[[6]])),fun=prod)

# making stacked raster and taking maximum
err_tf_all <- stack(c(err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6))
res_pred_max_err <- calc(err_tf_all,fun=max)
errs <- values(res_pred_max_err)
errs[which(values(res_pred_max) %in% NA)] <- NA
values(res_pred_max_err) <- errs

saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf1.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf2.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf3.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf4.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf5.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='err_tf6.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='tf_rast1.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='tf_rast2.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='tf_rast3.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='tf_rast4.tif')
saveRast(rast=err_tf1
         ,wd=wd_raster
         ,rastnm='tf_rast5.tif')
saveRast(rast=res_pred_max_err
         ,wd=wd_raster
         ,rastnm='res_pred_max_err.tif')

# thresholds
res_min <- 3*10^-5
res_max <- 301
res_thresh3 <- c(res_min,c(0.1,1.0),res_max)
res_thresh5 <- c(res_min,c(0.01,0.1,1.0,10),res_max)

# histogram
setwd(wd_image)
makeHist(rast=calc(res_pred_max,fun=log10)
         ,thresh3=log10(res_thresh3)
         ,thresh5=log10(res_thresh5)
         ,rev_sc=FALSE
         ,plotnm='re_hist.png'
         ,yloc=-0.15
         ,yshift=0.02
         ,title='log10 of Reservoir Productivity Index (L/MPa-s)')

# converting into the play fairway scheme
# three color
re_3_0_3_NA <- convRastPFRank(rast=calc(res_pred_max,fun=log10)
                              ,thresholds=log10(res_thresh3)
                              ,ignore=-9999
                              ,rev_scale=FALSE)
saveRast(rast=re_3_0_3_NA
         ,wd=wd_raster
         ,rastnm='re_3_0_3_x.tif')
makeMap(rast=re_3_0_3_NA
        ,plotnm='re_3_0_3_x.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

# five color-----
re_5_0_5_NA <- convRastPFRank(rast=calc(res_pred_max,fun=log10)
                              ,thresholds=log10(res_thresh5)
                              ,ignore=-9999
                              ,rev_scale=FALSE)
saveRast(rast=re_5_0_5_NA
         ,wd=wd_raster
         ,rastnm='re_5_0_5_x.tif')
makeMap (rast=re_5_0_5_NA
         ,plotnm='re_5_0_5_x.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)

# making uncertainty map
# reading-in tables for interpolated values
re_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/re_pfvar3.xlsx',1,header=FALSE))
re_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/re_pfvar5.xlsx',1,header=FALSE))

# values to interpolate variance
preds <- values(res_pred_max)
preds[values(res_pred_max) == -9999] <- NA
values(res_pred_max) <- preds
re_means <- log(values(res_pred_max))
re_uncer <- values(res_pred_max_err)

# values used in making the interpolation table
mean_re <- seq(-9.5,6.25,0.25) # range of means
uncer_re <- seq(0,2,by=0.1) # range of coefficient of variation

# substituting in values for NAs so interpolation algorithm will not crash
re_means[which(re_means %in% NA)] <- min(mean_re)
re_uncer[which(re_uncer %in% NA)] <- min(uncer_re)

re_means[which(re_means %in% 0)] <- min(mean_re)
re_uncer[which(re_uncer %in% 0)] <- min(uncer_re)


# interpolating for the 3 color scheme
revecPFvar3 <- interp2(x=uncer_re
                       ,y=mean_re
                       ,Z=re_interp_tab3
                       ,xp=re_uncer
                       ,yp=re_means
                       ,method='linear'
)

# interpolating for the 5 color scheme
revecPFvar5 <- interp2(x=uncer_re
                       ,y=mean_re
                       ,Z=re_interp_tab5
                       ,xp=re_uncer
                       ,yp=re_means
                       ,method='linear'
)

# setting values back to NAs
revecPFvar3[which(values(res_pred_max) %in% NA)] <- NA
revecPFvar5[which(values(res_pred_max) %in% NA)] <- NA

revecPFvar3[which(values(res_pred_max) %in% 0)] <- 0
revecPFvar5[which(values(res_pred_max) %in% 0)] <- 0

# initializing raster for the stored values
re_pfa_var3 <- res_pred_max_err
re_pfa_var5 <- res_pred_max_err

# updating values of the raster
# note: setValues() did not work
values(re_pfa_var3) <- revecPFvar3
values(re_pfa_var5) <- revecPFvar5



# saving rasters and making maps
saveRast(rast=re_pfa_var3
         ,wd=wd_raster
         ,rastnm='re_pfa_var3.tif')
makeMap (rast=re_pfa_var3
         ,plotnm='re_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=re_pfa_var5
         ,wd=wd_raster
         ,rastnm='re_pfa_var5.tif')
makeMap (rast=re_pfa_var5
         ,plotnm='re_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)



# saving rasters and making maps
saveRast(rast=calc(re_pfa_var3,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='re_pfa_sd3.tif')
makeMap (rast=calc(re_pfa_var3,fun=sqrt)
         ,plotnm='re_pfa_sd3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(re_pfa_var5,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='re_pfa_sd5.tif')
makeMap (rast=calc(re_pfa_var5,fun=sqrt)
         ,plotnm='re_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

rm(res_err_1015,res_err_1520,res_err_2025,res_err_2530,res_err_3035,res_err_3540,res_err_l1000
   ,res_pred_1015,res_pred_1520,res_pred_2025,res_pred_2530,res_pred_3035,res_pred_3540,res_pred_l1000
   ,re_3_0_3_NA,re_5_0_5_NA,revecPFvar3,revecPFvar5
   ,preds,errs,preds1,preds2,preds3,preds4,preds5,preds6,preds7
   ,errs1,errs2,errs3,errs4,errs5,errs6,errs7
   ,err_tf_all,err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6
   ,re_pfa_var3,re_pfa_var5,re_uncer,res_err,res_pred
   ,res_pred_max,res_pred_max_err
   ,tf_rast1,tf_rast2,tf_rast3,tf_rast4,tf_rast5,tf_rast6
   ,tf_stack,tf_stacksum)
