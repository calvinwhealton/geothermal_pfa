# code to run the thermal combination of variables

## reading-in rasters for prediction and error
# depth to 80 deg C
therm_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/d80p2-')
therm_err  <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/d80e2-')

# predicted temperature at 3500 m
#therm_pred35km <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/t35kmp2-')
#therm_err35km  <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/t35kme2-')

# substituting NA for -9999
therm_pred[therm_pred < 0] <- NA
therm_err[therm_err < 0] <- NA

#therm_pred35km[therm_pred35km < 0] <- NA
#therm_err35km[therm_err35km < 0] <- NA

## setting thresholds
# depth to 80 deg C thresholds [min, thresholds, max]
therm_thresh3 <- rev(c(8750,3000,2000,500))
therm_thresh5 <- rev(c(8750,4000,3000,2300,1500,500))

# creating histogram for depth to 80 degC
setwd(wd_image)
makeHist(rast=therm_pred
         ,thresh3=therm_thresh3
         ,thresh5=therm_thresh5
         ,rev_sc=TRUE
         ,plotnm='th_hist.png'
         ,yloc=-1.2*10^-4
         ,yshift=1.5*10^-5
         ,title='Depth to 80 Deg C (m)')

# converting into the play fairway scheme, saving, and making map
# three color----
th_3_0_3_NA <- convRastPFRank(rast=therm_pred
                              ,thresholds=therm_thresh3
                              ,ignore=-9999
                              ,rev_scale=TRUE)
saveRast(rast=th_3_0_3_NA
         ,wd=wd_raster
         ,rastnm='th_3_0_3_NA.tif')
makeMap (rast=th_3_0_3_NA
         ,plotnm='th_3_0_3_NA.png'
         ,wd=wd_image
         ,numCol=3
         ,comTy=NA
         ,numRF=1)

# five color-----
th_5_0_5_NA <- convRastPFRank(rast=therm_pred
                              ,thresholds=therm_thresh5
                              ,ignore=-9999
                              ,rev_scale=TRUE)
saveRast(rast=th_5_0_5_NA
         ,wd=wd_raster
         ,rastnm='th_5_0_5_NA.tif')
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)

# calcualting uncertainty maps for the play fairway scheme
# reading-in tables for interpolated values
th_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/th_d80_pfvar3.xlsx',1,header=FALSE))
th_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/th_d80_pfvar5.xlsx',1,header=FALSE))

# values to interpolate variance
th_means <- therm_pred@data@values
th_ses <- therm_err@data@values

# values used in making the interpolation table
mean_thd80 <- seq(750,6350,by=200) # range of means
std_thd80 <- seq(40,1740,by=50) # range of standard deviations

# substituting in values for NAs so interpolation algorithm will not crash
th_means[which(th_means %in% NA)] <- max(mean_thd80)
th_ses[which(th_ses %in% NA)] <- min(std_thd80)

# interpolating for the 3 color scheme
thvecPFvar3 <- interp2(x=std_thd80
                    ,y=mean_thd80
                    ,Z=th_interp_tab3
                    ,xp=th_ses
                    ,yp=th_means
                    ,method='linear'
                    )

# interpolating for the 5 color scheme
thvecPFvar5 <- interp2(x=std_thd80
                       ,y=mean_thd80
                       ,Z=th_interp_tab5
                       ,xp=th_ses
                       ,yp=th_means
                       ,method='linear'
)

# setting values back to NAs
thvecPFvar3[which(therm_pred@data@values %in% NA)] <- NA
thvecPFvar5[which(therm_pred@data@values %in% NA)] <- NA

# initializing raster for the stored values
therm_pfa_var3 <- therm_err
therm_pfa_var5 <- therm_err

# updating values of the raster
# note: setValues() did not work
values(therm_pfa_var3) <- thvecPFvar3
values(therm_pfa_var5) <- thvecPFvar5

# saving rasters and making maps
saveRast(rast=therm_pfa_var3
         ,wd=wd_raster
         ,rastnm='th_pfa_var3.tif')
makeMap (rast=therm_pfa_var3
         ,plotnm='th_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=therm_pfa_var5
         ,wd=wd_raster
         ,rastnm='th_pfa_var5.tif')
makeMap (rast=therm_pfa_var5
         ,plotnm='th_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)