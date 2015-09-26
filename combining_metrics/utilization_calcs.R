# calcuations for utilization risk factor

# reading-in prediction raster
util_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Utilization/lch4000p-')

util_pred[(util_pred %in% -9999)] <- NA

util_err <- util_pred
util_err[(util_pred > 0)] <- 2

# thresholds
util_min <- 5
util_max <- 25
util_thresh3 <- c(util_min,c(13.5,16),util_max)
util_thresh5 <- c(util_min,c(12,13.5,16,20),util_max)

# creating buffer images
makeWeightBuf(dist=5
              ,wd=wd_image
              ,plotnm='ut_buf_5.png')
makeWeightBuf(dist=3
              ,wd=wd_image
              ,plotnm='ut_buf_3.png')
makeWeightBuf(dist=7
              ,wd=wd_image
              ,plotnm='ut_buf_7.png')

# histogram
setwd(wd_image)
makeHist(rast=util_pred[(util_pred<100)]
         ,thresh3=util_thresh3
         ,thresh5=util_thresh5
         ,rev_sc=TRUE
         ,plotnm='ut_hist.png'
         ,yloc=-0.025
         ,yshift=0.003
         ,title='Utilization Surface Cost ($/MMBTU) (only < 100 plotted)')

# converting into the play fairway scheme
# three color
ut0_3_0_3_NA <- convRastPFRank(rast=util_pred
                               ,thresholds=util_thresh3
                               ,ignore=-9999
                               ,rev_scale=TRUE)
saveRast(rast=ut0_3_0_3_NA 
         ,wd=wd_raster
         ,rastnm='ut0_3_0_3_NA.tif')
makeMap(rast=ut0_3_0_3_NA 
        ,plotnm='ut0_3_0_3_NA.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

# five color-----
ut0_5_0_5_NA <- convRastPFRank(rast=util_pred
                               ,thresholds=util_thresh5
                               ,ignore=-9999
                               ,rev_scale=TRUE)
saveRast(rast=ut0_5_0_5_NA 
         ,wd=wd_raster
         ,rastnm='ut0_5_0_5_NA.tif')
makeMap(rast=ut0_5_0_5_NA 
        ,plotnm='ut0_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

# buffering utilization (5 km)
# three color
ut5_3_0_3_NA <- focal(ut0_3_0_3_NA
                      ,w=makeUtilBufWeight(5)
                      ,fun=max
                      ,na.rm=TRUE
                      ,pad=TRUE)
saveRast(rast=ut5_3_0_3_NA 
         ,wd=wd_raster
         ,rastnm='ut5_3_0_3_NA.tif')
makeMap(rast=ut5_3_0_3_NA 
        ,plotnm='ut5_3_0_3_NA.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy=NA
        ,numRF=1)

# five color
ut5_5_0_5_NA <- focal(ut0_5_0_5_NA
                      ,w=makeUtilBufWeight(5)
                      ,fun=max
                      ,na.rm=TRUE
                      ,pad=TRUE)
saveRast(rast=ut5_5_0_5_NA 
         ,wd=wd_raster
         ,rastnm='ut5_5_0_5_NA.tif')
makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

