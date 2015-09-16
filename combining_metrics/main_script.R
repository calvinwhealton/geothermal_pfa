##################################################
# script for combining risk factor rasters (RFs)
# into the play fairway metric raster
#
# steps include:
# -1. importing functions/libraries
# 0. importing rasters
# 1. verifying that the rasters are in the same system
# 2. converting individual rasters into the play fairway rank
# 3. combining the RFs using a specified function
# 4. saving output

##### importing functions/libraries #####

# libraries
library(sp)
library(raster)
library(rgdal)
library(rasterVis)
library(maps)
library(maptools)

##### defining working directories #####
wd_raster <- '/Users/calvinwhealton/Dropbox/PFA_Rasters'
wd_image <- '/Users/calvinwhealton/Desktop/combining_rf'

##### loading-in state/county shapefiles #####
# states
States = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# counties
Counties = readOGR(dsn='/Users/calvinwhealton/GitHub/geothermal/combining_metrics/', layer="us_county_WGS84_prj", stringsAsFactors=FALSE)
NY_co = Counties[which(Counties$STATEFP == "36"),]
PA_co = Counties[which(Counties$STATEFP == "42"),]
WV_co = Counties[which(Counties$STATEFP == "54"),]

# converting to UTM 17N system
NY2 <- spTransform(NY,CRS("+init=epsg:31986"))
PA2 <- spTransform(PA,CRS("+init=epsg:31986"))
WV2 <- spTransform(WV,CRS("+init=epsg:31986"))

NY_co2 <- spTransform(NY_co,CRS("+init=epsg:31986"))
PA_co2 <- spTransform(PA_co,CRS("+init=epsg:31986"))
WV_co2 <- spTransform(WV_co,CRS("+init=epsg:31986"))

##### user-defined functions #####
setwd('/Users/calvinwhealton/GitHub/geothermal/combining_metrics')
source('checkSameProjCoords.R')
source('convertRasterPFAMetric.R')
source('combineRFs.R')
source('makeUtilBuf.R')
source('makeHist.R')
source('makeMap.R')
source('saveRast.R')
source('plotWeightBuf.R')

##### THERMAL ######
# importing rasters
therm_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/d80p2-')
therm_err  <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Thermal/d80e2-')

therm_pred[therm_pred < 0] <- NA
therm_err[therm_err < 0] <- NA

# setting thresholds
therm_thresh3 <- rev(c(8750,3000,2000,500))  # thermal (made up)
therm_thresh5 <- rev(c(8750,4000,3000,2300,1500,500)) # thermal (made up)

# creating histogram
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

##### RESERVOIR ######
res_pred_l1000 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_l1000p')
res_pred_1015 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1015p2')
res_pred_1520 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1520p')
res_pred_2025 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2025p')
res_pred_2530 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2530p')
res_pred_3035 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3035p')
res_pred_3540 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3540p')

# reservoir error
res_pred_l1000 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_l1000e')
res_err_1015 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1015e')
res_err_1520 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_1520e')
res_err_2025 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2025e')
res_err_2530 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_2530e')
res_err_3035 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3035e')
res_err_3540 <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Reservoirs/fff_3540e')

# stacking rasters and taking maximum, assuming going for highest quality reservoir
# ignoring reservoirs shallower than 1000 m
res_pred <- stack(c(res_pred_1015,res_pred_1520,res_pred_2025,res_pred_2530,res_pred_3035,res_pred_3540))
res_pred_max <- calc(res_pred,fun=max,na.rm=TRUE)

res_pred_max2 <- res_pred_max
res_pred_max2[res_pred_max <0] <- NA

# stacking reservoir errors
res_err <- stack(c(res_err_1015,res_err_1520,res_err_2025,res_err_2530,res_err_3035,res_err_3540))

# making reservoir error term
# raster of whether the max value is equal to the value in the layer
tf_rast1 <- (res_pred_max == res_pred[[1]])
tf_rast2 <- (res_pred_max == res_pred[[2]])
tf_rast3 <- (res_pred_max == res_pred[[3]])
tf_rast4 <- (res_pred_max == res_pred[[4]])
tf_rast5 <- (res_pred_max == res_pred[[5]])
tf_rast6 <- (res_pred_max == res_pred[[6]])

# stacking raster with error raster
err_tf1 <- calc(stack(c(tf_rast1,res_err[[1]])),fun=prod)
err_tf2 <- calc(stack(c(tf_rast2,res_err[[2]])),fun=prod)
err_tf3 <- calc(stack(c(tf_rast3,res_err[[3]])),fun=prod)
err_tf4 <- calc(stack(c(tf_rast4,res_err[[4]])),fun=prod)
err_tf5 <- calc(stack(c(tf_rast5,res_err[[5]])),fun=prod)
err_tf6 <- calc(stack(c(tf_rast6,res_err[[6]])),fun=prod)

# making stacked raster and taking maximum
res_pred_max_err <- calc(stack(c(err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6)),fun=max)
res_pred_max_err[(res_pred_max_err < 0)] <- NA

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
         ,yloc=-0.22
         ,yshift=0.02
         ,title='log10 of Reservoir Ideality')


# converting into the play fairway scheme
# three color
re_3_0_3_NA <- convRastPFRank(rast=calc(res_pred_max2,fun=log10)
                             ,thresholds=log10(res_thresh3)
                             ,ignore=-9999
                             ,rev_scale=FALSE)
saveRast(rast=re_3_0_3_NA
         ,wd=wd_raster
         ,rastnm='re_3_0_3_NA.tif')
makeMap(rast=re_3_0_3_NA
         ,plotnm='re_3_0_3_NA.png'
         ,wd=wd_image
         ,numCol=3
         ,comTy=NA
         ,numRF=1)

# five color-----
re_5_0_5_NA <- convRastPFRank(rast=calc(res_pred_max2,fun=log10)
                              ,thresholds=log10(res_thresh5)
                              ,ignore=-9999
                              ,rev_scale=FALSE)
saveRast(rast=re_5_0_5_NA
         ,wd=wd_raster
         ,rastnm='re_5_0_5_NA.tif')
makeMap (rast=re_5_0_5_NA
         ,plotnm='re_5_0_5_NA.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)


##### UTILIZATION ####
# utilization prediciton and error
util_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Utilization/lch4000p-')

util_pred[(util_pred %in% -9999)] <- NA

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

##### SEISMIC ######
seis_eq_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kmp-')
seis_eq_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kme-')

seis_stress_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2p-')
seis_stress_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2e-')

seis_eq_pred[(seis_eq_pred %in% -9999)] <- NA
seis_stress_pred[(seis_stress_pred %in% -9999)] <- NA

# thresholds
seis_stress_min <- 0.001
seis_stress_max <- 25
seis_stress_thresh3 <- c(seis_stress_min,c(8,16),seis_stress_max)
seis_stress_thresh5<- c(seis_stress_min,c(5,10,15,20),seis_stress_max)

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
         ,plotnm='se_eq_hist.png'
         ,yloc=-1.5*10^-5
         ,yshift=3*10^-6
         ,title='Seismic Risk for Proximity to Earthquake (m)')

setwd(wd_image)
makeHist(rast=seis_stress_pred[(seis_stress_pred<100)]
         ,thresh3=seis_stress_thresh3
         ,thresh5=seis_stress_thresh5
         ,rev_sc=FALSE
         ,plotnm='se_stress_hist.png'
         ,yloc=-0.017
         ,yshift=0.003
         ,title='Seismic Risk for Angle in Stress (diff degree)')


# combining into the play fairway scheme EARTHQUAKES-----
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

# combining into the play fairway scheme STRESS-----
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
makeMap(rast=seEq_5_0_5_NA
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

##### checking rasters are in the same system #####
check_therm_res <- checkSameProjCoords(therm_pred,res_pred)
check_util_therm <- checkSameProjCoords(util_pred,therm_pred)
check_util_res <- checkSameProjCoords(util_pred,res_pred)


##### making a quick and simple map to compare results #####
# combining maps
# creating a stacked raster
comb_pfa3 <- stack(c(re_3_0_3_NA,th_3_0_3_NA,ut5_3_0_3_NA,se_3_0_3_a))
comb_pfa5 <- stack(c(re_5_0_5_NA,th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))

# using sums
co_3_0_12_s <- calc(comb_pfa3,fun=sum,na.rm=FALSE)
co_5_0_20_s <- calc(comb_pfa5,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_12_s
         ,thresh3=c(0,1,2,3)*4
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_12_s_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Combined Play Fairway Metric: Sum')
saveRast(rast=co_3_0_12_s
         ,wd=wd_raster
         ,rastnm='co_3_0_12_s.tif')
makeMap(rast=co_3_0_12_s
        ,plotnm='co_3_0_12_s.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_20_s
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*4
         ,rev_sc=FALSE
         ,plotnm='co_5_0_20_s_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='Combined Play Fairway Metric: Sum')
saveRast(rast=co_5_0_20_s
         ,wd=wd_raster
         ,rastnm='co_5_0_20_s.tif')
makeMap(rast=co_5_0_20_s
        ,plotnm='co_5_0_20_s.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=4)

# using products
co_3_0_81_p <- calc(comb_pfa3,fun=prod,na.rm=FALSE)
co_5_0_625_p <- calc(comb_pfa5,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_81_p
         ,thresh3=c(0,1,2,3)^4
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_81_p_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Combined Play Fairway Metric: Product')
saveRast(rast=co_3_0_81_p
         ,wd=wd_raster
         ,rastnm='co_3_0_81_p.tif')
makeMap(rast=co_3_0_81_p
        ,plotnm='co_3_0_81_p.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_625_p
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^4
         ,rev_sc=FALSE
         ,plotnm='co_5_0_625_p_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='Combined Play Fairway Metric: Product')
saveRast(rast=co_5_0_625_p
         ,wd=wd_raster
         ,rastnm='co_5_0_625_p.tif')
makeMap(rast=co_5_0_625_p
        ,plotnm='co_5_0_625_p.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=4)

# using minimums
co_3_0_3_m <- calc(comb_pfa3,fun=min,na.rm=FALSE)
co_5_0_5_m <- calc(comb_pfa5,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Combined Play Fairway Metric: Minimum')
saveRast(rast=co_3_0_3_m
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m.tif')
makeMap(rast=co_3_0_3_m
        ,plotnm='co_3_0_3_m.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=4)


setwd(wd_image)
makeHist(rast=co_5_0_5_m
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='Combined Play Fairway Metric: Minimum')
saveRast(rast=co_5_0_5_m
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m.tif')
makeMap(rast=co_5_0_5_m
        ,plotnm='co_5_0_5_m.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=1)


# geologic only
comb_pfa3_geo <- stack(c(re_3_0_3_NA,th_3_0_3_NA,se_3_0_3_a))
comb_pfa5_geo <- stack(c(re_5_0_5_NA,th_5_0_5_NA,se_5_0_5_a))


co_3_0_9_s_geo <- calc(comb_pfa3_geo,fun=sum,na.rm=FALSE)
co_5_0_15_s_geo <- calc(comb_pfa5_geo,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_9_s_geo
         ,thresh3=c(0,1,2,3)*3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_9_s_geo_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Combined Play Fairway Metric: Sum')
saveRast(rast=co_3_0_9_s
         ,wd=wd_raster
         ,rastnm='co_3_0_9_s_geo.tif')
makeMap(rast=co_3_0_12_s
        ,plotnm='co_3_0_9_s_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_15_s_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_20_s_geo_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='Combined Play Fairway Metric: Sum')
saveRast(rast=co_5_0_15_s_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_15_s_geo.tif')
makeMap(rast=co_5_0_15_s_geo
        ,plotnm='co_5_0_15_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=3)