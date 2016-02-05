##################################################
# script for combining risk factor rasters (RFs)
# into the play fairway metric raster
# and completing additional analysis for some cities
#
# steps include:
# --importing functions/libraries
# --importing rasters for each risk factor
# --converting individual rasters into the play fairway rank (0-3 or 0-5)
# --combining the RFs using a specified function
# --saving output
# --completing more detailed comparisons for specific locations
#
# modifications when running on a different machine:
# -- change all working directories

##### importing functions/libraries #####

# libraries
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
#library(rasterVis)    # 
#library(maps)         #
#library(maptools)     #
library(xlsx)         # for reading in xlsx tables
library(rgeos)        # for buffering places of interest
library(RColorBrewer) # R color brewer palettes for parallel axis plot
library(pracma)       # for interpolation in tables

##### defining working directories #####
# need to be changed based on machine

# location to SAVE rasters
wd_raster_out <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/Rasters_out'

# location to READ rasters
wd_raster_in <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/Rasters_in'

# location to save images
wd_image <- '/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics/figs'

# location of code/scripts are stored
wd_code <- '/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics'

# location of input shapefiles
wd_shapefiles <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/'

##### loading-in state/county shapefiles #####
# states
States = readOGR(dsn=paste(wd_shapefiles,'us_state_WGS',sep=''), layer="us_state_WGS", stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# counties
Counties = readOGR(dsn=paste(wd_shapefiles,'us_county_WGS84_prj',sep=''), layer="us_county_WGS84_prj", stringsAsFactors=FALSE)
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

# importing city locations, the US Census Places
cities <- readOGR(dsn=paste(wd_shapefiles,'usCensusPlaces',sep=''), layer="usCensusPlaces", stringsAsFactors=FALSE)
cities2 <- spTransform(cities,CRS("+init=epsg:31986"))

# removing the untransformed shapefiles
rm(NY,PA,WV,States,Counties,NY_co,PA_co,WV_co)

## importing places of interest
poi <- readOGR(dsn=paste(wd_shapefiles,'placesofinterest2',sep=''), layer="placesofinterest2", stringsAsFactors=FALSE)

# converting to UTM 17 N system
poi2 <- spTransform(poi,CRS("+init=epsg:31986"))

# buffering places of interest by 5000 and 10000 m
poi2Buf5 <- gBuffer(poi2,width=5000)
poi2Buf10 <- gBuffer(poi2,width=10000)

# deleting old file
rm(poi)

##### user-defined functions #####
# functions read-in from other R scripts
setwd(wd_code)
source('convertRasterPFAMetric.R')
source('combineRFs.R')
source('makeUtilBuf.R')
source('makeHist.R')
source('makeMap.R')
source('saveRast.R')
source('plotWeightBuf.R')

##### THERMAL ######
# importing rasters, depth to 80 DegC and standard error of prediction
therm_pred <- raster(paste(wd_raster_in,'/Thermal/d80p2-',sep=''))
therm_err  <- raster(paste(wd_raster_in,'/Thermal/d80e2-',sep=''))

# values of -9999 are no data, so replaced with NA
therm_pred[therm_pred < 0] <- NA
therm_err[therm_err < 0] <- NA

# setting thresholds
# deep values are less favorable
therm_thresh3 <- rev(c(8750,3000,2000,500))
therm_thresh5 <- rev(c(8750,4000,3000,2300,1500,500))

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
          ,wd=wd_raster_out
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
         ,wd=wd_raster_out
         ,rastnm='th_5_0_5_NA.tif')
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)

# making thermal uncertainty maps
th_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics/th_d80_pfvar3.xlsx',1,header=FALSE))
th_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics/th_d80_pfvar5.xlsx',1,header=FALSE))

# values to interpolate variance
th_means <- values(therm_pred)
th_ses <- values(therm_err)

# values used in making the interpolation table
# must check values from make_interp_table.R
mean_thd80 <- seq(750,6350,by=200) # range of means
std_thd80 <- seq(40,1740,by=50) # range of standard deviations

# substituting in values for NAs so interpolation algorithm will not crash
th_means[which(th_means %in% NA)] <- min(mean_thd80)
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
thvecPFvar3[which(values(therm_pred) %in% NA)] <- NA
thvecPFvar5[which(values(therm_pred) %in% NA)] <- NA

# initializing raster for the stored values
th_pfa_var3 <- therm_err
th_pfa_var5 <- therm_err

# updating values of the raster
# note: setValues() did not work
values(th_pfa_var3) <- thvecPFvar3
values(th_pfa_var5) <- thvecPFvar5

# saving rasters and making maps of variance
saveRast(rast=th_pfa_var3
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var3.tif')
makeMap (rast=th_pfa_var3
         ,plotnm='th_pfa_var3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=th_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var5.tif')
makeMap (rast=th_pfa_var5
         ,plotnm='th_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# saving rasters and making maps of standard deviation
saveRast(rast=calc(th_pfa_var3,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='th_pfa_sd3.tif')
makeMap (rast=calc(th_pfa_var3,fun=sqrt)
         ,plotnm='th_pfa_sd3.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(th_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='th_pfa_sd5.tif')
makeMap (rast=calc(th_pfa_var5,fun=sqrt)
         ,plotnm='th_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# deleting unneeded variables
rm(th_ses,th_means)
rm(thvecPFvar3,thvecPFvar5)
rm(therm_thresh3,therm_thresh5)
rm(std_thd80,mean_thd80)

##### RESERVOIR ######
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

rm(tf_rast1,tf_rast2,tf_rast3,tf_rast4,tf_rast5,tf_rast6)
rm(err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6)

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

# making uncertainty map
# reading-in tables for interpolated values
re_interp_tab3 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/re_pfvar3.xlsx',1,header=FALSE))
re_interp_tab5 <- as.matrix(read.xlsx('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/re_pfvar5.xlsx',1,header=FALSE))

# making the uncertainty maps
re_means <- log(values(res_pred_max2))
re_uncer <- values(res_pred_max_err)

# ranges for mean and cv
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

revecPFvar3_2 <- revecPFvar3
revecPFvar5_2 <- revecPFvar5

# setting values back to NAs
revecPFvar3[which(values(res_pred_max2) %in% NA)] <- NA
revecPFvar5[which(values(res_pred_max2) %in% NA)] <- NA

revecPFvar3[which(re_means %in% -Inf)] <- 0
revecPFvar5[which(re_means %in% -Inf)] <- 0

# initializing raster for the stored values
re_pfa_var3 <- res_pred_max_err
re_pfa_var5 <- res_pred_max_err

# updating values of the raster
# note: setValues() did not work
values(re_pfa_var3) <- revecPFvar3
values(re_pfa_var5) <- revecPFvar5

# rm(res_pred,res_err,res_pred_max,res_pred_max2
#    ,res_pred_1015,res_pred_1520,res_pred_2025,res_pred_2530,res_pred_3035,res_pred_3540
#    ,res_err_1015,res_err_1520,res_err_2025,res_err_2530,res_err_3035,res_err_3540
#    ,revecPFvar3,revecPFvar5,re_means,re_uncer)
# rm(res_pred_l1000,res_pred_max_err,res_max,res_min,res_thresh3,res_thresh5)
rm(revecPFvar3_2,revecPFvar5_2)
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
rm(uncer_re,mean_re)

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


rm(util_max,util_min,util_thresh3,util_thresh5)

##### SEISMIC ######
seis_eq_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kmp-')
seis_eq_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/EarthquakeBased/eqrisk_2kme-')

seis_stress_pred <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2p-')
seis_stress_err <- raster('/Users/calvinwhealton/GitHub/geothermal/combining_metrics/Seismic/StressFieldBased - Use 2km/stressrisk2e-')

seis_eq_pred[(seis_eq_pred %in% -9999)] <- NA
seis_stress_pred[(seis_stress_pred %in% -9999)] <- NA

seis_eq_err[(seis_eq_err %in% -9999)] <- NA
seis_stress_err[(seis_stress_err %in% -9999)] <- NA

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

se_eq_sds[which(se_eq_sds%in% 3000)] <-  2550
se_eq_sds[which(se_eq_sds %in% 2750)] <-  2515



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

saveRast(rast=calc(seEq_pfa_var5,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='seEq_pfa_sd5.tif')
makeMap (rast=calc(seEq_pfa_var5,fun=sqrt)
         ,plotnm='seEq_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

rm(seis_eq_max,seis_eq_min,seis_eq_thresh5,seis_eq_thresh3)
rm(seis_stress_max,seis_stress_min,seis_stress_thresh5,seis_stress_thresh3)
#rm(seis_eq_err,seis_eq_pred,seis_stress_err,seis_stress_pred)
rm(seEqvecPFvar3,seEqvecPFvar5,seStvecPFvar3,seStvecPFvar5)
rm(se_stress_sds,se_stress_means,se_eq_means,se_eq_sds)
rm(mean_seEq,mean_seSt,std_seEq,std_seSt)

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
         ,title='')
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
         ,title='')
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
         ,title='')
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
         ,title='')
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
         ,yloc=-1.2
         ,yshift=0
         ,title='')
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
         ,yloc=-0.8
         ,yshift=0
         ,title='')
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

# sums
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
         ,title='')
saveRast(rast=co_3_0_9_s_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_9_s_geo.tif')
makeMap(rast=co_3_0_9_s_geo
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
         ,plotnm='co_5_0_15_s_geo_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_15_s_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_15_s_geo.tif')
makeMap(rast=co_5_0_15_s_geo
        ,plotnm='co_5_0_15_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=3)

# product
co_3_0_27_p_geo <- calc(comb_pfa3_geo,fun=prod,na.rm=FALSE)
co_5_0_125_p_geo <- calc(comb_pfa5_geo,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_27_p_geo
         ,thresh3=c(0,1,2,3)^3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_27_p_geo_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_27_p_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_27_p_geo.tif')
makeMap(rast=co_3_0_27_p_geo
        ,plotnm='co_3_0_27_p_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_125_p_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_125_p_geo_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_125_p_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_125_p_geo.tif')
makeMap(rast=co_5_0_125_p_geo
        ,plotnm='co_5_0_125_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=3)

# minimum
co_3_0_3_m_geo <- calc(comb_pfa3_geo,fun=min,na.rm=FALSE)
co_5_0_5_m_geo <- calc(comb_pfa5_geo,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m_geo
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_geo_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_3_m_geo
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m_geo.tif')
makeMap(rast=co_3_0_3_m_geo
        ,plotnm='co_3_0_3_m_geo.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_5_m_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_geo_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_5_m_geo
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m_geo.tif')
makeMap(rast=co_5_0_5_m_geo
        ,plotnm='co_5_0_5_m_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3)

# uncertainty map for geology
co_uncer_geo_s3 <- stack(re_pfa_var3,th_pfa_var3,se_pfa_var3)
co_uncer_geo_s5 <- stack(re_pfa_var5,th_pfa_var5,se_pfa_var5)

co_pfa_var3_s_geo <- calc(co_uncer_geo_s3,fun=sum)
saveRast(rast=co_pfa_var3_s_geo
         ,wd=wd_raster
         ,rastnm='co_pfa_var3_s_geo.tif')
makeMap(rast=co_pfa_var3_s_geo
        ,plotnm='co_pfa_var3_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)
saveRast(rast=calc(co_pfa_var3_s_geo,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_pfa_sd3_s_geo.tif')
makeMap(rast=calc(co_pfa_var3_s_geo,fun=sqrt)
        ,plotnm='co_pfa_sd3_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)


co_pfa_var5_s_geo <- calc(co_uncer_geo_s5,fun=sum)
saveRast(rast=co_pfa_var5_s_geo
         ,wd=wd_raster
         ,rastnm='co_pfa_var5_s_geo.tif')
makeMap(rast=co_pfa_var5_s_geo
        ,plotnm='co_pfa_var5_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)
saveRast(rast=calc(co_pfa_var5_s_geo,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_pfa_sd5_s_geo.tif')
makeMap(rast=calc(co_pfa_var5_s_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_s_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

co_uncer_geo_pt3_1 <- calc(stack(c(re_3_0_3_NA,re_3_0_3_NA,th_3_0_3_NA,th_3_0_3_NA,se_pfa_var3)),fun=prod)
co_uncer_geo_pt3_2 <- calc(stack(c(re_3_0_3_NA,re_3_0_3_NA,th_pfa_var3,se_3_0_3_a,se_3_0_3_a)), fun=prod)
co_uncer_geo_pt3_3 <- calc(stack(c(re_pfa_var3,th_3_0_3_NA,th_3_0_3_NA,se_3_0_3_a,se_3_0_3_a)), fun=prod)

co_pfa_var3_p_geo <- calc(stack(c(co_uncer_geo_pt3_1,co_uncer_geo_pt3_2,co_uncer_geo_pt3_3)),fun=sum)

co_uncer_geo_pt5_1 <- calc(stack(c(re_5_0_5_NA,re_5_0_5_NA,th_5_0_5_NA,th_5_0_5_NA,se_pfa_var5)),fun=prod)
co_uncer_geo_pt5_2 <- calc(stack(c(re_5_0_5_NA,re_5_0_5_NA,th_pfa_var5,se_5_0_5_a,se_5_0_5_a)), fun=prod)
co_uncer_geo_pt5_3 <- calc(stack(c(re_pfa_var5,th_5_0_5_NA,th_5_0_5_NA,se_5_0_5_a,se_5_0_5_a)), fun=prod)

co_pfa_var5_p_geo <- calc(stack(c(co_uncer_geo_pt5_1,co_uncer_geo_pt5_2,co_uncer_geo_pt5_3)),fun=sum)

saveRast(rast=co_pfa_var3_p_geo
         ,wd=wd_raster
         ,rastnm='co_pfa_var3_p_geo.tif')
makeMap(rast=co_pfa_var3_p_geo
        ,plotnm='co_pfa_var3_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)
saveRast(rast=calc(co_pfa_var3_p_geo,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_pfa_sd3_p_geo.tif')
makeMap(rast=calc(co_pfa_var3_p_geo,fun=sqrt)
        ,plotnm='co_pfa_sd3_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)


saveRast(rast=co_pfa_var5_p_geo
         ,wd=wd_raster
         ,rastnm='co_pfa_var5_p_geo.tif')
makeMap(rast=co_pfa_var5_p_geo
        ,plotnm='co_pfa_var5_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)
saveRast(rast=calc(co_pfa_var5_p_geo,fun=sqrt)
         ,wd=wd_raster
         ,rastnm='co_pfa_sd5_p_geo.tif')
makeMap(rast=calc(co_pfa_var5_p_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_p_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

# no reservoirs
comb_pfa3_egs <- stack(c(ut5_3_0_3_NA,th_3_0_3_NA,se_3_0_3_a))
comb_pfa5_egs <- stack(c(ut5_5_0_5_NA,th_5_0_5_NA,se_5_0_5_a))

# sums
co_3_0_9_s_egs <- calc(comb_pfa3_egs,fun=sum,na.rm=FALSE)
co_5_0_15_s_egs <- calc(comb_pfa5_egs,fun=sum,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_9_s_egs
         ,thresh3=c(0,1,2,3)*3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_9_s_egs_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_9_s_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_9_s_egs.tif')
makeMap(rast=co_3_0_9_s_egs
        ,plotnm='co_3_0_9_s_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='sum'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_15_s_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)*3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_15_s_egs_hist.png'
         ,yloc=-0.04
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_15_s_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_15_s_egs.tif')
makeMap(rast=co_5_0_15_s_egs
        ,plotnm='co_5_0_15_s_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='sum'
        ,numRF=3)

# product
co_3_0_27_p_egs <- calc(comb_pfa3_egs,fun=prod,na.rm=FALSE)
co_5_0_125_p_egs <- calc(comb_pfa5_egs,fun=prod,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_27_p_egs
         ,thresh3=c(0,1,2,3)^3
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_27_p_egs_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_27_p_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_27_p_egs.tif')
makeMap(rast=co_3_0_27_p_egs
        ,plotnm='co_3_0_27_p_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='prod'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_125_p_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)^3
         ,rev_sc=FALSE
         ,plotnm='co_5_0_125_p_egs_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_125_p_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_125_p_egs.tif')
makeMap(rast=co_5_0_125_p_egs
        ,plotnm='co_5_0_125_p_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='prod'
        ,numRF=3)

# minimum
co_3_0_3_m_egs <- calc(comb_pfa3_egs,fun=min,na.rm=FALSE)
co_5_0_5_m_egs <- calc(comb_pfa5_egs,fun=min,na.rm=FALSE)

setwd(wd_image)
makeHist(rast=co_3_0_3_m_egs
         ,thresh3=c(0,1,2,3)
         ,thresh5=c()
         ,rev_sc=FALSE
         ,plotnm='co_3_0_3_m_egs_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='')
saveRast(rast=co_3_0_3_m_egs
         ,wd=wd_raster
         ,rastnm='co_3_0_3_m_egs.tif')
makeMap(rast=co_3_0_3_m_egs
        ,plotnm='co_3_0_3_m_egs.png'
        ,wd=wd_image
        ,numCol=3
        ,comTy='min'
        ,numRF=3)


setwd(wd_image)
makeHist(rast=co_5_0_5_m_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_egs_hist.png'
         ,yloc=-0.01
         ,yshift=0
         ,title='')
saveRast(rast=co_5_0_5_m_egs
         ,wd=wd_raster
         ,rastnm='co_5_0_5_m_egs.tif')
makeMap(rast=co_5_0_5_m_geo
        ,plotnm='co_5_0_5_m_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3)


########################
# extracting values of layers for cities
reservoir <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/re_5_0_5_NA.tif')
reservoir[reservoir < 0] <- NA

thermal <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/th_5_0_5_NA.tif')
thermal[thermal < 0] <- NA

seismic <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/se_5_0_5_a.tif')
seismic[seismic  < 0] <- NA

utilization <- raster('/Users/calvinwhealton/Dropbox/PFA_Rasters/ut5_5_0_5_NA.tif')
utilization[utilization  < 0] <- NA

comb_extract <- stack(c(reservoir,thermal,seismic,utilization
                        ,re_pfa_var5,th_pfa_var5,se_pfa_var5
                        ,res_pred_max2,res_pred_max_err
                        ,therm_pred,therm_err
                        ,seis_eq_pred,seis_eq_err,seis_stress_pred,seis_stress_err
                        ,util_pred
                        ,co_5_0_20_s,co_5_0_625_p,co_5_0_5_m
                        ,co_5_0_15_s_geo,co_5_0_125_p_geo,co_5_0_5_m_geo))

comb_names <- c('reservoir','thermal','seismic','utilization'
                ,'re_pfa_var5','th_pfa_var5','se_pfa_var5'
                ,'res_pred_max2','res_pred_max_err'
                ,'therm_pred','therm_err'
                ,'seis_eq_pred','seis_eq_err','seis_stress_pred','seis_stress_err'
                ,'util_pred'
                ,'co_5_0_20_s','co_5_0_625_p','co_5_0_5_m'
                ,'co_5_0_15_s_geo','co_5_0_125_p_geo','co_5_0_5_m_geo')

trialRun <- rasterToPoints(comb_extract,spatial=TRUE)
trialRun_df <- as.data.frame(trialRun)

# renaming the columns
names(trialRun_df) <- c('x','y',comb_names)

# dropping points that are NA for geology or all four risk factors
trialRun2_df <- trialRun_df[setdiff(seq(1,nrow(trialRun_df),1),intersect(which(trialRun_df$co_5_0_125_p %in% NA),which(trialRun_df$co_5_0_625_p %in% NA))),]

# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2

# calculating distance to each of the key points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  trialRun2_df[nm] <- sqrt((trialRun2_df$x-points$x[i])^2 + (trialRun2_df$y-points$y[i])^2)
}

points[c('x2','y2',comb_names)] <- NA

# extracting the balues corresponding to the maximum
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(trialRun2_df[nm] < 10000)
  
  if(length(inds) != 0){
    ind_max <- which(trialRun2_df$co_5_0_20_s[inds] %in% max(trialRun2_df$co_5_0_20_s[inds]))
    
    points[i,c('x','y',comb_names)] <- trialRun2_df[inds[ind_max],seq(1,length(c('x','y',comb_names)),1)]
    
  }
  
}

# making parallel axis plot
# setting exporting parameters
cols <- brewer.pal(10,'Set1')

setwd(wd_image)
png('parallel_axis.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,3)
     ,ylim=c(0,5)
     ,xaxs='i'
     ,yaxs='i'
     ,ylab='Play Fairway Metric'
     ,xlab=''
     ,xaxt='n')

lines(c(1,1)
     ,c(0,5)
     ,col='black')
lines(c(2,2)
      ,c(0,5)
      ,col='black')
j=1
names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

for(i in 1:nrow(points)){
  
  check <- c(points$reservoir[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[j])
    j <- j + 1
    names <- c(names,points$names[i])
  }  
}

par(xpd=TRUE)
text(c(0,1,2,3)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic','Utilization')
     ,adj=0.5)
text(-0.35,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.35,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=3.2,y=5
       ,legend=names
       ,col=cols
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()

## making boxplots/violin plots for distriutions from MC
ind_use <- c(1,2,4,5,6,7,8,9,10)

dists <- matrix(0,10000,9)
dist_vars <- matrix(0,4,9)

for(i in 1:length(ind_use)){
  
  set.seed(10)
  
  mat_mc <- matrix(0,10000,4)
  
  # reservoir MC
  pfm5 <- rep(0,10000)
  
  res_mean <- log(points$res_pred_max2[ind_use[i]])
  res_cv <- points$res_pred_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(10000,mu,sqrt(sigma2))

  pfm5[rand < log(3*10^-5)] <-5
  pfm5[rand > log(301)] <- 0
  pfm5[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.01)))] <- 5- (rand[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.01)))] - log(3*10^-5))/(log(0.01)-log(3*10^-5))
  pfm5[intersect(which(rand >= log(0.01)),which(rand < log(0.1)))] <- 4- (rand[intersect(which(rand >= log(0.01)),which(rand < log(0.1)))] - log(0.01))/(log(0.1)-log(0.01))
  pfm5[intersect(which(rand >= log(0.1)),which(rand < log(1)))] <- 3- (rand[intersect(which(rand >= log(0.1)),which(rand < log(1)))] - log(0.1))/(log(1)-log(0.1))
  pfm5[intersect(which(rand >= log(1)),which(rand < log(10)))] <- 2- (rand[intersect(which(rand >= log(1)),which(rand < log(10)))] - log(1))/(log(10)-log(1))
  pfm5[intersect(which(rand >= log(10)),which(rand < log(301)))] <- 1- (rand[intersect(which(rand >= log(10)),which(rand < log(301)))] - log(10))/(log(301)-log(10))
  
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal MC
  pfm5 <- rep(0,10000)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(10000,th_mean,th_se)
  
  pfm5[rand < 500] <- 5
  pfm5[rand > 8750] <- 0
  pfm5[intersect(which(rand >= 500),which(rand < 1500))] <- 5 - (rand[intersect(which(rand >= 500),which(rand < 1500))] - 500)/(1500-500)
  pfm5[intersect(which(rand >= 1500),which(rand < 2300))] <- 4- (rand[intersect(which(rand >= 1500),which(rand < 2300))] - 1500)/(2300-1500)
  pfm5[intersect(which(rand >= 2300),which(rand < 3000))] <- 3- (rand[intersect(which(rand >= 2300),which(rand < 3000))] - 2300)/(3000-2300)
  pfm5[intersect(which(rand >= 3000),which(rand < 4000))] <- 2- (rand[intersect(which(rand >= 3000),which(rand < 4000))] - 3000)/(4000-3000)
  pfm5[intersect(which(rand >= 4000),which(rand < 8750))] <- 1- (rand[intersect(which(rand >= 4000),which(rand < 8750))] - 4000)/(8750-3000)
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,10000)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(10000,se_eq_mean,se_eq_se)
  
  pfm5[rand < 1] <-5
  pfm5[rand > 25000] <-0
  pfm5[intersect(which(rand >= 1),which(rand < 5000))] <- 5- (rand[intersect(which(rand >= 1),which(rand < 5000))] - 1)/(5000-1)
  pfm5[intersect(which(rand >= 5000),which(rand < 10000))] <- 4- (rand[intersect(which(rand >= 5000),which(rand < 10000))] - 5000)/5000
  pfm5[intersect(which(rand >= 10000),which(rand < 15000))] <- 3- (rand[intersect(which(rand >= 10000),which(rand < 15000))] - 10000)/5000
  pfm5[intersect(which(rand >= 15000),which(rand < 20000))] <- 2- (rand[intersect(which(rand >= 15000),which(rand < 20000))] - 15000)/5000
  pfm5[intersect(which(rand >= 20000),which(rand < 25000))] <- 1- (rand[intersect(which(rand >= 20000),which(rand < 25000))] - 20000)/5000
  
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,10000)
  
  se_st_mean <- points$seis_stress_pred[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- abs(rnorm(10000,se_st_mean,se_st_se))
  
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  pfm5[rand < 0.01] <-5
  pfm5[intersect(which(rand >= 0.01),which(rand < 5))] <- 5- (rand[intersect(which(rand >= 0.01),which(rand < 5))] - 0.01)/(5-0.01)
  pfm5[intersect(which(rand >= 5),which(rand < 10))] <- 4- (rand[intersect(which(rand >= 5),which(rand < 10))] - 5)/5
  pfm5[intersect(which(rand >= 10),which(rand < 15))] <- 3- (rand[intersect(which(rand >= 10),which(rand < 15))] - 10)/5
  pfm5[intersect(which(rand >= 15),which(rand < 20))] <- 2- (rand[intersect(which(rand >= 15),which(rand < 20))] - 15)/5
  pfm5[intersect(which(rand >= 20),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 20),which(rand < 25))] - 20)/5
  
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  dists[,i] <- mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]
  
  print(mean(mat_mc[,1]))
  print(mean(mat_mc[,2]))
  print(mean(mat_mc[,3]))
  print(mean(mat_mc[,4]))
  
}

# making boxplots
setwd(wd_image)
png('boxplot.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,15)
        ,ylab='Combined Metric: Sum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(dists)){
  
  boxplot(dists[,i]
          ,at=i
          ,col='skyblue'
          ,add=TRUE
          ,pars=list(outcex=0.5,outpch=19)
  )
  
}
par(xpd=TRUE)
text(seq(1,9,1)
     ,rep(-0.5,9)
     ,names
     ,adj=0
     ,srt=270)
text(-1.2,15
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-1.2,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
par(xpd=FALSE)
dev.off()


# making violin plots
library(vioplot)
dists2 <- as.data.frame(dists)

setwd(wd_image)
png('violin.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
  
  vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
         ,col='skyblue'
         #,xlim=c(0,10)
         #,xlab=''
         #,xaxt='n'
         #,ylab='Combined Metric: Sum'
         ,names=rep('',9)
         ,ylim=c(0,15)
  )

par(xpd=TRUE)
rect(-1,-1,0.12,17
     ,col='white'
     ,border='white')
rect(0,-2,10,-0.64
     ,col='white'
     ,border='white')
text(seq(1,9,1)
     ,rep(-0.7,9)
     ,names
     ,adj=0
     ,srt=270)

text(-1
     ,7.5
     ,"Combined Metric:Sum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,14,2))
par(xpd=FALSE)


dev.off()



# extracting values of layers for cities
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA,se_5_0_5_a,ut5_5_0_5_NA
                                ,co_5_0_20_s,co_5_0_625_p,co_5_0_5_m))

combined_all_rasters[combined_all_rasters<0] <- NA
extract_pts <- extract(x=combined_all_rasters
                       ,y=cities2
                       ,sp=TRUE
                       ,nl=7
                       ,df=TRUE
                       ,na.rm=TRUE
                       ,method='simple'
                       ,buffer=2000
                       ,fun=mean
)

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[15] <- "th_5"
names(extract_pts2)[16] <- "re_5"
names(extract_pts2)[17] <- "se_5"
names(extract_pts2)[18] <- "ut_5"
names(extract_pts2)[19] <- "co_s_5"
names(extract_pts2)[20] <- "co_p_5"
names(extract_pts2)[21] <- "co_m_5"


extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]

ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5 > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5 > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
png('scatter_p_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]
     ,extract_pts3$co_p_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Sum'
     ,ylab='Combined Product'
)
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Sum'
     ,ylab='Combined Minimum'
)
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_p.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_p_5[inds_gtr0]
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlab='Combined Product'
     ,ylab='Combined Minimum'
)
dev.off()

