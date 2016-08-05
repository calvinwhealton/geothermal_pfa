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
# --saving output for the combined output
# --completing more detailed comparisons for specific locations
# --parallel axis plot
# --boxplot with Monte Carlo distributions
# --violin plot with Monte Carlo distributions
# --scatterplots of different metrics
# saving all output

# modifications when running on a different machine:
# -- change all working directories

##### importing functions/libraries #####

# libraries
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
library(xlsx)         # for reading in xlsx tables
library(rgeos)        # for buffering places of interest
library(RColorBrewer) # R color brewer palettes for parallel axis plot
library(pracma)       # for interpolation in tables
library(vioplot)      # for violin plots
#library(GISTools)     # adding stuff to plots
library(prettymapr)
library(maps)
##### defining working directories #####
# need to be changed based on machine

# location to SAVE rasters
wd_raster_out <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/Rasters_out'

# location to READ rasters
wd_raster_in <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/Rasters_in'

# location to save images
wd_image <- '/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics/figs_paper'

# location of code/scripts are stored
wd_code <- '/Users/calvinwhealton/GitHub/geothermal_pfa/combining_metrics'

# location of input shapefiles
wd_shapefiles <- '/Users/calvinwhealton/GitHub/geothermal_pfa/GISData/'

# location of error interoplation tables
wd_error_interp <- '/Users/calvinwhealton/GitHub/geothermal_pfa/error_interp_tabs'

# location to save workspace
wd_workspace <- '/Users/calvinwhealton/GitHub/geothermal_pfa/data'

##### loading-in state/county shapefiles #####
# states
States = readOGR(dsn=paste(wd_shapefiles,'us_state_WGS',sep='')
                 , layer="us_state_WGS"
                 , stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# counties
Counties = readOGR(dsn=paste(wd_shapefiles,'us_county_WGS84_prj',sep='')
                   , layer="us_county_WGS84_prj"
                   , stringsAsFactors=FALSE)
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
cities <- readOGR(dsn=paste(wd_shapefiles,'usCensusPlaces',sep='')
                  , layer="usCensusPlaces"
                  , stringsAsFactors=FALSE)
cities2 <- spTransform(cities,CRS("+init=epsg:31986"))

# removing the untransformed shapefiles
rm(NY,PA,WV,States,Counties,NY_co,PA_co,WV_co)

## importing places of interest
poi <- readOGR(dsn=paste(wd_shapefiles,'placesofinterest2',sep='')
               , layer="placesofinterest2"
               , stringsAsFactors=FALSE)

# converting to UTM 17 N system
poi2 <- spTransform(poi,CRS("+init=epsg:31986"))

# buffering places of interest by 5000 and 10000 m
poi2Buf5 <- gBuffer(poi2,width=5000)
poi2Buf10 <- gBuffer(poi2,width=10000)

# shapefile of the US Census places merged
pla <- readOGR(dsn='/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/Utilization',layer='mergerIhopeThisWorksZachFile')
plaBuf5 <- gBuffer(pla,width=5000)
plaBuf10 <- gBuffer(pla,width=10000)
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
source('rw_functions.R')

##### THERMAL ######
# importing rasters, depth to 80 DegC and standard error of prediction
therm_pred <- raster(paste(wd_raster_in,'/Thermal/d80p2-',sep=''))
therm_err  <- raster(paste(wd_raster_in,'/Thermal/d80e2-',sep=''))

# values of -9999 are no data, so replaced with NA
therm_pred[therm_pred < 0] <- NA
therm_err[therm_err < 0] <- NA

# setting thresholds
# deep values are less favorable
#therm_thresh3 <- rev(c(8750,3000,2000,500))
therm_thresh5 <- rev(c(8750,4000,3000,2300,1500,500))

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
th_interp_tab5 <- as.matrix(read.xlsx(paste(wd_error_interp,'/th_d80_pfvar5.xlsx',sep=''),1,header=FALSE))
th_interp_tab5_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/th_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

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

# interpolating for the 5 color scheme
thvecPFvar5 <- interp2(x=std_thd80
                       ,y=mean_thd80
                       ,Z=th_interp_tab5
                       ,xp=th_ses
                       ,yp=th_means
                       ,method='linear')

thvecPFvar5_ls <- interp2(x=std_thd80
                       ,y=mean_thd80
                       ,Z=th_interp_tab5_ls
                       ,xp=th_ses
                       ,yp=th_means
                       ,method='linear')

# setting values back to NAs
#thvecPFvar3[which(values(therm_pred) %in% NA)] <- NA
thvecPFvar5[which(values(therm_pred) %in% NA)] <- NA
thvecPFvar5_ls[which(values(therm_pred) %in% NA)] <- NA

# initializing raster for the stored values
th_pfa_var5 <- therm_err
th_pfa_var5_ls <- therm_err

# updating values of the raster
# note: setValues() did not work
values(th_pfa_var5) <- thvecPFvar5
values(th_pfa_var5_ls) <- thvecPFvar5_ls

# saving rasters and making maps of variance
saveRast(rast=th_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var5.tif')
makeMap(rast=th_pfa_var5
         ,plotnm='th_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=th_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var5_ls.tif')
makeMap(rast=th_pfa_var5_ls
        ,plotnm='th_pfa_var5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# saving rasters and making maps of standard deviation
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

# saving rasters and making maps of standard deviation
saveRast(rast=calc(th_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='th_pfa_sd5_ls.tif')
makeMap (rast=calc(th_pfa_var5_ls,fun=sqrt)
         ,plotnm='th_pfa_sd5_ls.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)


# deleting unneeded variables
rm(th_ses,th_means,thvecPFvar5,std_thd80,mean_thd80)

##### RESERVOIR ######
# l1000 = less than 1000 m
# 1015 = 1000 to 1500 m
# 1520 = 1500 to 2000 m, etc.
res_pred_l1000 <- raster(paste(wd_raster_in,'/Reservoirs/fff_l1000p',sep=''))
res_pred_1015 <- raster(paste(wd_raster_in,'/Reservoirs/fff_1015p2',sep=''))
res_pred_1520 <- raster(paste(wd_raster_in,'/Reservoirs/fff_1520p',sep=''))
res_pred_2025 <- raster(paste(wd_raster_in,'/Reservoirs/fff_2025p',sep=''))
res_pred_2530 <- raster(paste(wd_raster_in,'/Reservoirs/fff_2530p',sep=''))
res_pred_3035 <- raster(paste(wd_raster_in,'/Reservoirs/fff_3035p',sep=''))
res_pred_3540 <- raster(paste(wd_raster_in,'/Reservoirs/fff_3540p',sep=''))

# reservoir error
res_err_l1000 <- raster(paste(wd_raster_in,'/Reservoirs/fff_l1000e',sep=''))
res_err_1015 <- raster(paste(wd_raster_in,'/Reservoirs/fff_1015e',sep=''))
res_err_1520 <- raster(paste(wd_raster_in,'/Reservoirs/fff_1520e',sep=''))
res_err_2025 <- raster(paste(wd_raster_in,'/Reservoirs/fff_2025e',sep=''))
res_err_2530 <- raster(paste(wd_raster_in,'/Reservoirs/fff_2530e',sep=''))
res_err_3035 <- raster(paste(wd_raster_in,'/Reservoirs/fff_3035e',sep=''))
res_err_3540 <- raster(paste(wd_raster_in,'/Reservoirs/fff_3540e',sep=''))

# stacking rasters and taking maximum, assuming going for highest quality reservoir
# ignoring reservoirs shallower than 1000 m (the res_pred_l1000 raster)
res_pred <- stack(c(res_pred_1015,res_pred_1520,res_pred_2025,res_pred_2530,res_pred_3035,res_pred_3540))
res_pred_max <- calc(res_pred,fun=max,na.rm=TRUE)

# setting values less than zero (-9999) to NA because they have no data
res_pred_max2 <- res_pred_max
res_pred_max2[res_pred_max <0] <- NA

# stacking reservoir errors
res_err <- stack(c(res_err_1015,res_err_1520,res_err_2025,res_err_2530,res_err_3035,res_err_3540))

## making reservoir error term
# step 1: raster of whether the max value is equal to the value in the layer
tf_rast1 <- (res_pred_max == res_pred[[1]])
tf_rast2 <- (res_pred_max == res_pred[[2]])
tf_rast3 <- (res_pred_max == res_pred[[3]])
tf_rast4 <- (res_pred_max == res_pred[[4]])
tf_rast5 <- (res_pred_max == res_pred[[5]])
tf_rast6 <- (res_pred_max == res_pred[[6]])

# step 2: stacking raster with error raster
# product means only values where the max was matched will be assigned a positive value
# otherwise, it will be 0
err_tf1 <- calc(stack(c(tf_rast1,res_err[[1]])),fun=prod)
err_tf2 <- calc(stack(c(tf_rast2,res_err[[2]])),fun=prod)
err_tf3 <- calc(stack(c(tf_rast3,res_err[[3]])),fun=prod)
err_tf4 <- calc(stack(c(tf_rast4,res_err[[4]])),fun=prod)
err_tf5 <- calc(stack(c(tf_rast5,res_err[[5]])),fun=prod)
err_tf6 <- calc(stack(c(tf_rast6,res_err[[6]])),fun=prod)

# step 3: making stacked raster and taking maximum
# only values where the raster matched the maximum of the rasters will be used
res_pred_max_err <- calc(stack(c(err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6)),fun=max)
res_pred_max_err[(res_pred_max_err < 0)] <- NA

# thresholds
res_min <- 3*10^-5 # minimum for reservoir (both 3-color and 5-color)
res_max <- 301     # maximum for reservoir (both 3-color and 5-color)
#res_thresh3 <- c(res_min,c(0.1,1.0),res_max)
res_thresh5 <- c(res_min,c(0.01,0.1,1.0,10),res_max)

# deleting unneeded rasters
rm(tf_rast1,tf_rast2,tf_rast3,tf_rast4,tf_rast5,tf_rast6)
rm(err_tf1,err_tf2,err_tf3,err_tf4,err_tf5,err_tf6)

# converting into the play fairway scheme
# five color-----
re_5_0_5_NA <- convRastPFRank(rast=calc(res_pred_max2,fun=log10)
                              ,thresholds=log10(res_thresh5)
                              ,ignore=-9999
                              ,rev_scale=FALSE)
saveRast(rast=re_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='re_5_0_5_NA.tif')
makeMap(rast=re_5_0_5_NA
         ,plotnm='re_5_0_5_NA.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)

# making uncertainty map
# reading-in tables for interpolated values
re_interp_tab5 <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_pfvar5.xlsx',sep=''),1,header=FALSE))
re_interp_tab5_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

# making the uncertainty maps
# values are in base-e (ln) to match calculations in make_interp_tabs
re_means <- log(values(res_pred_max2))
re_uncer <- values(res_pred_max_err)

# ranges for mean and cv
# must check values from make_interp_table.R
mean_re <- seq(-9.5,6.25,0.25) # range of means
uncer_re <- seq(0,2,by=0.1) # range of coefficient of variation

# substituting in values for NAs so interpolation algorithm will not crash
re_means[which(re_means %in% NA)] <- min(mean_re)
re_uncer[which(re_uncer %in% NA)] <- min(uncer_re)

re_means[which(re_means %in% 0)] <- min(mean_re)
re_uncer[which(re_uncer %in% 0)] <- min(uncer_re)

# interpolating for the 5 color scheme
revecPFvar5 <- interp2(x=uncer_re
                       ,y=mean_re
                       ,Z=re_interp_tab5
                       ,xp=re_uncer
                       ,yp=re_means
                       ,method='linear')

revecPFvar5_ls <- interp2(x=uncer_re
                       ,y=mean_re
                       ,Z=re_interp_tab5_ls
                       ,xp=re_uncer
                       ,yp=re_means
                       ,method='linear')

# dummy variables
revecPFvar5_2 <- revecPFvar5
revecPFvar5_2_ls <- revecPFvar5_ls

# setting values back to NAs
revecPFvar5[which(values(res_pred_max2) %in% NA)] <- NA
revecPFvar5_ls[which(values(res_pred_max2) %in% NA)] <- NA

# converting any values that were -Inf to zero
# -Inf result from zero mean value reservoirs
revecPFvar5[which(re_means %in% -Inf)] <- 0
revecPFvar5_ls[which(re_means %in% -Inf)] <- 0

# initializing raster for the stored values
re_pfa_var5 <- res_pred_max_err
re_pfa_var5_ls <- res_pred_max_err

# updating values of the raster
# note: setValues() did not work
values(re_pfa_var5) <- revecPFvar5
values(re_pfa_var5_ls) <- revecPFvar5_ls

# saving rasters and making maps
saveRast(rast=re_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5.tif')
makeMap(rast=re_pfa_var5
         ,plotnm='re_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)


saveRast(rast=re_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_ls.tif')
makeMap(rast=re_pfa_var5_ls
        ,plotnm='re_pfa_var5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)


# saving rasters and making maps
saveRast(rast=calc(re_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5.tif')
makeMap(rast=calc(re_pfa_var5,fun=sqrt)
         ,plotnm='re_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(re_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_ls.tif')
makeMap(rast=calc(re_pfa_var5_ls,fun=sqrt)
        ,plotnm='re_pfa_sd5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# removing unneeded files
rm(uncer_re,mean_re)
rm(res_err_3540,res_err_3035,res_err_2530,res_err_1520,res_err_2025,res_err_1015)
rm(res_pred_l1000,res_pred_3540,res_pred_1520,res_pred_3035,res_pred_2530,res_pred_2025,res_pred_1015)
rm(res_pred_max,res_pred)
rm(revecPFvar5,revecPFvar5_2)
rm(res_err,res_err_l1000,re_means,re_uncer)

##### UTILIZATION ####
# utilization prediciton and error
util_pred <- raster(paste(wd_raster_in,'/Utilization/slcoh_f/slcoh_p4.tif',sep=''))

# replacing no data (-9999) with NA
util_pred[(util_pred %in% -9999)] <- NA
util_pred[(util_pred %in% 0)] <- NA

# making up a utilization error
util_err <- util_pred
util_err[is.na(util_err) == F] <- 5

  
# thresholds
util_min <- 5
util_max <- 25
#util_thresh3 <- c(util_min,c(13.5,16),util_max)
util_thresh5 <- c(util_min,c(12,13.5,16,20),util_max)

# creating buffer images
makeWeightBuf(dist=5
              ,wd=wd_image
              ,plotnm='ut_buf_5.png')

# converting into the play fairway scheme
# five color-----
ut0_5_0_5_NA <- convRastPFRank(rast=util_pred
                               ,thresholds=util_thresh5
                               ,ignore=-9999
                               ,rev_scale=TRUE)
saveRast(rast=ut0_5_0_5_NA 
         ,wd=wd_raster_out
         ,rastnm='ut0_5_0_5_NA.tif')
makeMap(rast=ut0_5_0_5_NA 
        ,plotnm='ut0_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F)

# buffering utilization (5 km)
# five color
ut5_5_0_5_NA <- focal(ut0_5_0_5_NA
                      ,w=makeUtilBufWeight(5)
                      ,fun=max
                      ,na.rm=TRUE
                      ,pad=TRUE)
saveRast(rast=ut5_5_0_5_NA 
         ,wd=wd_raster_out
         ,rastnm='ut5_5_0_5_NA.tif')
makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

# making uncertainty map
util_pred2 <- util_pred
util_pred2[which(values(util_pred) %in% NA)] <- 1000
util_pred2 <- calc(util_pred2,fun=function(x){2000-x})
utTemp <- focal(util_pred2
                ,w=makeUtilBufWeight(5)
                ,fun=max
                ,na.rm=TRUE
                ,pad=TRUE
                ,padValue=-2000)

util_means <- 2000 - values(utTemp) 
util_uncer <- rep(5,length(util_means))

# ranges for mean and % uncertainty
# must check values from make_interp_table.R
mean_util <- c(seq(5,65,by=2),900,1100) # range of means
uncer_util <- seq(1,10,by=0.5) # range of standard deviations

util_interp_tab5 <-  as.matrix(read.xlsx(paste(wd_error_interp,'/ut_slcoh_pfvar5.xlsx',sep=''),1,header=FALSE))
util_interp_tab5_ls <-  as.matrix(read.xlsx(paste(wd_error_interp,'/ut_slcoh_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

# interpolating for the 5 color scheme
utilvecPFvar5 <- interp2(x=uncer_util
                       ,y=mean_util
                       ,Z=util_interp_tab5
                       ,xp=util_uncer
                       ,yp=util_means
                       ,method='linear')

utilvecPFvar5_ls <- interp2(x=uncer_util
                         ,y=mean_util
                         ,Z=util_interp_tab5_ls
                         ,xp=util_uncer
                         ,yp=util_means
                         ,method='linear')

# setting values back to NAs
utilvecPFvar5[which(values(ut5_5_0_5_NA) %in% -Inf)] <- NA
utilvecPFvar5_ls[which(values(ut5_5_0_5_NA) %in% -Inf)] <- NA

# initializing raster for the stored values
# updating values of the raster
# note: setValues() did not work
util_pfa_var5 <- util_err
values(util_pfa_var5) <- utilvecPFvar5

util_pfa_var5_ls <- util_err
values(util_pfa_var5_ls) <- utilvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=util_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='util_pfa_var5.tif')
makeMap(rast=util_pfa_var5
         ,plotnm='util_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=util_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='util_pfa_var5_ls.tif')
makeMap(rast=util_pfa_var5_ls
         ,plotnm='util_pfa_var5_ls.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

# saving rasters and making maps
saveRast(rast=calc(util_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='util_pfa_sd5.tif')
makeMap(rast=calc(util_pfa_var5,fun=sqrt)
         ,plotnm='util_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(util_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='util_pfa_sd5_ls.tif')
makeMap(rast=calc(util_pfa_var5_ls,fun=sqrt)
        ,plotnm='util_pfa_sd5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# deleting unneeded files
rm(util_max,util_min,util_thresh5,ut0_5_0_5_NA)

##### SEISMIC ######
# reading-in the earthquake based risk rasters
seis_eq_pred <- raster(paste(wd_raster_in,'/Seismic/EarthquakeBased/eqrisk_2kmp-',sep=''))
seis_eq_err <- raster(paste(wd_raster_in,'/Seismic/EarthquakeBased/eqrisk_2kme-',sep=''))

# reading-in the stress field based risk rasters
seis_stress_pred <- raster(paste(wd_raster_in,'/Seismic/StressFieldBased - Use 2km/stressrisk2p-',sep=''))
seis_stress_err <- raster(paste(wd_raster_in,'/Seismic/StressFieldBased - Use 2km/stressrisk2e-',sep=''))

# setting no data (-9999) to NA
seis_eq_pred[(seis_eq_pred %in% -9999)] <- NA
seis_stress_pred[(seis_stress_pred %in% -9999)] <- NA

seis_eq_err[(seis_eq_err %in% -9999)] <- NA
seis_stress_err[(seis_stress_err %in% -9999)] <- NA

# thresholds for stress-based risk
seis_stress_min <- 0.001 # to avoid problems with numerically zero values
seis_stress_max <- 25
seis_stress_thresh5<- c(seis_stress_min,c(5,10,15,20),seis_stress_max)

# thresholds for earthquake-based risk
seis_eq_min <- 0.001 # to avoid problems with numerically zero values
seis_eq_max <- 25
seis_eq_thresh5<- 10^3*c(seis_eq_min,c(5,10,15,20),seis_eq_max)

# combining into the play fairway scheme EARTHQUAKES-----
# five color
seEq_5_0_5_NA <- convRastPFRank(rast=seis_eq_pred
                                ,thresholds=seis_eq_thresh5
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seEq_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='seEq_5_0_5_NA.tif')
makeMap(rast=seEq_5_0_5_NA
        ,plotnm='seEq_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

# combining into the play fairway scheme STRESS-----
# five color
seSt_5_0_5_NA <- convRastPFRank(rast=seis_stress_pred
                                ,thresholds=seis_stress_thresh5
                                ,ignore=-9999
                                ,rev_scale=FALSE)
saveRast(rast=seSt_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='seSt_5_0_5_NA.tif')
makeMap(rast=seSt_5_0_5_NA
        ,plotnm='seSt_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)


# combining stress and earthquakes into play fairway seismic map
se_5_0_5_a  <- calc(stack(c(seEq_5_0_5_NA,seSt_5_0_5_NA)),fun=mean)

# writing rasters
saveRast(rast=se_5_0_5_a
         ,wd=wd_raster_out
         ,rastnm='se_5_0_5_a.tif')
makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

## calcualting uncertainty maps for the play fairway scheme
# reading-in tables for interpolated values
se_stress_interp_tab5 <- as.matrix(read.xlsx(paste(wd_error_interp,'/se_stress_pfvar5.xlsx',sep=''),1,header=FALSE))
se_stress_interp_tab5_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/seSt_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

se_eq_interp_tab5 <-  as.matrix(read.xlsx(paste(wd_error_interp,'/se_eq_pfvar5.xlsx',sep=''),1,header=FALSE))
se_eq_interp_tab5_ls <-  as.matrix(read.xlsx(paste(wd_error_interp,'/seEq_pfvar5.xlsx',sep=''),1,header=FALSE))

## values to interpolate variance
se_stress_means <- values(seis_stress_pred)
se_stress_sds <- values(seis_stress_err)

se_stress_means[se_stress_means > 70] <- 70
se_stress_sds[se_stress_means == 70] <- 0

se_eq_means <- values(seis_eq_pred)
se_eq_sds <- values(seis_eq_err)

# correction because errors were added
# when they should have been squared, summed, and then square rooted
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
# interpolating for the 5 color scheme
seStvecPFvar5 <- interp2(x=std_seSt
                         ,y=mean_seSt
                         ,Z=se_stress_interp_tab5
                         ,xp=se_stress_sds
                         ,yp=se_stress_means
                         ,method='linear'
                         )

seStvecPFvar5_ls <- interp2(x=std_seSt
                         ,y=mean_seSt
                         ,Z=se_stress_interp_tab5_ls
                         ,xp=se_stress_sds
                         ,yp=se_stress_means
                         ,method='linear'
                         )

# setting values back to NAs
seStvecPFvar5[which(values(seis_stress_pred) %in% NA)] <- NA
seStvecPFvar5_ls[which(values(seis_stress_pred) %in% NA)] <- NA

seStvecPFvar5[which(values(seis_stress_pred) %in% 70)] <- 0
seStvecPFvar5_ls[which(values(seis_stress_pred) %in% 70)] <- 0

# initializing raster for the stored values
seSt_pfa_var5 <- seis_stress_err
seSt_pfa_var5_ls <- seis_stress_err

# updating values of the raster
# note: setValues() did not work
values(seSt_pfa_var5) <- seStvecPFvar5
values(seSt_pfa_var5_ls) <- seStvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=seSt_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='seSt_pfa_var5.tif')
makeMap(rast=seSt_pfa_var5
         ,plotnm='seSt_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=seSt_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='seSt_pfa_var5_ls.tif')
makeMap(rast=seSt_pfa_var5_ls
        ,plotnm='seSt_pfa_var5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# saving rasters and making maps
saveRast(rast=calc(seSt_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='seSt_pfa_sd5.tif')
makeMap(rast=calc(seSt_pfa_var5,fun=sqrt)
         ,plotnm='seSt_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(seSt_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='seSt_pfa_sd5_ls.tif')
makeMap(rast=calc(seSt_pfa_var5_ls,fun=sqrt)
        ,plotnm='seSt_pfa_sd5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# for EARTHQUAKE
# interpolating for the 5 color scheme
seEqvecPFvar5 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab5
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)

seEqvecPFvar5_ls <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab5_ls
                         ,xp=se_eq_sds
                         ,yp=se_eq_means
                         ,method='linear'
)

# setting values back to NAs
seEqvecPFvar5[which(values(seis_eq_pred) %in% NA)] <- NA
seEqvecPFvar5_ls[which(values(seis_eq_pred) %in% NA)] <- NA

# initializing raster for the stored values
seEq_pfa_var5 <- seis_eq_err
seEq_pfa_var5_ls <- seis_eq_err

# updating values of the raster
# note: setValues() did not work
values(seEq_pfa_var5) <- seEqvecPFvar5
values(seEq_pfa_var5_ls) <- seEqvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=seEq_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_var5.tif')
makeMap(rast=seEq_pfa_var5
         ,plotnm='seEq_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=seEq_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_var5_ls.tif')
makeMap(rast=seEq_pfa_var5_ls
        ,plotnm='seEq_pfa_var5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# saving rasters and making maps
saveRast(rast=calc(seEq_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_sd5.tif')
makeMap(rast=calc(seEq_pfa_var5,fun=sqrt)
         ,plotnm='seEq_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(seEq_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_sd5_ls.tif')
makeMap(rast=calc(seEq_pfa_var5_ls,fun=sqrt)
        ,plotnm='seEq_pfa_sd5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# for COMBINED
se_pfa_var5s <- stack(c(seEq_pfa_var5,seSt_pfa_var5))
se_pfa_var5 <- calc(se_pfa_var5s,fun=sum)

# 0.25 because in the averaging of the two each raster is multiplied 
# by 0.5, so the variance will be multiplied by 0.25
values(se_pfa_var5) <- values(se_pfa_var5)*0.25

# calculating the variance in log-space from real-space Taylor series
se_pfa_var5_temp1 <- calc(se_pfa_var5s,fun=sum)
se_pfa_var5_temp2 <- calc(stack(seSt_5_0_5_NA,seEq_5_0_5_NA),fun=sum)
se_pfa_var5_temp3 <- calc(se_pfa_var5_temp2,fun=function(x){return(x^(-2))})
se_pfa_var5_temp4 <- calc(se_pfa_var5_temp3,fun=function(x){return(ifelse(x > 0.2^2,0.2,x))})
se_pfa_var5_ls <- calc(stack(se_pfa_var5_temp4,se_pfa_var5_temp1),fun=prod)

# saving rasters and making maps
saveRast(rast=se_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='se_pfa_var5.tif')
makeMap(rast=se_pfa_var5
         ,plotnm='se_pfa_var5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=se_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='se_pfa_var5_ls.tif')
makeMap(rast=se_pfa_var5_ls
        ,plotnm='se_pfa_var5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)


# saving rasters and making maps
saveRast(rast=calc(se_pfa_var5,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='se_pfa_sd5.tif')
makeMap(rast=calc(se_pfa_var5,fun=sqrt)
         ,plotnm='se_pfa_sd5.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)

saveRast(rast=calc(se_pfa_var5_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='se_pfa_sd5_ls.tif')
makeMap(rast=calc(se_pfa_var5_ls,fun=sqrt)
        ,plotnm='se_pfa_sd5_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# removing unneeded files
# rm(seis_eq_max,seis_eq_min,seis_eq_thresh5)
# rm(seis_stress_max,seis_stress_min,seis_stress_thresh5)
# rm(seEqvecPFvar5,seStvecPFvar5)
# rm(se_stress_sds,se_stress_means,se_eq_means,se_eq_sds)
# rm(mean_seEq,mean_seSt,std_seEq,std_seSt)

##### making a quick and simple map to compare results #####
#### combining maps, all variables ###
# creating a stacked raster
#comb_pfa3 <- stack(c(re_3_0_3_NA,th_3_0_3_NA,ut5_3_0_3_NA,se_3_0_3_a))
re_5_0_5_NA <- raster(paste(wd_raster_out,'re_5_0_5_NA.tif',sep='/'))
th_5_0_5_NA <- raster(paste(wd_raster_out,'th_5_0_5_NA.tif',sep='/'))
ut5_5_0_5_NA <- raster(paste(wd_raster_out,'ut5_5_0_5_NA.tif',sep='/'))
se_5_0_5_a <- raster(paste(wd_raster_out,'se_5_0_5_a.tif',sep='/'))

re_pfa_var5 <- raster(paste(wd_raster_out,'re_pfa_var5.tif',sep='/'))
th_pfa_var5  <- raster(paste(wd_raster_out,'th_pfa_var5.tif',sep='/'))
util_pfa_var5  <- raster(paste(wd_raster_out,'util_pfa_var5.tif',sep='/'))
se_pfa_var5  <- raster(paste(wd_raster_out,'se_pfa_var5.tif',sep='/'))

re_pfa_var5_ls <- raster(paste(wd_raster_out,'re_pfa_var5_ls.tif',sep='/'))
th_pfa_var5_ls  <- raster(paste(wd_raster_out,'th_pfa_var5_ls.tif',sep='/'))
util_pfa_var5_ls  <- raster(paste(wd_raster_out,'util_pfa_var5_ls.tif',sep='/'))
se_pfa_var5_ls  <- raster(paste(wd_raster_out,'se_pfa_var5_ls.tif',sep='/'))

comb_pfa5 <- stack(c(re_5_0_5_NA,th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))
comb_check <- calc(comb_pfa5,fun=min)

comb_pfa5[comb_pfa5 < 0] <- NA

## using sums
co_5_0_20_s <- calc(comb_pfa5,fun=mean,na.rm=FALSE)

setwd(wd_image)
saveRast(rast=co_5_0_20_s
         ,wd=wd_raster_out
         ,rastnm='co_5_0_20_s.tif')
makeMap(rast=co_5_0_20_s
        ,plotnm='co_5_0_20_s.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4)

## using products
co_5_0_625_p2 <- calc(comb_pfa5,fun=prod,na.rm=FALSE)
co_5_0_625_p <- calc(co_5_0_625_p2,fun=function(x){return(x^0.25)})

setwd(wd_image)
saveRast(rast=co_5_0_625_p
         ,wd=wd_raster_out
         ,rastnm='co_5_0_625_p.tif')
makeMap(rast=co_5_0_625_p
        ,plotnm='co_5_0_625_p.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4)

## using minimums
co_5_0_5_m <- calc(comb_pfa5,fun=min,na.rm=FALSE)

setwd(wd_image)
saveRast(rast=co_5_0_5_m
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m.tif')
makeMap(rast=co_5_0_5_m
        ,plotnm='co_5_0_5_m.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=1)

# for average
uncer_temp <- stack(re_pfa_var5,th_pfa_var5,se_pfa_var5,util_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg <- calc(co_uncer_avg1,fun=function(x){return(x/16)})

# for geometric mean
uncer_temp <- stack(re_pfa_var5_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/16)})
co_pfa_var5_geomean <- calc(stack(co_5_0_625_p,co_5_0_625_p,co_uncer_geomean2),fun=prod)

# saving rasters and making maps
saveRast(rast=co_pfa_var5_avg
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg.tif')
makeMap(rast=co_pfa_var5_avg
        ,plotnm='co_pfa_var5_avg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(co_pfa_var5_avg,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg.tif')
makeMap(rast=calc(co_pfa_var5_avg,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

# saving rasters and making maps
saveRast(rast=co_pfa_var5_geomean
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_geomean.tif')
makeMap(rast=co_pfa_var5_geomean
        ,plotnm='co_pfa_var5_geomean.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(co_pfa_var5_geomean,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean.tif')
makeMap(rast=calc(co_pfa_var5_geomean,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

########################
# extracting values of layers for cities
# reading-in the scaled rasters
# note that scaled rasters are in the wd_raster_out
# because they were output from the first part of the code
reservoir <- raster(paste(wd_raster_out,'/re_5_0_5_NA.tif',sep=''))
reservoir[reservoir < 0] <- NA

re_pfa_var5 <- raster(paste(wd_raster_out,'/re_pfa_var5.tif',sep=''))
re_pfa_var5[re_pfa_var5 < 0] <- NA

re_pfa_var5_ls <- raster(paste(wd_raster_out,'/re_pfa_var5_ls.tif',sep=''))
re_pfa_var5_ls[re_pfa_var5_ls < 0] <- NA


thermal <- raster(paste(wd_raster_out,'/th_5_0_5_NA.tif',sep=''))
thermal[thermal < 0] <- NA

th_pfa_var5 <- raster(paste(wd_raster_out,'/th_pfa_var5.tif',sep=''))
th_pfa_var5[th_pfa_var5 < 0] <- NA

th_pfa_var5_ls <- raster(paste(wd_raster_out,'/th_pfa_var5_ls.tif',sep=''))
th_pfa_var5_ls[th_pfa_var5_ls < 0] <- NA

seismic <- raster(paste(wd_raster_out,'/se_5_0_5_a.tif',sep=''))
seismic[seismic  < 0] <- NA

se_pfa_var5 <- raster(paste(wd_raster_out,'/se_pfa_var5.tif',sep=''))
se_pfa_var5[se_pfa_var5 < 0] <- NA

se_pfa_var5_ls <- raster(paste(wd_raster_out,'/se_pfa_var5_ls.tif',sep=''))
se_pfa_var5_ls[se_pfa_var5_ls < 0] <- NA

utilization <- raster(paste(wd_raster_out,'/ut5_5_0_5_NA.tif',sep=''))
utilization[utilization  < 0] <- NA

util_pfa_var5 <- raster(paste(wd_raster_out,'/util_pfa_var5.tif',sep=''))
util_pfa_var5[util_pfa_var5 < 0] <- NA

util_pfa_var5_ls <- raster(paste(wd_raster_out,'/util_pfa_var5_ls.tif',sep=''))
util_pfa_var5_ls[util_pfa_var5_ls < 0] <- NA

util_pred2 <- calc(util_pred,fun=function(x){return(1000-x)})

util_pred3 <- focal(util_pred2
                   ,w=makeUtilBufWeight(5)
                   ,fun=max
                   ,na.rm=T
                   ,padValue=-2000
                   )

util_pred <- calc(util_pred3,fun=function(x){return(1000-x)})
util_pred <- calc(util_pred,fun=function(x){return(ifelse(x==Inf,NA,x))})


co_5_0_20_s <- raster(paste(wd_raster_out,'/co_5_0_20_s.tif',sep=''))
co_5_0_625_p <- raster(paste(wd_raster_out,'/co_5_0_625_p.tif',sep=''))
co_5_0_5_m <- raster(paste(wd_raster_out,'/co_5_0_5_m.tif',sep=''))

# stacking rasters that will have their values
# extracted at points
comb_extract <- stack(c(reservoir,thermal,seismic,utilization
                        ,re_pfa_var5,th_pfa_var5,se_pfa_var5,util_pfa_var5
                        ,re_pfa_var5_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls
                        ,res_pred_max2,res_pred_max_err
                        ,therm_pred,therm_err
                        ,seis_eq_pred,seis_eq_err,seis_stress_pred,seis_stress_err
                        ,util_pred,util_err
                        ,co_5_0_20_s,co_5_0_625_p,co_5_0_5_m))
                        #,co_5_0_15_s_geo,co_5_0_125_p_geo,co_5_0_5_m_geo))

# names for each layer in the raster stack
comb_names <- c('reservoir','thermal','seismic','utilization'
                ,'re_pfa_var5','th_pfa_var5','se_pfa_var5','util_pfa_var5'
                ,'re_pfa_var5_ls','th_pfa_var5_ls','se_pfa_var5_ls','util_pfa_var5_ls'
                ,'res_pred_max2','res_pred_max_err'
                ,'therm_pred','therm_err'
                ,'seis_eq_pred','seis_eq_err','seis_stress_pred','seis_stress_err'
                ,'util_pred','util_err'
                ,'co_5_0_20_s','co_5_0_625_p','co_5_0_5_m')
                #,'co_5_0_15_s_geo','co_5_0_125_p_geo','co_5_0_5_m_geo')

extracted <- rasterToPoints(comb_extract,spatial=TRUE)
extracted_df <- as.data.frame(extracted)


# renaming the columns
names(extracted_df) <- c(comb_names,'x','y')

# dropping points that are NA for geology or all four risk factors
extracted2_df <- extracted_df[setdiff(seq(1,nrow(extracted_df),1),intersect(which(extracted_df$co_5_0_125_p_geo%in% NA),which(extracted_df$co_5_0_625_p %in% NA))),]

# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
points[comb_names] <- NA

# calculating distance to each of the key points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  extracted2_df[nm] <- sqrt((extracted2_df$x-points$x[i])^2 + (extracted2_df$y-points$y[i])^2)
}


# extracting the values corresponding to the maximum
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1 ){
      ind_max <- which(extracted2_df$co_5_0_625_p[inds] %in% max(extracted2_df$co_5_0_625_p[inds]))
    }else{
      ind_max <- which(extracted2_df$co_5_0_625_p[inds] %in% max(extracted2_df$co_5_0_625_p[inds]))[1]
    }
    points[i,c(comb_names,'x','y')] <- extracted2_df[inds[ind_max],seq(1,length(c(comb_names,'x','y')),1)]
  }
}

points$util_err <- 5


## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))

# setting working directory and name
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
     ,ylab='Favorability'
     ,xlab=''
     ,xaxt='n')
# making vertical lines for the middle objectives
lines(c(1,1)
     ,c(0,5)
     ,col='black')
lines(c(2,2)
      ,c(0,5)
      ,col='black')

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

# initializing counter variable
j  <- 1

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # condition for no NAs
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
# ignoring 3 because it corresponds to Ithaca, which was not defined
ind_use <- c(1,2,4,5,6,7,8,9,10)

# number of replicates
rps <- 10000

# initialzing matrices
distsavg <- matrix(0,rps,9)
distsgeomean<- matrix(0,rps,9)
distsmin <- matrix(0,rps,9)
dist_vars <- matrix(0,4,9)

distsavg_g <- matrix(0,rps,9)
distsgeomean_g <- matrix(0,rps,9)
distsmin_g <- matrix(0,rps,9)
dist_vars_g <- matrix(0,4,9)

minRe <- rep(0,9)
minTh <- rep(0,9)
minSe <- rep(0,9)
minUt <- rep(0,9)

minRe_g <- rep(0,9)
minTh_g <- rep(0,9)
minSe_g <-rep(0,9)

reMean <- rep(0,9)
thMean <- rep(0,9)
seMean <- rep(0,9)
utMean <- rep(0,9)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5)
  
  # reservoir MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_max2[ind_use[i]])
  res_cv <- points$res_pred_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))

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
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
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
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
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
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$seis_stress_pred[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- abs(rnorm(rps,se_st_mean,se_st_se))
  
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
  
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100)

  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < 5] <- 5
  pfm5[rand > 25] <- 0
  pfm5[intersect(which(rand >= 5),which(rand < 12))] <- 5 - (rand[intersect(which(rand >= 5),which(rand < 12))] - 5)/(12-5)
  pfm5[intersect(which(rand >= 12),which(rand < 13.5))] <- 4- (rand[intersect(which(rand >= 12),which(rand < 13.5))] - 12)/(13.5-12)
  pfm5[intersect(which(rand >= 13.5),which(rand < 16))] <- 3- (rand[intersect(which(rand >= 13.5),which(rand < 16))] - 13.5)/(16-13.5)
  pfm5[intersect(which(rand >= 16),which(rand < 20))] <- 2- (rand[intersect(which(rand >= 16),which(rand < 20))] - 16)/(20-16)
  pfm5[intersect(which(rand >= 20),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 20),which(rand < 25))] - 20)/(25-20)
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  for(j in 1:rps){
    if(distsmin[j,i] == mat_mc[j,1]){
      minRe[i] <- minRe[i] + 1
    }
    if(distsmin[j,i] == mat_mc[j,2]){
      minTh[i] <- minTh[i] + 1
    }
    if(distsmin[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
      minSe[i] <- minSe[i] + 1
    }
    if(distsmin[j,i] == mat_mc[j,5]){
      minUt[i] <- minUt[i] + 1
    }
  }

  for(j in 1:rps){
    if(distsmin_g[j,i] == mat_mc[j,1]){
      minRe_g[i] <- minRe_g[i] + 1
    }
    if(distsmin_g[j,i] == mat_mc[j,2]){
      minTh_g[i] <- minTh_g[i] + 1
    }
    if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
      minSe_g[i] <- minSe_g[i] + 1
    }
  }
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
#   print(mean(mat_mc[,1]))
#   print(mean(mat_mc[,2]))
#   print(mean(mat_mc[,3]))
#   print(mean(mat_mc[,4]))
#   print(mean(mat_mc[,5]))
  
  setwd(wd_image)
  par(xpd=T)
  png(paste('scdist',i,'.png',sep='')
      ,height=6
      ,width=6
      ,units='in'
      ,res=300
  )
  
  par(mfrow=c(2,2)
      ,oma=c(0.5,0.5,0.5,0.5)+0.1
      ,mar=c(5,5,3,0)+0.1)

  hist(mat_mc[,1]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Reservoir"
       ,xlab="SFF"
       )
  lines(seq(0,5,0.01)
        ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
        ,col='royalblue'
        ,lwd=2
        )
  points(points$reservoir[ind_use[i]],0
         ,pch=19
         )
  
  hist(mat_mc[,2]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Thermal"
       ,xlab="SFF"
  )
  lines(seq(0,5,0.01)
        ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
        ,col='royalblue'
        ,lwd=2
        )
  points(points$thermal[ind_use[i]],0
         ,pch=19
  )
  
  hist(0.5*(mat_mc[,3]+mat_mc[,4])
       ,freq=F
       ,xlim=c(0,5)
       ,main="Seismic"
       ,xlab="SFF"
  )
  lines(seq(0,5,0.01)
        ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
        ,col='royalblue'
        ,lwd=2
        )
  points(points$seismic[ind_use[i]],0
         ,pch=19
  )
  
  hist(mat_mc[,5]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Utilization"
       ,xlab="SFF"
  )
  lines(seq(0,5,0.01)
        ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
        ,col='royalblue'
        ,lwd=2
  )
  points(points$utilization[ind_use[i]],0
         ,pch=19
  )
  
  par(xpd=T)
  dev.off()
}

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2)
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir[i],points$re_pfa_var5[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3)
roots2 <- matrix(NA,10,3)
converge1 <- matrix(NA,10,3)
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,]
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)]
  
  ps <- 0.05
  roots2[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge2[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots2[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge2[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots2[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge2[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
}

###################################
# making plot comparing analytical and monte carlo results
# all variables

setwd(wd_image)
#par(mar=c(1,1,1,1)*0.1)
png('three_panel2_all.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)

par(mfrow=c(1,5)
    ,oma=c(1,0,0,3)+0.1
    ,mar=c(5,0,0,0)+0.1)
dshift1 <- 0.4
dshift2 <- 0.2

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')
par(xpd=T)
for(i in 1:length(ind_use)){
  text(x=1,y=i-0.5
       ,labels=points$names[ind_use[i]]
       ,adj=1
       ,cex=1.2)
  
}
par(xpd=F)

plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Minimum'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  res_mean <- points$reservoir[ind_use[i]]
  #res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  #therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  #seis_var <- points$se_pfa_var5[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  
  # calculating empirical quantiles
  empQuant <- quantile(distsmin[,i],c(0.05,0.95))
  
  # calculated values
  if(max(is.na(dataParams[ind_use[i],c(7,8)])) == 0){
    lines(roots1[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
          ,lwd=2
          ,col='black'
    )
  }
  
  points(min(res_mean,therm_mean,seis_mean,util_mean),i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsmin[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}



plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Geometric Mean'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  # calculating empirical quantiles
  empQuant <- quantile(distsgeomean[,i],c(0.05,0.95))
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_ls[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var_ls <- points$th_pfa_var5_ls[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var_ls <- points$se_pfa_var5_ls[ind_use[i]]
  if(is.na(seis_var_ls)){
    seis_var_ls <- 0
  }
  
  util_mean <- points$utilization[ind_use[i]]
  util_var_ls <- points$util_pfa_var5_ls[ind_use[i]]
  
  var_geomean_ls <- (1/16)*(res_var_ls + therm_var_ls + seis_var_ls + util_var_ls)
  mean_geomean_ls <- 0.25*(log(res_mean) + log(therm_mean) + log(seis_mean) + log(util_mean))
  
  theoQuants <- exp(qnorm(c(0.05,0.95),mean_geomean_ls,sqrt(var_geomean_ls)))
  
  # calculated values
  lines(theoQuants,c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='black'
  )
  
  points(exp(mean_geomean_ls),i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsgeomean[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}


plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Average'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  # calculating empirical quantiles
  empQuant <- quantile(distsavg[,i],c(0.05,0.95))
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  util_var <- points$util_pfa_var5[ind_use[i]]
  
  var_avg <- (1/16)*(seis_var+res_var+therm_var+util_var)
  mean_avg <- (res_mean+therm_mean +seis_mean+util_mean)/4
  # calculated valuesseis_mean
  lines(qnorm(c(0.05,0.95),mean_avg,sqrt(var_avg)),c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='black'
  )
  
  points(mean_avg,i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsavg[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')

legend('center'
       ,legend=c('Favorability \n Metric Value','Monte Carlo \n Mean','Approximate \n 90% PI','Monte Carlo \n 90% PI')
       ,pch=c(4,4,NA,NA)
       ,col=c('black','gray48','black','gray48')
       ,lwd=c(NA,NA,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()


setwd(wd_image)
#par(mar=c(1,1,1,1)*0.1)
png('three_panel2_geo.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)

par(mfrow=c(1,5)
    ,oma=c(1,0,0,3)+0.1
    ,mar=c(5,0,0,0)+0.1)
dshift1 <- 0.4
dshift2 <- 0.2

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')
par(xpd=T)
for(i in 1:length(ind_use)){
  text(x=1,y=i-0.5
       ,labels=points$names[ind_use[i]]
       ,adj=1
       ,cex=1.2)
  
}
par(xpd=F)

plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Minimum'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  res_mean <- points$reservoir[ind_use[i]]
  #res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  #therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  #seis_var <- points$se_pfa_var5[ind_use[i]]
  
  #util_mean <- points$utilization[ind_use[i]]
  
  # calculating empirical quantiles
  empQuant <- quantile(distsmin_g[,i],c(0.05,0.95))
  
  # calculated values
  lines(roots2[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='black'
  )
  
  points(min(res_mean,therm_mean,seis_mean),i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsmin_g[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}



plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Geometric Mean'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  # calculating empirical quantiles
  empQuant <- quantile(distsgeomean_g[,i],c(0.05,0.95))
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_ls[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var_ls <- points$th_pfa_var5_ls[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var_ls <- points$se_pfa_var5_ls[ind_use[i]]
  if(is.na(seis_var_ls)){
    seis_var_ls <- 0
  }
  
  var_geomean_ls <- (1/9)*(res_var_ls + therm_var_ls + seis_var_ls)
  mean_geomean_ls <- (log(res_mean) + log(therm_mean) + log(seis_mean))/3
  
  theoQuants <- exp(qnorm(c(0.05,0.95),mean_geomean_ls,sqrt(var_geomean_ls)))
  
  # calculated values
  lines(theoQuants,c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='black'
  )
  
  points(exp(mean_geomean_ls),i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsgeomean_g[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}


plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Average'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  
  # calculating empirical quantiles
  empQuant <- quantile(distsavg_g[,i],c(0.05,0.95))
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]
  
  var_avg <- (1/9)*(seis_var+res_var+therm_var)
  mean_avg <- (res_mean+therm_mean +seis_mean)/3
  # calculated valuesseis_mean
  lines(qnorm(c(0.05,0.95),mean_avg,sqrt(var_avg)),c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='black'
  )
  
  points(mean_avg,i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # empirical values
  lines(c(as.numeric(empQuant)),c(i-1+dshift1,i-1+dshift1)
        ,lwd=2
        ,col='gray48'
  )
  
  points(mean(distsavg_g[,i]),i-1+dshift1
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
}

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')

legend('center'
       ,legend=c('Favorability \n Metric Value','Monte Carlo \n Mean','Approximate \n 90% PI','Monte Carlo \n 90% PI')
       ,pch=c(4,4,NA,NA)
       ,col=c('black','gray48','black','gray48')
       ,lwd=c(NA,NA,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()

###################################
# making plot of different metrics for all favorability factors

setwd(wd_image)
#par(mar=c(1,1,1,1)*0.1)
png('single_panel2_all.png'
    ,height=6
    ,width=7
    ,units='in'
    ,res=300
)
layout(mat=matrix(c(1,2,3),1,3)
       ,widths=c(1,3,1.5))
par(oma=c(1,0,0,3)+0.1
    ,mar=c(5,0,0,0)+0.1)

dshift1 <- 0.3
dshift2 <- 0.2

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')
par(xpd=T)
for(i in 1:length(ind_use)){
  text(x=1,y=i-0.5
       ,labels=points$names[ind_use[i]]
       ,adj=1
       ,cex=1.2)
  
}
par(xpd=F)

plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Favorability Metric'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  res_mean <- points$reservoir[ind_use[i]]

  therm_mean <- points$thermal[ind_use[i]]

  seis_mean <- points$seismic[ind_use[i]]

  util_mean <- points$utilization[ind_use[i]]
  
  # calculated values
  if(max(is.na(dataParams[ind_use[i],c(7,8)])) == 0){
    lines(roots1[ind_use[i],c(1,3)],c(i-1+dshift1,i-1+dshift1)
          ,lwd=2
          ,col='black'
    )
  }
  
  points(min(res_mean,therm_mean,seis_mean,util_mean),i-1+dshift1
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  # calculating empirical quantiles
  empQuant <- quantile(distsgeomean[,i],c(0.05,0.95))
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_ls[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var_ls <- points$th_pfa_var5_ls[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var_ls <- points$se_pfa_var5_ls[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  util_var_ls <- points$util_pfa_var5_ls[ind_use[i]]
  
  var_geomean_ls <- (1/16)*(res_var_ls + therm_var_ls + seis_var_ls + util_var_ls)
  mean_geomean_ls <- 0.25*(log(res_mean) + log(therm_mean) + log(seis_mean) + log(util_mean))
  
  theoQuants <- exp(qnorm(c(0.05,0.95),mean_geomean_ls,sqrt(var_geomean_ls)))
  
  # calculated values
  lines(theoQuants,c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='gray48'
  )
  
  points(exp(mean_geomean_ls),i-1+dshift1+dshift2
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  util_var <- points$util_pfa_var5[ind_use[i]]
  
  var_avg <- (1/16)*(seis_var+res_var+therm_var+util_var)
  mean_avg <- (res_mean+therm_mean +seis_mean+util_mean)/4
  # calculated valuesseis_mean
  lines(qnorm(c(0.05,0.95),mean_avg,sqrt(var_avg)),c(i-1+dshift1+2*dshift2,i-1+dshift1+2*dshift2)
        ,lwd=2
        ,col='royalblue'
  )
  
  points(mean_avg,i-1+dshift1+2*dshift2
         ,pch=4
         ,col='royalblue'
         ,cex=1.5
  )
}

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')

legend('center'
       ,legend=c(rev(c('FM Minimum','FM Geometric \n  Mean','FM Average')),rev(c('Approximate 90% \n PI Minimum','Approximate 90% \n PI Geometric Mean','Approximate 90% \n PI Average')))
       ,pch=c(4,4,4,NA,NA,NA)
       ,col=rev(c('black','gray48','royalblue'))
       ,lwd=c(NA,NA,NA,2,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()



###################################
# making plot of different metrics for all favorability factors

setwd(wd_image)
#par(mar=c(1,1,1,1)*0.1)
png('single_panel2_geo.png'
    ,height=6
    ,width=7
    ,units='in'
    ,res=300
)
layout(mat=matrix(c(1,2,3),1,3)
       ,widths=c(1,3,1.5))
par(oma=c(1,0,0,3)+0.1
    ,mar=c(5,0,0,0)+0.1)

dshift1 <- 0.3
dshift2 <- 0.2

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')
par(xpd=T)
for(i in 1:length(ind_use)){
  text(x=1,y=i-0.5
       ,labels=points$names[ind_use[i]]
       ,adj=1
       ,cex=1.2)
  
}
par(xpd=F)

plot(NA,NA
     ,xlim=c(0,5)
     ,ylim=c(0,length(ind_use))
     ,main=''
     ,xlab='Favorability Metric'
     ,yaxt='n'
     ,ylab=''
     ,cex.lab=1.2)
grid(ny=NA)
for(i in 1:length(ind_use)){
  res_mean <- points$reservoir[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  
  # calculated values
    lines(roots2[ind_use[i],c(1,3)],c(i-1+dshift1,i-1+dshift1)
          ,lwd=2
          ,col='black'
    )
  
  points(min(res_mean,therm_mean,seis_mean),i-1+dshift1
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_ls[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var_ls <- points$th_pfa_var5_ls[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var_ls <- points$se_pfa_var5_ls[ind_use[i]]
  
  var_geomean_ls <- (1/9)*(res_var_ls + therm_var_ls + seis_var_ls)
  mean_geomean_ls <- (1/3)*(log(res_mean) + log(therm_mean) + log(seis_mean))
  
  theoQuants <- exp(qnorm(c(0.05,0.95),mean_geomean_ls,sqrt(var_geomean_ls)))
  
  # calculated values
  lines(theoQuants,c(i-1+dshift1+dshift2,i-1+dshift1+dshift2)
        ,lwd=2
        ,col='gray48'
  )
  
  points(exp(mean_geomean_ls),i-1+dshift1+dshift2
         ,pch=4
         ,col='gray48'
         ,cex=1.5
  )
  
  #extracting values value
  res_mean <- points$reservoir[ind_use[i]]
  res_var <- points$re_pfa_var5[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  util_var <- points$util_pfa_var5[ind_use[i]]
  
  var_avg <- (1/9)*(seis_var+res_var+therm_var)
  mean_avg <- (res_mean+therm_mean +seis_mean)/3
  # calculated valuesseis_mean
  lines(qnorm(c(0.05,0.95),mean_avg,sqrt(var_avg)),c(i-1+dshift1+2*dshift2,i-1+dshift1+2*dshift2)
        ,lwd=2
        ,col='royalblue'
  )
  
  points(mean_avg,i-1+dshift1+2*dshift2
         ,pch=4
         ,col='royalblue'
         ,cex=1.5
  )
}

plot(NA,NA
     ,xlim=c(0,1)
     ,ylim=c(0,length(ind_use))
     ,xaxt='n'
     ,yaxt='n'
     ,ylab=''
     ,xlab=''
     ,bty='n')

legend('center'
       ,legend=c(rev(c('FM Minimum','FM Geometric \n  Mean','FM Average')),rev(c('Approximate 90% \n PI Minimum','Approximate 90% \n PI Geometric Mean','Approximate 90% \n PI Average')))
       ,pch=c(4,4,4,NA,NA,NA)
       ,col=rev(c('black','gray48','royalblue'))
       ,lwd=c(NA,NA,NA,2,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()




##############################
# extracting values of layers for cities
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA,se_5_0_5_a,ut5_5_0_5_NA))
                                #,co_5_0_20_s,co_5_0_625_p,co_5_0_5_m))

co_min <- stack(c(th_5_0_5_NA,re_5_0_5_NA,se_5_0_5_a,ut5_5_0_5_NA))
combined_all_rasters[which(values(combined_all_rasters) %in% -9999)] <- NA

extract_pts <- extract(x=combined_all_rasters
                       ,y=cities2
                       ,sp=TRUE
                       ,nl=7
                       ,df=TRUE
                       ,na.rm=TRUE
                       ,method='simple'
                       ,buffer=2000
                       ,fun=max
)

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[13] <- "th_5"
names(extract_pts2)[14] <- "re_5"
names(extract_pts2)[15] <- "se_5"
names(extract_pts2)[16] <- "ut_5"
#names(extract_pts2)[17] <- "co_s_5"
#names(extract_pts2)[18] <- "co_p_5"
#names(extract_pts2)[19] <- "co_m_5"


extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]
extract_pts3$co_p_5 <- extract_pts3$th_5*extract_pts3$re_5*extract_pts3$se_5*extract_pts3$ut_5
extract_pts3$co_s_5 <- extract_pts3$th_5+extract_pts3$re_5+extract_pts3$se_5+extract_pts3$ut_5
extract_pts3$co_m_5 <- apply(cbind(extract_pts3$th_5,extract_pts3$re_5,extract_pts3$se_5,extract_pts3$ut_5),FUN=min,MARGIN=1)



ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5  > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5  > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
png('scatter_p_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]/4
     ,extract_pts3$co_p_5[inds_gtr0]^0.25
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Average'
     ,ylab='Geometric Mean'
     ,xaxt='n'
     ,yaxt='n'
     ,cex=0.8
)
grid()
points(extract_pts3$co_s_5[inds_gtr0]/4
     ,extract_pts3$co_p_5[inds_gtr0]^0.25
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Average'
     ,ylab='Geometric Mean'
     ,xaxt='n'
     ,yaxt='n'
     ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_s.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_s_5[inds_gtr0]/4
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Average'
     ,ylab='Minimum'
     ,xaxt='n'
     ,yaxt='n'
     ,cex=0.8
)
grid()
points(extract_pts3$co_s_5[inds_gtr0]/4
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Average'
     ,ylab='Minimum'
     ,xaxt='n'
     ,yaxt='n'
     ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
axis(2,
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_m_p.png'
    ,height=6
    ,width=6
    ,units='in'
    ,res=300
)
plot(extract_pts3$co_p_5[inds_gtr0]^0.25
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Geometric Mean'
     ,ylab='Minimum'
     ,xaxt='n'
     ,yaxt='n'
     ,cex=0.8
)
grid()
points(extract_pts3$co_p_5[inds_gtr0]^0.25
     ,extract_pts3$co_m_5[inds_gtr0]
     ,pch=19
     ,xlim=c(0,5)
     ,ylim=c(0,5)
     ,xlab='Geometric Mean'
     ,ylab='Minimum'
     ,xaxt='n'
     ,yaxt='n'
    ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# saving workspace
setwd(wd_workspace)
filename <- paste('geotherma_pfa_analysis_',Sys.Date(),'.RData',sep='')
save.image(file=filename)