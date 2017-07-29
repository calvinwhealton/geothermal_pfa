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
library(xlsx)         # for reading in xlsx tables. Need 64 bit Java on computer to use this package.
library(rgeos)        # for buffering places of interest
library(RColorBrewer) # R color brewer palettes for parallel axis plot
library(pracma)       # for interpolation in tables
library(vioplot)      # for violin plots
#Edit the vioplot function
#Change the col parameter to col[i], and add a cex=0.75 to the points for median
vioplot = edit(vioplot)
library(prettymapr)
library(maps)
library(rootSolve)    #For root solving of Weibull and other distributions
library(Hmisc)        #For minor tick marks on plots
library(circular)     #For the von Mises distribution
library(ExtDist)      #For Generalized Beta distribution
library(msm)          #For doubly truncated normal distribution
library(VGAM)         #For folded normal distribution

##### defining working directories #####
# need to be changed based on machine

setwd('C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa')

# location to SAVE rasters
wd_raster_out <- paste(getwd(), '/Rasters/Rasters_out/CorrMeans', sep='')

# location to READ rasters
wd_raster_in <- paste(getwd(), '/Rasters/Rasters_in', sep='')

# location to save images
wd_image <- paste(getwd(), '/Results/CorrMeans/PaperGraphics', sep='')

# location of code/scripts are stored
wd_code <- paste(getwd(), '/combining_metrics', sep='')

# location of input shapefiles
wd_shapefiles <- paste(getwd(), '/Shapefiles', sep='')

# location of error interoplation tables
wd_error_interp <- paste(getwd(), '/Rasters/Rasters_in/error_interp_tabs', sep='') #Output from make_interp_tables.R

# location to save workspace
wd_workspace <- paste(getwd(), '/Results', sep='')

##### loading-in state/county shapefiles #####
# states
States = readOGR(dsn="C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_state_500k", 
                 layer="us_state_WGS", 
                 stringsAsFactors=FALSE)
NY = States[which(States$STATEFP == "36"),]
PA = States[which(States$STATEFP == "42"),]
WV = States[which(States$STATEFP == "54"),]

# counties
Counties = readOGR(dsn="C:/Users/Jared/Documents/Cornell/Research/Masters - Spatial Assessment/GIS/Population Density/2013_us_countys_500k", 
                   layer="us_county_WGS", 
                   stringsAsFactors=FALSE) 
NY_co = Counties[which(Counties$STATEFP == "36"),]
PA_co = Counties[which(Counties$STATEFP == "42"),]
WV_co = Counties[which(Counties$STATEFP == "54"),]

# converting to UTM 17N system
NY2 <- spTransform(NY,CRS("+init=epsg:26917"))
PA2 <- spTransform(PA,CRS("+init=epsg:26917"))
WV2 <- spTransform(WV,CRS("+init=epsg:26917"))

NY_co2 <- spTransform(NY_co,CRS("+init=epsg:26917"))
PA_co2 <- spTransform(PA_co,CRS("+init=epsg:26917"))
WV_co2 <- spTransform(WV_co,CRS("+init=epsg:26917"))

# importing city locations, the US Census Places
cities <- readOGR(dsn=wd_shapefiles
                  , layer="usCensusPlaces"
                  , stringsAsFactors=FALSE)
cities2 <- spTransform(cities,CRS("+init=epsg:26917"))

# removing the untransformed shapefiles
rm(NY,PA,WV,States,Counties,NY_co,PA_co,WV_co,cities)

## importing places of interest
poi <- readOGR(dsn=wd_shapefiles
               , layer="placesofinterest2"
               , stringsAsFactors=FALSE)

# converting to UTM 17 N system
poi2 <- spTransform(poi,CRS("+init=epsg:26917"))

# buffering places of interest by 5000 and 10000 m
poi2Buf5 <- gBuffer(poi2,width=5000, quadsegs=1000) #Circles created with this function are not the best. I increased quadsegs to make them rounder.

# shapefile of the US Census places merged
pla <- readOGR(dsn=wd_shapefiles,layer='mergerIhopeThisWorksZachFile', stringsAsFactors = FALSE)

# converting to UTM 17 N system
pla <- spTransform(pla,CRS("+init=epsg:26917"))

# buffering places of interest by 5000 and 10000 m
plaBuf5 <- gBuffer(pla,width=5000, quadsegs=1000)

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
source('rw_functions_lb.R') #Weibull and beta dists. Note that the new functions have not been added to all variables yet.

#### Set Thresholds ####
# Reservoirs
res_rfc_min <- 0.001 # minimum for reservoir (both 3-color and 5-color)
res_rfc_max <- 10000 # maximum for reservoir (both 3-color and 5-color)
res_rfc_thresh5 <- c(res_rfc_min,c(1,10,100,1000),res_rfc_max)

res_RPIw_min <- 0.0001 # minimum for reservoir (both 3-color and 5-color)
res_RPIw_max <- 10    # maximum for reservoir (both 3-color and 5-color)
res_RPIw_thresh5 <- c(res_RPIw_min,c(0.001,0.01,0.1,1.0),res_RPIw_max)

res_RPIg_min <- 0.0001 # minimum for reservoir (both 3-color and 5-color)
res_RPIg_max <- 10     # maximum for reservoir (both 3-color and 5-color)
res_RPIg_thresh5 <- c(res_RPIg_min,c(0.001,0.01,0.1,1.0),res_RPIg_max)

#Seismic Stress
seis_stress_min <- 0.001 # to avoid problems with numerically zero values
seis_stress_max <- 25
seis_stress_thresh5<- c(seis_stress_min,c(5,10,15,20),seis_stress_max)

# set critical angles:
critical_ang1 = 65.2
critical_ang2 = 114.8

#Seismic angle
seis_eq_min <- 0.001 # to avoid problems with numerically zero values
seis_eq_max <- 25
seis_eq_thresh5<- 10^3*c(seis_eq_min,c(5,10,15,20),seis_eq_max)

#Thermal
therm_thresh5 <- rev(c(5000,4000,3000,2500,2000,1000)) #3000 is about $6.4 mil/well and 23 C/km. 2500 is about $4.8 mil/well and 28C/km.

#Utilization
util_thresh5 <- c(5,12,13.5,16,20,25)


##### THERMAL ######
# importing rasters, depth to 80 DegC and standard error of prediction
therm_pred <- raster(paste(wd_raster_in,'/Thermal/d80cp-',sep=''))
therm_err  <- raster(paste(wd_raster_in,'/Thermal/d80ce-',sep=''))

# values of -9999 are no data, so replaced with NA
therm_pred[therm_pred < 0] <- NA
therm_err[therm_err < 0] <- NA

#Load interpolation table for the mean
th_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/th_pfmean5.xlsx',sep=''),1,header=FALSE))

#Set the range of means and variances of the interpolation tables
# must match values from make_interp_table.R
mean_thd80 <- seq(1690,6490,by=200) # range of means
std_thd80 <- seq(40,1840,by=50) # range of standard deviations

# values to interpolate the scaled mean
th_means <- values(therm_pred)
th_ses <- values(therm_err)

# substituting in values for NAs so interpolation algorithm will not crash
# These are set back to NA after the interpolation
th_means[which(th_means %in% NA)] <- min(mean_thd80)
th_ses[which(th_ses %in% NA)] <- min(std_thd80)

# five color
thvecPFmean5 <- interp2(x=std_thd80
                       ,y=mean_thd80
                       ,Z=th_interp_tab5_mean
                       ,xp=th_ses
                       ,yp=th_means
                       ,method='linear')

# setting values back to NAs
thvecPFmean5[which(values(therm_pred) %in% NA)] <- NA

# initializing raster for the stored values
th_5_0_5_NA <- therm_pred

# updating values of the raster
# note: setValues() did not work
values(th_5_0_5_NA) <- thvecPFmean5

saveRast(rast=th_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='th_5_0_5_NA.tif')
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1, leg2 = T, grey = F, County = F, RawThreshVals = therm_thresh5/1000, Unit = 'km', dpi = 300, FigFun = 'tiff')
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA_Grey.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1, leg2 = T, grey = T, County = F, RawThreshVals = therm_thresh5/1000, Unit = 'km', dpi = 600, FigFun = 'tiff')

# making thermal uncertainty maps
th_interp_tab5 <- as.matrix(read.xlsx(paste(wd_error_interp,'/th_pfvar5.xlsx',sep=''),1,header=FALSE))
th_interp_tab5_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/th_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

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
thvecPFvar5[which(values(therm_pred) %in% NA)] <- NA
thvecPFvar5_ls[which(values(therm_pred) %in% NA)] <- NA

# initializing raster for the stored values
th_pfa_var5 <- therm_err
th_pfa_var5_ls <- therm_err

# updating values of the raster
values(th_pfa_var5) <- thvecPFvar5
values(th_pfa_var5_ls) <- thvecPFvar5_ls

# saving rasters and making maps of variance
# NAs are set to -9999 in the raster.
saveRast(rast=th_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var5.tif')
saveRast(rast=th_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='th_pfa_var5_ls.tif')

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
makeMap (rast=calc(th_pfa_var5,fun=sqrt)
         ,plotnm='th_pfa_sd5.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap (rast=calc(th_pfa_var5,fun=sqrt)
         ,plotnm='th_pfa_sd5_Grey.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE, grey = T, County = F, dpi = 600, FigFun = 'tiff')

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
makeMap (rast=calc(th_pfa_var5_ls,fun=sqrt)
         ,plotnm='th_pfa_sd5_ls.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap (rast=calc(th_pfa_var5_ls,fun=sqrt)
         ,plotnm='th_pfa_sd5_ls_Grey.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE, grey = T, County = F, dpi = 600, FigFun = 'tiff')

# deleting unneeded variables
rm(th_ses,th_means,thvecPFvar5,std_thd80,mean_thd80)
rm(th_interp_tab5, th_interp_tab5_ls, th_interp_tab5_mean)
rm(thvecPFmean5, thvecPFvar5_ls)

##### RESERVOIR ######
# l1000 = less than 1000 m
# 1015 = 1000 to 1500 m
# 1520 = 1500 to 2000 m, etc.

# reservoir RFC pred
res_pred_l1000_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_l1000p-',sep=''))
res_pred_1015_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_1015p-',sep=''))
res_pred_1520_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_1520p-',sep=''))
res_pred_2025_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_2025p-',sep=''))
res_pred_2530_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_2530p-',sep=''))
res_pred_3035_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_3035p-',sep=''))
res_pred_3540_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_3540p-',sep=''))

# reservoir RFC error - Corr Overlay
res_err_l1000_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_l1000e-',sep=''))
res_err_1015_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_1015e-',sep=''))
res_err_1520_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_1520e-',sep=''))
res_err_2025_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_2025e-',sep=''))
res_err_2530_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_2530e-',sep=''))
res_err_3035_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_3035e-',sep=''))
res_err_3540_rfc <- raster(paste(wd_raster_in,'/Reservoirs/RFC/CorrOverlay/fff_3540e-',sep=''))

# reservoir RPIw pred - Corr Overlay
res_pred_l1000_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_l1000p-',sep=''))
res_pred_1015_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_1015p-',sep=''))
res_pred_1520_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_1520p-',sep=''))
res_pred_2025_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_2025p-',sep=''))
res_pred_2530_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_2530p-',sep=''))
res_pred_3035_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_3035p-',sep=''))
res_pred_3540_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_3540p-',sep=''))

# reservoir RPIw error - Corr Overlay
res_err_l1000_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_l1000e-',sep=''))
res_err_1015_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_1015e-',sep=''))
res_err_1520_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_1520e-',sep=''))
res_err_2025_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_2025e-',sep=''))
res_err_2530_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_2530e-',sep=''))
res_err_3035_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_3035e-',sep=''))
res_err_3540_RPIw <- raster(paste(wd_raster_in,'/Reservoirs/RPIw/CorrOverlay/fff_3540e-',sep=''))

# reservoir RPIg pred
res_pred_l1000_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_l1000p-',sep=''))
res_pred_1015_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_1015p-',sep=''))
res_pred_1520_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_1520p-',sep=''))
res_pred_2025_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_2025p-',sep=''))
res_pred_2530_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_2530p-',sep=''))
res_pred_3035_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_3035p-',sep=''))
res_pred_3540_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_3540p-',sep=''))

# reservoir RPIg error
res_err_l1000_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_l1000e-',sep=''))
res_err_1015_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_1015e-',sep=''))
res_err_1520_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_1520e-',sep=''))
res_err_2025_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_2025e-',sep=''))
res_err_2530_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_2530e-',sep=''))
res_err_3035_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_3035e-',sep=''))
res_err_3540_RPIg <- raster(paste(wd_raster_in,'/Reservoirs/RPIg/CorrOverlay/fff_3540e-',sep=''))

# stacking rasters and taking maximum, assuming going for highest quality reservoir
# ignoring reservoirs shallower than 1000 m (the res_pred_l1000 raster)
res_pred_rfc <- stack(c(res_pred_1015_rfc,res_pred_1520_rfc,res_pred_2025_rfc,res_pred_2530_rfc,res_pred_3035_rfc,res_pred_3540_rfc))
res_pred_rfc_max <- calc(res_pred_rfc,fun=max,na.rm=TRUE)

res_pred_RPIw <- stack(c(res_pred_1015_RPIw,res_pred_1520_RPIw,res_pred_2025_RPIw,res_pred_2530_RPIw,res_pred_3035_RPIw,res_pred_3540_RPIw))
res_pred_RPIw_max <- calc(res_pred_RPIw,fun=max,na.rm=TRUE)

res_pred_RPIg <- stack(c(res_pred_1015_RPIg,res_pred_1520_RPIg,res_pred_2025_RPIg,res_pred_2530_RPIg,res_pred_3035_RPIg,res_pred_3540_RPIg))
res_pred_RPIg_max <- calc(res_pred_RPIg,fun=max,na.rm=TRUE)

# setting values less than zero (-9999) to NA because they have no data
res_pred_rfc_max2 <- res_pred_rfc_max
res_pred_rfc_max2[res_pred_rfc_max < 0] <- NA

res_pred_RPIw_max2 <- res_pred_RPIw_max
res_pred_RPIw_max2[res_pred_RPIw_max < 0] <- NA

res_pred_RPIg_max2 <- res_pred_RPIg_max
res_pred_RPIg_max2[res_pred_RPIg_max < 0] <- NA

# stacking reservoir errors
res_err_rfc <- stack(c(res_err_1015_rfc,res_err_1520_rfc,res_err_2025_rfc,res_err_2530_rfc,res_err_3035_rfc,res_err_3540_rfc))

res_err_RPIw <- stack(c(res_err_1015_RPIw,res_err_1520_RPIw,res_err_2025_RPIw,res_err_2530_RPIw,res_err_3035_RPIw,res_err_3540_RPIw))

res_err_RPIg <- stack(c(res_err_1015_RPIg,res_err_1520_RPIg,res_err_2025_RPIg,res_err_2530_RPIg,res_err_3035_RPIg,res_err_3540_RPIg))

## making reservoir error term
# step 1: raster of whether the max value is equal to the value in the layer
tf_rast1_rfc <- (res_pred_rfc_max == res_pred_rfc[[1]])
tf_rast2_rfc <- (res_pred_rfc_max == res_pred_rfc[[2]])
tf_rast3_rfc <- (res_pred_rfc_max == res_pred_rfc[[3]])
tf_rast4_rfc <- (res_pred_rfc_max == res_pred_rfc[[4]])
tf_rast5_rfc <- (res_pred_rfc_max == res_pred_rfc[[5]])
tf_rast6_rfc <- (res_pred_rfc_max == res_pred_rfc[[6]])

tf_rast1_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[1]])
tf_rast2_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[2]])
tf_rast3_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[3]])
tf_rast4_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[4]])
tf_rast5_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[5]])
tf_rast6_RPIw <- (res_pred_RPIw_max == res_pred_RPIw[[6]])

tf_rast1_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[1]])
tf_rast2_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[2]])
tf_rast3_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[3]])
tf_rast4_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[4]])
tf_rast5_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[5]])
tf_rast6_RPIg <- (res_pred_RPIg_max == res_pred_RPIg[[6]])

# step 2: stacking raster with error raster
# product means only values where the max was matched will be assigned a positive value
# otherwise, it will be 0
err_tf1_rfc <- calc(stack(c(tf_rast1_rfc,res_err_rfc[[1]])),fun=prod)
err_tf2_rfc <- calc(stack(c(tf_rast2_rfc,res_err_rfc[[2]])),fun=prod)
err_tf3_rfc <- calc(stack(c(tf_rast3_rfc,res_err_rfc[[3]])),fun=prod)
err_tf4_rfc <- calc(stack(c(tf_rast4_rfc,res_err_rfc[[4]])),fun=prod)
err_tf5_rfc <- calc(stack(c(tf_rast5_rfc,res_err_rfc[[5]])),fun=prod)
err_tf6_rfc <- calc(stack(c(tf_rast6_rfc,res_err_rfc[[6]])),fun=prod)

err_tf1_RPIw <- calc(stack(c(tf_rast1_RPIw,res_err_RPIw[[1]])),fun=prod)
err_tf2_RPIw <- calc(stack(c(tf_rast2_RPIw,res_err_RPIw[[2]])),fun=prod)
err_tf3_RPIw <- calc(stack(c(tf_rast3_RPIw,res_err_RPIw[[3]])),fun=prod)
err_tf4_RPIw <- calc(stack(c(tf_rast4_RPIw,res_err_RPIw[[4]])),fun=prod)
err_tf5_RPIw <- calc(stack(c(tf_rast5_RPIw,res_err_RPIw[[5]])),fun=prod)
err_tf6_RPIw <- calc(stack(c(tf_rast6_RPIw,res_err_RPIw[[6]])),fun=prod)

err_tf1_RPIg <- calc(stack(c(tf_rast1_RPIg,res_err_RPIg[[1]])),fun=prod)
err_tf2_RPIg <- calc(stack(c(tf_rast2_RPIg,res_err_RPIg[[2]])),fun=prod)
err_tf3_RPIg <- calc(stack(c(tf_rast3_RPIg,res_err_RPIg[[3]])),fun=prod)
err_tf4_RPIg <- calc(stack(c(tf_rast4_RPIg,res_err_RPIg[[4]])),fun=prod)
err_tf5_RPIg <- calc(stack(c(tf_rast5_RPIg,res_err_RPIg[[5]])),fun=prod)
err_tf6_RPIg <- calc(stack(c(tf_rast6_RPIg,res_err_RPIg[[6]])),fun=prod)

# step 3: making stacked raster and taking maximum
# only values where the raster matched the maximum of the rasters will be used
res_pred_rfc_max_err <- calc(stack(c(err_tf1_rfc,err_tf2_rfc,err_tf3_rfc,err_tf4_rfc,err_tf5_rfc,err_tf6_rfc)),fun=max)
res_pred_rfc_max_err[(res_pred_rfc_max_err < 0)] <- NA

res_pred_RPIw_max_err <- calc(stack(c(err_tf1_RPIw,err_tf2_RPIw,err_tf3_RPIw,err_tf4_RPIw,err_tf5_RPIw,err_tf6_RPIw)),fun=max)
res_pred_RPIw_max_err[(res_pred_RPIw_max_err < 0)] <- NA

res_pred_RPIg_max_err <- calc(stack(c(err_tf1_RPIg,err_tf2_RPIg,err_tf3_RPIg,err_tf4_RPIg,err_tf5_RPIg,err_tf6_RPIg)),fun=max)
res_pred_RPIg_max_err[(res_pred_RPIg_max_err < 0)] <- NA

# deleting unneeded rasters
rm(tf_rast1_rfc,tf_rast2_rfc,tf_rast3_rfc,tf_rast4_rfc,tf_rast5_rfc,tf_rast6_rfc)
rm(tf_rast1_RPIw,tf_rast2_RPIw,tf_rast3_RPIw,tf_rast4_RPIw,tf_rast5_RPIw,tf_rast6_RPIw)
rm(tf_rast1_RPIg,tf_rast2_RPIg,tf_rast3_RPIg,tf_rast4_RPIg,tf_rast5_RPIg,tf_rast6_RPIg)
rm(err_tf1_rfc,err_tf2_rfc,err_tf3_rfc,err_tf4_rfc,err_tf5_rfc,err_tf6_rfc)
rm(err_tf1_RPIw,err_tf2_RPIw,err_tf3_RPIw,err_tf4_RPIw,err_tf5_RPIw,err_tf6_RPIw)
rm(err_tf1_RPIg,err_tf2_RPIg,err_tf3_RPIg,err_tf4_RPIg,err_tf5_RPIg,err_tf6_RPIg)

#Load interpolation table for the mean
re_rfc_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_rfc_pfmean5.xlsx',sep=''),1,header=FALSE))
re_RPIw_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIw_pfmean5_NewThresh.xlsx',sep=''),1,header=FALSE))
re_RPIg_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIg_pfmean5_NewThresh.xlsx',sep=''),1,header=FALSE))

#Set the range of means and CVs of the interpolation tables
mean_re_rfc <- seq(-7,9.75,0.25) # range of means in log base e space
mean_re_RPIw <- seq(-14,3,0.25) # range of means in log base e space
mean_re_RPIg <- seq(-14,6.5,0.25) # range of means in log base e space
uncer_re <- seq(0,0.4,by=0.05) # range of real space coefficient of variations

# values are in base-e (ln) to match calculations in make_interp_tables
re_means_rfc <- log(values(res_pred_rfc_max2))
re_uncer_rfc <- values(res_pred_rfc_max_err)
re_means_RPIw <- log(values(res_pred_RPIw_max2))
re_uncer_RPIw <- values(res_pred_RPIw_max_err)
re_means_RPIg <- log(values(res_pred_RPIg_max2))
re_uncer_RPIg <- values(res_pred_RPIg_max_err)

# substituting in values for NAs so interpolation algorithm will not crash
re_means_rfc[which(re_means_rfc %in% NA)] <- min(mean_re_rfc)
re_uncer_rfc[which(re_uncer_rfc %in% NA)] <- min(uncer_re)

re_means_RPIw[which(re_means_RPIw %in% NA)] <- min(mean_re_RPIw)
re_uncer_RPIw[which(re_uncer_RPIw %in% NA)] <- min(uncer_re)

re_means_RPIg[which(re_means_RPIg %in% NA)] <- min(mean_re_RPIg)
re_uncer_RPIg[which(re_uncer_RPIg %in% NA)] <- min(uncer_re)

#Only needed when values are 0 because 0s are not valid in log space. Not the case in this most recent dataset.
#re_means[which(re_means %in% 0)] <- min(mean_re)
#re_uncer[which(re_uncer %in% 0)] <- min(uncer_re)

# five color
revecPFmean5_rfc <- interp2(x=uncer_re
                        ,y=mean_re_rfc
                        ,Z=re_rfc_interp_tab5_mean
                        ,xp=re_uncer_rfc
                        ,yp=re_means_rfc
                        ,method='linear')
revecPFmean5_RPIw <- interp2(x=uncer_re
                            ,y=mean_re_RPIw
                            ,Z=re_RPIw_interp_tab5_mean
                            ,xp=re_uncer_RPIw
                            ,yp=re_means_RPIw
                            ,method='linear')
revecPFmean5_RPIg <- interp2(x=uncer_re
                            ,y=mean_re_RPIg
                            ,Z=re_RPIg_interp_tab5_mean
                            ,xp=re_uncer_RPIg
                            ,yp=re_means_RPIg
                            ,method='linear')

# setting values back to NAs
revecPFmean5_rfc[which(values(res_pred_rfc_max2) %in% NA)] <- NA
revecPFmean5_RPIw[which(values(res_pred_RPIw_max2) %in% NA)] <- NA
revecPFmean5_RPIg[which(values(res_pred_RPIg_max2) %in% NA)] <- NA

# initializing raster for the stored values
re_5_0_5_NA_rfc <- res_pred_rfc_max2
re_5_0_5_NA_RPIw <- res_pred_RPIw_max2
re_5_0_5_NA_RPIg <- res_pred_RPIg_max2

# updating values of the raster
values(re_5_0_5_NA_rfc) <- revecPFmean5_rfc
values(re_5_0_5_NA_RPIw) <- revecPFmean5_RPIw
values(re_5_0_5_NA_RPIg) <- revecPFmean5_RPIg

saveRast(rast=re_5_0_5_NA_rfc
         ,wd=wd_raster_out
         ,rastnm='re_5_0_5_NA_rfc.tif')
saveRast(rast=re_5_0_5_NA_RPIw
         ,wd=wd_raster_out
         ,rastnm='re_5_0_5_NA_RPIw.tif')
saveRast(rast=re_5_0_5_NA_RPIg
         ,wd=wd_raster_out
         ,rastnm='re_5_0_5_NA_RPIg.tif')

makeMap(rast=re_5_0_5_NA_rfc
         ,plotnm='re_5_0_5_NA_rfc.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1)
makeMap(rast=re_5_0_5_NA_RPIw
        ,plotnm='re_5_0_5_NA_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)
makeMap(rast=re_5_0_5_NA_RPIg
        ,plotnm='re_5_0_5_NA_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)

makeMap(rast=re_5_0_5_NA_rfc
        ,plotnm='re_5_0_5_NA_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,leg2 = T
        ,grey = F
        ,County = F
        ,Unit = 'mD-m'
        ,dpi = 300
        ,FigFun = 'tiff'
        ,RawThreshVals = c(bquote(10^.(log10(res_rfc_thresh5)[1])),bquote(10^.(log10(res_rfc_thresh5)[2])),bquote(10^.(log10(res_rfc_thresh5)[3])),bquote(10^.(log10(res_rfc_thresh5)[4])),bquote(10^.(log10(res_rfc_thresh5)[5])),expression(10^4)))
makeMap(rast=re_5_0_5_NA_rfc
        ,plotnm='re_5_0_5_NA_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,leg2 = T, grey = T, County = F, Unit = 'mD-m', dpi = 300, FigFun = 'tiff',
        RawThreshVals = c(bquote(10^.(log10(res_rfc_thresh5)[1])),bquote(10^.(log10(res_rfc_thresh5)[2])),bquote(10^.(log10(res_rfc_thresh5)[3])),bquote(10^.(log10(res_rfc_thresh5)[4])),bquote(10^.(log10(res_rfc_thresh5)[5])),expression(10^4)))

# making uncertainty map
# reading-in tables for interpolated values
re_interp_tab5_rfc <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_rfc_pfvar5.xlsx',sep=''),1,header=FALSE))
re_interp_tab5_RPIw <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIw_pfvar5_NewThresh.xlsx',sep=''),1,header=FALSE))
re_interp_tab5_RPIg <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIg_pfvar5_NewThresh.xlsx',sep=''),1,header=FALSE))

re_interp_tab5_rfc_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_rfc_pfvar5_ls.xlsx',sep=''),1,header=FALSE))
re_interp_tab5_RPIw_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIw_pfvar5_ls_NewThresh.xlsx',sep=''),1,header=FALSE))
re_interp_tab5_RPIg_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/re_RPIg_pfvar5_ls_NewThresh.xlsx',sep=''),1,header=FALSE))

# making the uncertainty maps
# interpolating for the 5 color scheme
revecPFvar5_rfc <- interp2(x=uncer_re
                       ,y=mean_re_rfc
                       ,Z=re_interp_tab5_rfc
                       ,xp=re_uncer_rfc
                       ,yp=re_means_rfc
                       ,method='linear')

revecPFvar5_RPIw <- interp2(x=uncer_re
                           ,y=mean_re_RPIw
                           ,Z=re_interp_tab5_RPIw
                           ,xp=re_uncer_RPIw
                           ,yp=re_means_RPIw
                           ,method='linear')

revecPFvar5_RPIg <- interp2(x=uncer_re
                           ,y=mean_re_RPIg
                           ,Z=re_interp_tab5_RPIg
                           ,xp=re_uncer_RPIg
                           ,yp=re_means_RPIg
                           ,method='linear')

revecPFvar5_rfc_ls <- interp2(x=uncer_re
                       ,y=mean_re_rfc
                       ,Z=re_interp_tab5_rfc_ls
                       ,xp=re_uncer_rfc
                       ,yp=re_means_rfc
                       ,method='linear')

revecPFvar5_RPIw_ls <- interp2(x=uncer_re
                              ,y=mean_re_RPIw
                              ,Z=re_interp_tab5_RPIw_ls
                              ,xp=re_uncer_RPIw
                              ,yp=re_means_RPIw
                              ,method='linear')

revecPFvar5_RPIg_ls <- interp2(x=uncer_re
                              ,y=mean_re_RPIg
                              ,Z=re_interp_tab5_RPIg_ls
                              ,xp=re_uncer_RPIg
                              ,yp=re_means_RPIg
                              ,method='linear')

# setting values back to NAs
revecPFvar5_rfc[which(values(res_pred_rfc_max2) %in% NA)] <- NA
revecPFvar5_RPIw[which(values(res_pred_RPIw_max2) %in% NA)] <- NA
revecPFvar5_RPIg[which(values(res_pred_RPIg_max2) %in% NA)] <- NA
revecPFvar5_rfc_ls[which(values(res_pred_rfc_max2) %in% NA)] <- NA
revecPFvar5_RPIw_ls[which(values(res_pred_RPIw_max2) %in% NA)] <- NA
revecPFvar5_RPIg_ls[which(values(res_pred_RPIg_max2) %in% NA)] <- NA

# converting any values that were -Inf to zero
# -Inf result from zero mean value reservoirs. none in this dataset
# revecPFvar5_rfc[which(re_means_rfc %in% -Inf)] <- 0
# revecPFvar5_RPIw[which(re_means_RPIw %in% -Inf)] <- 0
# revecPFvar5_RPIg[which(re_means_RPIg %in% -Inf)] <- 0
# revecPFvar5_rfc_ls[which(re_means_rfc %in% -Inf)] <- 0
# revecPFvar5_RPIw_ls[which(re_means_RPIw %in% -Inf)] <- 0
# revecPFvar5_RPIg_ls[which(re_means_RPIg %in% -Inf)] <- 0

# initializing raster for the stored values
re_pfa_var5_rfc <- res_pred_rfc_max_err
re_pfa_var5_RPIw <- res_pred_RPIw_max_err
re_pfa_var5_RPIg <- res_pred_RPIg_max_err
re_pfa_var5_rfc_ls <- res_pred_rfc_max_err
re_pfa_var5_RPIw_ls <- res_pred_RPIw_max_err
re_pfa_var5_RPIg_ls <- res_pred_RPIg_max_err

# updating values of the raster
# note: setValues() did not work
values(re_pfa_var5_rfc) <- revecPFvar5_rfc
values(re_pfa_var5_RPIw) <- revecPFvar5_RPIw
values(re_pfa_var5_RPIg) <- revecPFvar5_RPIg
values(re_pfa_var5_rfc_ls) <- revecPFvar5_rfc_ls
values(re_pfa_var5_RPIw_ls) <- revecPFvar5_RPIw_ls
values(re_pfa_var5_RPIg_ls) <- revecPFvar5_RPIg_ls

# saving rasters and making maps
saveRast(rast=re_pfa_var5_rfc
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_rfc.tif')

saveRast(rast=re_pfa_var5_RPIw
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_RPIw.tif')

saveRast(rast=re_pfa_var5_RPIg
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_RPIg.tif')

saveRast(rast=re_pfa_var5_rfc_ls
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_rfc_ls.tif')

saveRast(rast=re_pfa_var5_RPIw_ls
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_RPIw_ls.tif')

saveRast(rast=re_pfa_var5_RPIg_ls
         ,wd=wd_raster_out
         ,rastnm='re_pfa_var5_RPIg_ls.tif')

# saving rasters and making maps
saveRast(rast=calc(re_pfa_var5_rfc,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_rfc.tif')
makeMap(rast=calc(re_pfa_var5_rfc,fun=sqrt)
         ,plotnm='re_pfa_sd5_rfc.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)
makeMap(rast=calc(re_pfa_var5_rfc,fun=sqrt)
        ,plotnm='re_pfa_sd5_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(re_pfa_var5_rfc,fun=sqrt)
        ,plotnm='re_pfa_sd5_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE
        ,grey = T, County = F, dpi = 600, FigFun = 'tiff')

saveRast(rast=calc(re_pfa_var5_RPIw,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_RPIw.tif')
makeMap(rast=calc(re_pfa_var5_RPIw,fun=sqrt)
        ,plotnm='re_pfa_sd5_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

saveRast(rast=calc(re_pfa_var5_RPIg,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_RPIg.tif')
makeMap(rast=calc(re_pfa_var5_RPIg,fun=sqrt)
        ,plotnm='re_pfa_sd5_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

saveRast(rast=calc(re_pfa_var5_rfc_ls,fun=sqrt)
          ,wd=wd_raster_out
          ,rastnm='re_pfa_sd5_rfc_ls.tif')
makeMap(rast=calc(re_pfa_var5_rfc_ls,fun=sqrt)
         ,plotnm='re_pfa_sd5_rfc_ls.png'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1
         ,sdMap=TRUE)
makeMap(rast=calc(re_pfa_var5_rfc_ls,fun=sqrt)
        ,plotnm='re_pfa_sd5_rfc_ls.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(re_pfa_var5_rfc_ls,fun=sqrt)
        ,plotnm='re_pfa_sd5_rfc_ls_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE
        ,grey = T, County = F, dpi = 600, FigFun = 'tiff')

saveRast(rast=calc(re_pfa_var5_RPIw_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_RPIw_ls.tif')
makeMap(rast=calc(re_pfa_var5_RPIw_ls,fun=sqrt)
        ,plotnm='re_pfa_sd5_RPIw_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

saveRast(rast=calc(re_pfa_var5_RPIg_ls,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='re_pfa_sd5_RPIg_ls.tif')
makeMap(rast=calc(re_pfa_var5_RPIg_ls,fun=sqrt)
        ,plotnm='re_pfa_sd5_RPIg_ls.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE)

# removing unneeded files
rm(uncer_re,mean_re_rfc,mean_re_RPIw,mean_re_RPIg)
rm(res_err_3540_rfc,res_err_3035_rfc,res_err_2530_rfc,res_err_1520_rfc,res_err_2025_rfc,res_err_1015_rfc)
rm(res_err_3540_RPIw,res_err_3035_RPIw,res_err_2530_RPIw,res_err_1520_RPIw,res_err_2025_RPIw,res_err_1015_RPIw)
rm(res_err_3540_RPIg,res_err_3035_RPIg,res_err_2530_RPIg,res_err_1520_RPIg,res_err_2025_RPIg,res_err_1015_RPIg)
rm(res_pred_l1000_rfc,res_pred_3540_rfc,res_pred_1520_rfc,res_pred_3035_rfc,res_pred_2530_rfc,res_pred_2025_rfc,res_pred_1015_rfc)
rm(res_pred_l1000_RPIw,res_pred_3540_RPIw,res_pred_1520_RPIw,res_pred_3035_RPIw,res_pred_2530_RPIw,res_pred_2025_RPIw,res_pred_1015_RPIw)
rm(res_pred_l1000_RPIg,res_pred_3540_RPIg,res_pred_1520_RPIg,res_pred_3035_RPIg,res_pred_2530_RPIg,res_pred_2025_RPIg,res_pred_1015_RPIg)
rm(res_pred_rfc_max,res_pred_rfc)
rm(res_pred_RPIw_max,res_pred_RPIw)
rm(res_pred_RPIg_max,res_pred_RPIg)
rm(revecPFvar5_rfc)
rm(revecPFvar5_RPIw)
rm(revecPFvar5_RPIg)
rm(revecPFvar5_rfc_ls)
rm(revecPFvar5_RPIw_ls)
rm(revecPFvar5_RPIg_ls)
rm(res_err_rfc,res_err_l1000_rfc,re_means_rfc,re_uncer_rfc)
rm(res_err_RPIw,res_err_l1000_RPIw,re_means_RPIw,re_uncer_RPIw)
rm(res_err_RPIg,res_err_l1000_RPIg,re_means_RPIg,re_uncer_RPIg)
rm(re_interp_tab5_rfc, re_interp_tab5_rfc_ls, re_interp_tab5_RPIw, re_interp_tab5_RPIw_ls, re_interp_tab5_RPIg, re_interp_tab5_RPIg_ls)
rm(re_rfc_interp_tab5_mean, re_RPIw_interp_tab5_mean, re_RPIg_interp_tab5_mean)
rm(revecPFmean5_rfc, revecPFmean5_RPIg, revecPFmean5_RPIw)

##### UTILIZATION ####
# utilization prediciton and error
util_pred <- raster(paste(wd_raster_in,'/Utilization/slcoh_p4.tif',sep=''))

# replacing no data (-9999) with NA
util_pred[(util_pred %in% -9999)] <- NA
#util_pred[(util_pred %in% 0)] <- NA #None of them are 0

# making up a utilization error of 5%
util_err <- util_pred
util_err[is.na(util_err) == FALSE] <- 5

# creating buffer images
makeWeightBuf(dist=5 #pixels
              ,wd=wd_image
              ,plotnm='ut_buf_5.png')

#Load interpolation table for the mean
ut_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/ut_slcoh_pfmean5.xlsx',sep=''),1,header=FALSE))

#Set the range of means and variances of the interpolation tables
mean_util <- c(seq(5,65,by=2),900,1100) # range of means
std_util_pct <- seq(1,10,by=0.5) # range of standard deviations

# values to interpolate the scaled mean
util_pred2 <- util_pred
util_pred2[which(values(util_pred) %in% NA)] <- 1100

ut_means <- values(util_pred2) 
ut_ses <- rep(5,length(ut_means))

# five color
utvecPFmean5 <- interp2(x=std_util_pct
                        ,y=mean_util
                        ,Z=ut_interp_tab5_mean
                        ,xp=ut_ses
                        ,yp=ut_means
                        ,method='linear')

# setting values back to NAs
utvecPFmean5[which(values(util_pred2) == 1100)] <- NA

#remove unneeded data
rm(util_pred2, ut_means)

# initializing raster for the stored values
ut0_5_0_5_NA <- util_pred

# updating values of the raster
values(ut0_5_0_5_NA) <- utvecPFmean5

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
makeMap(rast=ut0_5_0_5_NA 
        ,plotnm='ut0_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F, leg2 = T, dpi = 300, FigFun = 'tiff', County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU')
makeMap(rast=ut0_5_0_5_NA 
        ,plotnm='ut0_5_0_5_NA_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F
        ,leg2 = T, grey = T, County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU', dpi = 600, FigFun = 'tiff')

# buffering utilization (5 km)
# five color
# Note that warnings will appear that are a result of all NA locations. This is intended.
ut5_5_0_5_NA <- focal(ut0_5_0_5_NA
                      ,w=makeUtilBufWeight(5)
                      ,fun=max
                      ,na.rm=TRUE
                      ,pad=TRUE
                      ,padValue=NA)

# setting values to NAs. The -Inf result from the focal at NA locations.
ut5_5_0_5_NA[which(values(ut5_5_0_5_NA) %in% -Inf)] <- NA

saveRast(rast=ut5_5_0_5_NA 
         ,wd=wd_raster_out
         ,rastnm='ut5_5_0_5_NA.tif')
makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)
makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F, leg2 = T, dpi = 300, FigFun = 'tiff', County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU')
makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,leg2 = T, grey = T, County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU', dpi = 600, FigFun = 'tiff')

# making uncertainty map
util_interp_tab5 <-  as.matrix(read.xlsx(paste(wd_error_interp,'/ut_slcoh_pfvar5.xlsx',sep=''),1,header=FALSE))
util_interp_tab5_ls <-  as.matrix(read.xlsx(paste(wd_error_interp,'/ut_slcoh_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

#Need to obtain a new mean vector for the buffered data
util_pred2 <- util_pred
util_pred2[which(values(util_pred) %in% NA)] <- 1100 #Use a number to increase speed of computation.
util_pred2 <- calc(util_pred2,fun=function(x){2000-x})
utTemp <- focal(util_pred2
                ,w=makeUtilBufWeight(5)
                ,fun=max
                ,na.rm=TRUE #Buffer has NAs and the pad values are NAs
                ,pad=TRUE
                ,padValue=NA)

ut_means5km <- 2000 - values(utTemp) 
ut_ses <- rep(5,length(ut_means5km))

# interpolating for the 5 color scheme
utilvecPFvar5 <- interp2(x=std_util_pct
                       ,y=mean_util
                       ,Z=util_interp_tab5
                       ,xp=ut_ses
                       ,yp=ut_means5km
                       ,method='linear')

utilvecPFvar5_ls <- interp2(x=std_util_pct
                        ,y=mean_util
                        ,Z=util_interp_tab5_ls
                        ,xp=ut_ses
                        ,yp=ut_means5km
                        ,method='linear')

# setting values back to NAs. The NA values are set above.
utilvecPFvar5[which(values(utTemp) %in% 900)] <- NA
utilvecPFvar5_ls[which(values(utTemp) %in% 900)] <- NA

# initializing raster for the stored values
# updating values of the raster
util_pfa_var5 <- util_err
values(util_pfa_var5) <- utilvecPFvar5

util_pfa_var5_ls <- util_err
values(util_pfa_var5_ls) <- utilvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=util_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='util_pfa_var5.tif')

saveRast(rast=util_pfa_var5_ls
          ,wd=wd_raster_out
          ,rastnm='util_pfa_var5_ls.tif')

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
makeMap(rast=calc(util_pfa_var5,fun=sqrt)
        ,plotnm='util_pfa_sd5.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(util_pfa_var5,fun=sqrt)
        ,plotnm='util_pfa_sd5_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

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
makeMap(rast=calc(util_pfa_var5_ls,fun=sqrt)
        ,plotnm='util_pfa_sd5_ls.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(util_pfa_var5_ls,fun=sqrt)
        ,plotnm='util_pfa_sd5_ls_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

# deleting unneeded files
rm(ut0_5_0_5_NA, utvecPFmean5, utilvecPFvar5, utilvecPFvar5_ls)
rm(mean_util, std_util_pct, ut_ses, ut_means5km)
rm(ut_interp_tab5_mean, util_interp_tab5, util_interp_tab5_ls)
rm(utTemp, util_pred, util_pred2)

##### SEISMIC ######
# reading-in the earthquake based risk rasters
seis_eq_pred <- raster(paste(wd_raster_in,'/Seismic/EQ/eqrisk_2kmp-',sep=''))
seis_eq_err <- raster(paste(wd_raster_in,'/Seismic/EQ/eqrisk_2kme-',sep=''))

# reading-in the stress field based risk rasters
seis_stress_pred <- raster(paste(wd_raster_in,'/Seismic/Stress2km/stressrisk2p-',sep=''))
seis_stress_err <- raster(paste(wd_raster_in,'/Seismic/Stress2km/stressrisk2e-',sep=''))

#Reading in normal angles for the stress:
NormAngs = raster(paste(wd_raster_in,'/Seismic/Stress2km/normang',sep=''))

# setting no data (-9999) to NA
seis_eq_pred[(seis_eq_pred %in% -9999)] <- NA
seis_stress_pred[(seis_stress_pred %in% -9999)] <- NA
NormAngs[(NormAngs %in% -9999)] <- NA

seis_eq_err[(seis_eq_err %in% -9999)] <- NA
seis_stress_err[(seis_stress_err %in% -9999)] <- NA

#Load interpolation table for the mean
seEq_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/seEq_pfmean5.xlsx',sep=''),1,header=FALSE))
seSt_interp_tab5_mean <- as.matrix(read.xlsx(paste(wd_error_interp,'/seSt_pfmean5.xlsx',sep=''),1,header=FALSE))

#Set the range of means and variances of the interpolation tables
mean_seEq <- seq(0,26000,by=200) # range of means
std_seEq <- seq(0,2600,by=100) # range of standard deviations
mean_seSt <- seq(0,180,by=2) # range of means.
std_seSt <- seq(0,310,by=10) # range of standard deviations in degrees

# values to interpolate the scaled mean
seEq_means <- values(seis_eq_pred)
seEq_ses <- values(seis_eq_err)
seSt_means <- values(NormAngs)
seSt_ses <- values(seis_stress_err)

# correction because errors were added when they should have been 
# squared, summed, and then square rooted
seEq_ses[which(seEq_ses%in% 3000)] <-  2550
seEq_ses[which(seEq_ses %in% 2750)] <-  2515

# substituting in values for NAs so interpolation algorithm will not crash
seEq_means[which(seEq_means %in% NA)] <- min(mean_seEq)
seEq_ses[which(seEq_ses %in% NA)] <- min(std_seEq)
seSt_means[which(seSt_means %in% NA)] <- min(mean_seSt)
seSt_ses[which(seSt_ses %in% NA)] <- min(std_seSt)

#Substituting in values for the no-worm locations so the algorithm will not crash
seEq_means[which(seEq_means > 1234566)] <- max(mean_seEq) #results in a scaled mean of 5.
seEq_ses[which(seEq_means > 1234566)] <- 0

# five color
seEqvecPFmean5 <- interp2(x=std_seEq
                        ,y=mean_seEq
                        ,Z=seEq_interp_tab5_mean
                        ,xp=seEq_ses
                        ,yp=seEq_means
                        ,method='linear')

seStvecPFmean5 <- interp2(x=std_seSt
                          ,y=mean_seSt
                          ,Z=seSt_interp_tab5_mean
                          ,xp=seSt_ses
                          ,yp=seSt_means
                          ,method='linear')

# setting values back to NAs
seEqvecPFmean5[which(values(seis_eq_pred) %in% NA)] <- NA
seStvecPFmean5[which(values(NormAngs) %in% NA)] <- NA

# initializing raster for the stored values
seEq_5_0_5_NA <- seis_eq_pred
seSt_5_0_5_NA <- seis_stress_pred

# updating values of the raster
values(seEq_5_0_5_NA) <- seEqvecPFmean5
values(seSt_5_0_5_NA) <- seStvecPFmean5

saveRast(rast=seEq_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='seEq_5_0_5_NA.tif')
makeMap(rast=seEq_5_0_5_NA
        ,plotnm='seEq_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)
makeMap(rast=seEq_5_0_5_NA
        ,plotnm='seEq_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 300, FigFun = 'tiff')
makeMap(rast=seEq_5_0_5_NA
        ,plotnm='seEq_5_0_5_NA_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 600, FigFun = 'tiff', grey = T)

saveRast(rast=seSt_5_0_5_NA
         ,wd=wd_raster_out
         ,rastnm='seSt_5_0_5_NA.tif')
makeMap(rast=seSt_5_0_5_NA
        ,plotnm='seSt_5_0_5_NA.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1)
makeMap(rast=seSt_5_0_5_NA
        ,plotnm='seSt_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 300, FigFun = 'tiff')
makeMap(rast=seSt_5_0_5_NA
        ,plotnm='seSt_5_0_5_NA_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 600, FigFun = 'tiff', grey = T)


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
makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 300, FigFun = 'tiff')
makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, grey = T, dpi = 600, FigFun = 'tiff')

## calcualting uncertainty maps for the play fairway scheme
# reading-in tables for interpolated values
se_stress_interp_tab5 <- as.matrix(read.xlsx(paste(wd_error_interp,'/seSt_pfvar5.xlsx',sep=''),1,header=FALSE))
se_stress_interp_tab5_ls <- as.matrix(read.xlsx(paste(wd_error_interp,'/seSt_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

se_eq_interp_tab5 <-  as.matrix(read.xlsx(paste(wd_error_interp,'/seEq_pfvar5.xlsx',sep=''),1,header=FALSE))
se_eq_interp_tab5_ls <-  as.matrix(read.xlsx(paste(wd_error_interp,'/seEq_pfvar5_ls.xlsx',sep=''),1,header=FALSE))

# for STRESS
# interpolating for the 5 color scheme
seStvecPFvar5 <- interp2(x=std_seSt
                         ,y=mean_seSt
                         ,Z=se_stress_interp_tab5
                         ,xp=seSt_ses
                         ,yp=seSt_means
                         ,method='linear'
                         )

seStvecPFvar5_ls <- interp2(x=std_seSt
                          ,y=mean_seSt
                          ,Z=se_stress_interp_tab5_ls
                          ,xp=seSt_ses
                          ,yp=seSt_means
                          ,method='linear'
                          )

# setting values back to NAs
seStvecPFvar5[which(values(NormAngs) %in% NA)] <- NA
seStvecPFvar5_ls[which(values(NormAngs) %in% NA)] <- NA

# initializing raster for the stored values
seSt_pfa_var5 <- seis_stress_err
seSt_pfa_var5_ls <- seis_stress_err

# updating values of the raster
values(seSt_pfa_var5) <- seStvecPFvar5
values(seSt_pfa_var5_ls) <- seStvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=seSt_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='seSt_pfa_var5.tif')

saveRast(rast=seSt_pfa_var5_ls
          ,wd=wd_raster_out
          ,rastnm='seSt_pfa_var5_ls.tif')

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
makeMap(rast=calc(seSt_pfa_var5,fun=sqrt)
        ,plotnm='seSt_pfa_sd5.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(seSt_pfa_var5,fun=sqrt)
        ,plotnm='seSt_pfa_sd5_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

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
makeMap(rast=calc(seSt_pfa_var5_ls,fun=sqrt)
        ,plotnm='seSt_pfa_sd5_ls.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(seSt_pfa_var5_ls,fun=sqrt)
        ,plotnm='seSt_pfa_sd5_ls_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

# for EARTHQUAKE
# interpolating for the 5 color scheme
seEqvecPFvar5 <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab5
                         ,xp=seEq_ses
                         ,yp=seEq_means
                         ,method='linear'
                         )

seEqvecPFvar5_ls <- interp2(x=std_seEq
                         ,y=mean_seEq
                         ,Z=se_eq_interp_tab5_ls
                         ,xp=seEq_ses
                         ,yp=seEq_means
                         ,method='linear'
                         )

# setting values back to NAs
seEqvecPFvar5[which(values(seis_eq_pred) %in% NA)] <- NA
seEqvecPFvar5_ls[which(values(seis_eq_pred) %in% NA)] <- NA

# initializing raster for the stored values
seEq_pfa_var5 <- seis_eq_err
seEq_pfa_var5_ls <- seis_eq_err

# updating values of the raster
values(seEq_pfa_var5) <- seEqvecPFvar5
values(seEq_pfa_var5_ls) <- seEqvecPFvar5_ls

# saving rasters and making maps
saveRast(rast=seEq_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_var5.tif')

saveRast(rast=seEq_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='seEq_pfa_var5_ls.tif')

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
makeMap(rast=calc(seEq_pfa_var5,fun=sqrt)
        ,plotnm='seEq_pfa_sd5.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(seEq_pfa_var5,fun=sqrt)
        ,plotnm='seEq_pfa_sd5_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

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
makeMap(rast=calc(seEq_pfa_var5_ls,fun=sqrt)
        ,plotnm='seEq_pfa_sd5_ls.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(seEq_pfa_var5_ls,fun=sqrt)
        ,plotnm='seEq_pfa_sd5_ls_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 600, FigFun = 'tiff', grey = T)

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
se_pfa_var5_temp4 <- calc(se_pfa_var5_temp3,fun=function(x){return(ifelse(x > 1/(0.4^2), 1/(0.4^2),x))})
se_pfa_var5_ls <- calc(stack(se_pfa_var5_temp4,se_pfa_var5_temp1),fun=prod)

# saving rasters and making maps
saveRast(rast=se_pfa_var5
         ,wd=wd_raster_out
         ,rastnm='se_pfa_var5.tif')

saveRast(rast=se_pfa_var5_ls
         ,wd=wd_raster_out
         ,rastnm='se_pfa_var5_ls.tif')

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
makeMap(rast=calc(se_pfa_var5,fun=sqrt)
        ,plotnm='se_pfa_sd5.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(se_pfa_var5,fun=sqrt)
        ,plotnm='se_pfa_sd5_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, grey = T, dpi = 600, FigFun = 'tiff')

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
makeMap(rast=calc(se_pfa_var5_ls,fun=sqrt)
        ,plotnm='se_pfa_sd5_ls.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(se_pfa_var5_ls,fun=sqrt)
        ,plotnm='se_pfa_sd5_ls_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=TRUE, grey = T, dpi = 600, FigFun = 'tiff')

# removing unneeded files
rm(seEqvecPFvar5,seStvecPFvar5)
rm(seEqvecPFmean5, seEqvecPFvar5_ls)
rm(seStvecPFmean5, seStvecPFvar5_ls)
rm(seSt_ses,seSt_means,seEq_means,seEq_ses)
rm(mean_seEq,mean_seSt,std_seEq,std_seSt)
rm(se_pfa_var5_temp1,se_pfa_var5_temp2,se_pfa_var5_temp3,se_pfa_var5_temp4)
rm(se_eq_interp_tab5, se_eq_interp_tab5_ls)
rm(se_stress_interp_tab5, se_stress_interp_tab5_ls)
rm(seEq_interp_tab5_mean, seSt_interp_tab5_mean)
rm(se_pfa_var5s)

#### Panel Map of All Risk Factors ####
tiff('All_5_0_5_NA.tiff'
     ,height=8
     ,width=8
     ,units='in'
     ,res=300
)
layout(rbind(c(1,2), c(3,4)))
par(mar = c(4,4,2,2))
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1, leg2 = T, grey = F, County = F, RawThreshVals = therm_thresh5/1000, Unit = 'km', dpi = 300, FigFun = NA
         ,xUnitLoc = 11.5e5, yUnitLine = -15
         ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
         ,ScalePos = 'bottomright')
title('Thermal Resource', cex.main = 1.5)

makeMap(rast=re_5_0_5_NA_rfc
        ,plotnm='re_5_0_5_NA_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,leg2 = T
        ,grey = F
        ,County = F
        ,Unit = 'mD-m'
        ,dpi = 300
        ,FigFun = NA
        ,RawThreshVals = c(bquote(10^.(log10(res_rfc_thresh5)[1])),bquote(10^.(log10(res_rfc_thresh5)[2])),bquote(10^.(log10(res_rfc_thresh5)[3])),bquote(10^.(log10(res_rfc_thresh5)[4])),bquote(10^.(log10(res_rfc_thresh5)[5])),expression(10^4))
        ,xUnitLoc = 11.2e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Natural Reservoirs', cex.main = 1.5)

makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F, leg2 = T, dpi = 300, FigFun = NA, County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU'
        ,xUnitLoc = 11.6e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Surface Utilization', cex.main = 1.5)

makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, dpi = 300, FigFun = NA
        ,xUnitLoc = 11.5e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Seismic Risk', cex.main = 1.5)
dev.off()


tiff('All_5_0_5_NA_Grey.tiff'
     ,height=8
     ,width=8
     ,units='in'
     ,res=600
)
layout(rbind(c(1,2), c(3,4)))
par(mar = c(4,4,2,2))
makeMap (rast=th_5_0_5_NA
         ,plotnm='th_5_0_5_NA.tiff'
         ,wd=wd_image
         ,numCol=5
         ,comTy=NA
         ,numRF=1, leg2 = T, grey = T, County = F, RawThreshVals = therm_thresh5/1000, Unit = 'km', dpi = 600, FigFun = NA
         ,xUnitLoc = 11.5e5, yUnitLine = -15
         ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
         ,ScalePos = 'bottomright')
title('Thermal Resource', cex.main = 1.5)

makeMap(rast=re_5_0_5_NA_rfc
        ,plotnm='re_5_0_5_NA_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,leg2 = T
        ,grey = T
        ,County = F
        ,Unit = 'mD-m'
        ,dpi = 600
        ,FigFun = NA
        ,RawThreshVals = c(bquote(10^.(log10(res_rfc_thresh5)[1])),bquote(10^.(log10(res_rfc_thresh5)[2])),bquote(10^.(log10(res_rfc_thresh5)[3])),bquote(10^.(log10(res_rfc_thresh5)[4])),bquote(10^.(log10(res_rfc_thresh5)[5])),expression(10^4))
        ,xUnitLoc = 11.2e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Natural Reservoirs', cex.main = 1.5)

makeMap(rast=ut5_5_0_5_NA 
        ,plotnm='ut5_5_0_5_NA.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1
        ,sdMap=F, leg2 = T, grey = T, dpi = 600, FigFun = NA, County = F, RawThreshVals = util_thresh5, Unit = '$/MMBTU'
        ,xUnitLoc = 11.6e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Surface Utilization', cex.main = 1.5)

makeMap(rast=se_5_0_5_a
        ,plotnm='se_5_0_5_a.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy=NA
        ,numRF=1, grey = T, dpi = 600, FigFun = NA
        ,xUnitLoc = 11.5e5, yUnitLine = -15
        ,ScaleLegmgp = c(-3,-1,0), UnscaleLegmgp = c(3,0.75,0)
        ,ScalePos = 'bottomright')
title('Seismic Risk', cex.main = 1.5)
dev.off()

##### combining maps, all variables ###----
# creating a stacked raster
re_5_0_5_NA_rfc <- raster(paste(wd_raster_out,'re_5_0_5_NA_rfc.tif',sep='/'))
re_5_0_5_NA_RPIw <- raster(paste(wd_raster_out,'re_5_0_5_NA_RPIw.tif',sep='/'))
re_5_0_5_NA_RPIg <- raster(paste(wd_raster_out,'re_5_0_5_NA_RPIg.tif',sep='/'))
th_5_0_5_NA <- raster(paste(wd_raster_out,'th_5_0_5_NA.tif',sep='/'))
ut5_5_0_5_NA <- raster(paste(wd_raster_out,'ut5_5_0_5_NA.tif',sep='/'))
se_5_0_5_a <- raster(paste(wd_raster_out,'se_5_0_5_a.tif',sep='/'))

re_pfa_var5_rfc <- raster(paste(wd_raster_out,'re_pfa_var5_rfc.tif',sep='/'))
re_pfa_var5_RPIw <- raster(paste(wd_raster_out,'re_pfa_var5_RPIw.tif',sep='/'))
re_pfa_var5_RPIg <- raster(paste(wd_raster_out,'re_pfa_var5_RPIg.tif',sep='/'))
th_pfa_var5  <- raster(paste(wd_raster_out,'th_pfa_var5.tif',sep='/'))
util_pfa_var5  <- raster(paste(wd_raster_out,'util_pfa_var5.tif',sep='/'))
se_pfa_var5  <- raster(paste(wd_raster_out,'se_pfa_var5.tif',sep='/'))

re_pfa_var5_rfc_ls <- raster(paste(wd_raster_out,'re_pfa_var5_rfc_ls.tif',sep='/'))
re_pfa_var5_RPIw_ls <- raster(paste(wd_raster_out,'re_pfa_var5_RPIw_ls.tif',sep='/'))
re_pfa_var5_RPIg_ls <- raster(paste(wd_raster_out,'re_pfa_var5_RPIg_ls.tif',sep='/'))
th_pfa_var5_ls  <- raster(paste(wd_raster_out,'th_pfa_var5_ls.tif',sep='/'))
util_pfa_var5_ls  <- raster(paste(wd_raster_out,'util_pfa_var5_ls.tif',sep='/'))
se_pfa_var5_ls  <- raster(paste(wd_raster_out,'se_pfa_var5_ls.tif',sep='/'))

comb_pfa5_rfc <- stack(c(re_5_0_5_NA_rfc,th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))
comb_pfa5_RPIw <- stack(c(re_5_0_5_NA_RPIw,th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))
comb_pfa5_RPIg <- stack(c(re_5_0_5_NA_RPIg,th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))

comb_pfa5_rfc_geo <- stack(c(re_5_0_5_NA_rfc,th_5_0_5_NA,se_5_0_5_a))
comb_pfa5_RPIw_geo <- stack(c(re_5_0_5_NA_RPIw,th_5_0_5_NA,se_5_0_5_a))
comb_pfa5_RPIg_geo <- stack(c(re_5_0_5_NA_RPIg,th_5_0_5_NA,se_5_0_5_a))

comb_pfa5_egs <- stack(c(th_5_0_5_NA,ut5_5_0_5_NA,se_5_0_5_a))


#All values less than 0 should be set to NA
comb_pfa5_rfc[comb_pfa5_rfc < 0] <- NA
comb_pfa5_RPIw[comb_pfa5_RPIw < 0] <- NA
comb_pfa5_RPIg[comb_pfa5_RPIg < 0] <- NA

comb_pfa5_rfc_geo[comb_pfa5_rfc_geo < 0] <- NA
comb_pfa5_RPIw_geo[comb_pfa5_RPIw_geo < 0] <- NA
comb_pfa5_RPIg_geo[comb_pfa5_RPIg_geo < 0] <- NA

comb_pfa5_egs[comb_pfa5_egs < 0] <- NA

## using mean
co_5_0_20_s_rfc <- calc(comb_pfa5_rfc,fun=mean,na.rm=FALSE)
co_5_0_20_s_RPIw <- calc(comb_pfa5_RPIw,fun=mean,na.rm=FALSE)
co_5_0_20_s_RPIg <- calc(comb_pfa5_RPIg,fun=mean,na.rm=FALSE)

co_5_0_16_s_rfc_geo <- calc(comb_pfa5_rfc_geo,fun=mean,na.rm=FALSE)
co_5_0_16_s_RPIw_geo <- calc(comb_pfa5_RPIw_geo,fun=mean,na.rm=FALSE)
co_5_0_16_s_RPIg_geo <- calc(comb_pfa5_RPIg_geo,fun=mean,na.rm=FALSE)

co_5_0_16_s_egs <- calc(comb_pfa5_egs,fun=mean,na.rm=FALSE)

#All RFs
saveRast(rast=co_5_0_20_s_rfc
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_rfc.tif')
makeMap(rast=co_5_0_20_s_rfc
        ,plotnm='co_5_0_5_a_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=4)
makeMap(rast=co_5_0_20_s_rfc
        ,plotnm='co_5_0_5_a_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = 'tiff')
makeMap(rast=co_5_0_20_s_rfc
        ,plotnm='co_5_0_5_a_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = T, dpi = 600, FigFun = 'tiff')

saveRast(rast=co_5_0_20_s_RPIw
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_RPIw.tif')
makeMap(rast=co_5_0_20_s_RPIw
        ,plotnm='co_5_0_5_a_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=4)

saveRast(rast=co_5_0_20_s_RPIg
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_RPIg.tif')
makeMap(rast=co_5_0_20_s_RPIg
        ,plotnm='co_5_0_5_a_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=4)

# Geology only
saveRast(rast=co_5_0_16_s_rfc_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_rfc_geo.tif')
makeMap(rast=co_5_0_16_s_rfc_geo
        ,plotnm='co_5_0_5_a_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=3)

saveRast(rast=co_5_0_16_s_RPIw_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_RPIw_geo.tif')
makeMap(rast=co_5_0_16_s_RPIw_geo
        ,plotnm='co_5_0_5_a_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=3)

saveRast(rast=co_5_0_16_s_RPIg_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_RPIg_geo.tif')
makeMap(rast=co_5_0_16_s_RPIg_geo
        ,plotnm='co_5_0_5_a_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=3)


# EGS
saveRast(rast=co_5_0_16_s_egs
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_a_egs.tif')
makeMap(rast=co_5_0_16_s_egs
        ,plotnm='co_5_0_5_a_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:3
        ,numRF=3)

## using geometric mean
co_5_0_625_p2_rfc <- calc(comb_pfa5_rfc,fun=prod,na.rm=FALSE)
co_5_0_625_p2_RPIw <- calc(comb_pfa5_RPIw,fun=prod,na.rm=FALSE)
co_5_0_625_p2_RPIg <- calc(comb_pfa5_RPIg,fun=prod,na.rm=FALSE)

co_5_0_125_p2_rfc_geo <- calc(comb_pfa5_rfc_geo,fun=prod,na.rm=FALSE)
co_5_0_125_p2_RPIw_geo <- calc(comb_pfa5_RPIw_geo,fun=prod,na.rm=FALSE)
co_5_0_125_p2_RPIg_geo <- calc(comb_pfa5_RPIg_geo,fun=prod,na.rm=FALSE)

co_5_0_125_p2_egs <- calc(comb_pfa5_egs,fun=prod,na.rm=FALSE)

co_5_0_625_p_rfc <- calc(co_5_0_625_p2_rfc,fun=function(x){return(x^0.25)})
co_5_0_625_p_RPIw <- calc(co_5_0_625_p2_RPIw,fun=function(x){return(x^0.25)})
co_5_0_625_p_RPIg <- calc(co_5_0_625_p2_RPIg,fun=function(x){return(x^0.25)})

co_5_0_125_p_rfc_geo <- calc(co_5_0_125_p2_rfc_geo,fun=function(x){return(x^(1/3))})
co_5_0_125_p_RPIw_geo <- calc(co_5_0_125_p2_RPIw_geo,fun=function(x){return(x^(1/3))})
co_5_0_125_p_RPIg_geo <- calc(co_5_0_125_p2_RPIg_geo,fun=function(x){return(x^(1/3))})

co_5_0_125_p_egs <- calc(co_5_0_125_p2_egs,fun=function(x){return(x^(1/3))})

rm(co_5_0_125_p2_egs, co_5_0_125_p2_rfc_geo, co_5_0_125_p2_RPIw_geo, co_5_0_125_p2_RPIg_geo)
rm(co_5_0_625_p2_rfc, co_5_0_625_p2_RPIw, co_5_0_625_p2_RPIg)

#All RFs
setwd(wd_image)
saveRast(rast=co_5_0_625_p_rfc
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_rfc.tif')
makeMap(rast=co_5_0_625_p_rfc
        ,plotnm='co_5_0_5_g_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4)
makeMap(rast=co_5_0_625_p_rfc
        ,plotnm='co_5_0_5_g_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = 'tiff')
makeMap(rast=co_5_0_625_p_rfc
        ,plotnm='co_5_0_5_g_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = T, dpi = 600, FigFun = 'tiff')

saveRast(rast=co_5_0_625_p_RPIw
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_RPIw.tif')
makeMap(rast=co_5_0_625_p_RPIw
        ,plotnm='co_5_0_5_g_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4)

saveRast(rast=co_5_0_625_p_RPIg
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_RPIg.tif')
makeMap(rast=co_5_0_625_p_RPIg
        ,plotnm='co_5_0_5_g_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4)

# Geology Only
saveRast(rast=co_5_0_125_p_rfc_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_rfc_geo.tif')
makeMap(rast=co_5_0_125_p_rfc_geo
        ,plotnm='co_5_0_5_g_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=3)

saveRast(rast=co_5_0_125_p_RPIw_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_RPIw_geo.tif')
makeMap(rast=co_5_0_125_p_RPIw_geo
        ,plotnm='co_5_0_5_g_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=3)

saveRast(rast=co_5_0_125_p_RPIg_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_RPIg_geo.tif')
makeMap(rast=co_5_0_125_p_RPIg_geo
        ,plotnm='co_5_0_5_g_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=3)

# EGS
saveRast(rast=co_5_0_125_p_egs
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_g_egs.tif')
makeMap(rast=co_5_0_125_p_egs
        ,plotnm='co_5_0_5_g_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:3
        ,numRF=3)

## using minimums
co_5_0_5_m_rfc <- calc(comb_pfa5_rfc,fun=min,na.rm=FALSE)
co_5_0_5_m_RPIw <- calc(comb_pfa5_RPIw,fun=min,na.rm=FALSE)
co_5_0_5_m_RPIg <- calc(comb_pfa5_RPIg,fun=min,na.rm=FALSE)

co_5_0_5_m_rfc_geo <- calc(comb_pfa5_rfc_geo,fun=min,na.rm=FALSE)
co_5_0_5_m_RPIw_geo <- calc(comb_pfa5_RPIw_geo,fun=min,na.rm=FALSE)
co_5_0_5_m_RPIg_geo <- calc(comb_pfa5_RPIg_geo,fun=min,na.rm=FALSE)

co_5_0_5_m_egs <- calc(comb_pfa5_egs,fun=min,na.rm=FALSE)

#All RFs
setwd(wd_image)
saveRast(rast=co_5_0_5_m_rfc
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_rfc.tif')
makeMap(rast=co_5_0_5_m_rfc
        ,plotnm='co_5_0_5_m_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=4)
makeMap(rast=co_5_0_5_m_rfc
        ,plotnm='co_5_0_5_m_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = 'tiff')
makeMap(rast=co_5_0_5_m_rfc
        ,plotnm='co_5_0_5_m_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = T, dpi = 600, FigFun = 'tiff')

saveRast(rast=co_5_0_5_m_RPIw
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIw.tif')
makeMap(rast=co_5_0_5_m_RPIw
        ,plotnm='co_5_0_5_m_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=4)

saveRast(rast=co_5_0_5_m_RPIg
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIg.tif')
makeMap(rast=co_5_0_5_m_RPIg
        ,plotnm='co_5_0_5_m_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=4)

# Geology Only
saveRast(rast=co_5_0_5_m_rfc_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_rfc_geo.tif')
makeMap(rast=co_5_0_5_m_rfc_geo
        ,plotnm='co_5_0_5_m_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=3)

saveRast(rast=co_5_0_5_m_RPIw_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIw_geo.tif')
makeMap(rast=co_5_0_5_m_RPIw_geo
        ,plotnm='co_5_0_5_m_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=3)

saveRast(rast=co_5_0_5_m_RPIg_geo
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIg_geo.tif')
makeMap(rast=co_5_0_5_m_RPIg_geo
        ,plotnm='co_5_0_5_m_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:5
        ,numRF=3)

# EGS
saveRast(rast=co_5_0_5_m_egs
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_egs.tif')
makeMap(rast=co_5_0_5_m_egs
        ,plotnm='co_5_0_5_m_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the min is on 0:3
        ,numRF=3)

#removing unneeded data
rm(comb_pfa5_egs, comb_pfa5_rfc, comb_pfa5_rfc_geo, comb_pfa5_RPIw,
   comb_pfa5_RPIw_geo, comb_pfa5_RPIg, comb_pfa5_RPIg_geo)


##### COMBINED UNCERTAINTY -----
# for average
uncer_temp <- stack(re_pfa_var5_rfc,th_pfa_var5,se_pfa_var5,util_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_rfc <- calc(co_uncer_avg1,fun=function(x){return(x/16)})

uncer_temp <- stack(re_pfa_var5_RPIw,th_pfa_var5,se_pfa_var5,util_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_RPIw <- calc(co_uncer_avg1,fun=function(x){return(x/16)})

uncer_temp <- stack(re_pfa_var5_RPIg,th_pfa_var5,se_pfa_var5,util_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_RPIg <- calc(co_uncer_avg1,fun=function(x){return(x/16)})

#Geology only
uncer_temp <- stack(re_pfa_var5_rfc,th_pfa_var5,se_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_rfc_geo <- calc(co_uncer_avg1,fun=function(x){return(x/9)})

uncer_temp <- stack(re_pfa_var5_RPIw,th_pfa_var5,se_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_RPIw_geo <- calc(co_uncer_avg1,fun=function(x){return(x/9)})

uncer_temp <- stack(re_pfa_var5_RPIg,th_pfa_var5,se_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_RPIg_geo <- calc(co_uncer_avg1,fun=function(x){return(x/9)})


#EGS
uncer_temp <- stack(th_pfa_var5,se_pfa_var5,util_pfa_var5)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_avg1 <- calc(uncer_temp,fun=sum)
co_pfa_var5_avg_egs <- calc(co_uncer_avg1,fun=function(x){return(x/9)}) #Divide by 9 because each RF is divided by 3

rm(co_uncer_avg1, uncer_temp)

# saving rasters and making maps
#All RFs
saveRast(rast=co_pfa_var5_avg_rfc
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_rfc.tif')
saveRast(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_rfc.tif')
makeMap(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)
makeMap(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = T, dpi = 600, FigFun = 'tiff')

saveRast(rast=co_pfa_var5_avg_RPIw
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_RPIw.tif')
saveRast(rast=calc(co_pfa_var5_avg_RPIw,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_RPIw.tif')
makeMap(rast=calc(co_pfa_var5_avg_RPIw,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

saveRast(rast=co_pfa_var5_avg_RPIg
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_RPIg.tif')
saveRast(rast=calc(co_pfa_var5_avg_RPIg,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_RPIg.tif')
makeMap(rast=calc(co_pfa_var5_avg_RPIg,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

# Geology Only
saveRast(rast=co_pfa_var5_avg_rfc_geo
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_rfc_geo.tif')
saveRast(rast=calc(co_pfa_var5_avg_rfc_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_rfc_geo.tif')
makeMap(rast=calc(co_pfa_var5_avg_rfc_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=co_pfa_var5_avg_RPIw_geo
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_RPIw_geo.tif')
saveRast(rast=calc(co_pfa_var5_avg_RPIw_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_RPIw_geo.tif')
makeMap(rast=calc(co_pfa_var5_avg_RPIw_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=co_pfa_var5_avg_RPIg_geo
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_RPIg_geo.tif')
saveRast(rast=calc(co_pfa_var5_avg_RPIg_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_RPIg_geo.tif')
makeMap(rast=calc(co_pfa_var5_avg_RPIg_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)


# EGS
saveRast(rast=co_pfa_var5_avg_egs
         ,wd=wd_raster_out
         ,rastnm='co_pfa_var5_avg_egs.tif')
saveRast(rast=calc(co_pfa_var5_avg_egs,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_avg_egs.tif')
makeMap(rast=calc(co_pfa_var5_avg_egs,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

# for geometric mean
# All RFs
uncer_temp <- stack(re_pfa_var5_rfc_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/16)})
co_pfa_var5_geomean_rfc <- calc(stack(co_5_0_625_p_rfc,co_5_0_625_p_rfc,co_uncer_geomean2),fun=prod) #Multiply by the squared scaled mean to transform to real space

uncer_temp <- stack(re_pfa_var5_RPIw_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/16)})
co_pfa_var5_geomean_RPIw <- calc(stack(co_5_0_625_p_RPIw,co_5_0_625_p_RPIw,co_uncer_geomean2),fun=prod)

uncer_temp <- stack(re_pfa_var5_RPIg_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/16)})
co_pfa_var5_geomean_RPIg <- calc(stack(co_5_0_625_p_RPIg,co_5_0_625_p_RPIg,co_uncer_geomean2),fun=prod)

# Geology Only
uncer_temp <- stack(re_pfa_var5_rfc_ls,th_pfa_var5_ls,se_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/9)})
co_pfa_var5_geomean_rfc_geo <- calc(stack(co_5_0_125_p_rfc_geo,co_5_0_125_p_rfc_geo,co_uncer_geomean2),fun=prod)

uncer_temp <- stack(re_pfa_var5_RPIw_ls,th_pfa_var5_ls,se_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/9)})
co_pfa_var5_geomean_RPIw_geo <- calc(stack(co_5_0_125_p_RPIw_geo,co_5_0_125_p_RPIw_geo,co_uncer_geomean2),fun=prod)

uncer_temp <- stack(re_pfa_var5_RPIg_ls,th_pfa_var5_ls,se_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/9)})
co_pfa_var5_geomean_RPIg_geo <- calc(stack(co_5_0_125_p_RPIg_geo,co_5_0_125_p_RPIg_geo,co_uncer_geomean2),fun=prod)

# EGS
uncer_temp <- stack(th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls)
uncer_temp[uncer_temp < 0] <- NA
co_uncer_geomean1 <- calc(uncer_temp,fun=sum)
co_uncer_geomean2 <- calc(co_uncer_geomean1,fun=function(x){return(x/9)})
co_pfa_var5_geomean_egs <- calc(stack(co_5_0_125_p_egs,co_5_0_125_p_egs,co_uncer_geomean2),fun=prod)

rm(co_uncer_geomean1, co_uncer_geomean2, uncer_temp)

# saving rasters and making maps. Note that variances are so small that they require modification to the makeMap.R function to plot.
# Variance maps are not plotted to avoid modifications to the plotting function.

#All RFs
saveRast(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_rfc.tif')
makeMap(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)
makeMap(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = 'tiff')
makeMap(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc_Grey.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = T, dpi = 600, FigFun = 'tiff')

saveRast(rast=calc(co_pfa_var5_geomean_RPIw,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_RPIw.tif')
makeMap(rast=calc(co_pfa_var5_geomean_RPIw,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

saveRast(rast=calc(co_pfa_var5_geomean_RPIg,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_RPIg.tif')
makeMap(rast=calc(co_pfa_var5_geomean_RPIg,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)


# Geology Only
saveRast(rast=calc(co_pfa_var5_geomean_rfc_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_rfc_geo.tif')
makeMap(rast=calc(co_pfa_var5_geomean_rfc_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(co_pfa_var5_geomean_RPIw_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_RPIw_geo.tif')
makeMap(rast=calc(co_pfa_var5_geomean_RPIw_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

saveRast(rast=calc(co_pfa_var5_geomean_RPIg_geo,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_RPIg_geo.tif')
makeMap(rast=calc(co_pfa_var5_geomean_RPIg_geo,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)


# EGS
saveRast(rast=calc(co_pfa_var5_geomean_egs,fun=sqrt)
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_geomean_egs.tif')
makeMap(rast=calc(co_pfa_var5_geomean_egs,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

#removing unneeded data
rm(co_pfa_var5_avg_rfc, co_pfa_var5_avg_rfc_geo, co_pfa_var5_avg_RPIw,
   co_pfa_var5_avg_RPIw_geo, co_pfa_var5_avg_RPIg, co_pfa_var5_avg_RPIg_geo)
rm(co_pfa_var5_geomean_rfc, co_pfa_var5_geomean_rfc_geo, co_pfa_var5_geomean_RPIw,
   co_pfa_var5_geomean_RPIw_geo, co_pfa_var5_geomean_RPIg, co_pfa_var5_geomean_RPIg_geo)


#### Panel plot of all FMs ####
tiff('co_pfa_all_rfc.tiff'
     ,height=8
     ,width=8
     ,units='in'
     ,res=300
)
layout(rbind(c(1,2), c(3,4)))
par(mar = c(4,4,2,2))
makeMap(rast=co_5_0_20_s_rfc
        ,plotnm='co_5_0_5_a_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Favorability Metric: FM'['avg'])), cex.main = 1.5)

makeMap(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Standard Deviation of FM'['avg'])), cex.main = 1.5)

makeMap(rast=co_5_0_625_p_rfc
        ,plotnm='co_5_0_5_g_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Favorability Metric: FM'['gm'])), cex.main = 1.5)

makeMap(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = F, dpi = 300, FigFun = NA        
        ,ScalePos = 'bottomright')
title(expression(bold('Standard Deviation of FM'['gm'])), cex.main = 1.5)
dev.off()


tiff('co_pfa_all_rfc_Grey.tiff'
     ,height=8
     ,width=8
     ,units='in'
     ,res=600
)
layout(rbind(c(1,2), c(3,4)))
par(mar = c(4,4,2,2))
makeMap(rast=co_5_0_20_s_rfc
        ,plotnm='co_5_0_5_a_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the mean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = T, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Favorability Metric: FM'['avg'])), cex.main = 1.5)

makeMap(rast=calc(co_pfa_var5_avg_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_avg_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = T, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Standard Deviation of FM'['avg'])), cex.main = 1.5)

makeMap(rast=co_5_0_625_p_rfc
        ,plotnm='co_5_0_5_g_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min' #Use min because the geomean is on 0:5
        ,numRF=4
        ,leg2 = F, County = F, grey = T, dpi = 300, FigFun = NA
        ,ScalePos = 'bottomright')
title(expression(bold('Favorability Metric: FM'['gm'])), cex.main = 1.5)

makeMap(rast=calc(co_pfa_var5_geomean_rfc,fun=sqrt)
        ,plotnm='co_pfa_sd5_geomean_rfc.tiff'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE
        ,leg2 = F, County = F, grey = T, dpi = 300, FigFun = NA        
        ,ScalePos = 'bottomright')
title(expression(bold('Standard Deviation of FM'['gm'])), cex.main = 1.5)
dev.off()

##### Extracting City Info----
# extracting values of layers for cities
# reading-in the scaled rasters
# note that scaled rasters are in the wd_raster_out
# because they were output from the first part of the code

#Five Color:
reservoir_rfc <- raster(paste(wd_raster_out,'/re_5_0_5_NA_rfc.tif',sep=''))
reservoir_rfc[reservoir_rfc < 0] <- NA

reservoir_RPIw <- raster(paste(wd_raster_out,'/re_5_0_5_NA_RPIw.tif',sep=''))
reservoir_RPIw[reservoir_RPIw < 0] <- NA

reservoir_RPIg <- raster(paste(wd_raster_out,'/re_5_0_5_NA_RPIg.tif',sep=''))
reservoir_RPIg[reservoir_RPIg < 0] <- NA

re_pfa_var5_rfc <- raster(paste(wd_raster_out,'/re_pfa_var5_rfc.tif',sep=''))
re_pfa_var5_rfc[re_pfa_var5_rfc < 0] <- NA

re_pfa_var5_RPIw <- raster(paste(wd_raster_out,'/re_pfa_var5_RPIw.tif',sep=''))
re_pfa_var5_RPIw[re_pfa_var5_RPIw < 0] <- NA

re_pfa_var5_RPIg <- raster(paste(wd_raster_out,'/re_pfa_var5_RPIg.tif',sep=''))
re_pfa_var5_RPIg[re_pfa_var5_RPIg < 0] <- NA

re_pfa_var5_rfc_ls <- raster(paste(wd_raster_out,'/re_pfa_var5_rfc_ls.tif',sep=''))
re_pfa_var5_rfc_ls[re_pfa_var5_rfc_ls < 0] <- NA

re_pfa_var5_RPIw_ls <- raster(paste(wd_raster_out,'/re_pfa_var5_RPIw_ls.tif',sep=''))
re_pfa_var5_RPIw_ls[re_pfa_var5_RPIw_ls < 0] <- NA

re_pfa_var5_RPIg_ls <- raster(paste(wd_raster_out,'/re_pfa_var5_RPIg_ls.tif',sep=''))
re_pfa_var5_RPIg_ls[re_pfa_var5_RPIg_ls < 0] <- NA

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

util_pred <- raster(paste(wd_raster_in,'/Utilization/slcoh_p4.tif',sep=''))

# replacing no data (-9999) with NA
util_pred[(util_pred %in% -9999)] <- NA
#util_pred[(util_pred %in% 0)] <- NA #None of them are 0

util_pred2 <- util_pred
util_pred2[which(values(util_pred) %in% NA)] <- 1100 #Use a number to increase speed of computation.
util_pred2 <- calc(util_pred2,fun=function(x){2000-x})
utTemp <- focal(util_pred2
                ,w=makeUtilBufWeight(5)
                ,fun=max
                ,na.rm=TRUE #Buffer has NAs and the pad values are NAs
                ,pad=TRUE
                ,padValue=NA)

util_pred <- calc(utTemp,fun=function(x){return(2000-x)})
util_pred <- calc(util_pred,fun=function(x){return(ifelse(x==1100,NA,x))})

#Five Color Combinations:
co_5_0_20_s_rfc <- raster(paste(wd_raster_out,'/co_5_0_5_a_rfc.tif',sep=''))
co_5_0_625_p_rfc <- raster(paste(wd_raster_out,'/co_5_0_5_g_rfc.tif',sep=''))
co_5_0_5_m_rfc <- raster(paste(wd_raster_out,'/co_5_0_5_m_rfc.tif',sep=''))

co_5_0_20_s_rfc[co_5_0_20_s_rfc < 0] <- NA
co_5_0_625_p_rfc[co_5_0_625_p_rfc < 0] <- NA
co_5_0_5_m_rfc[co_5_0_5_m_rfc < 0] <- NA

co_5_0_20_s_RPIw <- raster(paste(wd_raster_out,'/co_5_0_5_a_RPIw.tif',sep=''))
co_5_0_625_p_RPIw <- raster(paste(wd_raster_out,'/co_5_0_5_g_RPIw.tif',sep=''))
co_5_0_5_m_RPIw <- raster(paste(wd_raster_out,'/co_5_0_5_m_RPIw.tif',sep=''))

co_5_0_20_s_RPIw[co_5_0_20_s_RPIw < 0] <- NA
co_5_0_625_p_RPIw[co_5_0_625_p_RPIw < 0] <- NA
co_5_0_5_m_RPIw[co_5_0_5_m_RPIw < 0] <- NA

co_5_0_20_s_RPIg <- raster(paste(wd_raster_out,'/co_5_0_5_a_RPIg.tif',sep=''))
co_5_0_625_p_RPIg <- raster(paste(wd_raster_out,'/co_5_0_5_g_RPIg.tif',sep=''))
co_5_0_5_m_RPIg <- raster(paste(wd_raster_out,'/co_5_0_5_m_RPIg.tif',sep=''))

co_5_0_20_s_RPIg[co_5_0_20_s_RPIg < 0] <- NA
co_5_0_625_p_RPIg[co_5_0_625_p_RPIg < 0] <- NA
co_5_0_5_m_RPIg[co_5_0_5_m_RPIg < 0] <- NA

co_5_0_16_s_rfc_geo <- raster(paste(wd_raster_out,'/co_5_0_5_a_rfc_geo.tif',sep=''))
co_5_0_125_p_rfc_geo <- raster(paste(wd_raster_out,'/co_5_0_5_g_rfc_geo.tif',sep=''))
co_5_0_5_m_rfc_geo <- raster(paste(wd_raster_out,'/co_5_0_5_m_rfc_geo.tif',sep=''))

co_5_0_16_s_rfc_geo[co_5_0_16_s_rfc_geo < 0] <- NA
co_5_0_125_p_rfc_geo[co_5_0_125_p_rfc_geo < 0] <- NA
co_5_0_5_m_rfc_geo[co_5_0_5_m_rfc_geo < 0] <- NA

co_5_0_16_s_RPIw_geo <- raster(paste(wd_raster_out,'/co_5_0_5_a_RPIw_geo.tif',sep=''))
co_5_0_125_p_RPIw_geo <- raster(paste(wd_raster_out,'/co_5_0_5_g_RPIw_geo.tif',sep=''))
co_5_0_5_m_RPIw_geo <- raster(paste(wd_raster_out,'/co_5_0_5_m_RPIw_geo.tif',sep=''))

co_5_0_16_s_RPIw_geo[co_5_0_16_s_RPIw_geo < 0] <- NA
co_5_0_125_p_RPIw_geo[co_5_0_125_p_RPIw_geo < 0] <- NA
co_5_0_5_m_RPIw_geo[co_5_0_5_m_RPIw_geo < 0] <- NA

co_5_0_16_s_RPIg_geo <- raster(paste(wd_raster_out,'/co_5_0_5_a_RPIg_geo.tif',sep=''))
co_5_0_125_p_RPIg_geo <- raster(paste(wd_raster_out,'/co_5_0_5_g_RPIg_geo.tif',sep=''))
co_5_0_5_m_RPIg_geo <- raster(paste(wd_raster_out,'/co_5_0_5_m_RPIg_geo.tif',sep=''))

co_5_0_16_s_RPIg_geo[co_5_0_16_s_RPIg_geo < 0] <- NA
co_5_0_125_p_RPIg_geo[co_5_0_125_p_RPIg_geo < 0] <- NA
co_5_0_5_m_RPIg_geo[co_5_0_5_m_RPIg_geo < 0] <- NA

co_pfa_sd5_avg_rfc = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_rfc.tif',sep=''))
co_pfa_sd5_geomean_rfc = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_rfc.tif',sep=''))
co_pfa_sd5_avg_rfc_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_rfc_geo.tif',sep=''))
co_pfa_sd5_geomean_rfc_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_rfc_geo.tif',sep=''))

co_pfa_sd5_avg_rfc[co_pfa_sd5_avg_rfc < 0] <- NA
co_pfa_sd5_geomean_rfc[co_pfa_sd5_geomean_rfc < 0] <- NA
co_pfa_sd5_avg_rfc_geo[co_pfa_sd5_avg_rfc_geo < 0] <- NA
co_pfa_sd5_geomean_rfc_geo[co_pfa_sd5_geomean_rfc_geo < 0] <- NA

co_pfa_sd5_avg_RPIw = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_RPIw.tif',sep=''))
co_pfa_sd5_geomean_RPIw = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_RPIw.tif',sep=''))
co_pfa_sd5_avg_RPIw_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_RPIw_geo.tif',sep=''))
co_pfa_sd5_geomean_RPIw_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_RPIw_geo.tif',sep=''))

co_pfa_sd5_avg_RPIw[co_pfa_sd5_avg_RPIw < 0] <- NA
co_pfa_sd5_geomean_RPIw[co_pfa_sd5_geomean_RPIw < 0] <- NA
co_pfa_sd5_avg_RPIw_geo[co_pfa_sd5_avg_RPIw_geo < 0] <- NA
co_pfa_sd5_geomean_RPIw_geo[co_pfa_sd5_geomean_RPIw_geo < 0] <- NA

co_pfa_sd5_avg_RPIg = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_RPIg.tif',sep=''))
co_pfa_sd5_geomean_RPIg = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_RPIg.tif',sep=''))
co_pfa_sd5_avg_RPIg_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_avg_RPIg_geo.tif',sep=''))
co_pfa_sd5_geomean_RPIg_geo = raster(paste(wd_raster_out,'/co_pfa_sd5_geomean_RPIg_geo.tif',sep=''))

co_pfa_sd5_avg_RPIg[co_pfa_sd5_avg_RPIg < 0] <- NA
co_pfa_sd5_geomean_RPIg[co_pfa_sd5_geomean_RPIg < 0] <- NA
co_pfa_sd5_avg_RPIg_geo[co_pfa_sd5_avg_RPIg_geo < 0] <- NA
co_pfa_sd5_geomean_RPIg_geo[co_pfa_sd5_geomean_RPIg_geo < 0] <- NA

# stacking rasters that will have their values
# extracted at points
comb_extract_5_rfc <- stack(c(reservoir_rfc,thermal,seismic,utilization
                        ,re_pfa_var5_rfc,th_pfa_var5,se_pfa_var5,util_pfa_var5
                        ,re_pfa_var5_rfc_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls
                        ,res_pred_rfc_max2,res_pred_rfc_max_err
                        ,therm_pred,therm_err
                        ,seis_eq_pred,seis_eq_err,NormAngs,seis_stress_pred,seis_stress_err
                        ,util_pred,util_err
                        ,co_5_0_20_s_rfc,co_5_0_625_p_rfc,co_5_0_5_m_rfc
                        ,co_pfa_sd5_avg_rfc,co_pfa_sd5_geomean_rfc
                        ,co_5_0_16_s_rfc_geo,co_5_0_125_p_rfc_geo,co_5_0_5_m_rfc_geo
                        ,co_pfa_sd5_avg_rfc_geo,co_pfa_sd5_geomean_rfc_geo))

# names for each layer in the raster stack
comb_names_5_rfc <- c('reservoir_rfc','thermal','seismic','utilization'
                ,'re_pfa_var5_rfc','th_pfa_var5','se_pfa_var5','util_pfa_var5'
                ,'re_pfa_var5_rfc_ls','th_pfa_var5_ls','se_pfa_var5_ls','util_pfa_var5_ls'
                ,'res_pred_rfc_max2','res_pred_rfc_max_err'
                ,'therm_pred','therm_err'
                ,'seis_eq_pred','seis_eq_err','NormAngs','seis_stress_pred','seis_stress_err'
                ,'util_pred','util_err'
                ,'co_5_0_20_s_rfc','co_5_0_625_p_rfc','co_5_0_5_m_rfc'
                ,'co_pfa_sd5_avg_rfc','co_pfa_sd5_geomean_rfc'
                ,'co_5_0_16_s_rfc_geo','co_5_0_125_p_rfc_geo','co_5_0_5_m_rfc_geo'
                ,'co_pfa_sd5_avg_rfc_geo','co_pfa_sd5_geomean_rfc_geo')

comb_extract_5_RPIw <- stack(c(reservoir_RPIw,thermal,seismic,utilization
                              ,re_pfa_var5_RPIw,th_pfa_var5,se_pfa_var5,util_pfa_var5
                              ,re_pfa_var5_RPIw_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls
                              ,res_pred_RPIw_max2,res_pred_RPIw_max_err
                              ,therm_pred,therm_err
                              ,seis_eq_pred,seis_eq_err,NormAngs,seis_stress_pred,seis_stress_err
                              ,util_pred,util_err
                              ,co_5_0_20_s_RPIw,co_5_0_625_p_RPIw,co_5_0_5_m_RPIw
                              ,co_pfa_sd5_avg_RPIw,co_pfa_sd5_geomean_RPIw
                              ,co_5_0_16_s_RPIw_geo,co_5_0_125_p_RPIw_geo,co_5_0_5_m_RPIw_geo
                              ,co_pfa_sd5_avg_RPIw_geo,co_pfa_sd5_geomean_RPIw_geo))

# names for each layer in the raster stack
comb_names_5_RPIw <- c('reservoir_RPIw','thermal','seismic','utilization'
                      ,'re_pfa_var5_RPIw','th_pfa_var5','se_pfa_var5','util_pfa_var5'
                      ,'re_pfa_var5_RPIw_ls','th_pfa_var5_ls','se_pfa_var5_ls','util_pfa_var5_ls'
                      ,'res_pred_RPIw_max2','res_pred_RPIw_max_err'
                      ,'therm_pred','therm_err'
                      ,'seis_eq_pred','seis_eq_err','NormAngs','seis_stress_pred','seis_stress_err'
                      ,'util_pred','util_err'
                      ,'co_5_0_20_s_RPIw','co_5_0_625_p_RPIw','co_5_0_5_m_RPIw'
                      ,'co_pfa_sd5_avg_RPIw','co_pfa_sd5_geomean_RPIw'
                      ,'co_5_0_16_s_RPIw_geo','co_5_0_125_p_RPIw_geo','co_5_0_5_m_RPIw_geo'
                      ,'co_pfa_sd5_avg_RPIw_geo','co_pfa_sd5_geomean_RPIw_geo')

comb_extract_5_RPIg <- stack(c(reservoir_RPIg,thermal,seismic,utilization
                              ,re_pfa_var5_RPIg,th_pfa_var5,se_pfa_var5,util_pfa_var5
                              ,re_pfa_var5_RPIg_ls,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls
                              ,res_pred_RPIg_max2,res_pred_RPIg_max_err
                              ,therm_pred,therm_err
                              ,seis_eq_pred,seis_eq_err,NormAngs,seis_stress_pred,seis_stress_err
                              ,util_pred,util_err
                              ,co_5_0_20_s_RPIg,co_5_0_625_p_RPIg,co_5_0_5_m_RPIg
                              ,co_pfa_sd5_avg_RPIg,co_pfa_sd5_geomean_RPIg
                              ,co_5_0_16_s_RPIg_geo,co_5_0_125_p_RPIg_geo,co_5_0_5_m_RPIg_geo
                              ,co_pfa_sd5_avg_RPIg_geo,co_pfa_sd5_geomean_RPIg_geo))

# names for each layer in the raster stack
comb_names_5_RPIg <- c('reservoir_RPIg','thermal','seismic','utilization'
                      ,'re_pfa_var5_RPIg','th_pfa_var5','se_pfa_var5','util_pfa_var5'
                      ,'re_pfa_var5_RPIg_ls','th_pfa_var5_ls','se_pfa_var5_ls','util_pfa_var5_ls'
                      ,'res_pred_RPIg_max2','res_pred_RPIg_max_err'
                      ,'therm_pred','therm_err'
                      ,'seis_eq_pred','seis_eq_err','NormAngs','seis_stress_pred','seis_stress_err'
                      ,'util_pred','util_err'
                      ,'co_5_0_20_s_RPIg','co_5_0_625_p_RPIg','co_5_0_5_m_RPIg'
                      ,'co_pfa_sd5_avg_RPIg','co_pfa_sd5_geomean_RPIg'
                      ,'co_5_0_16_s_RPIg_geo','co_5_0_125_p_RPIg_geo','co_5_0_5_m_RPIg_geo'
                      ,'co_pfa_sd5_avg_RPIg_geo','co_pfa_sd5_geomean_RPIg_geo')

co_pfa_sd5_avg_egs = calc(co_pfa_var5_avg_egs,fun=sqrt)
co_pfa_sd5_geomean_egs = calc(co_pfa_var5_geomean_egs,fun=sqrt)
comb_extract_5_egs <- stack(c(thermal,seismic,utilization
                              ,th_pfa_var5,se_pfa_var5,util_pfa_var5
                              ,th_pfa_var5_ls,se_pfa_var5_ls,util_pfa_var5_ls
                              ,therm_pred,therm_err
                              ,seis_eq_pred,seis_eq_err,NormAngs,seis_stress_pred,seis_stress_err
                              ,util_pred,util_err
                              ,co_5_0_16_s_egs,co_5_0_125_p_egs,co_5_0_5_m_egs
                              ,co_pfa_sd5_avg_egs,co_pfa_sd5_geomean_egs))

comb_names_5_egs <- c('thermal','seismic','utilization'
                      ,'th_pfa_var5','se_pfa_var5','util_pfa_var5'
                      ,'th_pfa_var5_ls','se_pfa_var5_ls','util_pfa_var5_ls'
                      ,'therm_pred','therm_err'
                      ,'seis_eq_pred','seis_eq_err','NormAngs','seis_stress_pred','seis_stress_err'
                      ,'util_pred','util_err'
                      ,'co_5_0_16_s_egs','co_5_0_125_p_egs','co_5_0_5_m_egs'
                      ,'co_pfa_sd5_avg_egs','co_pfa_sd5_geomean_egs')

#Remove the individual files for space:
rm(co_5_0_5_m_rfc, co_5_0_5_m_RPIg, co_5_0_5_m_RPIw)
rm(co_5_0_5_m_rfc_geo, co_5_0_5_m_RPIg_geo, co_5_0_5_m_RPIw_geo)
rm(co_5_0_20_s_rfc, co_5_0_20_s_RPIg, co_5_0_20_s_RPIw)
rm(co_5_0_16_s_rfc_geo, co_5_0_16_s_RPIg_geo, co_5_0_16_s_RPIw_geo)
rm(co_5_0_625_p_rfc, co_5_0_625_p_RPIg, co_5_0_625_p_RPIw)
rm(co_5_0_125_p_rfc_geo, co_5_0_125_p_RPIg_geo, co_5_0_125_p_RPIw_geo)
rm(co_pfa_sd5_avg_rfc, co_pfa_sd5_avg_RPIw, co_pfa_sd5_avg_RPIg)
rm(co_pfa_sd5_avg_rfc_geo, co_pfa_sd5_avg_RPIw_geo, co_pfa_sd5_avg_RPIg_geo)
rm(co_pfa_sd5_geomean_rfc, co_pfa_sd5_geomean_RPIw, co_pfa_sd5_geomean_RPIg)
rm(co_pfa_sd5_geomean_rfc_geo, co_pfa_sd5_geomean_RPIw_geo, co_pfa_sd5_geomean_RPIg_geo)

extracted_5_rfc <- rasterToPoints(comb_extract_5_rfc,spatial=TRUE)
extracted_df_5_rfc <- as.data.frame(extracted_5_rfc)
rm(extracted_5_rfc)

extracted_5_rfc_geo <- rasterToPoints(comb_extract_5_rfc,spatial=TRUE)
extracted_df_5_rfc_geo <- as.data.frame(extracted_5_rfc_geo)
rm(extracted_5_rfc_geo)
rm(comb_extract_5_rfc)

extracted_5_RPIw <- rasterToPoints(comb_extract_5_RPIw,spatial=TRUE)
extracted_df_5_RPIw <- as.data.frame(extracted_5_RPIw)
rm(extracted_5_RPIw)

extracted_5_RPIw_geo <- rasterToPoints(comb_extract_5_RPIw,spatial=TRUE)
extracted_df_5_RPIw_geo <- as.data.frame(extracted_5_RPIw_geo)
rm(extracted_5_RPIw_geo)
rm(comb_extract_5_RPIw)

extracted_5_RPIg <- rasterToPoints(comb_extract_5_RPIg,spatial=TRUE)
extracted_df_5_RPIg <- as.data.frame(extracted_5_RPIg)
rm(extracted_5_RPIg)

extracted_5_RPIg_geo <- rasterToPoints(comb_extract_5_RPIg,spatial=TRUE)
extracted_df_5_RPIg_geo <- as.data.frame(extracted_5_RPIg_geo)
rm(extracted_5_RPIg_geo)
rm(comb_extract_5_RPIg)

extracted_5_egs <- rasterToPoints(comb_extract_5_egs,spatial=TRUE)
extracted_df_5_egs <- as.data.frame(extracted_5_egs)
#rm(extracted_5_egs)
rm(comb_extract_5_egs)

# renaming the columns
names(extracted_df_5_rfc) <- c(comb_names_5_rfc,'x','y')
names(extracted_df_5_RPIw) <- c(comb_names_5_RPIw,'x','y')
names(extracted_df_5_RPIg) <- c(comb_names_5_RPIg,'x','y')

names(extracted_df_5_rfc_geo) <- c(comb_names_5_rfc,'x','y')
names(extracted_df_5_RPIw_geo) <- c(comb_names_5_RPIw,'x','y')
names(extracted_df_5_RPIg_geo) <- c(comb_names_5_RPIg,'x','y')

names(extracted_df_5_egs) <- c(comb_names_5_egs,'x','y')

# dropping points that are NA for geology or all four risk factors
#extracted2_df_5_rfc <- extracted_df_5_rfc[setdiff(seq(1,nrow(extracted_df_5_rfc),1),intersect(which(extracted_df_5_rfc$co_5_0_16_s_rfc_geo %in% NA),which(extracted_df_5_rfc$co_5_0_20_s_rfc %in% NA))),]

#For only 1 at a time, use this:
extracted2_df_5_rfc_a <- extracted_df_5_rfc[setdiff(seq(1,nrow(extracted_df_5_rfc),1),which(extracted_df_5_rfc$co_5_0_20_s_rfc %in% NA)),]
extracted2_df_5_rfc_g <- extracted_df_5_rfc[setdiff(seq(1,nrow(extracted_df_5_rfc),1),which(extracted_df_5_rfc$co_5_0_625_p_rfc %in% NA)),]
extracted2_df_5_rfc_m <- extracted_df_5_rfc[setdiff(seq(1,nrow(extracted_df_5_rfc),1),which(extracted_df_5_rfc$co_5_0_5_m_rfc %in% NA)),]
rm(extracted_df_5_rfc)

extracted2_df_5_RPIw_a <- extracted_df_5_RPIw[setdiff(seq(1,nrow(extracted_df_5_RPIw),1),which(extracted_df_5_RPIw$co_5_0_20_s_RPIw %in% NA)),]
extracted2_df_5_RPIw_g <- extracted_df_5_RPIw[setdiff(seq(1,nrow(extracted_df_5_RPIw),1),which(extracted_df_5_RPIw$co_5_0_625_p_RPIw %in% NA)),]
extracted2_df_5_RPIw_m <- extracted_df_5_RPIw[setdiff(seq(1,nrow(extracted_df_5_RPIw),1),which(extracted_df_5_RPIw$co_5_0_5_m_RPIw %in% NA)),]
rm(extracted_df_5_RPIw)

extracted2_df_5_RPIg_a <- extracted_df_5_RPIg[setdiff(seq(1,nrow(extracted_df_5_RPIg),1),which(extracted_df_5_RPIg$co_5_0_20_s_RPIg %in% NA)),]
extracted2_df_5_RPIg_g <- extracted_df_5_RPIg[setdiff(seq(1,nrow(extracted_df_5_RPIg),1),which(extracted_df_5_RPIg$co_5_0_625_p_RPIg %in% NA)),]
extracted2_df_5_RPIg_m <- extracted_df_5_RPIg[setdiff(seq(1,nrow(extracted_df_5_RPIg),1),which(extracted_df_5_RPIg$co_5_0_5_m_RPIg %in% NA)),]
rm(extracted_df_5_RPIg)

#Geologic
extracted2_df_5_rfc_geo_a <- extracted_df_5_rfc_geo[setdiff(seq(1,nrow(extracted_df_5_rfc_geo),1),which(extracted_df_5_rfc_geo$co_5_0_16_s_rfc_geo %in% NA)),]
extracted2_df_5_rfc_geo_g <- extracted_df_5_rfc_geo[setdiff(seq(1,nrow(extracted_df_5_rfc_geo),1),which(extracted_df_5_rfc_geo$co_5_0_125_p_rfc_geo %in% NA)),]
extracted2_df_5_rfc_geo_m <- extracted_df_5_rfc_geo[setdiff(seq(1,nrow(extracted_df_5_rfc_geo),1),which(extracted_df_5_rfc_geo$co_5_0_5_m_rfc_geo %in% NA)),]
rm(extracted_df_5_rfc_geo)

extracted2_df_5_RPIw_geo_a <- extracted_df_5_RPIw_geo[setdiff(seq(1,nrow(extracted_df_5_RPIw_geo),1),which(extracted_df_5_RPIw_geo$co_5_0_16_s_RPIw_geo %in% NA)),]
extracted2_df_5_RPIw_geo_g <- extracted_df_5_RPIw_geo[setdiff(seq(1,nrow(extracted_df_5_RPIw_geo),1),which(extracted_df_5_RPIw_geo$co_5_0_125_p_RPIw_geo %in% NA)),]
extracted2_df_5_RPIw_geo_m <- extracted_df_5_RPIw_geo[setdiff(seq(1,nrow(extracted_df_5_RPIw_geo),1),which(extracted_df_5_RPIw_geo$co_5_0_5_m_RPIw_geo %in% NA)),]
rm(extracted_df_5_RPIw_geo)

extracted2_df_5_RPIg_geo_a <- extracted_df_5_RPIg_geo[setdiff(seq(1,nrow(extracted_df_5_RPIg_geo),1),which(extracted_df_5_RPIg_geo$co_5_0_16_s_RPIg_geo %in% NA)),]
extracted2_df_5_RPIg_geo_g <- extracted_df_5_RPIg_geo[setdiff(seq(1,nrow(extracted_df_5_RPIg_geo),1),which(extracted_df_5_RPIg_geo$co_5_0_125_p_RPIg_geo %in% NA)),]
extracted2_df_5_RPIg_geo_m <- extracted_df_5_RPIg_geo[setdiff(seq(1,nrow(extracted_df_5_RPIg_geo),1),which(extracted_df_5_RPIg_geo$co_5_0_5_m_RPIg_geo %in% NA)),]
rm(extracted_df_5_RPIg_geo)

#EGS
extracted2_df_5_egs_a <- extracted_df_5_egs[setdiff(seq(1,nrow(extracted_df_5_egs),1),which(extracted_df_5_egs$co_5_0_16_s_egs %in% NA)),]
extracted2_df_5_egs_g <- extracted_df_5_egs[setdiff(seq(1,nrow(extracted_df_5_egs),1),which(extracted_df_5_egs$co_5_0_125_p_egs %in% NA)),]
extracted2_df_5_egs_m <- extracted_df_5_egs[setdiff(seq(1,nrow(extracted_df_5_egs),1),which(extracted_df_5_egs$co_5_0_5_m_egs %in% NA)),]
rm(extracted_df_5_egs)

#### Setting up plotting parameters and points data----
# Setting parameters for boxplots/violin plots for distriutions after MC of the parameters in extracted2
# ignoring 3 because it corresponds to Ithaca, which did not have a reservoir within 10 km
ind_use <- c(1,2,4,5,6,7,8,9,10)

# number of MC replicates
rps <- 10000

# initialzing matrices
#All RFs
distsavg <- matrix(0,rps,length(ind_use))
distsgeomean <- matrix(0,rps,length(ind_use))
distsmin <- matrix(0,rps,length(ind_use))
dist_vars <- matrix(0,4,length(ind_use)) #4 is the number of RFs

#Geology only
distsavg_g <- matrix(0,rps,length(ind_use))
distsgeomean_g <- matrix(0,rps,length(ind_use))
distsmin_g <- matrix(0,rps,length(ind_use))
dist_vars_g <- matrix(0,3,length(ind_use)) #3 is the number of RFs.

#All RFs - may need to change the reservoir and theremal here.
minRe <- rep(0,length(ind_use))
minTh <- rep(0,length(ind_use))
minSe <- rep(0,length(ind_use))
minUt <- rep(0,length(ind_use))

#Geology only - may need to change the reservoir and theremal here.
minRe_g <- rep(0,length(ind_use))
minTh_g <- rep(0,length(ind_use))
minSe_g <-rep(0,length(ind_use))

#May need to change the reservoir and theremal here.
reMean <- rep(0,length(ind_use))
thMean <- rep(0,length(ind_use))
seMean <- rep(0,length(ind_use))
utMean <- rep(0,length(ind_use))

reMed <- rep(0,length(ind_use))
thMed <- rep(0,length(ind_use))
seMed <- rep(0,length(ind_use))
utMed <- rep(0,length(ind_use))

## Starting to collect information for poi2
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
           ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
           ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')
points$names <- names
# calculating distance to each of the key points to the cell centers. Does all extracted2 variables above at once.
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  extracted2_df_5_rfc_a[nm] <- sqrt((extracted2_df_5_rfc_a$x - points$x[i])^2 + (extracted2_df_5_rfc_a$y - points$y[i])^2)
  extracted2_df_5_rfc_g[nm] <- sqrt((extracted2_df_5_rfc_g$x - points$x[i])^2 + (extracted2_df_5_rfc_g$y - points$y[i])^2)
  extracted2_df_5_rfc_m[nm] <- sqrt((extracted2_df_5_rfc_m$x - points$x[i])^2 + (extracted2_df_5_rfc_m$y - points$y[i])^2)
  extracted2_df_5_RPIw_a[nm] <- sqrt((extracted2_df_5_RPIw_a$x - points$x[i])^2 + (extracted2_df_5_RPIw_a$y - points$y[i])^2)
  extracted2_df_5_RPIw_g[nm] <- sqrt((extracted2_df_5_RPIw_g$x - points$x[i])^2 + (extracted2_df_5_RPIw_g$y - points$y[i])^2)
  extracted2_df_5_RPIw_m[nm] <- sqrt((extracted2_df_5_RPIw_m$x - points$x[i])^2 + (extracted2_df_5_RPIw_m$y - points$y[i])^2)
  extracted2_df_5_RPIg_a[nm] <- sqrt((extracted2_df_5_RPIg_a$x - points$x[i])^2 + (extracted2_df_5_RPIg_a$y - points$y[i])^2)
  extracted2_df_5_RPIg_g[nm] <- sqrt((extracted2_df_5_RPIg_g$x - points$x[i])^2 + (extracted2_df_5_RPIg_g$y - points$y[i])^2)
  extracted2_df_5_RPIg_m[nm] <- sqrt((extracted2_df_5_RPIg_m$x - points$x[i])^2 + (extracted2_df_5_RPIg_m$y - points$y[i])^2)
}
rm(nm)

#### RFC Average ####
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for rfc
      ind_max <- which(extracted2_df_5_rfc_a$co_5_0_20_s_rfc[inds] %in% max(extracted2_df_5_rfc_a$co_5_0_20_s_rfc[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_a[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_a$co_5_0_20_s_rfc[inds] %in% max(extracted2_df_5_rfc_a$co_5_0_20_s_rfc[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_rfc_a$co_pfa_sd5_avg_rfc[inds][ind_max] %in% min(extracted2_df_5_rfc_a$co_pfa_sd5_avg_rfc[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_a[nm][inds[ind_max][ind_max2],]))[1]
        
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }

  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }

  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
          
          if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
            #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
            dataParams[i,6] <- 1.5
            k = dataParams[i,6]
            dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
            rm(k)
          }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
        
        if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
          #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          dataParams[i,6] <- 1.5
          k = dataParams[i,6]
          dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          rm(k)
        }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }

  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,4] = 0
    dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,4] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams[which(points$utilization == 0),7] <- NA
dataParams[which(points$utilization == 0),8] <- NA

roots1 <- matrix(NA,10,3) #Percentiles to find roots
converge1 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams[i,7]) == TRUE){
    roots1[i,1] = roots1[i,2] = roots1[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  }
}
rm(params, ps)

#Truncated Weibull
dataParams_trunc <- matrix(NA,10,4*2) #4 is for the number of RFs
right = 5.0
lb = 0
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams_trunc)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams_trunc[i,c(1,2)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams_trunc[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams_trunc[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_trunc[i,5] <- NA
      dataParams_trunc[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams_trunc[i,5] <- NA
          dataParams_trunc[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          
          # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
          #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          #   dataParams_trunc[i,6] <- 1.5
          #   k = dataParams_trunc[i,6]
          #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          #   rm(k)
          # }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        
        # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
        #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
        #   dataParams_trunc[i,6] <- 1.5
        #   k = dataParams_trunc[i,6]
        #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
        #   rm(k)
        # }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,4] = 0
    dataParams_trunc[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_trunc[which(points$utilization == 0),7] <- NA
dataParams_trunc[which(points$utilization == 0),8] <- NA

roots1_trunc <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_trunc[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_trunc[i,7]) == TRUE){
    roots1_trunc[i,1] = roots1_trunc[i,2] = roots1_trunc[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  }
}
rm(params, ps)

# solving for the Beta distribution parameters
dataParams_Beta <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
dpBeta <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Beta as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_Beta[i,c(1,2)] <- NA
    dpBeta[i,1] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_Beta[i,c(1,2)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,1] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }

  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,2] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }

  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_Beta[i,5] <- NA
      dataParams_Beta[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta, start=c(1,1), positive=TRUE)$root
      dpBeta[i,3] <-  multiroot(solveBeta,start=c(100,100), positive=TRUE)$iter

      if (dataParams_Beta[i,6] == 0){
        #Did not converge. Try a fixed high alpha low beta
        dataParams_Beta[i,c(5,6)] = c(100,0.001)
      }

      if (dataParams_Beta[i,5] == 0){
        #Did not converge. Try a fixed low alpha high beta
        dataParams_Beta[i,c(5,6)] = rev(c(100,0.001))
      }

      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,3] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }

  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    ub_mat[i,4] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,4] = 0
    ub_mat[i,4] = 5
    dataParams_Beta[i,c(7,8)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,4] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_Beta[which(points$utilization == 0),7] <- NA
dataParams_Beta[which(points$utilization == 0),8] <- NA

roots1_Beta <- matrix(NA,10,3) #Percentiles to find roots for
converge1_Beta <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams_Beta[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_Beta[i,7]) == TRUE){
    roots1_Beta[i,1] = roots1_Beta[i,2] = roots1_Beta[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
  }
}
rm(params, ps)

#Doubly Truncated Normal Distribution
# solving for the Beta distribution parameters
dataParams_dtn <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
#set the lower bound as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_dtn[i,c(1,2)] <- NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_dtn[i,c(1,2)] <- mom
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_dtn[i,c(3,4)] <- NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_dtn[i,c(3,4)] <- mom
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_dtn[i,5] <- NA
      dataParams_dtn[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_dtn[i,c(5,6)] <- mom
      
      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_dtn[i,c(5,6)] <- mom
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    ub_mat[i,4] = 5
    dataParams_dtn[i,c(3,4)] <- NA
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,4] = 0
    ub_mat[i,4] = 5
    dataParams_dtn[i,c(7,8)] <- mom
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_dtn[which(points$utilization == 0),7] <- NA
dataParams_dtn[which(points$utilization == 0),8] <- NA

roots1_dtn <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_dtn[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_dtn[i,7]) == TRUE){
    roots1_dtn[i,1] = roots1_dtn[i,2] = roots1_dtn[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_dtn[i,1] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_dtn[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_dtn[i,1] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_dtn[i,2] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_dtn[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_dtn[i,2] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_dtn[i,3] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_dtn[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_dtn[i,3] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_dtn[i,1] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    
    ps <- 0.5
    roots1_dtn[i,2] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
    
    ps <- 0.95
    roots1_dtn[i,3] <- uniroot(solveQuant_dtn,interval=c(min(lb_mat[i,]),5))$root
  }
}
rm(params, ps)

#Truncated Folded Normal Distribution
# solving for the Beta distribution parameters
dataParams_ftn <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
#set the lower bound as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_ftn[i,c(1,2)] <- NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_ftn[i,c(1,2)] <- mom
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_ftn[i,c(3,4)] <- NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_ftn[i,c(3,4)] <- mom
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_ftn[i,5] <- NA
      dataParams_ftn[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_ftn[i,c(5,6)] <- mom
      
      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_ftn[i,c(5,6)] <- mom
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    ub_mat[i,4] = 5
    dataParams_ftn[i,c(3,4)] <- NA
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,4] = 0
    ub_mat[i,4] = 5
    dataParams_ftn[i,c(7,8)] <- mom
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_ftn[which(points$utilization == 0),7] <- NA
dataParams_ftn[which(points$utilization == 0),8] <- NA

roots1_ftn <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_ftn[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_ftn[i,7]) == TRUE){
    roots1_ftn[i,1] = roots1_ftn[i,2] = roots1_ftn[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_ftn[i,1] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_ftn[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_ftn[i,1] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_ftn[i,2] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_ftn[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_ftn[i,2] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_ftn[i,3] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots1_ftn[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_ftn[i,3] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_ftn[i,1] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    
    ps <- 0.5
    roots1_ftn[i,2] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
    
    ps <- 0.95
    roots1_ftn[i,3] <- uniroot(solveQuant_ftn,interval=c(min(lb_mat[i,]),5))$root
  }
}
rm(params, ps)

# Seismic Doubly Truncated Normal when lb = 2.5, all else Weibull - Used for paper
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the upper bound as 5
ub = 5
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower and upper bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(5, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          #Use the doubly truncated normal
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- mom
          dp2[i,3] <-  1
          
          lb=0
        }
      }else{
        #Use the doubly truncated normal
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- mom
        dp2[i,3] <-  1
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,4] = 0
    dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,4] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams[which(points$utilization == 0),7] <- NA
dataParams[which(points$utilization == 0),8] <- NA

roots1 <- matrix(NA,10,3) #Percentiles to find roots
converge1 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams[i,7]) == TRUE){
    roots1[i,1] = roots1[i,2] = roots1[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge1[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge1[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge1[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
  }
} #Note that a warning message will pop up when Utilization is 0 (not fit by Weibull)
rm(params, ps)

# loop for Monte Carlo
right = Inf #Set to Inf for non-truncated Weibull plots
seis = TRUE #Set to TRUE for doubly truncated normal condition
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  reMed[i] <- median(mat_mc[,1])
  thMed[i] <-  median(mat_mc[,2])
  seMed[i] <-  median(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMed[i] <-  median(mat_mc[,5])
  
  setwd(wd_image)
  par(xpd=T)
  if (right != Inf){
    name = 'scdist_trunc'
  }else if (seis == TRUE){
    name = 'scdist_s'
  }else{
    name = 'scdist'
  }
  # png(paste(name,i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # tiff(paste(name,i,'.tiff',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=600
  # )
  setEPS()
  postscript(paste(name,i,'.eps',sep=''), height = 6, width = 6)
  
  par(mfrow=c(2,2)
      ,oma=c(0.5,0.5,0.5,0.5)+0.1
      ,mar=c(5,5,3,0)+0.1)
  
  ymax = ceil(max(hist(mat_mc[,1], plot = F, breaks = seq(0,5,0.1))$density, 
                  ifelse(right != Inf, 
                         max(dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])/pweibull(5,shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3]), na.rm=TRUE), 
                         max(dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3]), na.rm=TRUE)), na.rm=TRUE))
  hist(mat_mc[,1]
       ,freq=F
       ,xlim=c(0,5)
       ,ylim=c(0,ymax)
       ,main="Reservoir"
       ,xlab="SFF"
       ,breaks = seq(0,5,0.1)
  )
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])/pweibull(5,shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }else{
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }
  points(points$reservoir[ind_use[i]],0
         ,pch=19
         ,ylim=c(0,ymax)
  )
  
  ymax = ceil(max(hist(mat_mc[,2], plot = F, breaks = seq(0,5,0.1))$density, 
                  ifelse(right != Inf, 
                         max(dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])/pweibull(5,shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1]), na.rm=TRUE), 
                         max(dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1]), na.rm = TRUE)), na.rm=TRUE))
  hist(mat_mc[,2]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Thermal"
       ,xlab="SFF"
       ,ylim=c(0,ymax), breaks = seq(0,5,0.1)
  )
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])/pweibull(5,shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }else{
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }
  points(points$thermal[ind_use[i]],0
         ,pch=19
         ,ylim=c(0,ymax)
  )
  
  ymax = ceil(max(hist(0.5*(mat_mc[,3]+mat_mc[,4]), plot = F, breaks = seq(0,5,0.1))$density), digits = 1)
  hist(0.5*(mat_mc[,3]+mat_mc[,4])
       ,freq=F
       ,xlim=c(0,5)
       ,main="Seismic"
       ,xlab="SFF"
       ,ylim=c(0,ymax), breaks = seq(0,5,0.1)
  )
  if (right != Inf){
    #Use truncated Weibull pdf
    if (lb_mat[ind_use[i],3] != 0){
      lines(seq(lb_mat[ind_use[i],3],right,0.001)
            ,dweibull(seq(0,(right - lb_mat[ind_use[i],3]),0.001),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull((right - lb_mat[ind_use[i],3]), shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='black'
            ,lwd=2
            ,ylim=c(0,ymax)
      )
    }
    else{
      lines(seq(0,5,0.001)
            ,dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull(5,shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='black'
            ,lwd=2
            ,ylim=c(0,ymax)
      )
    }
  }else{
    if (seis != TRUE){
      if (lb_mat[ind_use[i],3] != 0){
        lines(seq(lb_mat[ind_use[i],3],5+lb_mat[ind_use[i],3],0.001)
              ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='black'
              ,lwd=2
              ,ylim=c(0,ymax)
        )
      }
      else{
        lines(seq(0,5,0.001)
              ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='black'
              ,lwd=2
              ,ylim=c(0,ymax)
        )
      }
    }else{
      if (lb_mat[ind_use[i],3] != 0){
        lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.001)
              ,dtnorm(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.001), dataParams[ind_use[i],5], sqrt(dataParams[ind_use[i],6]), lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
              ,col='black'
              ,lwd=2
              ,ylim=c(0,ymax)
        )
      }
      else{
        lines(seq(0,5,0.001)
              ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='black'
              ,lwd=2
              ,ylim=c(0,ymax)
        )
      }
    }
  }
  points(points$seismic[ind_use[i]],0
         ,pch=19
         ,ylim=c(0,ymax)
  )
  
  ymax = ceil(max(hist(mat_mc[,5], plot = F, breaks = seq(0,5,0.1))$density, 
                  ifelse(right != Inf, 
                         max(dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7])/pweibull(5,shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7]), na.rm=TRUE), 
                         max(dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7]), na.rm = TRUE)), na.rm=TRUE))
  hist(mat_mc[,5]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Utilization"
       ,xlab="SFF"
       ,ylim=c(0,ymax), breaks = seq(0,5,0.1)
  )
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7])/pweibull(5,shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }else{
    lines(seq(0,5,0.001)
          ,dweibull(seq(0,5,0.001),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
          ,col='black'
          ,lwd=2
          ,ylim=c(0,ymax)
    )
  }
  points(points$utilization[ind_use[i]],0
         ,pch=19
         ,ylim=c(0,ymax)
  )
  
  par(xpd=T)
  dev.off()
  
  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist_Beta',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],3],dataParams_Beta[ind_use[i],4], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],1],dataParams_Beta[ind_use[i],2], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
  #       ,dBeta_ab(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams_Beta[ind_use[i],5], dataParams_Beta[ind_use[i],6], lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],7],dataParams_Beta[ind_use[i],8], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
  
  # par(xpd=T)
  # png(paste('scdist_dtn',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dtnorm(seq(0,5,0.01),dataParams_dtn[ind_use[i],3],sqrt(dataParams_dtn[ind_use[i],4]), 0 , 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dtnorm(seq(0,5,0.01),dataParams_dtn[ind_use[i],1],sqrt(dataParams_dtn[ind_use[i],2]), 0 , 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
  #       ,dtnorm(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams_dtn[ind_use[i],5], sqrt(dataParams_dtn[ind_use[i],6]), lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dtnorm(seq(0,5,0.01),dataParams_dtn[ind_use[i],7],sqrt(dataParams_dtn[ind_use[i],8]), 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
  # 
  # par(xpd=T)
  # png(paste('scdist_ftn',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dfoldnorm(seq(0,5,0.01), dataParams_ftn[ind_use[i],3], sqrt(dataParams_ftn[ind_use[i],4]))/(pfoldnorm(ub_mat[ind_use[i],2], dataParams_ftn[ind_use[i],3], sqrt(dataParams_ftn[ind_use[i],4])) - pfoldnorm(lb_mat[ind_use[i],2], dataParams_ftn[ind_use[i],3], sqrt(dataParams_ftn[ind_use[i],4])))
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dfoldnorm(seq(0,5,0.01),dataParams_ftn[ind_use[i],1],sqrt(dataParams_ftn[ind_use[i],2]))/(pfoldnorm(ub_mat[ind_use[i],1], dataParams_ftn[ind_use[i],1], sqrt(dataParams_ftn[ind_use[i],2])) - pfoldnorm(lb_mat[ind_use[i],1], dataParams_ftn[ind_use[i],1], sqrt(dataParams_ftn[ind_use[i],2])))
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
  #       ,dfoldnorm(seq(0, (ub_mat[ind_use[i],3] - lb_mat[ind_use[i],3]),0.01), (dataParams_ftn[ind_use[i],5] - lb_mat[ind_use[i],3]), sqrt(dataParams_ftn[ind_use[i],6]))/(pfoldnorm((ub_mat[ind_use[i],3] - lb_mat[ind_use[i],3]), (dataParams_ftn[ind_use[i],5] - lb_mat[ind_use[i],3]), sqrt(dataParams_ftn[ind_use[i],6])) - pfoldnorm(0, (dataParams_ftn[ind_use[i],5] - lb_mat[ind_use[i],3]), sqrt(dataParams_ftn[ind_use[i],6])))
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dfoldnorm(seq(0,5,0.01),dataParams_ftn[ind_use[i],7],sqrt(dataParams_ftn[ind_use[i],8]))/(pfoldnorm(ub_mat[ind_use[i],4], dataParams_ftn[ind_use[i],7], sqrt(dataParams_ftn[ind_use[i],8])) - pfoldnorm(lb_mat[ind_use[i],4], dataParams_ftn[ind_use[i],7], sqrt(dataParams_ftn[ind_use[i],8])))
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg)){
  
  boxplot(distsavg[,i]
          ,at=i
          ,col=cols[i]
          ,add=TRUE
          ,pars=list(outcex=0.5,outpch=19)
  )
}
par(xpd=TRUE)
text(seq(1,9,1)
     ,rep(-0.5,9)
     ,names[ind_use]
     ,adj=0
     ,srt=270)
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg)

setwd(wd_image)
png('violin_5_rfc_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,names[ind_use]
     ,adj=0
     ,srt=270)

text(-1
     ,2.5
     ,"Combined Metric: Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
# png('parallel_axis_5_rfc_a.png'
#     ,height=5
#     ,width=7
#     ,units='in'
#     ,res=300
# )
cols = grey.colors(10)
# tiff('parallel_axis_5_rfc_a.tiff'
#     ,height=5
#     ,width=7
#     ,units='in'
#     ,res=600
# )
setEPS()
postscript('parallel_axis_5_rfc_a.eps', height=5, width=7)

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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))
namesNA = NULL

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    namesNA <- c(namesNA,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2,3)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic','Utilization')
     ,adj=0.5)
text(-0.35,5
     ,'Better'
     ,adj=1
     ,col='black'
     ,font=2)
text(-0.35,0
     ,'Worse'
     ,adj=1
     ,col='black'
     ,font=2)
legend(x=3.2,y=5
       ,legend=names[ind_use]
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RFC Geometric Mean ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
points$names = names
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for rfc
      ind_max <- which(extracted2_df_5_rfc_g$co_5_0_625_p_rfc[inds] %in% max(extracted2_df_5_rfc_g$co_5_0_625_p_rfc[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_g[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_g$co_5_0_625_p_rfc[inds] %in% max(extracted2_df_5_rfc_g$co_5_0_625_p_rfc[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_rfc_g$co_pfa_sd5_geomean_rfc[inds][ind_max] %in% min(extracted2_df_5_rfc_g$co_pfa_sd5_geomean_rfc[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
          
          if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
            #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
            dataParams[i,6] <- 1.5
            k = dataParams[i,6]
            dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
            rm(k)
          }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
        
        if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
          #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          dataParams[i,6] <- 1.5
          k = dataParams[i,6]
          dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          rm(k)
        }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,4] = 0
    dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,4] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams[which(points$utilization == 0),7] <- NA
dataParams[which(points$utilization == 0),8] <- NA

roots1 <- matrix(NA,10,3) #Percentiles to find roots
converge1 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams[i,7]) == TRUE){
    roots1[i,1] = roots1[i,2] = roots1[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  }
}
rm(params, ps)

#Truncated Weibull
dataParams_trunc <- matrix(NA,10,4*2) #4 is for the number of RFs
right = 5.0
lb = 0
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams_trunc)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams_trunc[i,c(1,2)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams_trunc[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams_trunc[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_trunc[i,5] <- NA
      dataParams_trunc[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams_trunc[i,5] <- NA
          dataParams_trunc[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          
          # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
          #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          #   dataParams_trunc[i,6] <- 1.5
          #   k = dataParams_trunc[i,6]
          #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          #   rm(k)
          # }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        
        # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
        #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
        #   dataParams_trunc[i,6] <- 1.5
        #   k = dataParams_trunc[i,6]
        #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
        #   rm(k)
        # }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,4] = 0
    dataParams_trunc[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_trunc[which(points$utilization == 0),7] <- NA
dataParams_trunc[which(points$utilization == 0),8] <- NA

roots1_trunc <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_trunc[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_trunc[i,7]) == TRUE){
    roots1_trunc[i,1] = roots1_trunc[i,2] = roots1_trunc[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_trunc[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_trunc[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1_trunc[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1_trunc[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
    converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  }
}
rm(params, ps)

# solving for the Beta distribution parameters
dataParams_Beta <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
dpBeta <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Beta as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_Beta[i,c(1,2)] <- NA
    dpBeta[i,1] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_Beta[i,c(1,2)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,1] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,2] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_Beta[i,5] <- NA
      dataParams_Beta[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta, start=c(1,1), positive=TRUE)$root
      dpBeta[i,3] <-  multiroot(solveBeta,start=c(100,100), positive=TRUE)$iter
      
      if (dataParams_Beta[i,6] == 0){
        #Did not converge. Try a fixed high alpha low beta
        dataParams_Beta[i,c(5,6)] = c(100,0.001)
      }
      
      if (dataParams_Beta[i,5] == 0){
        #Did not converge. Try a fixed low alpha high beta
        dataParams_Beta[i,c(5,6)] = rev(c(100,0.001))
      }
      
      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,3] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  #check if utilization is 5 with 0 uncertainty
  if ((points$utilization[i] == 5 & points$util_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,4] = 5
    ub_mat[i,4] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,4] = 0
    ub_mat[i,4] = 5
    dataParams_Beta[i,c(7,8)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,4] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
}
rm(mom)

# setting NAs for utilization 0 because this will have no uncertainty.
dataParams_Beta[which(points$utilization == 0),7] <- NA
dataParams_Beta[which(points$utilization == 0),8] <- NA

roots1_Beta <- matrix(NA,10,3) #Percentiles to find roots for
converge1_Beta <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams_Beta[i,]
  
  #Determine how to search for the root:
  #First check for utilization 0. All roots for this are 0.
  if (is.na(dataParams_Beta[i,7]) == TRUE){
    roots1_Beta[i,1] = roots1_Beta[i,2] = roots1_Beta[i,3] = 0.0
  }else if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots1_Beta[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.5
    roots1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.95
    roots1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge1_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
  }
}
rm(params, ps)


# loop for Monte Carlo
right = Inf #Set to Inf for non-truncated Weibull plots
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  reMed[i] <- median(mat_mc[,1])
  thMed[i] <-  median(mat_mc[,2])
  seMed[i] <-  median(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMed[i] <-  median(mat_mc[,5])
  
  # setwd(wd_image)
  # par(xpd=T)
  # if (right != Inf){
  #   name = 'scdist_trunc'
  # }else{
  #   name = 'scdist'
  # }
  # png(paste(name,i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # if (right != Inf){
  #   #Truncated Weibull pdf
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])/pweibull(5,shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }else{
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # if (right != Inf){
  #   #Truncated Weibull pdf
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])/pweibull(5,shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }else{
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # if (right != Inf){
  #   #Use truncated Weibull pdf
  #   if (lb_mat[ind_use[i],3] != 0){
  #     lines(seq(lb_mat[ind_use[i],3],right,0.01)
  #           ,dweibull(seq(0,(right - lb_mat[ind_use[i],3]),0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull((right - lb_mat[ind_use[i],3]), shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
  #           ,col='royalblue'
  #           ,lwd=2
  #     )
  #   }
  #   else{
  #     lines(seq(0,5,0.01)
  #           ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull(5,shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
  #           ,col='royalblue'
  #           ,lwd=2
  #     )
  #   }
  # }else{
  #   if (lb_mat[ind_use[i],3] != 0){
  #     lines(seq(lb_mat[ind_use[i],3],5+lb_mat[ind_use[i],3],0.01)
  #           ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #           ,col='royalblue'
  #           ,lwd=2
  #     )
  #   }
  #   else{
  #     lines(seq(0,5,0.01)
  #           ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #           ,col='royalblue'
  #           ,lwd=2
  #     )
  #   }
  # }
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # if (right != Inf){
  #   #Truncated Weibull pdf
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7])/pweibull(5,shape=dataParams_trunc[ind_use[i],8],scale=dataParams_trunc[ind_use[i],7])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }else{
  #   lines(seq(0,5,0.01)
  #         ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #         ,col='royalblue'
  #         ,lwd=2
  #   )
  # }
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
  # 
  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist_Beta',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],3],dataParams_Beta[ind_use[i],4], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],1],dataParams_Beta[ind_use[i],2], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
  #       ,dBeta_ab(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams_Beta[ind_use[i],5], dataParams_Beta[ind_use[i],6], lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],7],dataParams_Beta[ind_use[i],8], 0, 5)
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geometric Mean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean)){
  
  boxplot(distsgeomean[,i]
          ,at=i
          ,col=cols[i]
          ,add=TRUE
          ,pars=list(outcex=0.5,outpch=19)
  )
}
par(xpd=TRUE)
text(seq(1,9,1)
     ,rep(-0.5,9)
     ,names[ind_use]
     ,adj=0
     ,srt=270)
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean)
png('violin_5_rfc_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,names[ind_use]
     ,adj=0
     ,srt=270)

text(-1
     ,2.5
     ,"Combined Metric: Geometric Mean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_rfc_g.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))
namesNA <- NULL
# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    namesNA <- c(namesNA,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,legend=namesNA
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RFC Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for minimum for rfc
      ind_max <- which(extracted2_df_5_rfc_m$co_5_0_5_m_rfc[inds] %in% max(extracted2_df_5_rfc_m$co_5_0_5_m_rfc[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_m[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_m$co_5_0_5_m_rfc[inds] %in% max(extracted2_df_5_rfc_m$co_5_0_5_m_rfc[inds]))
      #Take the value with the lowest uncertainty. Use geomean for now.
      ind_max2 <- which(extracted2_df_5_rfc_m$co_pfa_sd5_geomean_rfc[inds][ind_max] %in% min(extracted2_df_5_rfc_m$co_pfa_sd5_geomean_rfc[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin)){
  
  boxplot(distsmin[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin)
png('violin_5_rfc_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_rfc_m.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIw Average ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_a$co_5_0_20_s_RPIw[inds] %in% max(extracted2_df_5_RPIw_a$co_5_0_20_s_RPIw[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_a[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_a$co_5_0_20_s_RPIw[inds] %in% max(extracted2_df_5_RPIw_a$co_5_0_20_s_RPIw[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_a$co_pfa_sd5_avg_RPIw[inds][ind_max] %in% min(extracted2_df_5_RPIw_a$co_pfa_sd5_avg_RPIw[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg)){
  
  boxplot(distsavg[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg)

setwd(wd_image)
png('violin_5_RPIw_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_a.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIw Geometric Mean ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_g$co_5_0_625_p_RPIw[inds] %in% max(extracted2_df_5_RPIw_g$co_5_0_625_p_RPIw[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_g[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_g$co_5_0_625_p_RPIw[inds] %in% max(extracted2_df_5_RPIw_g$co_5_0_625_p_RPIw[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_g$co_pfa_sd5_geomean_RPIw[inds][ind_max] %in% min(extracted2_df_5_RPIw_g$co_pfa_sd5_geomean_RPIw[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geometric Mean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean)){
  
  boxplot(distsgeomean[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean)
png('violin_5_RPIw_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: Geometric Mean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_g.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIw Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_m$co_5_0_5_m_RPIw[inds] %in% max(extracted2_df_5_RPIw_m$co_5_0_5_m_RPIw[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_m[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_m$co_5_0_5_m_RPIw[inds] %in% max(extracted2_df_5_RPIw_m$co_5_0_5_m_RPIw[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_m$co_pfa_sd5_avg_RPIw[inds][ind_max] %in% min(extracted2_df_5_RPIw_m$co_pfa_sd5_avg_RPIw[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin)){
  
  boxplot(distsmin[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin)
png('violin_5_RPIw_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_m.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Average ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_a$co_5_0_20_s_RPIg[inds] %in% max(extracted2_df_5_RPIg_a$co_5_0_20_s_RPIg[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_a[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_a$co_5_0_20_s_RPIg[inds] %in% max(extracted2_df_5_RPIg_a$co_5_0_20_s_RPIg[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_a$co_pfa_sd5_avg_RPIg[inds][ind_max] %in% min(extracted2_df_5_RPIg_a$co_pfa_sd5_avg_RPIg[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg)){
  
  boxplot(distsavg[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg)

setwd(wd_image)
png('violin_5_RPIg_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_a.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Geometric Mean ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_g$co_5_0_625_p_RPIg[inds] %in% max(extracted2_df_5_RPIg_g$co_5_0_625_p_RPIg[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_g[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_g$co_5_0_625_p_RPIg[inds] %in% max(extracted2_df_5_RPIg_g$co_5_0_625_p_RPIg[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_g$co_pfa_sd5_geomean_RPIg[inds][ind_max] %in% min(extracted2_df_5_RPIg_g$co_pfa_sd5_geomean_RPIg[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geometric Mean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean)){
  
  boxplot(distsgeomean[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean)
png('violin_5_RPIg_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: Geometric Mean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_g.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_m$co_5_0_5_m_RPIg[inds] %in% max(extracted2_df_5_RPIg_m$co_5_0_5_m_RPIg[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_m[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_m$co_5_0_5_m_RPIg[inds] %in% max(extracted2_df_5_RPIg_m$co_5_0_5_m_RPIg[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_m$co_pfa_sd5_geomean_RPIg[inds][ind_max] %in% min(extracted2_df_5_RPIg_m$co_pfa_sd5_geomean_RPIg[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$utilization[i],points$util_pfa_var5[i])
  dataParams[i,c(7,8)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,4] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

dataParams[which(dataParams[,7]<0),7] <- NA
dataParams[which(dataParams[,7]%in% NA),8] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  distsavg[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4]+mat_mc[,5])/4
  distsgeomean[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4])*mat_mc[,5])^(1/4)
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  #distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  #for(j in 1:rps){
  #  if(distsmin_g[j,i] == mat_mc[j,1]){
  #    minRe_g[i] <- minRe_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == mat_mc[j,2]){
  #    minTh_g[i] <- minTh_g[i] + 1
  #  }
  #  if(distsmin_g[j,i] == 0.5*(mat_mc[j,3]+mat_mc[j,4])){
  #    minSe_g[i] <- minSe_g[i] + 1
  #  }
  #}
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  utMean[i] <-  mean(mat_mc[,5])
  
  # print(mean(mat_mc[,1]))
  # print(mean(mat_mc[,2]))
  # print(mean(mat_mc[,3]))
  # print(mean(mat_mc[,4]))
  # print(mean(mat_mc[,5]))
  
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
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin)){
  
  boxplot(distsmin[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin)
png('violin_5_RPIg_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_m.png'
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i]
             ,points$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

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
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)


#### Geology Only ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
           ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
           ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')
points$names <- names

# calculating distance to each of the key points to the cell centers. Does all extracted2 variables above at once.
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  extracted2_df_5_rfc_geo_a[nm] <- sqrt((extracted2_df_5_rfc_geo_a$x - points$x[i])^2 + (extracted2_df_5_rfc_geo_a$y - points$y[i])^2)
  extracted2_df_5_rfc_geo_g[nm] <- sqrt((extracted2_df_5_rfc_geo_g$x - points$x[i])^2 + (extracted2_df_5_rfc_geo_g$y - points$y[i])^2)
  extracted2_df_5_rfc_geo_m[nm] <- sqrt((extracted2_df_5_rfc_geo_m$x - points$x[i])^2 + (extracted2_df_5_rfc_geo_m$y - points$y[i])^2)
  extracted2_df_5_RPIw_geo_a[nm] <- sqrt((extracted2_df_5_RPIw_geo_a$x - points$x[i])^2 + (extracted2_df_5_RPIw_geo_a$y - points$y[i])^2)
  extracted2_df_5_RPIw_geo_g[nm] <- sqrt((extracted2_df_5_RPIw_geo_g$x - points$x[i])^2 + (extracted2_df_5_RPIw_geo_g$y - points$y[i])^2)
  extracted2_df_5_RPIw_geo_m[nm] <- sqrt((extracted2_df_5_RPIw_geo_m$x - points$x[i])^2 + (extracted2_df_5_RPIw_geo_m$y - points$y[i])^2)
  extracted2_df_5_RPIg_geo_a[nm] <- sqrt((extracted2_df_5_RPIg_geo_a$x - points$x[i])^2 + (extracted2_df_5_RPIg_geo_a$y - points$y[i])^2)
  extracted2_df_5_RPIg_geo_g[nm] <- sqrt((extracted2_df_5_RPIg_geo_g$x - points$x[i])^2 + (extracted2_df_5_RPIg_geo_g$y - points$y[i])^2)
  extracted2_df_5_RPIg_geo_m[nm] <- sqrt((extracted2_df_5_RPIg_geo_m$x - points$x[i])^2 + (extracted2_df_5_RPIg_geo_m$y - points$y[i])^2)
}
rm(nm)

#### RFC Geologic Average ####
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_geo_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for rfc_geo
      ind_max <- which(extracted2_df_5_rfc_geo_a$co_5_0_16_s_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_a$co_5_0_16_s_rfc_geo[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_a[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_geo_a$co_5_0_16_s_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_a$co_5_0_16_s_rfc_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_rfc_geo_a$co_pfa_sd5_avg_rfc_geo[inds][ind_max] %in% min(extracted2_df_5_rfc_geo_a$co_pfa_sd5_avg_rfc_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_geo_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_geo_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,3*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,3) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
          
          if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
            #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
            dataParams[i,6] <- 1.5
            k = dataParams[i,6]
            dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
            rm(k)
          }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
        
        if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
          #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          dataParams[i,6] <- 1.5
          k = dataParams[i,6]
          dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          rm(k)
        }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

roots2 <- matrix(NA,10,3) #Percentiles to find roots
converge2 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
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
}
rm(params, ps)

# Seismic Doubly Truncated Normal when lb = 2.5, all else Weibull - Used for paper
dataParams <- matrix(NA,10,3*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,3) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- mom
          dp2[i,3] <-  1
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- mom
        dp2[i,3] <-  1
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

roots2 <- matrix(NA,10,3) #Percentiles to find roots
converge2 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge2[i,1] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
    
    ps <- 0.5
    roots2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge2[i,2] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
    
    ps <- 0.95
    roots2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$root
    converge2[i,3] <- uniroot(solveQuant_SpecialSeis,interval=c(0,5))$iter
  }
}
rm(params, ps)

#Truncated Weibull
dataParams_trunc <- matrix(NA,10,3*2) #4 is for the number of RFs
right = 5.0
lb = 0
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams_trunc)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams_trunc[i,c(1,2)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams_trunc[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams_trunc[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_trunc[i,5] <- NA
      dataParams_trunc[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams_trunc[i,5] <- NA
          dataParams_trunc[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          
          # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
          #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          #   dataParams_trunc[i,6] <- 1.5
          #   k = dataParams_trunc[i,6]
          #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          #   rm(k)
          # }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        
        # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
        #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
        #   dataParams_trunc[i,6] <- 1.5
        #   k = dataParams_trunc[i,6]
        #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
        #   rm(k)
        # }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
}
rm(mom)

roots2_trunc <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_trunc[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    
    ps <- 0.5
    roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    
    ps <- 0.95
    roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  }
}
rm(params, ps)

# solving for the Beta distribution parameters
dataParams_Beta <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
dpBeta <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Beta as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_Beta[i,c(1,2)] <- NA
    dpBeta[i,1] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_Beta[i,c(1,2)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,1] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,2] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_Beta[i,5] <- NA
      dataParams_Beta[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta, start=c(1,1), positive=TRUE)$root
      dpBeta[i,3] <-  multiroot(solveBeta,start=c(100,100), positive=TRUE)$iter
      
      if (dataParams_Beta[i,6] == 0){
        #Did not converge. Try a fixed high alpha low beta
        dataParams_Beta[i,c(5,6)] = c(100,0.001)
      }
      
      if (dataParams_Beta[i,5] == 0){
        #Did not converge. Try a fixed low alpha high beta
        dataParams_Beta[i,c(5,6)] = rev(c(100,0.001))
      }
      
      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,3] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
}
rm(mom)

roots2_Beta <- matrix(NA,10,3) #Percentiles to find roots for
converge2_Beta <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams_Beta[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.5
    roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.95
    roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
  }
}
rm(params, ps)

# loop for Monte Carlo
right = Inf #Set to Inf for non-truncated Weibull plots
seis = TRUE
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  
  reMed[i] <- median(mat_mc[,1])
  thMed[i] <-  median(mat_mc[,2])
  seMed[i] <-  median(0.5*(mat_mc[,3]+mat_mc[,4]))
  
  setwd(wd_image)
  par(xpd=T)
  if (right != Inf){
    name = 'scdist_trunc_g'
  }else if (seis == TRUE){
    name = 'scdist_s_g'
  }else{
    name = 'scdist_g'
  }
  png(paste(name,i,'.png',sep='')
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
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])/pweibull(5,shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])
          ,col='royalblue'
          ,lwd=2
    )
  }else{
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
          ,col='royalblue'
          ,lwd=2
    )
  }
  points(points$reservoir[ind_use[i]],0
         ,pch=19
  )
  
  hist(mat_mc[,2]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Thermal"
       ,xlab="SFF"
  )
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])/pweibull(5,shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])
          ,col='royalblue'
          ,lwd=2
    )
  }else{
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
          ,col='royalblue'
          ,lwd=2
    )
  }
  points(points$thermal[ind_use[i]],0
         ,pch=19
  )
  
  hist(0.5*(mat_mc[,3]+mat_mc[,4])
       ,freq=F
       ,xlim=c(0,5)
       ,main="Seismic"
       ,xlab="SFF"
  )
  if (right != Inf){
    #Use truncated Weibull pdf
    if (lb_mat[ind_use[i],3] != 0){
      lines(seq(lb_mat[ind_use[i],3],right,0.01)
            ,dweibull(seq(0,(right - lb_mat[ind_use[i],3]),0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull((right - lb_mat[ind_use[i],3]), shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
    else{
      lines(seq(0,5,0.01)
            ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull(5,shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
  }else{
    if (seis != TRUE){
      if (lb_mat[ind_use[i],3] != 0){
        lines(seq(lb_mat[ind_use[i],3],5+lb_mat[ind_use[i],3],0.01)
              ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='royalblue'
              ,lwd=2
        )
      }
      else{
        lines(seq(0,5,0.01)
              ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='royalblue'
              ,lwd=2
        )
      }
    }else{
      if (lb_mat[ind_use[i],3] != 0){
        lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
              ,dtnorm(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams[ind_use[i],5], sqrt(dataParams[ind_use[i],6]), lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
              ,col='royalblue'
              ,lwd=2
        )
      }
      else{
        lines(seq(0,5,0.01)
              ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
              ,col='royalblue'
              ,lwd=2
        )
      }
    }
    
  }
  points(points$seismic[ind_use[i]],0
         ,pch=19
  )
  
  par(xpd=T)
  dev.off()
  
  setwd(wd_image)
  par(xpd=T)
  png(paste('scdist_Beta_g',i,'.png',sep='')
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
        ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],3],dataParams_Beta[ind_use[i],4], 0, 5)
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
        ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],1],dataParams_Beta[ind_use[i],2], 0, 5)
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
  lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
        ,dBeta_ab(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams_Beta[ind_use[i],5], dataParams_Beta[ind_use[i],6], lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
        ,col='royalblue'
        ,lwd=2
  )
  points(points$seismic[ind_use[i]],0
         ,pch=19
  )
  
  par(xpd=T)
  dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg_g)){
  
  boxplot(distsavg_g[,i]
          ,at=i
          ,col=cols[i]
          ,add=TRUE
          ,pars=list(outcex=0.5,outpch=19)
  )
}
par(xpd=TRUE)
text(seq(1,9,1)
     ,rep(-0.5,9)
     ,names[ind_use]
     ,adj=0
     ,srt=270)
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg_g)

setwd(wd_image)
png('violin_5_rfc_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,names[ind_use]
     ,adj=0
     ,srt=270)

text(-1
     ,2.5
     ,"Combined Metric: Geologic Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_rfc_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))
namesNA <- NULL

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    namesNA <- c(namesNA,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=namesNA
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RFC Geologic Geometric Mean ####
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_geo_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for rfc
      ind_max <- which(extracted2_df_5_rfc_geo_g$co_5_0_125_p_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_g$co_5_0_125_p_rfc_geo[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_g[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_geo_g$co_5_0_125_p_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_g$co_5_0_125_p_rfc_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_rfc_geo_g$co_pfa_sd5_geomean_rfc_geo[inds][ind_max] %in% min(extracted2_df_5_rfc_geo_g$co_pfa_sd5_geomean_rfc_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_geo_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_geo_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,3*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,3) #Number of iterations
#set the lower bound of the Weibull as 0
lb = 0
#Set the right truncation point to Inf
right = Inf
#Record the value of the lower bound in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams[i,c(1,2)] <- NA
    dp2[i,1] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
    dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams[i,c(3,4)] <- NA
    dp2[i,2] <-  NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,2] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams[i,5] <- NA
      dataParams[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams[i,5] <- NA
          dataParams[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
          
          if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
            #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
            dataParams[i,6] <- 1.5
            k = dataParams[i,6]
            dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
            rm(k)
          }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
        
        if (dataParams[i,5] == 0 | dataParams[i,6] == 0){
          #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          dataParams[i,6] <- 1.5
          k = dataParams[i,6]
          dataParams[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          rm(k)
        }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
    dp2[i,3] <-  multiroot(solveWeibull,start=c(1,2), positive=TRUE)$iter
  }
}
rm(mom)

roots2 <- matrix(NA,10,3) #Percentiles to find roots
converge2 <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    converge2[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
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
}
rm(params, ps)

#Truncated Weibull
dataParams_trunc <- matrix(NA,10,3*2) #4 is for the number of RFs
right = 5.0
lb = 0
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams_trunc)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    dataParams_trunc[i,c(1,2)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,1] = 0
    dataParams_trunc[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1), positive=TRUE)$root
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    dataParams_trunc[i,c(3,4)] <- NA
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,2] = 0
    dataParams_trunc[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Use a 3 parameter Weibull. Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_trunc[i,5] <- NA
      dataParams_trunc[i,6] <- NA
    }
    else{
      if (right != Inf){
        if (points$seismic[i] == 5){
          #No matter what the variance is, this distribution by definition of being truncated at 5 cannot be fit.
          lb_mat[i,3] = 2.5
          dataParams_trunc[i,5] <- NA
          dataParams_trunc[i,6] <- NA
        }else{
          lb = 2.5
          lb_mat[i,3] = lb
          dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
          
          # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
          #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
          #   dataParams_trunc[i,6] <- 1.5
          #   k = dataParams_trunc[i,6]
          #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
          #   rm(k)
          # }
          
          lb=0
        }
      }else{
        lb = 2.5
        lb_mat[i,3] = lb
        dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull, start=c(1,2), positive=TRUE)$root
        
        # if (dataParams_trunc[i,5] == 0 | dataParams_trunc[i,6] == 0){
        #   #The lower bound is very close to 0, try fitting only 1 parameter with a fixed shape.
        #   dataParams_trunc[i,6] <- 1.5
        #   k = dataParams_trunc[i,6]
        #   dataParams_trunc[i,5] <- uniroot(solveWeibull_u, interval=c(0.001,1000))$root
        #   rm(k)
        # }
        
        lb=0
      }
    }
  }
  else{
    #Use a 2 parameter Weibull because the lower bound will be 0
    lb_mat[i,3] = 0
    dataParams_trunc[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,2), positive=TRUE)$root
  }
}
rm(mom)

roots2_trunc <- matrix(NA,10,3) #Percentiles to find roots for
for(i in ind_use){
  params <- dataParams_trunc[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(min(lb_mat[i,]),5))$root
    
    if (roots2_trunc[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots2_trunc[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
    
    ps <- 0.5
    roots2_trunc[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
    
    ps <- 0.95
    roots2_trunc[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  }
}
rm(params, ps)

# solving for the Beta distribution parameters
dataParams_Beta <- matrix(NA,10,4*2) #4 is for the number of RFs. 2 is number of parameters. 10 is number of places
dpBeta <- matrix(NA,10,4) #Number of iterations
#set the lower bound of the Beta as 0 and the upper bound as 5
lb = 0
ub = 5
#Record the value of the lower and upper bounds in a matrix
lb_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
ub_mat = matrix(0, nrow=nrow(points), ncol=ncol(dataParams)/2)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  #check if thermal is 5 with 0 uncertainty
  if ((points$thermal[i] == 5 & points$th_pfa_var5[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,1] = 5
    ub_mat[i,1] = 5
    #Set these parameters to NA because the roots are all 5
    dataParams_Beta[i,c(1,2)] <- NA
    dpBeta[i,1] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,1] = 0
    ub_mat[i,1] = 5
    dataParams_Beta[i,c(1,2)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,1] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  #check if reservoir rfc is 5 with 0 uncertainty
  if ((points$reservoir_rfc[i] == 5 & points$re_pfa_var5_rfc[i] == 0)){
    #Use a lower bound of 5
    lb_mat[i,2] = 5
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- NA
    dpBeta[i,2] <-  NA
  }
  else{
    #Solve for the roots with a lower bound of 0
    lb_mat[i,2] = 0
    ub_mat[i,2] = 5
    dataParams_Beta[i,c(3,4)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,2] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  #check if either seismic stress or earthquake is a 5 with no uncertainty
  if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) | (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
    #Need to find the lower bound.
    if ((points$seis_eq_pred[i] > 25000 & points$seis_eq_err[i] == 0) & (points$seis_stress_pred[i] > 25 & points$seis_stress_err[i] == 0)){
      lb_mat[i,3] = 5
      ub_mat[i,3] = 5
      #These have no need to be fit because the lower bound is a 5. Set parameters to NA.
      dataParams_Beta[i,5] <- NA
      dataParams_Beta[i,6] <- NA
    }
    else{
      lb = 2.5
      lb_mat[i,3] = lb
      ub_mat[i,3] = ub
      dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta, start=c(1,1), positive=TRUE)$root
      dpBeta[i,3] <-  multiroot(solveBeta,start=c(100,100), positive=TRUE)$iter
      
      if (dataParams_Beta[i,6] == 0){
        #Did not converge. Try a fixed high alpha low beta
        dataParams_Beta[i,c(5,6)] = c(100,0.001)
      }
      
      if (dataParams_Beta[i,5] == 0){
        #Did not converge. Try a fixed low alpha high beta
        dataParams_Beta[i,c(5,6)] = rev(c(100,0.001))
      }
      
      #Set the lower bound back to 0
      lb=0
    }
  }
  else{
    #Solve using a lower bound of 0
    lb_mat[i,3] = 0
    ub_mat[i,3] = 5
    dataParams_Beta[i,c(5,6)] <- multiroot(solveBeta,start=c(1,1), positive=TRUE)$root
    dpBeta[i,3] <-  multiroot(solveBeta,start=c(1,1), positive=TRUE)$iter
  }
}
rm(mom)

roots2_Beta <- matrix(NA,10,3) #Percentiles to find roots for
converge2_Beta <- matrix(NA,10,3) #Iterations to converge
for(i in ind_use){
  params <- dataParams_Beta[i,]
  
  #Determine how to search for the root:
  if (any(lb_mat[i,] != 0)){
    #The lower bound is non-zero for some terms. Solve for the value using only the smallest values.
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.05
    roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,1] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.50
    roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,2] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
    
    #Next Percentile
    ind_params = which(lb_mat[i,] == min(lb_mat[i,]))
    
    ps <- 0.95
    roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$root
    converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(min(lb_mat[i,]),5))$iter
    
    if (roots2_Beta[i,3] > sort(unique(lb_mat[i,]))[2]){
      #The value of this percentile should be checked for the lower bound.
      ind_params = which(lb_mat[i,] <= sort(unique(lb_mat[i,]))[2])
      
      roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$root
      converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(sort(unique(lb_mat[i,]))[2],5))$iter
    }
  }else{
    #The lower bound is 0 for all risk factors.
    ind_params = which(lb_mat[i,] == 0)
    
    ps <- 0.05
    roots2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,1] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.5
    roots2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,2] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
    
    ps <- 0.95
    roots2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$root
    converge2_Beta[i,3] <- uniroot(solveQuant_Beta,interval=c(0,5))$iter
  }
}
rm(params, ps)

# loop for Monte Carlo
right = Inf #Set to Inf for non-truncated Weibull plots
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  
  reMed[i] <- median(mat_mc[,1])
  thMed[i] <-  median(mat_mc[,2])
  seMed[i] <-  median(0.5*(mat_mc[,3]+mat_mc[,4]))
  
  setwd(wd_image)
  par(xpd=T)
  if (right != Inf){
    name = 'scdist_trunc_g'
  }else{
    name = 'scdist_g'
  }
  png(paste(name,i,'.png',sep='')
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
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])/pweibull(5,shape=dataParams_trunc[ind_use[i],4],scale=dataParams_trunc[ind_use[i],3])
          ,col='royalblue'
          ,lwd=2
    )
  }else{
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
          ,col='royalblue'
          ,lwd=2
    )
  }
  points(points$reservoir[ind_use[i]],0
         ,pch=19
  )
  
  hist(mat_mc[,2]
       ,freq=F
       ,xlim=c(0,5)
       ,main="Thermal"
       ,xlab="SFF"
  )
  if (right != Inf){
    #Truncated Weibull pdf
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])/pweibull(5,shape=dataParams_trunc[ind_use[i],2],scale=dataParams_trunc[ind_use[i],1])
          ,col='royalblue'
          ,lwd=2
    )
  }else{
    lines(seq(0,5,0.01)
          ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
          ,col='royalblue'
          ,lwd=2
    )
  }
  points(points$thermal[ind_use[i]],0
         ,pch=19
  )
  
  hist(0.5*(mat_mc[,3]+mat_mc[,4])
       ,freq=F
       ,xlim=c(0,5)
       ,main="Seismic"
       ,xlab="SFF"
  )
  if (right != Inf){
    #Use truncated Weibull pdf
    if (lb_mat[ind_use[i],3] != 0){
      lines(seq(lb_mat[ind_use[i],3],right,0.01)
            ,dweibull(seq(0,(right - lb_mat[ind_use[i],3]),0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull((right - lb_mat[ind_use[i],3]), shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
    else{
      lines(seq(0,5,0.01)
            ,dweibull(seq(0,5,0.01),shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])/pweibull(5,shape=dataParams_trunc[ind_use[i],6],scale=dataParams_trunc[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
  }else{
    if (lb_mat[ind_use[i],3] != 0){
      lines(seq(lb_mat[ind_use[i],3],5+lb_mat[ind_use[i],3],0.01)
            ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
    else{
      lines(seq(0,5,0.01)
            ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
            ,col='royalblue'
            ,lwd=2
      )
    }
  }
  points(points$seismic[ind_use[i]],0
         ,pch=19
  )
  
  par(xpd=T)
  dev.off()
  
  setwd(wd_image)
  par(xpd=T)
  png(paste('scdist_Beta_g',i,'.png',sep='')
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
        ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],3],dataParams_Beta[ind_use[i],4], 0, 5)
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
        ,dBeta_ab(seq(0,5,0.01),dataParams_Beta[ind_use[i],1],dataParams_Beta[ind_use[i],2], 0, 5)
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
  lines(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01)
        ,dBeta_ab(seq(lb_mat[ind_use[i],3], ub_mat[ind_use[i],3],0.01), dataParams_Beta[ind_use[i],5], dataParams_Beta[ind_use[i],6], lb_mat[ind_use[i],3], ub_mat[ind_use[i],3])
        ,col='royalblue'
        ,lwd=2
  )
  points(points$seismic[ind_use[i]],0
         ,pch=19
  )
  
  par(xpd=T)
  dev.off()
}
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Geomean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean_g)){
  
  boxplot(distsgeomean_g[,i]
          ,at=i
          ,col=cols[i]
          ,add=TRUE
          ,pars=list(outcex=0.5,outpch=19)
  )
}
par(xpd=TRUE)
text(seq(1,9,1)
     ,rep(-0.5,9)
     ,names[ind_use]
     ,adj=0
     ,srt=270)
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean_g)
png('violin_5_rfc_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,names[ind_use]
     ,adj=0
     ,srt=270)

text(-1
     ,2.5
     ,"Combined Metric: Geologic Geomean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_rfc_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))
namesNA <- NULL

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    namesNA <- c(namesNA,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=namesNA
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#With utilization
png('parallel_axis_5_rfc_geo_g_util.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,10))
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

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))
namesNA <- NULL

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,3,1)
          ,c(check, points$utilization[i])
          ,lwd=3
          ,col=cols[i])
    namesNA <- c(namesNA,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2,3)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic', 'Utilization')
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
legend(x=3.0,y=5
       ,legend=namesNA
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)


#### RFC Geologic Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_rfc] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_rfc_geo_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for minimum for rfc
      ind_max <- which(extracted2_df_5_rfc_geo_m$co_5_0_5_m_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_m$co_5_0_5_m_rfc_geo[inds]))
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_m[inds[ind_max],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_rfc_geo_m$co_5_0_5_m_rfc_geo[inds] %in% max(extracted2_df_5_rfc_geo_m$co_5_0_5_m_rfc_geo[inds]))
      #Take the value with the lowest uncertainty. Use geomean for now.
      ind_max2 <- which(extracted2_df_5_rfc_geo_m$co_pfa_sd5_geomean_rfc_geo[inds][ind_max] %in% min(extracted2_df_5_rfc_geo_m$co_pfa_sd5_geomean_rfc_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_rfc_geo_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_rfc_geo_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_rfc,'x','y')] <- extracted2_df_5_rfc_geo_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_rfc,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_rfc[i],points$re_pfa_var5_rfc[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir rfc MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[ind_use[i]])
  res_cv <- points$res_pred_rfc_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))
  
  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_rfc_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin_g)){
  
  boxplot(distsmin_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin_g)
png('violin_5_rfc_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_rfc_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_rfc[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIw Geologic Average ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_geo_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_geo_a$co_5_0_16_s_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_a$co_5_0_16_s_RPIw_geo[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_a[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_geo_a$co_5_0_16_s_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_a$co_5_0_16_s_RPIw_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_geo_a$co_pfa_sd5_avg_RPIw_geo[inds][ind_max] %in% min(extracted2_df_5_RPIw_geo_a$co_pfa_sd5_avg_RPIw_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_geo_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_geo_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg_g)){
  
  boxplot(distsavg_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg_g)

setwd(wd_image)
png('violin_5_RPIw_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)


#### RPIw Geologic Geomean ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_geo_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_geo_g$co_5_0_125_p_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_g$co_5_0_125_p_RPIw_geo[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_g[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_geo_g$co_5_0_125_p_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_g$co_5_0_125_p_RPIw_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_geo_g$co_pfa_sd5_geomean_RPIw_geo[inds][ind_max] %in% min(extracted2_df_5_RPIw_geo_g$co_pfa_sd5_geomean_RPIw_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_geo_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_geo_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Geomean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean_g)){
  
  boxplot(distsgeomean_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean_g)
png('violin_5_RPIw_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Geomean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIw Geologic Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIw] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIw_geo_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIw
      ind_max <- which(extracted2_df_5_RPIw_geo_m$co_5_0_5_m_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_m$co_5_0_5_m_RPIw_geo[inds]))
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_m[inds[ind_max],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIw_geo_m$co_5_0_5_m_RPIw_geo[inds] %in% max(extracted2_df_5_RPIw_geo_m$co_5_0_5_m_RPIw_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIw_geo_m$co_pfa_sd5_geomean_RPIw_geo[inds][ind_max] %in% min(extracted2_df_5_RPIw_geo_m$co_pfa_sd5_geomean_RPIw_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIw_geo_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIw_geo_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIw,'x','y')] <- extracted2_df_5_RPIw_geo_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIw,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIw[i],points$re_pfa_var5_RPIw[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIw_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIw_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin_g)){
  
  boxplot(distsmin_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin_g)
png('violin_5_RPIw_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIw_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIw[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Geologic Average ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_geo_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_geo_a$co_5_0_16_s_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_a$co_5_0_16_s_RPIg_geo[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_a[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_geo_a$co_5_0_16_s_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_a$co_5_0_16_s_RPIg_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_geo_a$co_pfa_sd5_avg_RPIg_geo[inds][ind_max] %in% min(extracted2_df_5_RPIg_geo_a$co_pfa_sd5_avg_RPIg_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_geo_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_geo_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg_g)){
  
  boxplot(distsavg_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
#dummy variable
dists2 <- as.data.frame(distsavg_g)

setwd(wd_image)
png('violin_5_RPIg_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_geo_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Geologic Geomean ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_geo_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_geo_g$co_5_0_125_p_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_g$co_5_0_125_p_RPIg_geo[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_g[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_geo_g$co_5_0_125_p_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_g$co_5_0_125_p_RPIg_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_geo_g$co_pfa_sd5_geomean_RPIg_geo[inds][ind_max] %in% min(extracted2_df_5_RPIg_geo_g$co_pfa_sd5_geomean_RPIg_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_geo_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_geo_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,4*2) #4 is for the number of RFs
dp2 <- matrix(NA,10,4)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
#Set all places with a 0 to NA. 
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Geomean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean_g)){
  
  boxplot(distsgeomean_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean_g)
png('violin_5_RPIg_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Geomean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_geo_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

#### RPIg Geologic Minimum ####
# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2
#Change for each different extracted2 variable.
points[comb_names_5_RPIg] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_RPIg_geo_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for geomean for RPIg
      ind_max <- which(extracted2_df_5_RPIg_geo_m$co_5_0_5_m_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_m$co_5_0_5_m_RPIg_geo[inds]))
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_m[inds[ind_max],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_RPIg_geo_m$co_5_0_5_m_RPIg_geo[inds] %in% max(extracted2_df_5_RPIg_geo_m$co_5_0_5_m_RPIg_geo[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_RPIg_geo_m$co_pfa_sd5_geomean_RPIg_geo[inds][ind_max] %in% min(extracted2_df_5_RPIg_geo_m$co_pfa_sd5_geomean_RPIg_geo[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_RPIg_geo_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_RPIg_geo_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_RPIg,'x','y')] <- extracted2_df_5_RPIg_geo_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_RPIg,'x','y')),1)]
    }
  }
}
points$util_err <- 5
rm(inds, ind_max, ind_max2, ind_max3, nm)

# solving for the weibull distribution parameters
dataParams <- matrix(NA,10,3*2) #3 is for the number of RFs
dp2 <- matrix(NA,10,3)
for(i in ind_use){
  
  mom <- c(points$thermal[i],points$th_pfa_var5[i])
  dataParams[i,c(1,2)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,1] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$reservoir_RPIg[i],points$re_pfa_var5_RPIg[i])
  dataParams[i,c(3,4)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,2] <-  multiroot(solveWeibull,start=c(1,1))$iter
  
  mom <- c(points$seismic[i],points$se_pfa_var5[i])
  dataParams[i,c(5,6)] <- multiroot(solveWeibull,start=c(1,1))$root
  dp2[i,3] <-  multiroot(solveWeibull,start=c(1,1))$iter
}
rm(mom)

# setting NAs for seismic 5 and utilization 0
dataParams[which(dataParams[,5]  > 4.9999),5] <- NA
dataParams[which(dataParams[,5] %in% NA),6] <- NA

roots1 <- matrix(NA,10,3) #All RFs
roots2 <- matrix(NA,10,3) #Geologic Only
converge1 <- matrix(NA,10,3) #Iterations to converge
converge2 <- matrix(NA,10,3)
for(i in ind_use){
  
  params <- dataParams[i,] #params is taken in solveQuant
  
  ps <- 0.05
  roots1[i,1] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,1] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.5
  roots1[i,2] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,2] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  ps <- 0.95
  roots1[i,3] <- uniroot(solveQuant,interval=c(0,5))$root
  converge1[i,3] <- uniroot(solveQuant,interval=c(0,5))$iter
  
  params <- dataParams[i,seq(1,6,1)] #params is taken in solveQuant
  
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
rm(params, ps)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[ind_use[i]])
  res_cv <- points$res_pred_RPIg_max_err[ind_use[i]]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Is this missing????
  #dist_vars[4,i] <- var(pfm5)
  
  # calcuating overall distribution
  #Geologic
  distsavg_g[,i] <- (mat_mc[,1] + mat_mc[,2] + 0.5*mat_mc[,3]+0.5*mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*mat_mc[,2]*(0.5*mat_mc[,3]+0.5*mat_mc[,4]))^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
  
  #Fraction of time each RF is the minimum for each city
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
  rm(j)
  
  reMean[i] <- mean(mat_mc[,1])
  thMean[i] <-  mean(mat_mc[,2])
  seMean[i] <-  mean(0.5*(mat_mc[,3]+mat_mc[,4]))

  # setwd(wd_image)
  # par(xpd=T)
  # png(paste('scdist',i,'.png',sep='')
  #     ,height=6
  #     ,width=6
  #     ,units='in'
  #     ,res=300
  # )
  # 
  # par(mfrow=c(2,2)
  #     ,oma=c(0.5,0.5,0.5,0.5)+0.1
  #     ,mar=c(5,5,3,0)+0.1)
  # 
  # hist(mat_mc[,1]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Reservoir"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],4],scale=dataParams[ind_use[i],3])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$reservoir[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,2]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Thermal"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],2],scale=dataParams[ind_use[i],1])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$thermal[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(0.5*(mat_mc[,3]+mat_mc[,4])
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Seismic"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],6],scale=dataParams[ind_use[i],5])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$seismic[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # hist(mat_mc[,5]
  #      ,freq=F
  #      ,xlim=c(0,5)
  #      ,main="Utilization"
  #      ,xlab="SFF"
  # )
  # lines(seq(0,5,0.01)
  #       ,dweibull(seq(0,5,0.01),shape=dataParams[ind_use[i],8],scale=dataParams[ind_use[i],7])
  #       ,col='royalblue'
  #       ,lwd=2
  # )
  # points(points$utilization[ind_use[i]],0
  #        ,pch=19
  # )
  # 
  # par(xpd=T)
  # dev.off()
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_RPIg_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: Geologic Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin_g)){
  
  boxplot(distsmin_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin_g)
png('violin_5_RPIg_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: Geologic Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()

## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_RPIg_geo_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points)){
  
  # check for NA values
  check <- c(points$reservoir_RPIg[i]
             ,points$thermal[i]
             ,points$seismic[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[i])
    names <- c(names,points$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Reservoir','Thermal','Seismic')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[-which(IndNA == 1)]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)



#### EGS ####
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  extracted2_df_5_egs_a[nm] <- sqrt((extracted2_df_5_egs_a$x - points$x[i])^2 + (extracted2_df_5_egs_a$y - points$y[i])^2)
  extracted2_df_5_egs_g[nm] <- sqrt((extracted2_df_5_egs_g$x - points$x[i])^2 + (extracted2_df_5_egs_g$y - points$y[i])^2)
  extracted2_df_5_egs_m[nm] <- sqrt((extracted2_df_5_egs_m$x - points$x[i])^2 + (extracted2_df_5_egs_m$y - points$y[i])^2)
}

# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2

#### EGS Minimum ####
#Change for each different extracted2 variable.
points[comb_names_5_egs] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_egs_m[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for egs minimum
      ind_max <- which(extracted2_df_5_egs_m$co_5_0_5_m_egs[inds] %in% max(extracted2_df_5_egs_m$co_5_0_5_m_egs[inds]))
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_m[inds[ind_max],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_egs_m$co_5_0_5_m_egs[inds] %in% max(extracted2_df_5_egs_m$co_5_0_5_m_egs[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_egs_m$co_pfa_sd5_geomean_egs[inds][ind_max] %in% min(extracted2_df_5_egs_m$co_pfa_sd5_geomean_egs[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_egs_m[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_egs_m[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_m[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }
  }
}
rm(inds, ind_max, ind_max2, ind_max3, nm)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,2] <- pfm5
  dist_vars[2,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,4] <- pfm5

  # calcuating overall distribution
  #EGS
  distsavg_g[,i] <- (mat_mc[,1] + 0.5*mat_mc[,2]+0.5*mat_mc[,3] + mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*(0.5*mat_mc[,2]+0.5*mat_mc[,3])*mat_mc[,4])^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],0.5*mat_mc[,2]+0.5*mat_mc[,3],mat_mc[,4]),1,min)
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_egs_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: EGS Minimum'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsmin_g)){
  
  boxplot(distsmin_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsmin_g)
png('violin_5_egs_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: EGS Minimum"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()


## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_egs_m.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points[ind_use])){
  
  # check for NA values
  check <- c(points[ind_use,]$thermal[i]
             ,points[ind_use,]$seismic[i]
             ,points[ind_use,]$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[ind_use][i])
    names <- c(names,points[ind_use,]$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Thermal','Seismic', 'Utilization')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[ind_use]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)


#### EGS Geometric Mean ####

# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2

#Change for each different extracted2 variable.
points[comb_names_5_egs] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_egs_g[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for egs geometric mean
      ind_max <- which(extracted2_df_5_egs_g$co_5_0_125_p_egs[inds] %in% max(extracted2_df_5_egs_g$co_5_0_125_p_egs[inds]))
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_g[inds[ind_max],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_egs_g$co_5_0_125_p_egs[inds] %in% max(extracted2_df_5_egs_g$co_5_0_125_p_egs[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_egs_g$co_pfa_sd5_geomean_egs[inds][ind_max] %in% min(extracted2_df_5_egs_g$co_pfa_sd5_geomean_egs[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_egs_g[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_egs_g[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_g[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }
  }
}
rm(inds, ind_max, ind_max2, ind_max3, nm)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,2] <- pfm5
  dist_vars[2,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  #EGS
  distsavg_g[,i] <- (mat_mc[,1] + 0.5*mat_mc[,2]+0.5*mat_mc[,3] + mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*(0.5*mat_mc[,2]+0.5*mat_mc[,3])*mat_mc[,4])^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],0.5*mat_mc[,2]+0.5*mat_mc[,3],mat_mc[,4]),1,min)
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_egs_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: EGS Geometric Mean'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsgeomean_g)){
  
  boxplot(distsgeomean_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsgeomean_g)
png('violin_5_egs_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

dists3=dists2
dists2 = as.data.frame(distsavg)
#Call a blank vioplot to get the dimensions of the plot correct. 
vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col='white'
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,rectCol='white'
        ,border='white'
        ,drawRect=FALSE
)
dists2=dists3
vioplot(dists2[,2],dists2[,3],dists2[,5],dists2[,6],dists2[,7],dists2[,8]
        ,col=cols[c(2,3,5,6,7,8)]
        ,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',6)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
        ,at=c(2,3,5,6,7,8)
        ,add=TRUE
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
     ,2.5
     ,"Combined Metric: EGS Geometric Mean"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()


## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_egs_g.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points[ind_use])){
  
  # check for NA values
  check <- c(points[ind_use,]$thermal[i]
             ,points[ind_use,]$seismic[i]
             ,points[ind_use,]$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[ind_use][i])
    names <- c(names,points[ind_use,]$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Thermal','Seismic', 'Utilization')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[ind_use]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)


#### EGS Average ####

# points of interest
points <- as.data.frame(poi2)
points$x <- points$coords.x1
points$y <- points$coords.x2

#Change for each different extracted2 variable.
points[comb_names_5_egs] <- NA 

# extracting the values corresponding to the maximum within 10 km for points
for(i in 1:nrow(points)){
  nm <- paste('dist',i,sep='')
  inds <- which(extracted2_df_5_egs_a[nm] < 10000)
  
  if(length(inds) != 0){
    if(length(inds) == 1){
      #One max for egs geometric mean
      ind_max <- which(extracted2_df_5_egs_a$co_5_0_16_s_egs[inds] %in% max(extracted2_df_5_egs_a$co_5_0_16_s_egs[inds]))
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_a[inds[ind_max],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }else{
      #Find all points
      ind_max <- which(extracted2_df_5_egs_a$co_5_0_16_s_egs[inds] %in% max(extracted2_df_5_egs_a$co_5_0_16_s_egs[inds]))
      #Take the value with the lowest uncertainty.
      ind_max2 <- which(extracted2_df_5_egs_a$co_pfa_sd5_avg_egs[inds][ind_max] %in% min(extracted2_df_5_egs_a$co_pfa_sd5_avg_egs[inds][ind_max]))
      #If there are more than 1, then take the one with the shortest distance to the point. If there are still more than 1, take the first one:
      ind_max3 <- which(extracted2_df_5_egs_a[nm][inds[ind_max][ind_max2],] %in% min(extracted2_df_5_egs_a[nm][inds[ind_max][ind_max2],]))[1]
      
      points[i,c(comb_names_5_egs,'x','y')] <- extracted2_df_5_egs_a[inds[ind_max][ind_max2][ind_max3],seq(1,length(c(comb_names_5_egs,'x','y')),1)]
    }
  }
}
rm(inds, ind_max, ind_max2, ind_max3, nm)

# loop for Monte Carlo
for(i in 1:length(ind_use)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[ind_use[i]]
  th_se <- points$therm_err[ind_use[i]]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,1] <- pfm5
  
  dist_vars[1,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[ind_use[i]]
  se_eq_se <- points$seis_eq_err[ind_use[i]]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,2] <- pfm5
  dist_vars[2,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[ind_use[i]]
  se_st_se <- points$seis_stress_err[ind_use[i]]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # Old Method - uses absolute value for rand
  #rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
  #rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[ind_use[i]],points$util_pred[ind_use[i]]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  #EGS
  distsavg_g[,i] <- (mat_mc[,1] + 0.5*mat_mc[,2]+0.5*mat_mc[,3] + mat_mc[,4])/3
  distsgeomean_g[,i] <- (mat_mc[,1]*(0.5*mat_mc[,2]+0.5*mat_mc[,3])*mat_mc[,4])^(1/3)
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],0.5*mat_mc[,2]+0.5*mat_mc[,3],mat_mc[,4]),1,min)
}
rm(rand, pfm5, i)

# making boxplot
cols <- c(brewer.pal(9,'Set1'))
setwd(wd_image)
png('boxplot_5_egs_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))
boxplot(NA
        ,xlim=c(0,10)
        ,ylim=c(0,5)
        ,ylab='Combined Metric: EGS Average'
        ,yaxs='i'
        ,yaxt='n'
)

for(i in 1:ncol(distsavg_g)){
  
  boxplot(distsavg_g[,i]
          ,at=i
          ,col=cols[i]
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
text(-1.2,5
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

# making violin plot
dists2 = as.data.frame(distsavg_g)
png('violin_5_egs_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(mar=c(8,5,2,2))

vioplot(dists2[,1],dists2[,2],dists2[,3],dists2[,4],dists2[,5],dists2[,6],dists2[,7],dists2[,8],dists2[,9]
        ,col=cols
        #,xlim=c(0,10)
        #,xlab=''
        #,xaxt='n'
        #,ylab='Combined Metric: Sum'
        ,names=rep('',9)
        ,ylim=c(0,5)
        ,h=0.04
        ,colMed='white'
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
     ,2.5
     ,"Combined Metric: EGS Average"
     ,adj=0.5
     ,srt=90)

axis(2,at=seq(0,5,1))
par(xpd=FALSE)
dev.off()


## making parallel axis plot
# setting exporting parameters
cols <- c(brewer.pal(9,'Set1'))
cols[3:10] = c('black', cols[3:9])

# setting working directory and name
setwd(wd_image)
png('parallel_axis_5_egs_a.png'
    ,height=5
    ,width=7
    ,units='in'
    ,res=300
)
par(xpd=FALSE)
par(mar=c(3,5,1,12))
plot(NA,NA
     ,xlim=c(0,2)
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

names <- NULL
points$names <- c('Corning, NY', 'Elmira, NY', 'Ithaca, NY', 'Jamestown, NY'
                  ,'Mansfield, PA', 'Meadville, PA','Sayre, PA'
                  ,'Charleston, WV', 'Morgantown, WV', 'Pineville, WV')

#Track NA indices:
IndNA = vector('numeric', length=nrow(points))

# loop to add lines to plots
for(i in 1:nrow(points[ind_use])){
  
  # check for NA values
  check <- c(points[ind_use,]$thermal[i]
             ,points[ind_use,]$seismic[i]
             ,points[ind_use,]$utilization[i])
  # plot lines only if there are no NAs
  if(sum(is.na(check)) == 0){
    
    lines(seq(0,2,1)
          ,check
          ,lwd=3
          ,col=cols[ind_use][i])
    names <- c(names,points[ind_use,]$names[i])
  }
  else{
    #Keep a record of the NA indexes
    IndNA[i] = 1
  }
}
rm(check)

par(xpd=TRUE)
text(c(0,1,2)
     ,c(-0.3)
     ,c('Thermal', 'Seismic', 'Utilization')
     ,adj=0.5)
text(-0.25,5
     ,'Better'
     ,adj=1
     ,col='darkgreen'
     ,font=2)
text(-0.25,0
     ,'Worse'
     ,adj=1
     ,col='firebrick'
     ,font=2)
legend(x=2.0,y=5
       ,legend=names
       ,col=cols[ind_use]
       ,lwd=3
       ,lty=1
       ,bty='n'
       ,pch=NA)
par(xpd=FALSE)
dev.off()
rm(IndNA)

##### Analytical vs Monte Carlo #####
# making plot comparing analytical and monte carlo results
# all variables
setwd(wd_image)
par(mar=c(1,1,1,1)*0.1)
png('three_panel2_all_med_beta_trunc_dtn_ftn_rfc.png'
     ,height=6
     ,width=6
     ,units='in'
     ,res=300
)
 
par(mfrow=c(1,5)
     ,oma=c(1,0,0,3)+0.1
     ,mar=c(5,0,0,0)+0.1)
dshift1 <- 0.15
dshift2 <- 0.15
dshift3 <- 0.15
dshift4 <- 0.15
dshift5 <- 0.15
dshift6 <- 0.15

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
        ,labels=names[ind_use[i]]
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
   
   res_mean <- points$reservoir_rfc[ind_use[i]]
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
           ,col='black')
     #Truncated folded normal
     lines(roots1_ftn[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2+dshift3+dshift4+dshift5+dshift6,i-1+dshift1+dshift2+dshift3+dshift4+dshift5+dshift6)
           ,lwd=2
           ,col='orange4')
     #Doubly truncated normal
     lines(roots1_dtn[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2+dshift3+dshift4+dshift5,i-1+dshift1+dshift2+dshift3+dshift4+dshift5)
           ,lwd=2
           ,col='purple4')
     #Beta
     lines(roots1_Beta[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2+dshift3+dshift4,i-1+dshift1+dshift2+dshift3+dshift4)
           ,lwd=2
           ,col='red4')
     #Truncated Weibull
     lines(roots1_trunc[ind_use[i],c(1,3)],c(i-1+dshift1+dshift2+dshift3,i-1+dshift1+dshift2+dshift3)
           ,lwd=2
           ,col='darkgreen')
     }
   
   points(min(res_mean,therm_mean,seis_mean,util_mean),i-1+dshift1+dshift2
          ,pch=4
          ,col='black'
          ,cex=1.5
          )
   
   points(roots1[ind_use[i],2],i-1+dshift1+dshift2
          ,pch=3
          ,col='blue'
          ,cex=1.5
   )
   points(roots1_trunc[ind_use[i],2],i-1+dshift1+dshift2+dshift3
          ,pch=3
          ,col='springgreen'
          ,cex=1.5
   )
   
   points(roots1_Beta[ind_use[i],2],i-1+dshift1+dshift2+dshift3+dshift4
          ,pch=3
          ,col='red'
          ,cex=1.5
   )
   
   points(roots1_dtn[ind_use[i],2],i-1+dshift1+dshift2+dshift3+dshift4+dshift5
          ,pch=3
          ,col='purple'
          ,cex=1.5
   )
   
   points(roots1_ftn[ind_use[i],2],i-1+dshift1+dshift2+dshift3+dshift4+dshift5+dshift6
          ,pch=3
          ,col='orange'
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
   
   points(median(distsmin[,i]),i-1+dshift1
          ,pch=3
          ,col='skyblue'
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
   res_mean <- points$reservoir_rfc[ind_use[i]]
   res_var_ls <- points$re_pfa_var5_rfc_ls[ind_use[i]]
   
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
   res_mean <- points$reservoir_rfc[ind_use[i]]
   res_var <- points$re_pfa_var5_rfc[ind_use[i]]
   
   therm_mean <- points$thermal[ind_use[i]]
   therm_var <- points$th_pfa_var5[ind_use[i]]
   
   seis_mean <- points$seismic[ind_use[i]]
   seis_var <- points$se_pfa_var5[ind_use[i]]
   
   util_mean <- points$utilization[ind_use[i]]
   util_var <- points$util_pfa_var5[ind_use[i]]
   
   var_avg <- (1/16)*(seis_var + res_var + therm_var + util_var)
   mean_avg <- (res_mean + therm_mean + seis_mean + util_mean)/4
   
   # calculated values seis_mean
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
 
# legend('center'
#         ,legend=c('Favorability \n Metric Value','Monte Carlo \n Mean','Approximate \n 90% PI','Monte Carlo \n 90% PI')
#         ,pch=c(4,4,NA,NA)
#         ,col=c('black','gray48','black','gray48')
#         ,lwd=c(NA,NA,2,2)
#         ,y.intersp = 1.5
# )

legend('right'
        ,legend=c('Approximate \n Mean', 'Monte Carlo \n Mean','3 Parameter \n Weibull Median','Truncated \n Weibull Median','Beta Median','Truncated \n Normal Median','Folded \n Normal Median','Monte Carlo \n Median','Approximate \n 90% PI','Truncated \n Weibull 90% PI','Beta \n 90% PI','Truncated \n Normal 90% PI','Folded \n Normal 90% PI','Monte Carlo \n 90% PI')
        ,pch=c(4,4,3,3,3,3,3,3,NA,NA,NA,NA,NA,NA)
        ,col=c('black','gray48','blue','springgreen','red','purple','orange','skyblue','black','darkgreen','red4','purple4','orange4','gray48')
        ,lwd=c(NA,NA,NA,NA,NA,NA,NA,NA,2,2,2,2,2,2)
        ,y.intersp = 1.5
)
par(xpd=F)
dev.off()


# par(mar=c(1,1,1,1)*0.1)
# png('three_panel2_all_med.png'
#     ,height=6
#     ,width=6
#     ,units='in'
#     ,res=300
# )
# tiff('three_panel2_all_med.tiff'
#     ,height=6
#     ,width=6
#     ,units='in'
#     ,res=600
# )

setEPS()
postscript(file = "three_panel2_all_med.eps", width = 6, height = 6)

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
       ,labels=names[ind_use[i]]
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
  
  res_mean <- points$reservoir_rfc[ind_use[i]]
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
          ,col='black')
  }
  
  points(min(res_mean,therm_mean,seis_mean,util_mean),i-1+dshift1+dshift2
         ,pch=4
         ,col='black'
         ,cex=1.5
  )
  
  points(roots1[ind_use[i],2],i-1+dshift1+dshift2
         ,pch="|"
         ,col='black'
         ,cex=1.2
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
  
  points(median(distsmin[,i]),i-1+dshift1
         ,pch="|"
         ,col='grey48'
         ,cex=1.2
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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_rfc_ls[ind_use[i]]
  
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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var <- points$re_pfa_var5_rfc[ind_use[i]]
  
  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]
  
  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]
  
  util_mean <- points$utilization[ind_use[i]]
  util_var <- points$util_pfa_var5[ind_use[i]]
  
  var_avg <- (1/16)*(seis_var + res_var + therm_var + util_var)
  mean_avg <- (res_mean + therm_mean + seis_mean + util_mean)/4
  
  # calculated values seis_mean
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

legend('right'
       ,legend=c('Approximate \n Mean', 'Monte Carlo \n Mean','3 Parameter \n Weibull Median', 'Monte Carlo \n Median','Approximate \n 90% PI','Monte Carlo \n 90% PI')
       ,pch=c("X","X","|","|",NA,NA)
       ,col=c('black','gray48','black','grey48','black','gray48')
       ,lwd=c(NA,NA,NA,NA,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()

 
#### Geologic Only #### 
setwd(wd_image)
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
dshift3 <- 0.2
dshift4 <- 0.2

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
       ,labels=names[i]
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

  res_mean <- points$reservoir_rfc[ind_use[i]]

  therm_mean <- points$thermal[ind_use[i]]

  seis_mean <- points$seismic[ind_use[i]]

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
  
  points(roots2[ind_use[i],2],i-1+dshift1+dshift2
         ,pch=3
         ,col='blue'
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
  points(median(distsmin_g[,i]),i-1+dshift1
         ,pch=3
         ,col='skyblue'
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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_rfc_ls[ind_use[i]]

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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var <- points$re_pfa_var5_rfc[ind_use[i]]

  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]

  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]

  var_avg <- (1/9)*(seis_var+res_var+therm_var)
  mean_avg <- (res_mean+therm_mean +seis_mean)/3
  
  # calculated values 
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
       ,legend=c('Approximate \n Mean','Monte Carlo \n Mean','3 Parameter \n Weibull Median','Monte Carlo \n Median','Approximate \n 90% PI','Monte Carlo \n 90% PI')
       ,pch=c(4,4,3,3,NA,NA)
       ,col=c('black','gray48','blue','skyblue','black','gray48')
       ,lwd=c(NA,NA,NA,NA,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()

##### Metric Comparison #####
# making plot of different metrics for all favorability factors
setwd(wd_image)
# png('single_panel2_all_rfc.png'
#     ,height=6
#     ,width=7
#     ,units='in'
#     ,res=300
# )

color = grey.colors(4)
setEPS()
postscript('single_panel2_all_rfc.eps', height=6, width=7)

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
       ,labels=names[ind_use[i]]
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
  res_mean <- points$reservoir_rfc[ind_use[i]]

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
  
  # points(roots1[ind_use[i],2],i-1+dshift1
  #        ,pch=5
  #        ,col='black'
  #        ,cex=1.5
  # )

  # calculating empirical quantiles
  empQuant <- quantile(distsgeomean[,i],c(0.05,0.95))

  #extracting values value
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_rfc_ls[ind_use[i]]

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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var <- points$re_pfa_var5_rfc[ind_use[i]]

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
        ,col=color[3]
  )

  points(mean_avg,i-1+dshift1+2*dshift2
         ,pch=4
         ,col=color[3]
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
       ,col=rev(c('black','gray48',color[3]))
       ,lwd=c(NA,NA,NA,2,2,2)
       ,y.intersp = 1.5
)
par(xpd=F)
dev.off()

#### Geology Only ####
# making plot of different metrics for geologic favorability factors
setwd(wd_image)
png('single_panel2_rfc_geo.png'
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
       ,labels=names[i]
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
  res_mean <- points$reservoir_rfc[ind_use[i]]

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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var_ls <- points$re_pfa_var5_rfc_ls[ind_use[i]]

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
  res_mean <- points$reservoir_rfc[ind_use[i]]
  res_var <- points$re_pfa_var5_rfc[ind_use[i]]

  therm_mean <- points$thermal[ind_use[i]]
  therm_var <- points$th_pfa_var5[ind_use[i]]

  seis_mean <- points$seismic[ind_use[i]]
  seis_var <- points$se_pfa_var5[ind_use[i]]

  util_mean <- points$utilization[ind_use[i]]
  util_var <- points$util_pfa_var5[ind_use[i]]

  var_avg <- (1/9)*(seis_var + res_var + therm_var)
  mean_avg <- (res_mean + therm_mean + seis_mean)/3
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

##### Scatter Plots of Metrics #####
# extracting values of layers for cities
th_5_0_5_NA <- raster(paste(wd_raster_out,'/th_5_0_5_NA.tif',sep=''))
th_5_0_5_NA[th_5_0_5_NA < 0] <- NA

re_5_0_5_NA_rfc <- raster(paste(wd_raster_out,'/re_5_0_5_NA_rfc.tif',sep=''))
re_5_0_5_NA_rfc[re_5_0_5_NA_rfc < 0] <- NA

re_5_0_5_NA_RPIw <- raster(paste(wd_raster_out,'/re_5_0_5_NA_RPIw.tif',sep=''))
re_5_0_5_NA_RPIw[re_5_0_5_NA_RPIw < 0] <- NA

re_5_0_5_NA_RPIg <- raster(paste(wd_raster_out,'/re_5_0_5_NA_RPIg.tif',sep=''))
re_5_0_5_NA_RPIg[re_5_0_5_NA_RPIg < 0] <- NA

se_5_0_5_a <- raster(paste(wd_raster_out,'/se_5_0_5_a.tif',sep=''))
se_5_0_5_a[se_5_0_5_a  < 0] <- NA

ut5_5_0_5_NA <- raster(paste(wd_raster_out,'/ut5_5_0_5_NA.tif',sep=''))
ut5_5_0_5_NA[ut5_5_0_5_NA  < 0] <- NA

#### RFC ####
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA_rfc,se_5_0_5_a,ut5_5_0_5_NA))

co_min <- stack(c(th_5_0_5_NA,re_5_0_5_NA_rfc,se_5_0_5_a,ut5_5_0_5_NA))
combined_all_rasters[which(values(combined_all_rasters) %in% -9999)] <- NA

extract_pts <- extract(x=combined_all_rasters
                       ,y=cities2
                       ,sp=TRUE
                       ,nl=4
                       ,df=TRUE
                       ,na.rm=TRUE
                       ,method='simple'
                       ,buffer=2000
                       ,fun=max
)

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[13] <- "th_5"
names(extract_pts2)[14] <- "rfc_5"
names(extract_pts2)[15] <- "se_5"
names(extract_pts2)[16] <- "ut_5"

extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]
extract_pts3$co_p_5 <- extract_pts3$th_5*extract_pts3$rfc_5*extract_pts3$se_5*extract_pts3$ut_5
extract_pts3$co_s_5 <- extract_pts3$th_5+extract_pts3$rfc_5+extract_pts3$se_5+extract_pts3$ut_5
extract_pts3$co_m_5 <- apply(cbind(extract_pts3$th_5,extract_pts3$rfc_5,extract_pts3$se_5,extract_pts3$ut_5),FUN=min,MARGIN=1)

ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5  > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5  > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
# png('scatter_rfc_p_s.png'
#     ,height=6
#     ,width=6
#     ,units='in'
#     ,res=300
# )
setEPS()
postscript('scatter_rfc_p_s.eps'
    ,height=6
    ,width=6)
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
# png('scatter_rfc_m_s.png'
#     ,height=6
#     ,width=6
#     ,units='in'
#     ,res=300
# )
setEPS()
postscript('scatter_rfc_m_s.eps'
           ,height=6
           ,width=6)
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
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# setting export criteria for plot
setwd(wd_image)
# png('scatter_rfc_m_p.png'
#     ,height=6
#     ,width=6
#     ,units='in'
#     ,res=300
# )
setEPS()
postscript('scatter_rfc_m_p.eps'
           ,height=6
           ,width=6)
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

#Correlations
cor(extract_pts3$co_p_5[inds_gtr0]^.25, extract_pts3$co_s_5[inds_gtr0]/4)
cor(extract_pts3$co_p_5[inds_gtr0]^.25, extract_pts3$co_m_5[inds_gtr0])
cor(extract_pts3$co_m_5[inds_gtr0], extract_pts3$co_s_5[inds_gtr0]/4)


#All 3 Together
setwd(wd_image)
setEPS()
postscript('scatter_rfc_all.eps'
           ,height=6
           ,width=6)
par(mar = c(5,5,0.5,0.5))
layout(rbind(c(1,1,2,2), c(4,3,3,5)), widths = c(0.25,0.25,0.25,0.25))
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
     ,cex.lab = 1.5
)
grid()
points(extract_pts3$co_s_5[inds_gtr0]/4
       ,extract_pts3$co_p_5[inds_gtr0]^0.25
       ,pch=19
       ,xlim=c(0,5)
       ,ylim=c(0,5)
       ,xlab=''
       ,ylab=''
       ,xaxt='n'
       ,yaxt='n'
       ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
text('A', x = 4.75, y = 4.75, cex=1.5)
rect(xleft = 4.5, ybottom = 4.5, xright = 5, ytop = 5)

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
     ,cex.lab = 1.5
)
grid()
points(extract_pts3$co_s_5[inds_gtr0]/4
       ,extract_pts3$co_m_5[inds_gtr0]
       ,pch=19
       ,xlim=c(0,5)
       ,ylim=c(0,5)
       ,xlab=''
       ,ylab=''
       ,xaxt='n'
       ,yaxt='n'
       ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
text('B', x = 4.75, y = 4.75, cex=1.5)
rect(xleft = 4.5, ybottom = 4.5, xright = 5, ytop = 5)

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
     ,cex.lab = 1.5
)
grid()
points(extract_pts3$co_p_5[inds_gtr0]^0.25
       ,extract_pts3$co_m_5[inds_gtr0]
       ,pch=19
       ,xlim=c(0,5)
       ,ylim=c(0,5)
       ,xlab=''
       ,ylab=''
       ,xaxt='n'
       ,yaxt='n'
       ,cex=0.8
)
axis(1
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5), cex.axis=1.5)
text('C', x = 4.75, y = 4.75, cex=1.5)
rect(xleft = 4.5, ybottom = 4.5, xright = 5, ytop = 5)

dev.off()


#### RPIw ####
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA_RPIw,se_5_0_5_a,ut5_5_0_5_NA))

co_min <- stack(c(th_5_0_5_NA,re_5_0_5_NA_RPIw,se_5_0_5_a,ut5_5_0_5_NA))
combined_all_rasters[which(values(combined_all_rasters) %in% -9999)] <- NA

extract_pts <- extract(x=combined_all_rasters
                       ,y=cities2
                       ,sp=TRUE
                       ,nl=4
                       ,df=TRUE
                       ,na.rm=TRUE
                       ,method='simple'
                       ,buffer=2000
                       ,fun=max
)

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[13] <- "th_5"
names(extract_pts2)[14] <- "RPIw_5"
names(extract_pts2)[15] <- "se_5"
names(extract_pts2)[16] <- "ut_5"

extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]
extract_pts3$co_p_5 <- extract_pts3$th_5*extract_pts3$RPIw_5*extract_pts3$se_5*extract_pts3$ut_5
extract_pts3$co_s_5 <- extract_pts3$th_5+extract_pts3$RPIw_5+extract_pts3$se_5+extract_pts3$ut_5
extract_pts3$co_m_5 <- apply(cbind(extract_pts3$th_5,extract_pts3$RPIw_5,extract_pts3$se_5,extract_pts3$ut_5),FUN=min,MARGIN=1)

ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5  > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5  > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
png('scatter_RPIw_p_s.png'
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
png('scatter_RPIw_m_s.png'
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
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_RPIw_m_p.png'
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

#### RPIg ####
combined_all_rasters <- stack(c(th_5_0_5_NA,re_5_0_5_NA_RPIg,se_5_0_5_a,ut5_5_0_5_NA))

co_min <- stack(c(th_5_0_5_NA,re_5_0_5_NA_RPIg,se_5_0_5_a,ut5_5_0_5_NA))
combined_all_rasters[which(values(combined_all_rasters) %in% -9999)] <- NA

extract_pts <- extract(x=combined_all_rasters
                       ,y=cities2
                       ,sp=TRUE
                       ,nl=4
                       ,df=TRUE
                       ,na.rm=TRUE
                       ,method='simple'
                       ,buffer=2000
                       ,fun=max
)

extract_pts2 <- as.data.frame(extract_pts)
names(extract_pts2)[13] <- "th_5"
names(extract_pts2)[14] <- "RPIg_5"
names(extract_pts2)[15] <- "se_5"
names(extract_pts2)[16] <- "ut_5"

extract_pts3 <- extract_pts2[complete.cases(extract_pts2),]
extract_pts3$co_p_5 <- extract_pts3$th_5*extract_pts3$RPIg_5*extract_pts3$se_5*extract_pts3$ut_5
extract_pts3$co_s_5 <- extract_pts3$th_5+extract_pts3$RPIg_5+extract_pts3$se_5+extract_pts3$ut_5
extract_pts3$co_m_5 <- apply(cbind(extract_pts3$th_5,extract_pts3$RPIg_5,extract_pts3$se_5,extract_pts3$ut_5),FUN=min,MARGIN=1)

ny_inds <- intersect(which(extract_pts3$USPS == 'NY'),which(extract_pts3$co_p_5 > 0))
pa_inds <- intersect(which(extract_pts3$USPS == 'PA'),which(extract_pts3$co_p_5  > 0))
wv_inds <- intersect(which(extract_pts3$USPS == 'WV'),which(extract_pts3$co_p_5  > 0))

inds_gtr0 <- c(ny_inds,pa_inds,wv_inds)
# setting export criteria for plot
setwd(wd_image)
png('scatter_RPIg_p_s.png'
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
png('scatter_RPIg_m_s.png'
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
axis(2
     ,at=c(0,1,2,3,4,5)
     ,labels=c(0,1,2,3,4,5))
dev.off()

# setting export criteria for plot
setwd(wd_image)
png('scatter_RPIg_m_p.png'
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

#remove unneeded files
rm(extract_pts3, extract_pts, extract_pts2, combined_all_rasters, co_min)

##### Monte Carlo Minimum Uncertainty Map #####
#Trying to make an uncertainty map for the minimum using Monte Carlo
# number of MC replicates
rps <- 10000

# initialzing matrices
#All RFs
distsmin <- matrix(0,rps,nrow(extracted2_df_5_RPIg_m))

# RPIg Minimum Monte Carlo
#Change for each different extracted2 variable.
points = extracted2_df_5_RPIg_m
points$util_err <- 5

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIg MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[i])
  res_cv <- points$res_pred_RPIg_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  #dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  #dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  #dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
}
rm(rand, pfm5, i)

#Add the mean and the standard deviation of the minimum to the dataset
extracted2_df_5_RPIg_m$Min = apply(distsmin, 2, mean)
extracted2_df_5_RPIg_m$MinSd = apply(distsmin, 2, sd)

#Make it a spatial dataset
coordinates(extracted2_df_5_RPIg_m) = c('x','y')
proj4string(extracted2_df_5_RPIg_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIg_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIg_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

test = merge(dummy, extracted2_df_5_RPIg_m[,(ncol(extracted2_df_5_RPIg_m)-2):ncol(extracted2_df_5_RPIg_m)], by = 'key')
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')

saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_RPIg_m_test2.tif')
makeMap(rast=Min_test
        ,plotnm='co_RPIg_m_test2.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_RPIg_m_test2.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_RPIg_m_test2.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)


# RPIw Minimum Monte Carlo All Risk Factors
#Change for each different extracted2 variable.
points = extracted2_df_5_RPIw_m
points$util_err <- 5

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[i])
  res_cv <- points$res_pred_RPIw_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  #dist_vars[1,i] <- var(pfm5)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  #dist_vars[2,i] <- var(pfm5)
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  #dist_vars[3,i] <- var(pfm5)
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  #distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
}
rm(rand, pfm5, i)

#Add the mean and the standard deviation of the minimum to the dataset
extracted2_df_5_RPIw_m$Min = apply(distsmin, 2, mean)
extracted2_df_5_RPIw_m$MinSd = apply(distsmin, 2, sd)

#Make it a spatial dataset
coordinates(extracted2_df_5_RPIw_m) = c('x','y')
proj4string(extracted2_df_5_RPIw_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIw_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIw_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

test = merge(dummy, extracted2_df_5_RPIw_m[,(ncol(extracted2_df_5_RPIw_m)-2):ncol(extracted2_df_5_RPIw_m)], by = 'key')
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')

saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_RPIw_m.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_RPIw_m.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_RPIw.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)


# RPIw Minimum Monte Carlo Geologic Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_RPIw_geo_m)

#Geology only
distsmin_g <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[i])
  res_cv <- points$res_pred_RPIw_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  #distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin_g, 2, mean)
sd_test = apply(distsmin_g, 2, sd)
rm(distsmin_g)
extracted2_df_5_RPIw_geo_m$Min = mean_test
extracted2_df_5_RPIw_geo_m$MinSd = sd_test


#Make it a spatial dataset
coordinates(extracted2_df_5_RPIw_geo_m) = c('x','y')
proj4string(extracted2_df_5_RPIw_geo_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIw_geo_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIw_geo_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_RPIw_geo_m[,(ncol(extracted2_df_5_RPIw_geo_m)-2):ncol(extracted2_df_5_RPIw_geo_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_RPIw_geo_m = as.data.frame(extracted2_df_5_RPIw_geo_m[,(ncol(extracted2_df_5_RPIw_geo_m)-2):ncol(extracted2_df_5_RPIw_geo_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIw_geo.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_RPIw_geo.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_RPIw_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)




# RPIg Minimum Monte Carlo Geologic Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_RPIg_geo_m)

#Geology only
distsmin_g <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[i])
  res_cv <- points$res_pred_RPIg_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  #distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin_g, 2, mean)
sd_test = apply(distsmin_g, 2, sd)
rm(distsmin_g)
extracted2_df_5_RPIg_geo_m$Min = mean_test
extracted2_df_5_RPIg_geo_m$MinSd = sd_test
rm(mean_test, sd_test)

#Make it a spatial dataset
coordinates(extracted2_df_5_RPIg_geo_m) = c('x','y')
proj4string(extracted2_df_5_RPIg_geo_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIg_geo_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIg_geo_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_RPIg_geo_m[,(ncol(extracted2_df_5_RPIg_geo_m)-2):ncol(extracted2_df_5_RPIg_geo_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_RPIg_geo_m = as.data.frame(extracted2_df_5_RPIg_geo_m[,(ncol(extracted2_df_5_RPIg_geo_m)-2):ncol(extracted2_df_5_RPIg_geo_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIg_geo.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_RPIg_geo.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_RPIg_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)


# RFC Minimum Monte Carlo Geologic Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_rfc_geo_m)

#Geology only
distsmin_g <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[i])
  res_cv <- points$res_pred_rfc_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  # calcuating overall distribution
  #distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
  
  #Geologic
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin_g, 2, mean)
sd_test = apply(distsmin_g, 2, sd)
rm(distsmin_g)
extracted2_df_5_rfc_geo_m$Min = mean_test
extracted2_df_5_rfc_geo_m$MinSd = sd_test
rm(mean_test, sd_test)

#Make it a spatial dataset
coordinates(extracted2_df_5_rfc_geo_m) = c('x','y')
proj4string(extracted2_df_5_rfc_geo_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_rfc_geo_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_rfc_geo_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_rfc_geo_m[,(ncol(extracted2_df_5_rfc_geo_m)-2):ncol(extracted2_df_5_rfc_geo_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_rfc_geo_m = as.data.frame(extracted2_df_5_rfc_geo_m[,(ncol(extracted2_df_5_rfc_geo_m)-2):ncol(extracted2_df_5_rfc_geo_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_rfc_geo.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_rfc_geo.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_rfc_geo.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)


# RPIw Minimum Monte Carlo All Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_RPIw_m)

#All RFs
distsmin <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIw_max2[i])
  res_cv <- points$res_pred_RPIw_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIw_thresh5[1])] <-5
  pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin, 2, mean)
sd_test = apply(distsmin, 2, sd)
rm(distsmin)
extracted2_df_5_RPIw_m$Min = mean_test
extracted2_df_5_RPIw_m$MinSd = sd_test


#Make it a spatial dataset
coordinates(extracted2_df_5_RPIw_m) = c('x','y')
proj4string(extracted2_df_5_RPIw_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIw_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIw_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_RPIw_m[,(ncol(extracted2_df_5_RPIw_m)-2):ncol(extracted2_df_5_RPIw_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_RPIw_m = as.data.frame(extracted2_df_5_RPIw_m[,(ncol(extracted2_df_5_RPIw_m)-2):ncol(extracted2_df_5_RPIw_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIw.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_RPIw.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_RPIw.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)




# RPIg Minimum Monte Carlo All Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_RPIg_m)

#All RFs
distsmin <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_RPIg_max2[i])
  res_cv <- points$res_pred_RPIg_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_RPIg_thresh5[1])] <-5
  pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin, 2, mean)
sd_test = apply(distsmin, 2, sd)
rm(distsmin)
extracted2_df_5_RPIg_m$Min = mean_test
extracted2_df_5_RPIg_m$MinSd = sd_test


#Make it a spatial dataset
coordinates(extracted2_df_5_RPIg_m) = c('x','y')
proj4string(extracted2_df_5_RPIg_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_RPIg_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_RPIg_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_RPIg_m[,(ncol(extracted2_df_5_RPIg_m)-2):ncol(extracted2_df_5_RPIg_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_RPIg_m = as.data.frame(extracted2_df_5_RPIg_m[,(ncol(extracted2_df_5_RPIg_m)-2):ncol(extracted2_df_5_RPIg_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_RPIg.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_RPIg.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_RPIg.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)


# RFC Minimum Monte Carlo All Risk Factors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_rfc_m)

#All RFs
distsmin <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,5) #5 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # reservoir RPIw MC
  pfm5 <- rep(0,rps)
  
  res_mean <- log(points$res_pred_rfc_max2[i])
  res_cv <- points$res_pred_rfc_max_err[i]
  
  sigma2 <- log(res_cv^2 + 1)
  mu <- res_mean- sigma2/2
  rand <- rnorm(rps,mu,sqrt(sigma2))
  rm(mu)
  
  pfm5[rand < log(res_rfc_thresh5[1])] <-5
  pfm5[rand > log(res_rfc_thresh5[6])] <- 0
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
  pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
  
  #Because the values are reversed, switch them.
  pfm5 <- 5-pfm5
  
  mat_mc[,1] <- pfm5
  
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,2] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,4] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,5] <- pfm5
  
  # calcuating overall distribution
  distsmin[,i] <- apply(cbind(mat_mc[,1],mat_mc[,2],0.5*mat_mc[,3]+0.5*mat_mc[,4],mat_mc[,5]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin, 2, mean)
sd_test = apply(distsmin, 2, sd)
rm(distsmin)
extracted2_df_5_rfc_m$Min = mean_test
extracted2_df_5_rfc_m$MinSd = sd_test


#Make it a spatial dataset
coordinates(extracted2_df_5_rfc_m) = c('x','y')
proj4string(extracted2_df_5_rfc_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_RPIg
#Add the rownames to each database for joining purposes
extracted2_df_5_rfc_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_rfc_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_rfc_m[,(ncol(extracted2_df_5_rfc_m)-2):ncol(extracted2_df_5_rfc_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_rfc_m = as.data.frame(extracted2_df_5_rfc_m[,(ncol(extracted2_df_5_rfc_m)-2):ncol(extracted2_df_5_rfc_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_rfc.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_rfc.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_rfc.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=4
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)


#EGS

# EGS Minimum Monte Carlo 5 Colors
#Change for each different extracted2 variable.
points = as.data.frame(extracted2_df_5_egs_m)

#EGS
distsmin_g <- matrix(0,rps,nrow(points))

# loop for Monte Carlo
for(i in 1:nrow(points)){
  
  # setting seed
  set.seed(10)
  
  # initializing matrix
  mat_mc <- matrix(0,rps,4) #4 is for the number of RFs, considering that seismic has 2 (stress and angle)
  
  # thermal new thresholds MC
  pfm5 <- rep(0,rps)
  
  th_mean <- points$therm_pred[i]
  th_se <- points$therm_err[i]
  
  rand <- rnorm(rps,th_mean,th_se)
  
  pfm5[rand < therm_thresh5[1]] <- 5
  pfm5[rand > therm_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
  pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
  pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
  pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
  pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])
  
  mat_mc[,1] <- pfm5
  
  
  # seismic earthquake MC
  pfm5 <- rep(0,rps)
  
  se_eq_mean <- points$seis_eq_pred[i]
  se_eq_se <- points$seis_eq_err[i]
  
  rand <- rnorm(rps,se_eq_mean,se_eq_se)
  
  pfm5[rand < seis_eq_thresh5[1]] <-5
  pfm5[rand > seis_eq_thresh5[6]] <-0
  pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
  pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
  pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
  pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
  pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
  
  #Because the values are reversed, switch them back here.
  pfm5 <- 5-pfm5
  
  mat_mc[,2] <- pfm5
  
  
  # seismic stress angle MC
  pfm5 <- rep(0,rps)
  
  se_st_mean <- points$NormAngs[i]
  se_st_se <- points$seis_stress_err[i]
  
  rand <- rnorm(rps,se_st_mean,se_st_se)
  
  # New Method:
  # Find the values that are negative and add 360 degrees to them until they are all positive.
  while (length(rand[which(rand < 0)]) != 0){
    rand[which(rand < 0)] = rand[which(rand < 0)] + 360
  }
  # Now convert all values to [0,180]
  while (length(rand[which(rand > 180)]) != 0){
    rand[which(rand > 180)] = rand[which(rand > 180)] - 180
  }
  # Now take all angles and convert them to a risk angle
  a1 = abs(rand - critical_ang1)
  a2 = abs(rand - critical_ang2)
  # Bind them
  b = rbind(a1, a2)
  # Assign the minimum value (most risky) to those points
  rand = apply(b, 2, min)
  rm(b,a1,a2)
  
  pfm5[rand < seis_stress_thresh5[1]] <- 5
  pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
  pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
  pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
  pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
  pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
  
  #Because the values are reversed, switch them back
  pfm5 <- 5-pfm5
  
  mat_mc[,3] <- pfm5
  
  #Utilization:
  # generating random values
  rand <- rnorm(rps,points$util_pred[i],points$util_pred[i]*5/100) #5/100 is the 5% uncertainty that is assumed for utilization.
  
  # play fairway 5
  pfm5 <- rep(0,rps)
  pfm5[rand < util_thresh5[1]] <- 5
  pfm5[rand > util_thresh5[6]] <- 0
  pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
  pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
  pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
  pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
  pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
  
  mat_mc[,4] <- pfm5
  
  #Geologic
  distsmin_g[,i] <- apply(cbind(mat_mc[,1],0.5*mat_mc[,2]+0.5*mat_mc[,3],mat_mc[,4]),1,min)
}
rm(rand, pfm5, i, th_mean, th_se, sigma2, res_mean, res_cv, se_eq_mean, se_eq_se, se_st_se, se_st_mean, mat_mc, points)

#Add the mean and the standard deviation of the minimum to the dataset
mean_test = apply(distsmin_g, 2, mean)
sd_test = apply(distsmin_g, 2, sd)
rm(distsmin_g)
extracted2_df_5_egs_m$Min = mean_test
extracted2_df_5_egs_m$MinSd = sd_test
rm(mean_test, sd_test)

#Make it a spatial dataset
coordinates(extracted2_df_5_egs_m) = c('x','y')
proj4string(extracted2_df_5_egs_m) = CRS("+init=epsg:26917")

#Add data back to the full region grid, which is stored in extracted_df_5_RPIg for all variables.
dummy = extracted_df_5_egs
#Add the rownames to each database for joining purposes
extracted2_df_5_egs_m$key = as.numeric(rownames(as.data.frame(extracted2_df_5_egs_m)))
dummy$key = as.numeric(rownames(as.data.frame(dummy)))

#Merge the data based on the key. Only bring over the Min, MinSd, and key columns.
test = merge(dummy, extracted2_df_5_egs_m[,(ncol(extracted2_df_5_egs_m)-2):ncol(extracted2_df_5_egs_m)], by = 'key')
rm(dummy)

#Save only the results in case it needs to be remade.
extracted2_df_5_egs_m = as.data.frame(extracted2_df_5_egs_m[,(ncol(extracted2_df_5_egs_m)-2):ncol(extracted2_df_5_egs_m)])

#Make a raster
Min_test = raster(test, layer = 'Min')
MinSd_test = raster(test, layer = 'MinSd')
rm(test)

#Save the raster
saveRast(rast=Min_test
         ,wd=wd_raster_out
         ,rastnm='co_5_0_5_m_egs2.tif')
makeMap(rast=Min_test
        ,plotnm='co_5_0_5_m_egs2.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=FALSE)

saveRast(rast=MinSd_test
         ,wd=wd_raster_out
         ,rastnm='co_pfa_sd5_min_egs.tif')
makeMap(rast=MinSd_test
        ,plotnm='co_pfa_sd5_min_egs.png'
        ,wd=wd_image
        ,numCol=5
        ,comTy='min'
        ,numRF=3
        ,sdMap=TRUE)

#remove files for space
rm(Min_test, MinSd_test)

##### Histograms of Metrics #####
#Making histograms for the metrics
# Note that the individual RFs are approximmate because they do not use the
# tables for interpolating the scaled means.
setwd(wd_image)
makeHist(rast=therm_pred
         ,thresh3=c()
         ,thresh5=therm_thresh5
         ,rev_sc=TRUE
         ,plotnm='th_hist.png'
         ,yloc=-1.2*10^-4
         ,yshift=1.5*10^-5
         ,title='Depth to 80 Deg C (m)')

makeHist(rast=calc(res_pred_RPIg_max2,fun=log10)
         ,thresh3=c()
         ,thresh5=log10(res_RPIg_thresh5)
         ,rev_sc=FALSE
         ,plotnm='re_RPIg_hist.png'
         ,yloc=-0.13
         ,yshift=0.02
         ,title='log10 of RPIg (L/MPa-s)')

makeHist(rast=calc(res_pred_RPIw_max2,fun=log10)
         ,thresh3=c()
         ,thresh5=log10(res_RPIw_thresh5)
         ,rev_sc=FALSE
         ,plotnm='re_RPIw_hist.png'
         ,yloc=-0.15
         ,yshift=0.02
         ,title='log10 of RPIw (L/MPa-s)')

makeHist(rast=calc(res_pred_rfc_max2,fun=log10)
         ,thresh3=c()
         ,thresh5=log10(res_rfc_thresh5)
         ,rev_sc=FALSE
         ,plotnm='re_rfc_hist.png'
         ,yloc=-0.13
         ,yshift=0.02
         ,title='log10 of RFC (L/MPa-s)')

makeHist(rast=util_pred[(util_pred<100)]
         ,thresh3=c()
         ,thresh5=util_thresh5
         ,rev_sc=TRUE
         ,plotnm='ut_hist.png'
         ,yloc=-0.025
         ,yshift=0.003
         ,title='Utilization Surface Cost ($/MMBTU) \n (only < 100 $/MMBTU plotted)')

makeHist(rast=seis_eq_pred[(seis_eq_pred<10^5)]
         ,thresh3=c()
         ,thresh5=seis_eq_thresh5
         ,rev_sc=FALSE
         ,plotnm='seEq_hist.png'
         ,yloc=-1.5*10^-5
         ,yshift=3*10^-6
         ,title='Seismic Risk for Proximity to Earthquake (m)')

makeHist(rast=seis_stress_pred[(seis_stress_pred<100)]
         ,thresh3=c()
         ,thresh5=seis_stress_thresh5
         ,rev_sc=FALSE
         ,plotnm='seSt_hist.png'
         ,yloc=-0.017
         ,yshift=0.003
         ,title=expression(paste('Seismic Risk for Angle in Stress (',degree,')')))


#Combined Risk Factors
setwd(wd_image)
makeHist(rast=co_5_0_20_s_rfc
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_rfc_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of All Risk Factors with RFC for Reservoirs')

makeHist(rast=co_5_0_20_s_RPIg
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_RPIg_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of All Risk Factors with RPIg for Reservoirs')

makeHist(rast=co_5_0_20_s_RPIw
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_RPIw_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of All Risk Factors with RPIw for Reservoirs')

makeHist(rast=co_5_0_625_p_rfc
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_rfc_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Geometric Mean of All Risk Factors \n with RFC for Reservoirs')

makeHist(rast=co_5_0_625_p_RPIw
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_RPIw_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Geometric Mean of All Risk Factors \n with RPIw for Reservoirs')

makeHist(rast=co_5_0_625_p_RPIg
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_RPIg_hist.png'
         ,yloc=-0.1
         ,yshift=0
         ,title='Geometric Mean of All Risk Factors \n with RPIg for Reservoirs')

makeHist(rast=co_5_0_5_m_rfc
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_rfc_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='Minimum of All Risk Factors \n with RFC for Reservoirs')

makeHist(rast=co_5_0_5_m_RPIw
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_RPIw_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='Minimum of All Risk Factors \n with RPIw for Reservoirs')

makeHist(rast=co_5_0_5_m_RPIg
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_RPIg_hist.png'
         ,yloc=-0.8
         ,yshift=0
         ,title='Minimum of All Risk Factors \n with RPIg for Reservoirs')

#Geologic Only
makeHist(rast=co_5_0_16_s_rfc_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_rfc_geo_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of Geologic Risk Factors \n with RFC for Reservoirs')

makeHist(rast=co_5_0_16_s_RPIg_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_RPIg_geo_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of Geologic Risk Factors \n with RPIg for Reservoirs')

makeHist(rast=co_5_0_16_s_RPIw_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_RPIw_geo_hist.png'
         ,yloc=-0.05
         ,yshift=0
         ,title='Average of Geologic Risk Factors \n with RPIw for Reservoirs')

makeHist(rast=co_5_0_125_p_rfc_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_rfc_geo_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Geometric Mean of Geologic Risk Factors \n with RFC for Reservoirs')

makeHist(rast=co_5_0_125_p_RPIw_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_RPIw_geo_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Geometric Mean of Geologic Risk Factors \n with RPIw for Reservoirs')

makeHist(rast=co_5_0_125_p_RPIg_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_RPIg_geo_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Geometric Mean of Geologic Risk Factors \n with RPIg for Reservoirs')

makeHist(rast=co_5_0_5_m_rfc_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_rfc_geo_hist.png'
         ,yloc=-0.3
         ,yshift=0
         ,title='Minimum of Geologic Risk Factors \n with RFC for Reservoirs')

makeHist(rast=co_5_0_5_m_RPIw_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_RPIw_geo_hist.png'
         ,yloc=-0.3
         ,yshift=0
         ,title='Minimum of Geologic Risk Factors \n with RPIw for Reservoirs')

makeHist(rast=co_5_0_5_m_RPIg_geo
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_RPIg_geo_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Minimum of Geologic Risk Factors \n with RPIg for Reservoirs')


#EGS
makeHist(rast=co_5_0_5_m_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_m_egs_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Minimum of Thermal, Seismic, and Utilization Risk Factors')

makeHist(rast=co_5_0_125_p_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_g_egs_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Geometric Mean of Thermal, Seismic, \n and Utilization Risk Factors')

makeHist(rast=co_5_0_16_s_egs
         ,thresh3=c()
         ,thresh5=c(0,1,2,3,4,5)
         ,rev_sc=FALSE
         ,plotnm='co_5_0_5_a_egs_hist.png'
         ,yloc=-0.15
         ,yshift=0
         ,title='Average of Thermal, Seismic, and Utilization Risk Factors')

##### saving workspace #####
setwd(wd_workspace)
filename <- paste('geothermal_pfa_analysis_FinalRun_',Sys.Date(),'.RData',sep='')
save.image(file=filename)
