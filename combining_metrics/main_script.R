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
library(raster)

# user-defined functions
setwd('/Users/calvinwhealton/GitHub/geothermal/combining_metrics')
source('checkSameProjCoords.R')
source('convertRasterPFAMetric.R')
source('combineRFs.R')

##### importing rasters #####
# thermal prediction and error
therm_pred <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Thermal/heatflowp-')
therm_err <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Thermal/heatflowe-')

# reservoir prediction and error
res_pred <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Reservoirs/fff_1015p')
res_err  <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Reservoirs/fff_1015e')

# utilization prediciton and error
util_pred_ny <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Utilization/nyjoinedp-')
util_pred_pa <- raster('/Users/calvinwhealton/Documents/Education/Cornell/Play Fairway/Final_Rasters/Utilization/pajoinedp-')

util_pred_ny[(util_pred_ny %in% -9999)] <- NA
util_pred_pa[(util_pred_pa %in% -9999)] <- NA

util_pred <- merge(util_pred_ny,util_pred_pa)
plot(util_pred)
##### checking rasters are in the same system #####
check_therm_res <- checkSameProjCoords(therm_pred,res_pred)

##### converting into play fairway scheme #####

# setting thresholds [min, intermediates, max]
therm_thresh <- c(30,35,40,60)  # thermal
res_thresh <- 10^seq(-3,0,1)    # reservoir

# converting into play fairway scheme
therm_pfa <- convRastPFRank(therm_pred,therm_thresh,ignore=-9999)
writeRaster(therm_pfa,'therm_pfa.grd')

##### making a quick and simple map to compare results #####
# maps should look the same if the thresholds were properly specified

# info on manually manipulating the color bar at
# http://gis.stackexchange.com/questions/17339/raster-legend-in-r-how-to-colour-specific-values


# play fairway map
breaks <- c(-9999,-0.01,1,2,5) # -0.01 used because using 0 makes 0 map to white
cols <- c('white','red','yellow','green')
plot(therm_pfa,breaks=breaks,col=cols,legend=FALSE)

# original raster map
breaks2 <- c(-9999,-0.01,therm_thresh,10^10) # from play fairway threhsolds with highest value higher than max
cols2 <- c('white','red','red','yellow','green','green') # red used twice for below minimums and minimum to first threshold (both 0)
plot(therm_pred,breaks=breaks2,col=cols2,legend=FALSE)

##### making a quick and simple map to compare results #####
# combining maps


