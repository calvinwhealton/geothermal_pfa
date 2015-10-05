# script to load libraries and user defined functions

# libraries
library(sp) # for spatial points and converting coordinate systems
library(raster) # raster calcuations
library(rgdal) # reading in rasters
library(xlsx) # for reading in xls and xlsx files
library(pracma) # for 2-d interpolation
library(RColorBrewer) # for colors
library(rgeos) # for buffering places of interest

# user defined functions
setwd('/Users/calvinwhealton/GitHub/geothermal/combining_metrics')
source('convertRasterPFAMetric.R') # converting into PF scheme
source('makeUtilBuf.R') # make utilization buffer weighting matrix
source('plotWeightBuf.R') # plot utilization weighting matrix 
source('makeHist.R') # make histogram of distribution and thresholds
source('makeMap.R') # make map of raster
source('saveRast.R') # save raster
