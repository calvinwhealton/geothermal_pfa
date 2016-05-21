#########################
# script to combine data for US Census places



# libraries
library(sp)           # for transforming coordinates
library(raster)       # for raster calcuations
library(rgdal)        # reading in data, readOGR() and writeOGR()
library(spdep)        # other spatial processing
library(rgeos)        # for buffering places of interest and finding intersections
library(maptools)

# setting working directory
setwd('/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper')

# loading in US Census Place Data
places <- readOGR(dsn='/Users/calvinwhealton/Documents/GIS Projects/COSUNA_Paper/cb_2015_us_state_500k_copy_clipAppSt'
                  ,layer='clipped_appBas_st')

# loading in US Economic Census Data
economic <- read.table("EC1200A1.dat"
                       ,header=TRUE
                       ,sep='|')
