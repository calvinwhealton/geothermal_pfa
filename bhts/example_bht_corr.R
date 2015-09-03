# code to run function "BHT_NY_PA_code2.R" for an example dataset
# example dataset designed to include bad observations possible
# all error codes should be produced in the output data frame

# set working directory, will need to change depending on user
setwd("/Users/calvinwhealton/GitHub/geothermal/bhts")

# loading function from same directory
source("func_BHT_NY_PA_WV_corr.R")

# test dataset imported from .csv file
test_data <- read.table('bht_test_data.csv', header=TRUE, sep=',')

# running function
corrected_BHTS <- NY_PA_WV_BHT2(test_data)

corrected_BHTS$bht_corr <- corrected_BHTS$corr_bht_c - corrected_BHTS$bht_c


plot(corrected_BHTS$calc_depth_m[corrected_BHTS$calc_depth_m<6000]
     , corrected_BHTS$corr_bht_c[corrected_BHTS$calc_depth_m<6000]-corrected_BHTS$bht_c[corrected_BHTS$calc_depth_m<6000]
     , pch = corrected_BHTS$reg[corrected_BHTS$calc_depth_m<6000]+15)
grid(lwd=3)