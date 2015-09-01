# function for converting raster into play fairway ranking
# inputs are:
#     rast        = raster object with values
#     thresholds  = play fairway thresholds with [min, threshold1,...,threshold2/4, max]
#     ignore      = numerical value to ignore
#     log_scale   = whether conversion should be done on a log-scale (thresholds and values on real-scale)
#     rev_scale   = whether the scale should be reversed (low is good, high is bad)

convRastPFRank <- function(rast               # raster
                           ,thresholds        # thresholds including min and max values
                           ,ignore = -9999    # values to ignore
                           ,log_scale=FALSE   # whether to convert using log scale
                           ,rev_scale=FALSE   # whether scale is reversed
                           ){
  
  # initializing raster and thresholds
  if(log_scale == TRUE){ # case for when log-scale
    
    # initialize raster/thresholds as log10 of input raster/thresholds
    rast_pfa <- calc(rast,fun=function(x){log10(x)}) 
    thresh <- log10(thresholds)
  
  }else if(log_scale == FALSE){ # case for not log-scale
    
    # initialize raster/thresholds as input raster/thresholds
    rast_pfa <- rast
    thresh <- thresholds
  
  }else{ # error message
    
    print('Not a valid selection for log_scale')
  }
  
  # setting number scheme based on number of thresholds
  # thresh includes min, max, and category separation values
  num_cat <- seq(0,length(thresh)-1,1)
  
  # assigning 0 to values < minimum criteria
  # only assigning values with raster cell is not ignored
  rast_pfa[(rast != ignore)&(rast < thresh[1])] <- num_cat[1]
  
  # assigning highest to values >= maximum criteria
  # only assigning values with raster cell is not ignored
  rast_pfa[(rast != ignore)&(rast >= thresh[length(thresh)])] <- num_cat[length(num_cat)]
  
  # looping over intermediate values (between min and max)
  for(i in 1:length(num_cat)){
    
    # values between the two thresholds are linearly scaled
    rast_pfa[(rast >= thresh[i])&(rast < thresh[i+1])] <- (rast[(rast >= thresh[i])&(rast < thresh[i+1])] - thresh[i])/(thresh[i+1]-thresh[i]) + i - 1
  }
  
  # returning calculated raster
  if(rev_scale == TRUE){ # case for reversing the scale
    
    # reversing scale by subtracting raster from the maximum value
    return(calc(rast_pfa,fun=fun(x){num_cat[length(num_cat)]-x}))
    
  }else if(rev_scale == FALSE){ # case for scale not being reversed
    
    # no adjustment of the calculated raster
    return(rast_pfa)
    
  }else{ # error message
    
    print('Not a valid selection for rev_scale')
  }
}