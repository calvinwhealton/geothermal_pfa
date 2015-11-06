# function for converting raster into play fairway ranking
# inputs are:
#     rast        = raster object with values
#     thresholds  = play fairway thresholds with [min, threshold1,...,threshold2/4, max]
#     ignore      = numerical value to ignore
#     rev_scale   = whether the scale should be reversed (low is good, high is bad)

convRastPFRank <- function(rast               # raster
                           ,thresholds        # thresholds including min and max values
                           ,ignore = -9999    # values to ignore
                           ,rev_scale=FALSE   # whether scale is reversed
                           ,log_scale=FALSE   # whether scale is logarithmic
                           ){
  
  # initialize raster/thresholds as input raster/thresholds
  rast_pfa <- rast
  thresh <- thresholds
  
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
    
    if(log_scale==FALSE){
      
      # values between the two thresholds are linearly scaled
      rast_pfa[(rast >= thresh[i])&(rast < thresh[i+1])] <- (rast[(rast >= thresh[i])&(rast < thresh[i+1])] - thresh[i])/(thresh[i+1]-thresh[i]) + i - 1
      
    }else{
      
      # case for logarithmic scalse
      rast_pfa[(rast >= thresh[i])&(rast < thresh[i+1])] <- (log10(rast[(rast >= thresh[i])&(rast < thresh[i+1])]) - log10(thresh[i]))/(log10(thresh[i+1])-log10(thresh[i])) + i - 1
    }
  }
  
  # returning calculated raster
  if(rev_scale == TRUE){ # case for reversing the scale
    
    # reversing scale by subtracting raster from the maximum value
    return(calc(rast_pfa,fun=function(x){num_cat[length(num_cat)]-x}))
    
  }else if(rev_scale == FALSE){ # case for scale not being reversed
    
    # no adjustment of the calculated raster
    return(rast_pfa)
    
  }else{ # error message
    
    print('Not a valid selection for rev_scale')
  }
}