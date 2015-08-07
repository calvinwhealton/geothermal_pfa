# function for converting raster into play fairway ranking
# inputs are:
#     rast        = raster object with values
#     thresholds  = play fairway thresholds with [min, threshold1,...,threshold2/4, max]
#     ignore      = numerical value to ignore

convRastPFRank <- function(rast         # raster
                           ,thresholds  # thresholds including min and max values
                           ,ignore = -9999
                           ){
  
  # initializing raster as the original value
  rast_pfa <- rast
  
  # setting number scheme based on number of thresholds
  num_cat <- seq(0,length(thresholds)-1,1)
  
  # assigning 0 to values < minimum criteria
  # only assigning values with raster cell is not ignored
  rast_pfa[intersect((rast != ignore),(rast < thresholds[1]))] <- num_cat[1]
  
  # assigning highest to values >= maximum criteria
  # only assigning values with raster cell is not ignored
  rast_pfa[intersect((rast != ignore),(rast >= thresholds[length(thresholds)]))] <- num_cat[length(num_cat)]
  
  # looping over intermediate values (between min and max)
  for(i in 1:(length(thresholds)-1)){
    
    rast_pfa[(rast >= thresholds[i])&(rast < thresholds[i+1])] <- (rast[(rast >= thresholds[i])&(rast < thresholds[i+1])] - thresholds[i])/(thresholds[i+1]-thresholds[i]) + i - 1
    print(i)
  }
}