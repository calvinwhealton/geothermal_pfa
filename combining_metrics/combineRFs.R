# functions for combining the risk maps
# and creating the play fairway map

# function for combining three risk factors
comb3rf <- function(stack_rast      # stacked raster
                    ,method = 'sum' # method used to calculate the combined metric
                    ,ignore = -9999 # cells to be ignored
                    ,min_eval=0     # minimum value used to evaluate combined, if less than minimum value will return minimum of cell
                    ){
  
  # combined PFA metric by summing all maps
  if(method == 'sum'){
    rast_pfa <- calc(stack_rast,sum)
  }
  
  # combined PFA metric by taking product of values
  else if (method == 'product'){
    rast_pfa <- calc(stack_rast,prod)
  }
  
  # combined PFA metric by  minimum of values
  else if (method == 'minimum'){
    rast_pfa <- calc(stack_rast,min)
  }
  
  # if values are below the specified minimum, the minimum cell value will be used
  rast_pfa[(calc(stack_rast,min) < min_eval)] <- calc(stack_rast,min)[(calc(stack_rast,min) < min_eval)]
  
  # setting values for ignored cells to ignored value
  for(i in 1:length(raster@layers)){
    rast_pfa[(stack_rast[[i]] %in% ignore)] <- ignore
  }
    
  # returning raster
  return(rast_pfa)
}
