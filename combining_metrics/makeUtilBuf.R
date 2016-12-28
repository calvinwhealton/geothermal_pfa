# function to make weighting matrix for utilization

makeUtilBufWeight <- function(dist # distance used in pumping
                              ){
  
  if((dist %% 2) == 0){
    print('The distance must be an odd number')    
  }else{
    # initializing matrix of ones
    weightMat <- matrix(1,2*dist+1,2*dist+1)
   
    # looping over values to find cells where center is within
    # the specified distance
    for(i in 1:nrow(weightMat)){
      for(j in 1:ncol(weightMat)){
          
        # center of weighting matrix has coordinates (dist+1,dist+1)
        cell_dist <- sqrt((i - dist - 1)^2 + (j - dist - 1)^2)
        
        # condition for setting value equal to zero
        # when center distance of cells is greater than
        # the specified distance
        if(cell_dist > dist){
          weightMat[i,j] <- NA
        }
        
      }
    }
  
    # returning the weighting matrix
    return(weightMat)
  
  }
}