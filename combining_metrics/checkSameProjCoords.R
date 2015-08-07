# function for checking if two rasters 
# have the same dimensions & coordinates

checkSameProjCoords <-function(rast1    # raster 1
                               ,rast2   # raster 2
                               ){
  
  # initializing vector to hold TRUE/FALSE
  tf_vec <- rep(NA,7)
  
  # checking if projection is the same
  tf_vec[1] <- (rast1@crs@projargs == rast2@crs@projargs)
  
  # checking if number of rows and columns is the same
  tf_vec[2] <- (rast1@nrows == rast2@nrows)
  tf_vec[3] <- (rast1@ncols == rast2@ncols)
  
  # checking if min and max values of coordinates are the same
  tf_vec[4] <- (rast1@extent@xmin == rast2@extent@xmin)
  tf_vec[5] <- (rast1@extent@ymin == rast2@extent@ymin)
  tf_vec[6] <- (rast1@extent@xmax == rast2@extent@xmax)
  tf_vec[7] <- (rast1@extent@ymax == rast2@extent@ymax)
  
  # printing warning messages
  if(tf_vec[1] == FALSE){
    print('Projections do not match.')
  }
  
  if(FALSE %in% tf_vec[c(2,3)]){
    print('Number of rows and/or columns do not match.')
  }
  
  if(FALSE %in% tf_vec[c(4,5,6,7)]){
    print('Minimum and/or maximum do not match.')
  }
  
  # returning value
  return(tf_vec)
}