# function to save raster

saveRast <- function(rast # raster
                     ,wd  # working directory
                     ,rastnm # name for raster
                     ){
  
  # initializing exported raster as input raster
  rast2 <- rast
  
  # substituting -9999 for NAs
  rast2[(is.na(rast)==TRUE)] <- -9999
  
  # substituting -9999 for -Inf
  rast2[(rast==-Inf)] <- -9999
  
  # setting working directory
  setwd(wd)
  
  # converting into spatial gird
  X <- as(rast2 ,'SpatialGridDataFrame')
  
  # writing file
  writeGDAL(X,rastnm,drivername='GTiff')
}