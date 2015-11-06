#This function is used to set up the dataframe for outlier analysis, and then perform the outlier analysis in outlier_iden. 
#All options for outlier identification Are available in this function.
#There is an option to set a threshold, below which all values will be dropped before the outlier identification is run. The default is to drop all negative values.

WhichOne = function(Data,            # Database to be tested for outliers
                    ThermalVar,      # Desired column name for the tested variable in the output file
                    NameInData,      # Field name of the thermal variable in Data. May be the same as ThermalVar.
                    X_coord,         # Field name of the UTM longitude in Data.
                    Y_coord,         # Field name of the UTM latitude in Data. 
                    Threshold = 0.0  # Points in field NameInData below Threshold will be dropped before outlier analysis. Default is 0.0 (negative values)
                    algo = 1,        # local detection algorithm
                    outcri = 1,      # outlier criterion
                    pt_eval = 25,    # minumum number of points to evaluate
                    rad_eval = 16,   # radius used for local outlier identification (km)
                    box_size = 32,   # edge length of the box (km)
                    pt_min = 25,     # minimum number of points to evaluate local outliers
                    rad_max = 16,    # maximum radius to search for local points (km)
                    k_glob = 3,      # constant for quartile-median range in global analysis
                    k_loc = 3,       # constant for quartile-median range in local analysis
                    type = 7,        # type of quantile interpolation
                    ){
  
  #Rename columns based on the necessary inputs to the outlier identification function.
  colnames(Data@data)[which(colnames(Data@data)==NameInData)] = "test"
  colnames(Data@data)[which(colnames(Data@data)==X_coord)] = "x_coord"
  colnames(Data@data)[which(colnames(Data@data)==Y_coord)] = "y_coord"
  
  #Remove all negative values. They are physically impossible and will affect the outlier identification.
  if (length(which(Data@data[, "test"] <= Threshold)) > 0) {
    Data = Data[-which(Data@data["test"] <= Threshold),] 
  }
  
  #Run outlier identification
  Outs = outlier_iden(Data
                      , algo = algo          # local detection algorithm
                      , outcri = outcri      # outlier criterion
                      , pt_eval = pt_eval    # minumum number of points to evaluate
                      , rad_eval = rad_eval  # radius used for local outlier identification (km)
                      , box_size = box_size  # edge length of the box (km)
                      , pt_min = pt_min      # minimum number of points to evaluate local outliers
                      , rad_max = rad_max    # maximum radius to search for local points (km)
                      , k_glob = k_glob      # constant for quartile-median range in global analysis
                      , k_loc = k_loc        # constant for quartile-median range in local analysis
                      , type = type)         # type of quantile interpolation
  
  #Change the column name of the UTM coordinates back to the original name.
  colnames(Outs@data)[which(colnames(Outs@data)=="x_coord")] = X_coord
  colnames(Outs@data)[which(colnames(Outs@data)=="y_coord")] = Y_coord
  
  #Rename the thermal variable column according to the ThermalVar
  colnames(Outs@data)[which(colnames(Outs@data)=="test")] = ThermalVar
  
  #Split the data into Outliers and NotOutliers
  NotOutliers = Outs[-which(Outs@data$outs != 0),]
  Outliers = Outs[which(Outs@data$outs != 0),]
  
  #List for returning data
  TestedOutliers = list("NotOutliers" = NotOutliers, "Outliers" = Outliers)
  
  return(TestedOutliers)
}