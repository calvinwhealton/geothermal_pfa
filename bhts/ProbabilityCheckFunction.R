#Function for identifying the drilling fluid weights for nearest neighbor wells to use in the Allegheny Plateau BHT correction

ProbabilityCheck = function(Data,             #Full database 
                            DrillFluidData,   #Drilling fluid data. Not necessarily the same as the database being used.
                            reg_field_dfd,    #Field name for BHT region in DrillFluidData
                            reg_field_d,      #Field name for BHT region in Data
                            ROI_num,          #Number in reg_field that corresponds to the BHT region of interest (ROI). Same for Data and Drilling Fluid Data.
                            df_field,         #Field name for drilling fluid field in Data and DrillFluidData
                            pt_eval = 25,     #Number of nearest neighbor points used to calculate the percentage of air or mud drilled wells
                            rad_max = 50000)  #Maximum distance to search for wells (m)
{   
  
  #sp package needed for the distance calculation.
  require(sp)
  
  #printing error if any necessary variable is not defined properly
  if (length(DrillFluidData@data[df_field]) == 0){
    print("'df_field' variable not in DrillingFluidData")
  }
  else if (sum(is.na(DrillFluidData@data[df_field])+1-1) != 0){ # +1-1 transforms TRUE/FALSE into 1 or 0 so that it can be summed.
    print("Test not performed because of NAs in 'df_field'.")
    print("Remove NAs from DrillFluidData and call function again.")
  }
  else if (length(Data@data[df_field]) == 0){
    print("'df_field' variable not in Data")
  }
  else if(length(Data@data[reg_field_d]) == 0){
    print("No 'reg_field_d' variable in Data")
  }
  else if(length(DrillFluidData@data[reg_field_dfd])==0){
    print("No 'reg_field_dfd' variable in DrillFluidData")
  }
  else if (Data@proj4string@projargs != "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
    print("Data not in correct coordinate system, convert to: +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  }
  else if (DrillFluidData@proj4string@projargs != "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
    print("Data not in correct coordinate system, convert to: +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  }
  else { #Coordinates defined properly, drilling fluid field exists, BHT region exists. Proceed with calculation.
    
    #Add columns to the dataframe
    Data@data$max_pt_rad = 0
    Data@data$max_pt_error = 0
    Data@data$PctAir = 0
    Data@data$PctMud = 0
    
    
    #Check for the nearest neighbors
    for(i in 1:nrow(Data)){ #loop over all observations in Data
      
      if (Data@data[reg_field_d][i,] == ROI_num) { #These are wells in the region of interest. They may need drilling fluid assigned.
        
        if (is.na(Data@data[df_field][i,]) == FALSE) {
          #This well has a drilling fluid. Use it. Assign Probability accordingly
          if (Data@data[df_field][i,] == 'all_agfs') { #The well was drilled using air
            Data@data$PctAir[i] = 1.0
            Data@data$PctMud[i] = 0.0
          }
          else { #The well was drilled using mud
            Data@data$PctAir[i] = 0.0
            Data@data$PctMud[i] = 1.0
          } 
        }
        else {
          #Initializing matrix with one column for the index and the other for distance to points
          dist_vec <- matrix(0,nrow(DrillFluidData[-which(DrillFluidData@data[reg_field_dfd] == ROI_num),]),1) 
          
          #Calculating Euclidian distance for all of the points in the region of interest
          dist_vec[,1] <- spDistsN1(DrillFluidData[-which(DrillFluidData@data[reg_field_dfd] == ROI_num),], Data[i,], longlat=FALSE)
          
          #Finding the distance to the point that is the (pt_eval)-th farthest from the point of interest.
          dist_cutoff <- sort(dist_vec)[pt_eval]
          
          #Storing the distance of the nearest pt_eval points
          Data@data$max_pt_rad[i] <- dist_cutoff
          
          #Evaluating if radius is too large to capture the points locally in rad_max
          if(Data@data$max_pt_rad[i] > rad_max){ # case for pt_eval points not within rad_max
            
            #Assigning code for not enough points within rad_max to evaluate
            Data@data$max_pt_error[i] <- 1
            
            #Assign the percentage of the regional averages of drilling fluid because there's not enough local data.
            Data@data$PctAir[i] = length(which(DrillFluidData@data[reg_field_dfd][which(DrillFluidData@data[df_field] == 'all_agfs'),] == ROI_num))/length(which(DrillFluidData@data[reg_field_dfd] == ROI_num))
            Data@data$PctMud[i] = 1.0 - Data@data$PctAir[i]
            
          }
          else{ # case for at least pt_eval points within rad_max
            
            #Use the neighboring points to determine the BHT correction weights
            # finding the indices of the closest points
            inds_within_cutoff <- which(dist_vec <= dist_cutoff)
            
            # finding percentage of air drilled wells from reference database
            PercentAir <- length(which(DrillFluidData@data[df_field][inds_within_cutoff,] == 'all_agfs'))/length(inds_within_cutoff)
            PercentMud <- 1.0 - PercentAir
            
            # storing the percentage of air and gas
            Data@data$PctAir[i] = PercentAir
            Data@data$PctMud[i] = PercentMud
          }
        }
      }
      else {
        #These are not Allegheny Plateau, assign NA to the added columns
        Data@data$PctAir[i] = NA
        Data@data$PctMud[i] = NA
      }
    } # end for loop
    
    return(Data) 
  }  
}
