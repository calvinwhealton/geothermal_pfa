#Function for identifying the drilling fluid weights for nearest neighbor wells to use in the Allegheny Plateau BHT correction

ProbabilityCheck = function(Data,             #Full database as SpatialPointsDataFrame
                            DrillFluidData,   #Drilling fluid data. Not necessarily the same as the database being used.
                            pt_eval = 25,     #Number of nearest neighbor points used to calculate the percentage of air or mud drilled wells
                            rad_max = 50000)  #Maximum distance to search for wells (m)
{   
  
  #sp package needed for the distance calculation.
  require(sp)
  
  # printing error if any necessary variable is not defined properly
  if (length(DrillFluidData@data$gen_name3)==0){
    print("No 'gen_name3' variable in DrillingFluidData")
  }
  else if (sum(is.na(DrillFluidData@data$gen_name3)+1-1)!=0){ # +1-1 transforms TRUE/FALSE into 1 or 0. condition for any NA present
    print("Test not performed because of NAs in 'gen_name3'.")
    print("Remove NAs from DrillFluidData and call function again.")
  }
  else if (length(Data@data$gen_name3)==0){
    print("No 'gen_name3' variable in Data")
  }
  else if(length(Data@data$reg)==0){
    print("No BHT region variable in Data")
  }
  else if (Data@proj4string@projargs != "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
    print("Data not in correct coordinate system, convert to: +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  }
  else if (DrillFluidData@proj4string@projargs != "+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
    print("Data not in correct coordinate system, convert to: +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  }
  else{ #Coordinates defined properly. Drilling fluid is present. Proceed with calculation.
    
    #Add columns to the dataframe
    Data@data$max_pt_rad = 0
    Data@data$max_pt_error = 0
    Data@data$PctAir = 0
    Data@data$PctMud = 0
    
    
    #Check for the nearest neighbors
    for(i in 1:nrow(Data)){ #loop over all observations in Data that need drilling fluid assigned
      
      if (Data@data$reg[i] == 1) { #These are Allegheny Plateau Wells
        # initializing matrix with one column for index and the other for distance
        dist_vec <- matrix(0,nrow(DrillFluidData),1) 
        
        # calculating weighted Euclidian distance for all of the points
        dist_vec[,1] <- spDistsN1(DrillFluidData, Data[i,], longlat=FALSE)
        
        # finding the empirical quantile for removing, pt_eval is an input
        dist_cutoff <- sort(dist_vec)[pt_eval]
        
        # storing the maximum distance of the nearest pt_eval points
        Data@data$max_pt_rad[i] <- dist_cutoff
        
        # evaluating if radius is too large to capture the points
        if(Data@data$max_pt_rad[i] > rad_max){ # case for points criterion not satisfied within rad_max
          
          # assigning error for not enough points to evaluate
          Data@data$max_pt_error[i] <- 1
          
        }
        else{ # case for points criterion satisfied within rad_max
          
          if (is.na(Data@data$gen_name3[i]) == FALSE) {
            #This well has a drilling fluid. Use it. Assign Probability accordingly
            if (Data@data$gen_name3[i] == 'all_agfs') { #The well was drilled using air
              Data@data$PctAir[i] = 1
              Data@data$PctMud[i] = 0
            }
            else { #The well was drilled using mud
              Data@data$PctAir[i] = 0
              Data@data$PctMud[i] = 1
            } 
          }
          else { #The well doesn't have a drilling fluid, use the neighboring points to determine the BHT correction weights
            
            # finding the indices of the closest points
            inds_within_cutoff <- which(dist_vec <= dist_cutoff)
            
            # finding percentage of air drilled wells from reference database
            PercentAir <- length(which(DrillFluidData@data$gen_name3[inds_within_cutoff] == 'all_agfs'))/length(inds_within_cutoff)
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
