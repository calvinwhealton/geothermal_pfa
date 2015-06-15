#Function for identifying the drilling fluid weights for nearest neighbor wells to use in the Allegheny Plateau BHT correction

ProbabilityCheck = function(Data,             #Full database 
                            DrillFluidData,   #Drilling fluid data. Not necessarily the same as the database being used.
                            pt_eval = 25,     #Number of nearest neighbor points used to calculate the percentage of air or mud drilled wells
                            rad_max = 50000)  #Maximum distance to search for wells (m)
{   
  
  # checking that x and y coordinates in Data are properly defined
  # printing error if the coordinates are not defined
  if(length(Data$POINT_X)==0){           
    print("No x-coordinate variable")
  }
  else if(sum(is.na(Data$POINT_X)+1-1)!=0){
    print("Test not performed beause of NAs in the 'POINT_X' of Data")
    print("Remove NAs from 'POINT_X' and call function again.")
  }
  else if(length(Data$POINT_Y)==0){
    print("No y-coordinate variable in Data")
  }
  else if(sum(is.na(Data$POINT_Y)+1-1)!=0){
    print("Test not performed beause of NAs in the 'POINT_Y' of Data.")
    print("Remove NAs from 'POINT_Y' and call function again.")
  }
  else if(length(DrillFluidData$POINT_X)==0){
    print("No x-coordinate variable in DrillFluidData")
  }
  else if(sum(is.na(DrillFluidData$POINT_X)+1-1)!=0){
    print("Test not performed beause of NAs in the 'POINT_X' of DrillFluidData.")
    print("Remove NAs from 'POINT_X' and call function again.")
  }else if(length(DrillFluidData$POINT_Y)==0){
    print("No y-coordinate variable")
  }
  else if(sum(is.na(DrillFluidData$POINT_Y)+1-1)!=0){
    print("Test not performed beause of NAs in the 'POINT_Y' of DrillFluidData.")
    print("Remove NAs from 'POINT_Y' and call function again.")
  }else if (length(DrillFluidData$gen_name3)==0){
    print("No 'gen_name3' variable in DrillingFluidData")
  }
  else if (sum(is.na(DrillFluidData$gen_name3)+1-1)!=0){ # +1-1 transforms TRUE/FALSE into 1 or 0. condition for any NA present
    print("Test not performed because of NAs in 'gen_name3'.")
    print("Remove NAs from DrillFluidData and call function again.")
  }else if (length(Data$gen_name3)==0){
    print("No 'gen_name3' variable in Data")
  }
  else{ #Coordinates defined properly. Drilling fluid is present. Proceed with calculation.
    
    #Check for the nearest neighbors
    for(i in 1:nrow(Data)){ #loop over all observations in Data that need drilling fluid assigned
      
      if (Data@data$reg[i] == 1) { #These are Allegheny Plateau Wells
        # initializing matrix with one column for index and the other for distance
        dist_vec <- matrix(0,nrow(Data),1) 
        
        # calculating weighted Euclidian distance for all of the points
        #Fixme: This can be done using spDistsN1 with longlat = FALSE. We may want to use the Great Circle ellipsoidal distance, though (TRUE).
        dist_vec[,1] <- sqrt((DrillFluidData@data$POINT_X-Data@data$POINT_X[i])^2 + (DrillFluidData@data$POINT_Y-Data@data$POINT_Y[i])^2)
        
        # finding the empirical quantile for removing, pt_eval is an input
        dist_cutoff <- sort(dist_vec)[pt_eval]
        
        # storing the maximum distance of the nearest pt_eval points
        Data$max_pt_rad[i] <- dist_cutoff
        
        # evaluating if radius is too large to capture the points
        if(Data$max_pt_rad[i] > rad_max){ # case for points criterion not satisfied within rad_max
          
          # assigning error for not enough points to evaluate
          Data$max_pt_error[i] <- 1
          
        }
        else{ # case for points criterion satisfied within rad_max
          
          if (is.na(Data$gen_name3[i]) == FALSE) {
            #This well has a drilling fluid. Use it. Assign Probability accordingly
            if (Data$gen_name3 == 'all_agfs') { #The well was drilled using air
              Data$PctAir[i] = 1
              Data$PctMud[i] = 0
            }
            else { #The well was drilled using mud
              Data$PctAir[i] = 0
              Data$PctMud[i] = 1
            } 
          }
          else { #The well doesn't have a drilling fluid, use the neighboring points to determine the BHT correction weights
            
            # finding the indices of the closest points
            inds_within_cutoff <- which(dist_vec <= dist_cutoff)
            
            # finding percentage of air drilled wells from reference database
            PercentAir <- length(which(Data@data$gen_name3[inds_within_cutoff] == 'all_agfs'))/length(inds_within_cutoff)
            PercentMud <- 1.0 - PercentAir
            
            # storing the percentage of air and gas
            Data$PctAir[i] = PercentAir
            Data$PctAir[i] = PercentMud
          } 
        }
      }
      else {
        #These are not Allegheny Plateau, assign NA to the added columns
        Data$PctAir[i] = NA
        Data$PctAir[i] = NA
      } 
    } # end for loop
    
    return(Data) 
  }  
}
