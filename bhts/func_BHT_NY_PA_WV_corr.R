# function to calculate NY, PA, and WV bottom-hole temperature corrections
# code developed for use in DOE geothermal play fairway analysis grant
# studying the Appalachian Basin

# function coded by Calvin Whealton (caw324@cornell.edu)
# code writting under R version 2.15.1

# function takes an R data frame with the following named variables
#     bht_c         [BHT in degrees C]
#     calc_depth_m  [calculated depth in m]
#     reg           [region coded by 
#                       1 = Allegheny Plateau 
#                       0 = Rome Trough and PA points south east
#                       2 = West Virginia
#                       3 = Modified Harrison]
# variable names are designed to match NGDS standard names when possible

# output matrix add columns for 
#     corr_bht_c    [BHT corrected in degrees C]
#     corr_error    [error codes described below]

# function returns error codes when data is missing or seems outside the normal range
# 0: no errors in the calculation
# 20: depth likely too deep
# 21: depth negative
# 22: depth is missing for Allegheny Plateau or West Virginia data
# 30: categorical variable defined but not 0, 1, 2
# 32: missing categorical variable
# 42: BHT is missing

NY_PA_BHT2 <- function(X){
  
  # initializing data frame to hold all initial data 
  # and new columns for the corrected values and errors
  BHT_corrected <- X
  BHT_corrected$corr_bht_c <- 0
  BHT_corrected$corr_error <- 0
  
  for(i in 1:nrow(X)){ # loop over all observations
    
    # case for missing BHT data
    if(is.na(X$bht_c[i])){
      BHT_corrected$corr_bht_c[i] <- NA # obviously bad value
      BHT_corrected$corr_error[i] <- 42 # error for missing bottom-hole temperature
    }
    
    # case for BHT present
    else{
    
      # case for missing categorical variable
      if(is.na(X$reg[i])){
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
        BHT_corrected$corr_error[i] <- 32 # no categorical variable error
      }
    
      # checking if observation is in Allegheny Plateau
      else if(X$reg[i] == 1){ 
      
        ## chekcing the depth variable & computing correction when appropriate
        # checking missing depth value
        if(is.na(X$calc_depth_m[i])==TRUE){ 
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
          BHT_corrected$corr_error[i] <- 22 # error for missing depth
        }
      
        # general case with no problems (0,6500)
        else if ((X$calc_depth_m[i] < 6500) && (X$calc_depth_m[i] > 0)){
        
          if(X$calc_depth_m[i] > 4000){
            BHT_corrected$corr_bht_c[i] <- X$bht_c[i] - 23.48 + 0.01791*4000 # correction at 4000 m used for interval deeper than 4000 m, avoids extrapolation
          }
          else{
            BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + max(0, - 23.48 + 0.01791*X$calc_depth_m[i])
          }
        }
      
        # checking for excessively deep measurements (6500, Inf)
        else if(X$calc_depth_m[i]>6500){ 
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] - 23.48 + 0.01791*4000 # correction at 4000 m used for interval deeper than 4000 m, avoids extrapolation
          BHT_corrected$corr_error[i] <- 20 # error for measurement too deep to be normal
        }
      
        # checking for negative depths (-Inf, 0)
        else{
          BHT_corrected$corr_bht_c[i]<- X$bht_c[i] # negative depth has correction 0
          BHT_corrected$corr_error[i] <- 21 # error for negative depth
        }
      
      } # end Allegheny Plateau calculation
    
      # checking if observation is in Rome Trough and PA points SE
      else if(X$reg[i] == 0){ 
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i]  # no temperature correction
      }
      
      
    # West Virginia correction  
    else if(X$reg[i] == 2){
        
      # checking for given depth
      if(is.na(X$calc_depth_m[i])){ 
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
        BHT_corrected$corr_error[i] <- 22 # error for missing depth
      }
        
      # correction uses maximum value within data
      else if(X$calc_depth_m[i] >= 6500){
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + 15
        BHT_corrected$corr_error[i] <- 20 # depth likely too deep
      }
        
      # correction for portion of data within bounds, becomes positive at 467 m which is shallower than minimum depth used
      else{
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + min(15, -3.562 + 0.00763*X$calc_depth_m[i])
      }
        
    }  
      
    # Modifed Harrison Correction  
    else if(X$reg[i] == 3){
      
      # checking for given depth
      if(is.na(X$calc_depth_m[i])){ 
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
        BHT_corrected$corr_error[i] <- 22 # error for missing depth
      }
      
      # correction uses maximum value for depths deeper than peak
      else if(X$calc_depth_m[i] >= 6500){
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + 19.07
        BHT_corrected$corr_error[i] <- 20 # depth likely too deep
      }
      
      else if(X$calc_depth_m[i] >= 3860){
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + 19.07
      }
      
      # correction for positive portion of correction less than peak correction
      else if(X$calc_depth_m[i] >= 1000){
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i] -16.51 + 0.0183*X$calc_depth_m[i] - (2.34*10^(-6))*((X$calc_depth_m[i])^2)
      }
     
      # no correction for data shallower than 1000 m
      else{
        BHT_corrected$corr_bht_c[i] <- X$bht_c[i]
      }
      
    }
    
    # categorical variable defined but not 0, 1, or 2
    else{
      BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
      BHT_corrected$corr_error[i] <- 30
    }
      
    } # end case of BHT defined
  } # end for loop
  
  return(BHT_corrected)
  
} # end function