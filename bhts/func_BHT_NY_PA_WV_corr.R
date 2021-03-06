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
#                       3 = Allegheny Plateau with Drilling Fluid]
#     PctMud        [estimated probability of mud, 0 < prob < 1]
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

NY_PA_WV_BHT2 <- function(X){
  
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
        
        # correction for portion of data within bounds, becomes positive at 305 m which is shallower than minimum depth used
        else if(X$calc_depth_m[i] > 305){
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + min(15, -1.99 + 0.00652*X$calc_depth_m[i])
        }
      
        # correction for portion greater than 0 < 305 m
        else if(X$calc_depth_m[i] > 0){
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i]
        }
      
        # no correction for negative depths and return an error
        else{
  
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i]
          BHT_corrected$corr_error[i] <- 21 # error for negative depth
        }
      }  
    
      # Allegheny Plateau with drilling fluid  correction  
      else if(X$reg[i] == 3){
        
        # checking for given depth
        if(is.na(X$calc_depth_m[i])){ 
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
          BHT_corrected$corr_error[i] <- 22 # error for missing depth
        }
      
        # correction for depth deeper than 4000 m
        else if(X$calc_depth_m[i] > 4000){
          #                                             mud portion             air portion
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + 37.8*X$PctMud[i] + 15.4*(1-X$PctMud[i])
        }
        # correction for depth > 2500m, < 4000 m
        else if(X$calc_depth_m[i] > 2500){
          #                                             mud portion             air portion
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + (0.0155*((1650^3 + X$calc_depth_m[i]^3)^(1/3) - 1650))*X$PctMud[i] + 15.4*(1-X$PctMud[i])
        }
        # correction for depth <  2500m, > 0
        else if(X$calc_depth_m[i] >= 0){
          #                                             mud portion + air portion
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i] + (0.0155*((1650^3 + X$calc_depth_m[i]^3)^(1/3) - 1650))*X$PctMud[i]  + (0.0104*((1090^3 + X$calc_depth_m[i]^3)^(1/3) - 1090))*(1-X$PctMud[i])
                                                      
                                                    
        }
        else{
          BHT_corrected$corr_bht_c[i] <- X$bht_c[i]
          BHT_corrected$corr_error[i] <- 21 # error for negative depth
        }
      }
    
    # categorical variable defined but not 0, 2, or 3
    else{
      BHT_corrected$corr_bht_c[i] <- X$bht_c[i] # providing no adjustment
      BHT_corrected$corr_error[i] <- 30
    }
      
    } # end case of BHT defined
  } # end for loop
  
  return(BHT_corrected)
  
} # end function