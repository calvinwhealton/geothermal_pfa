# code to run outlier identification

# function inputs are:
# X         dataframe of values including
#             x_coord    x coordinate
#             y_coord    y coordinate
#             test       value to be tested
# algo      algorithm for local testing the data
#             1    finding nearest pt_eval points to the outlier, must be within rad_max
#             2    finding all points within rad_eval from the outlier, must satisfy min_pt
#             3    for gridding data
# outcri    criteria used to determine outliers
#             1   only local outliers removed
#             2   only local and global outliers removed
#             3   only global outliers removed
# pt_eval   number of points to evaulate local outliers
# rad_eval  radius used for local evaluation
# box_size  size of box for gridded local search
# pt_min    minimum number of points to evaluate local outliers
# rad_max   maximum radus to use in search
# k_glob    constant used for outlier identification
# k_loc     constant used for local outlier identification
# wts       weights for distance in each direction
# type      type of quantile estimation algorithm
# rad_eval, rad_max, and box_size are assumed to be units of km, but can be specified for any coordinate system

# function output is:
# X         dataframe of all input with additional columns for
#             outs = 1 for outlier, outs = 0 for non-outlier
#             

# error codes
# 0   no errors
# 1   some of pt_eval points outside of rad_max
# 2   fewer than pt_min within rad_eval
outlier_iden <- function(X                # input dataframe
                        , algo = 1        # local detection algorithm
                        , outcri = 1      # outlier criterion
                        , pt_eval = 25    # minumum number of points to evaluate
                        , rad_eval = 16   # radius used for local outlier identification
                        , box_size = 32   # edge length of the box
                        , pt_min = 25     # minimum number of points to evaluate local outliers
                        , rad_max = 16    # maximum radius to search for local points
                        , k_glob = 3      # constant for quartile-median range in global analysis
                        , k_loc = 3       # constant for quartile-median range in local analysis
                        , type = 7        # type of quantile interpolation
                        ){
  
  # checking that x and y coordinates are properly defined
  # printing error if the coordinates are not defined
  if(length(X$x_coord)==0){           
    print("No x-coordinate variable")
  }else if(sum(is.na(X$x_coord)+1-1)!=0){
    print("Test not performed beause of NAs in the 'x_coord'.")
    print("Remove NAs from 'x_coord' and call function again.")
  }else if(length(X$y_coord)==0){
    print("No y-coordinate variable")
  }else if(sum(is.na(X$y_coord)+1-1)!=0){
    print("Test not performed beause of NAs in the 'y_coord'.")
    print("Remove NAs from 'y_coord' and call function again.")
  }else if (length(X$test)==0){
    print("No 'test' variable in the dataset")
  }else if (sum(is.na(X$test)+1-1)!=0){ # +1-1 transforms TRUE/FALSE into 1/0, conditoin for any NA present
    print("Test not performed because of NAs in 'test'.")
    print("Remove NAs from 'test' and call function again.")
  }else{
    
    if(outcri == 1){ # case for local outlier identification only
      
      if(algo == 1){ # case for finding nearest pt_eval points within rad_max
        X <- outlier_loc_pts(X,pt_eval,rad_max,k_loc,type)
      }else if (algo == 2){  # case for finding points within rad_eval but must have points min_pts
        X <- outlier_loc_rad(X,rad_eval,pt_min,k_loc,type)
      }else if (algo == 3){ # case for algorithm using fixed-grid boxes
       X <- outlier_log_grid(X,box_size,pt_min,type)
      }else{
       print("Not a valid local outlier detection algorithm")
      }
     
     # outlier identified as either low or high local outlier
     X$outs <- X$out_loc_lo + X$out_loc_hi
     
    }else if(outcri == 2){ # case for local and global outlier identification
      
      # finding local outliers
      if(algo == 1){ # case for finding nearest pt_eval points with-in rad_max
        X <- outlier_loc_pts(X,pt_eval,rad_max,k_loc,type)
      }else if (algo == 2){  # case for finding points within rad_eval but must have points min_pts
        X <- outlier_loc_rad(X,rad_eval,pt_min,k_loc,type)
      }else if (algo == 3){ # case for algorithm using fixed-grid boxes
        X <- outlier_log_grid(X,box_size,pt_min,type)
      }else{
        print("Not a valid local outlier detection algorithm")
      }
      
      # finding global outliers
      X <- outlier_glob(X,k_glob,type)
      
      # outliers identified must be local and global
      X$outs <- (X$out_loc_lo + X$out_loc_hi)*(X$out_glob_lo + X$out_glob_hi)
      
    }else if(outcri == 3){ # case for global outlier identification only
      
      # finding global outliers
      X <- outlier_glob(X,k_glob,type)
      
      # outliers identified must be local and global
      X$outs <- X$out_glob_lo + X$out_glob_hi
      
    }else{ # case for not a properly defined outlier criterion
      print("Not a valid outlier criterion")
    } 
  }
  
  return(X)
}

# function for local outlier identification
# finds pt_eval points that are within rad_max from tested point
# does not calculate outliers when pt_eval has points outside rad_max
outlier_loc_pts <- function(X,pt_eval,rad_max,k_loc,type){
  
  # initializing varaibles to hold flags and outputs
  X$out_loc_lo <- 0 # low outlier
  X$out_loc_hi <- 0 # high outlier
  X$out_loc_lq <- NA # lower quartile
  X$out_loc_uq <- NA # upper quartile
  X$out_loc_mq <- NA # middle quartile
  X$out_loc_rad <- 0 # radius of maximum distance
  X$out_loc_lb <- NA # lower bound 
  X$out_loc_ub <- NA # upper bound
  X$out_loc_error <- 0 # error variable
  
  for(i in 1:nrow(X)){ # loop over all observations
    
    # initializing matrix with one column for index and the other for distance
    dist_vec <- matrix(0,nrow(X),1) 
    
    # calculating weighted Euclidian distance for the points
    dist_vec[,1] <- sqrt((X$x_coord-X$x_coord[i])^2 + (X$y_coord-X$y_coord[i])^2)
    
    # finding the empirical quantile for removing, pt_eval is an input
    dist_cutoff <- sort(dist_vec)[pt_eval]
    
    # storing the maximum distance of the nearest 25 points
    X$out_loc_rad[i] <- dist_cutoff
    
    # evaluating if radius is too large to capture the points
    if(X$out_loc_rad[i] > rad_max){ # case for points criterion not satisfied within rad_max
      
      # assigning error for not enough points to evaluate
      X$out_loc_error[i] <- 1
      
    }else{ # case for points criterion satisfied within rad_max
      
      # finding the indices of the closest points, including the point tested
      inds_within_cutoff <- which(dist_vec <= dist_cutoff)
      
      # finding quartiles
      quarts <- as.numeric(quantile(X$test[inds_within_cutoff],probs=c(0.25,0.5,0.75),type=type))
      
      # storing quartile and bound values
      X$out_loc_lq[i] <- quarts[1]
      X$out_loc_mq[i] <- quarts[2]
      X$out_loc_uq[i] <- quarts[3]
      X$out_loc_lb[i] <- X$out_loc_lq[i] - k_loc*(X$out_loc_mq[i] - X$out_loc_lq[i])
      X$out_loc_ub[i] <- X$out_loc_uq[i] + k_loc*(X$out_loc_uq[i] - X$out_loc_mq[i])
      
      # testing for local outliers
      if(X$test[i] > X$out_loc_ub[i]){ # testing upper outliers
        
        # flagging as a high outlier
        X$out_loc_hi[i] <- 1
        
      }else if(X$test[i] < X$out_loc_lb[i]){ # testing low outliers
        
        # flagging as a low outlier
        X$out_loc_lo[i] <- 1
      }
    }
    
  } # end for loop
  
  return(X)
}

# function for local outlier identification
# finds points within rad_eval
# does not calculate outliers when fewer than pt_min points within rad_max
outlier_loc_rad <- function(X,rad_eval,pt_min,k_loc,type){
  
  # initializing varaibles to hold flags and outputs
  X$out_loc_lo <- 0 # low outlier
  X$out_loc_hi <- 0 # high outlier
  X$out_loc_lq <- NA # lower quartile
  X$out_loc_mq <- NA # middle quartile
  X$out_loc_uq <- NA # upper quartile
  X$out_loc_lb <- NA # lower bound
  X$out_loc_ub <- NA # upper bound
  X$out_loc_pts <- 0 # number of points
  X$out_loc_error <- 0 # error
  
  for(i in 1:nrow(X)){
    
    # initializing matrix with one column for index and the other for distance
    dist_vec <- matrix(0,nrow(X),1) 
    
    # calculating weighted Euclidian distance for the points
    dist_vec[,1] <- sqrt((X$x_coord-X$x_coord[i])^2 + (X$y_coord-X$y_coord[i])^2)
    
    # finding the indices of the closest points, including the point tested
    inds_within_cutoff <- which(dist_vec < rad_eval)
    
    X$out_loc_pts[i] <- as.numeric(length(inds_within_cutoff))
    
    # evaluating if enough points are within the radius
    if(X$out_loc_pts[i] < pt_min){ # case for not enough points
      
      # assigning error for not enough points to evaluate
      X$out_loc_error[i] <- 1
      
    }else{ # case for points criterion satisfied within rad_eval
      
      # finding quartiles
      quarts <- as.numeric(quantile(X$test[inds_within_cutoff],probs=c(0.25,0.5,0.75),type=type))
      
      # storing calculated values
      X$out_loc_lq[i] <- quarts[1]
      X$out_loc_mq[i] <- quarts[2]
      X$out_loc_uq[i] <- quarts[3]
      X$out_loc_lb[i] <- X$out_loc_lq[i] - k_loc*(X$out_loc_mq[i] - X$out_loc_lq[i])
      X$out_loc_ub[i] <- X$out_loc_uq[i] + k_loc*(X$out_loc_uq[i] - X$out_loc_mq[i])
      
      # testing for local outliers
      if(X$test[i] > X$out_loc_ub[i]){ # testing upper outliers
        
        # flagging as a high outlier
        X$out_loc_hi[i] <- 1
        
      }else if(X$test[i] < X$out_loc_lb[i]){ # testing low outliers
        
        # flagging as a low outlier
        X$out_loc_lo[i] <- 1
      }
    }
  } # end for loop
  
  return(X)
}

# function for local outlier identification
# finds points within a fixed grid system with box edge size box_size
# does not calculate outliers when fewer than pt_min points in a box
outlier_loc_grid <- function(X,box_size,pt_min,k_loc,type){
  
  # initializing varaibles to hold flags and outputs
  X$out_loc_lo <- 0 # low outlier
  X$out_loc_hi <- 0 # high outlier
  X$out_loc_lq <- NA # lower quartile
  X$out_loc_uq <- NA # upper quartile
  X$out_loc_mq <- NA # middle quartile
  X$out_loc_pts <- 0 # points in grid cell
  X$out_loc_lb <- NA # lower bound 
  X$out_loc_ub <- NA # upper bound
  X$out_loc_error <- 0 # error variable
  
  # finding minimum values of coordinates
  x_min <- min(X$x_coord)
  y_min <- min(X$y_coord)
  
  # calculating the box number
  X$x_box <- floor((X$x_coord - x_min)/box_size)
  X$y_box <- floor((X$y_coord - y_min)/box_size)
  
  # fining maximum boxes
  max_box_x <- max(X$x_box)
  max_box_y <- max(X$y_box)
  
  for(i in 1:nrow(X)){
    
    # finding indices for the points boxes
    inds_box <- intersect(which(X$x_box == X$x_box[i]), which(X$y_box == X$y_box[i]))
    
    # storing number of points
    X$out_loc_pts[i] <- as.numeric(length(inds_box))
    
    # evaluating if enough points are within the radius
    if(X$out_loc_pts[i] < pt_min){ # case for not enough points
      
      # assigning error for not enough points to evaluate
      X$out_loc_error[inds_box] <- 1
      
    }else{ # case for points criterion satisfied
      
      # finding quartiles
      quarts <- as.numeric(quantile(X$test[inds_box],probs=c(0.25,0.5,0.75),type=type))
      
      # storing quartile and bound values
      X$out_loc_lq[i] <- quarts[1]
      X$out_loc_mq[i] <- quarts[2]
      X$out_loc_uq[i] <- quarts[3]
      X$out_loc_lb[i] <- X$out_loc_lq[i] - k_loc*(X$out_loc_mq[i] - X$out_loc_lq[i])
      X$out_loc_ub[i] <- X$out_loc_uq[i] + k_loc*(X$out_loc_uq[i] - X$out_loc_mq[i])
      
      # testing for local outliers
      if(X$test[i] > X$out_loc_ub[i]){ # testing upper outliers
        
        # flagging as a high outlier
        X$out_loc_hi[i] <- 1
        
      }else if(X$test[i] < X$out_loc_lb[i]){ # testing low outliers
        
        # flagging as a low outlier
        X$out_loc_lo[i] <- 1
      }
    }
  } # end for loop

  return(X)
}


# function for global outlier identification
outlier_glob <- function(X,k_glob,type){
  
  X$out_glob_lo <- 0
  X$out_glob_hi <- 0
  
  # finding quartiles
  quants <- as.numeric(quantile(X$test,probs=c(0.25,0.5,0.75),type=type))
  
  # calculating bounds
  lb <- quants[1] - k_glob*(quants[2]-quants[1])
  ub <- quants[3] + k_glob*(quants[3]-quants[2])
  
  # indices for points outside bounds
  lb_inds <- which(X$test < lb)
  X$out_glob_lo[lb_inds] <- 1
  
  ub_inds <- which(X$test > ub)
  X$out_glob_hi[ub_inds] <- 1
  
  
  return(X)
}
