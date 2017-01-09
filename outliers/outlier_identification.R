# code to run outlier identification

# outlier_iden function inputs are:
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
# rad_eval  radius used for local evaluation (m)
# box_size  size of box for gridded local search (m)
# pt_min    minimum number of points to evaluate local outliers
# rad_max   maximum radius to use in search (m)
# k_glob    constant used for outlier identification
# k_loc     constant used for local outlier identification
# wts       weights for distance in each direction
# type      type of quantile estimation algorithm

# rad_eval, rad_max, and box_size are assumed to be units of m, but can be specified for any coordinate system

# outlier_iden function output is:
# X         dataframe of all input with additional columns for
#             outs = 1 for outlier, outs = 0 for non-outlier
#             

# error codes
# 0   no errors
# 1   some of pt_eval points outside of rad_max
# 2   fewer than pt_min within rad_eval

# Additional inputs are available in the select_out_algo function, which calls outlier_iden:
# Threshold   Exclude all points with values less than Threshold
# OutVarName  Desired column name for the tested variable in the output file
# InVarName   Column name of the thermal variable in Data. May be the same as OutVarName.
# X_coordName Column name of the UTM longitude in Data.
# Y_coordName Column name of the UTM latitude in Data. 

# select_out_algo function output is:
# TestedOutliers    List containing a dataframe of the outliers and
#                   a dataframe of the not outliers.

#### Functions ####
# Function used to set up the dataframe column names and run the outlier analysis. 
# All options for outlier identification Are available in this function.
# There is an option to set a threshold, below which all values will be dropped before the outlier identification. 
# The default is to drop all negative values.
select_out_algo <- function(Data,            # Dataframe to be tested for outliers
                            OutVarName,      # Desired column name for the tested variable in the output file
                            InVarName,       # Column name of the thermal variable in Data. May be the same as OutVarName.
                            X_coordName,     # Column name of the UTM longitude in Data.
                            Y_coordName,     # Column name of the UTM latitude in Data. 
                            Threshold = 0.0, # Points in column InVarName less than Threshold will be dropped before outlier analysis. Default is 0.0.
                            algo = 1,        # local detection algorithm
                            outcri = 1,      # outlier criterion
                            pt_eval = 25,    # minumum number of points to evaluate
                            rad_eval = 16000,# radius used for local outlier identification (m)
                            box_size = 32000,# edge length of the box (m)
                            pt_min = 25,     # minimum number of points to evaluate local outliers
                            rad_max = 16000, # maximum radius to search for local points (m)
                            k_glob = 3,      # constant for quartile-median range in global analysis
                            k_loc = 3,       # constant for quartile-median range in local analysis
                            type = 7         # type of quantile interpolation
){
  
  #Rename columns based on the necessary inputs to the outlier identification function.
  colnames(Data@data)[which(colnames(Data@data)==InVarName)] = "test"
  colnames(Data@data)[which(colnames(Data@data)==X_coordName)] = "x_coord"
  colnames(Data@data)[which(colnames(Data@data)==Y_coordName)] = "y_coord"
  
  #Remove all values less than Threshold.
  if (length(which(Data@data[, "test"] <= Threshold)) > 0) {
    Data = Data[-which(Data@data["test"] <= Threshold),] 
  }
  
  #Run outlier identification
  Outs = outlier_iden(Data
                      , algo = algo          # local detection algorithm
                      , outcri = outcri      # outlier criterion
                      , pt_eval = pt_eval    # minumum number of points to evaluate
                      , rad_eval = rad_eval  # radius used for local outlier identification (m)
                      , box_size = box_size  # edge length of the box (m)
                      , pt_min = pt_min      # minimum number of points to evaluate local outliers
                      , rad_max = rad_max    # maximum radius to search for local points (m)
                      , k_glob = k_glob      # constant for quartile-median range in global analysis
                      , k_loc = k_loc        # constant for quartile-median range in local analysis
                      , type = type)         # type of quantile interpolation
  
  #Change the column name of the UTM coordinates back to the original name.
  colnames(Outs@data)[which(colnames(Outs@data)=="x_coord")] = X_coordName
  colnames(Outs@data)[which(colnames(Outs@data)=="y_coord")] = Y_coordName
  
  #Rename the thermal variable column according to the OutVarName
  colnames(Outs@data)[which(colnames(Outs@data)=="test")] = OutVarName
  
  #Split the data into Outliers and NotOutliers
  NotOutliers = Outs[which(Outs@data$outs == 0),]
  Outliers = Outs[which(Outs@data$outs != 0),]
  
  #List for returning data
  TestedOutliers = list("NotOutliers" = NotOutliers, "Outliers" = Outliers)
  
  return(TestedOutliers)
}

# Wrapper function that calls the outlier identification algorithm functions.
outlier_iden <- function(X                # input dataframe
                         , algo = 1        # local detection algorithm
                         , outcri = 1      # outlier criterion
                         , pt_eval = 25    # minumum number of points to evaluate
                         , rad_eval = 16000# radius used for local outlier identification
                         , box_size = 32000# edge length of the box
                         , pt_min = 25     # minimum number of points to evaluate local outliers
                         , rad_max = 16000 # maximum radius to search for local points
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

#### Sensitivity Analysis for Local Points Algorithm ####
library(rgdal)
#Load data
DataTest = readOGR(dsn = "C:\\Users\\Jared\\Documents\\Cornell\\Research\\Masters - Spatial Assessment\\Figures", layer = "DataForOutlierTesting", stringsAsFactors = FALSE)

pts_sens <- c(10,25,50,100,200) # number of points for local neighborhood
rad_sens <- c(4000,8000,16000,32000,64000) # maximum size of radius

outs_iden <- matrix(0,length(pts_sens),length(rad_sens)) # matrix to hold number of outliers
outs_sparse <- matrix(0,length(pts_sens),length(rad_sens)) # number of points in sparse areas

#Calculate Outliers
for(i in 1:length(pts_sens)){
  for(j in 1:length(rad_sens)){
    sens_data <- DataTest
    sens_data2 <- select_out_algo(Data = sens_data,
                                  OutVarName = "Qs",
                                  InVarName = "Qs",
                                  X_coordName = "POINT_X",
                                  Y_coordName = "POINT_Y",
                                  Threshold = 0.0,
                                  algo = 1,
                                  outcri = 1,
                                  pt_eval = pts_sens[i],
                                  rad_eval = 16000,
                                  box_size = 32000,
                                  pt_min = 25,
                                  rad_max = rad_sens[j],
                                  k_glob = 3,
                                  k_loc = 3,
                                  type = 7)
    
    outs_iden[i,j] <- nrow(sens_data2$Outliers)
    outs_sparse[i,j] <- sum(sens_data2$NotOutliers$out_loc_error)
  }
}

# creating plot
setwd("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\ESDA")
setEPS()
postscript(file = "outlier_sens.eps", title = "Sensitivity Outliers Local Points", width = 5, height = 5)

#Make color ramp
Pal = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
cols <- Pal(max(outs_iden)+1)
par(mar =c(3,3,0,9)+0.1) 

data <- matrix(0,length(pts_sens)*length(rad_sens),4)

for(i in 1:length(rad_sens)){
  inds <- seq((i-1)*length(rad_sens)+1,i*length(rad_sens), by=1)
  
  data[c(inds),1] <- outs_iden[,i]
  data[c(inds),2] <- outs_sparse[,i]
  data[c(inds),3] <- pts_sens
  data[c(inds),4] <- rad_sens[i]/1000
}

dataplot <- data.frame(data)

colnames(dataplot) <- c("iden", "sparse", "pts", "rad")

#Changing the plotting location for pts = 10 to 12.5 for equal spacing in log base 2
dataplot$pts[dataplot$pts == 10] = 12.5

#Assigning color and size of points
dataplot$cols <- cols[dataplot$iden+1]
dataplot$cex <- sqrt(nrow(DataTest) - dataplot$sparse)/30

plot(dataplot$pts
     , log(dataplot$rad, base=2)
     , log='x'
     , cex = 1.02*sqrt(nrow(DataTest))/30
     , col = "black"
     , pch = 19
     , ylab = "Max Radius (km)"
     , xlab = "Points to Evaluate"
     , xaxt = "n"
     , yaxt = "n"
     , ylim = c(1.9,6.1)
     , xlim = c(10.85,230)
     , line = 2
)
points(dataplot$pts
       , log(dataplot$rad, base=2)
       , cex = 0.98*sqrt(nrow(DataTest))/30
       , col = "white"
       , pch = 19
)
points(dataplot$pts # x value
       , log(dataplot$rad, base=2)   # y value
       , cex = dataplot$cex
       , col = dataplot$cols
       , pch = 19
)

axis(1, at=c(12.5, pts_sens[-1]), labels=pts_sens, padj = -0.5)
axis(2, at=seq(2,6,1), labels=rad_sens/1000, padj = 0.5)

par(xpd = TRUE)
legend(x = 300
       , y = 6 # location
       , legend=c("# Outs/ # All (%)", seq(0,8,1)) # legend entries
       , pch = c(NA, rep(19,9))
       , col = c(NA, cols[round(max(outs_iden)/(max(outs_iden)/nrow(DataTest))*seq(0,0.08,0.01),0) + 1])
       , ncol = 1
       )
text(x=300
     , y = 0.75 + c(2.5,2.1, 1.7, 1.1)
     , labels=c("Point size is ", "  proportional to %", "  of data tested.", "Black circle is 100%.")
     , adj = c(0,0)
     )
dev.off()
