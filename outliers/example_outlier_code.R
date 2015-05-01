# code to generate dataset to test outlier algorithm

# packages/libraries----
library(rgdal)
library(data.table)
library(fields) # for colorbar
library(aqfig)

# set working directory, will need to change depending on user
setwd("/Users/calvinwhealton/GitHub/Geothermal_Codes")

# loading functions from a script
source("outlier_identification.R")

# importing data----
cornell_data <- read.table('cornell_data.csv', header=TRUE, sep=',')
grid_test <- read.table('out_test_grid.csv', header=TRUE, sep=',')
rad_test <- read.table('out_test_rad.csv', header=TRUE, sep=',')
pt_test <- read.table('out_test_pt.csv', header=TRUE, sep=',')

# will be used for testing with 'test' is not defined
cornell_datanotest <- cornell_data 

# transforming into UTM coordinate system----

# making matrix for the coordinates from input data
conv_coord <- matrix(0,nrow(cornell_data),2)
conv_coord[,1] <- cornell_data$decimal_lo
conv_coord[,2] <- cornell_data$decimal_la
conv_data <- data.frame(conv_coord)
colnames(conv_data) <- c("lat", "long")

# converting coordinates
coordinates(conv_data) <- c("long", "lat")
(conv_data.crs<-CRS("+proj=longlat +datum=WGS84")) # spacing is important
proj4string(conv_data)<-conv_data.crs
conv_data <- project(conv_coord, "+proj=utm +zone=18 ellps=WGS84")

# converting into km as distance units
cornell_data$x_coord <- conv_data[,1]*10^(-3)
cornell_data$y_coord <- conv_data[,2]*10^(-3)

# calculating the harrison gradient----
cornell_data$harr_grad <- 1000*(cornell_data$bht_c - 16.512 + 0.018268*cornell_data$calc_depth - 
                             (2.3449*10^(-6))*cornell_data$calc_depth^2 - 9)/cornell_data$calc_depth

cornell_data$test <- cornell_data$harr_grad

# creating datasets with NAs for error check----
cornell_dataNAx <- cornell_data
cornell_dataNAy <- cornell_data
cornell_dataNAtest <- cornell_data

cornell_dataNAx$x_coord[c(4,5,6)] <- NA
cornell_dataNAy$y_coord[c(1,2,3)] <- NA
cornell_dataNAtest$test[c(10,20,30)] <- NA

## testing the individual functions
# testing global function----
cornell_data1 <- outlier_glob(X=cornell_data, k_glob=3, type=7)

glob_outs <- sum(cornell_data1$out_glob_lo) + sum(cornell_data1$out_glob_hi)

# calculating by hand
glob_quant <- as.numeric(quantile(cornell_data$harr_grad,c(0.25,0.5,0.75)))
glob_lb <- glob_quant[1] - 3*(glob_quant[2]-glob_quant[1])
glob_ub <- glob_quant[3] + 3*(glob_quant[3]-glob_quant[2])

glob_lo <- as.numeric(length(which(cornell_data$harr_grad < glob_lb)))
glob_hi <- as.numeric(length(which(cornell_data$harr_grad > glob_ub)))

# testing gridded local function----
grid_test2 <- outlier_loc_grid(X=grid_test,box_size=32,pt_min=25,k_loc=3,type=7)

# testing radius-based local function-----
rad_test2 <- outlier_loc_rad(X=rad_test,rad_eval=16,pt_min=25,k_loc=3,type=7)

# testing point-based local function-----
pt_test2 <- outlier_loc_pts(X=pt_test,pt_eval=25,rad_max=32,k_loc=3,type=7)

## testing whole algorithm
# testing NA errors in data ----
test_xNA <- outlier_iden(X=cornell_dataNAx
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16
                         , box_size = 32
                         , pt_min = 25
                         , rad_max = 16
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)

test_yNA <- outlier_iden(X=cornell_dataNAy
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16
                         , box_size = 32
                         , pt_min = 25
                         , rad_max = 16
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)

test_NAtest <- outlier_iden(X=cornell_dataNAtest
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16
                         , box_size = 32
                         , pt_min = 25
                         , rad_max = 16
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)
# testing improperly defined variables ----

# deleting the 'test' column
cornell2 <- cornell_data
cornell2$test <- NULL

test_noname1 <- outlier_iden(X=cornell2
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16
                         , box_size = 32
                         , pt_min = 25
                         , rad_max = 16
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)


# deleting the 'x_coord' column
cornell2 <- cornell_data
cornell2$x_coord <- NULL

test_noname2 <- outlier_iden(X=cornell2
                             , algo = 1
                             , outcri = 1
                             , pt_eval = 25
                             , rad_eval = 16
                             , box_size = 32
                             , pt_min = 25
                             , rad_max = 16
                             , k_glob = 3 
                             , k_loc = 3  
                             , type = 7)

# deleting the 'y_coord' column
cornell2 <- cornell_data
cornell2$y_coord <- NULL

test_noname3 <- outlier_iden(X=cornell2
                             , algo = 1
                             , outcri = 1
                             , pt_eval = 25
                             , rad_eval = 16
                             , box_size = 32
                             , pt_min = 25
                             , rad_max = 16
                             , k_glob = 3 
                             , k_loc = 3
                             , type = 7)

# running global outlier function----
cd_glob <- cornell_data
cd_glob <- outlier_iden(X=cd_glob
                             , algo = 1
                             , outcri = 3
                             , pt_eval = 25
                             , rad_eval = 16
                             , box_size = 32
                             , pt_min = 25
                             , rad_max = 16
                             , k_glob = 3 
                             , k_loc = 3  
                             , type = 7)

# running local outlier functions----
# gridded
cd_loc_grid <- cornell_data
cd_loc_grid <- outlier_iden(X=cd_loc_grid
                        , algo = 1
                        , outcri = 1
                        , pt_eval = 25
                        , rad_eval = 16
                        , box_size = 32
                        , pt_min = 25
                        , rad_max = 16
                        , k_glob = 3 
                        , k_loc = 3  
                        , type = 7)

cd_loc_rad <- cornell_data
cd_loc_rad <- outlier_iden(X=cd_loc_rad
                            , algo = 2
                            , outcri = 1
                            , pt_eval = 25
                            , rad_eval = 16
                            , box_size = 32
                            , pt_min = 25
                            , rad_max = 16
                            , k_glob = 3 
                            , k_loc = 3  
                            , type = 7)

cd_loc_pt <- cornell_data
cd_loc_pt <- outlier_iden(X=cd_loc_pt
                           , algo = 2
                           , outcri = 1
                           , pt_eval = 25
                           , rad_eval = 16
                           , box_size = 32
                           , pt_min = 25
                           , rad_max = 16
                           , k_glob = 3 
                           , k_loc = 3  
                           , type = 7)

# running local and global outlier functions----
cd_loc_glob_grid <- cornell_data
cd_loc_glob_grid <- outlier_iden(X=cd_loc_glob_grid
                            , algo = 1
                            , outcri = 2
                            , pt_eval = 25
                            , rad_eval = 16
                            , box_size = 32
                            , pt_min = 25
                            , rad_max = 16
                            , k_glob = 3 
                            , k_loc = 3  
                            , type = 7)

cd_loc_glob_rad <- cornell_data
cd_loc_glob_rad <- outlier_iden(X=cd_loc_glob_rad
                           , algo = 2
                           , outcri = 2
                           , pt_eval = 25
                           , rad_eval = 16
                           , box_size = 32
                           , pt_min = 25
                           , rad_max = 16
                           , k_glob = 3 
                           , k_loc = 3  
                           , type = 7)

cd_loc_glob_pt <- cornell_data
cd_loc_glob_pt <- outlier_iden(X=cd_loc_glob_pt
                          , algo = 2
                          , outcri = 2
                          , pt_eval = 25
                          , rad_eval = 16
                          , box_size = 32
                          , pt_min = 25
                          , rad_max = 16
                          , k_glob = 3 
                          , k_loc = 3  
                          , type = 7)

# sensitivity analysis for points algorithm----
pts_sens <- c(10,25,50,100,200) # number of points for local neighborhood
rad_sens <- c(4,8,16,32,64) # maximum size of radius

outs_iden <- matrix(0,length(pts_sens),length(rad_sens)) # matrix to hold number of outliers
outs_sparse <- matrix(0,length(pts_sens),length(rad_sens)) # number of points in sparse areas


for(i in 1:length(pts_sens)){
  for(j in 1:length(rad_sens)){
    sens_data <- cornell_data
    
    sens_data2 <- outlier_iden(X=sens_data
                              , algo = 1
                              , outcri = 1
                              , pt_eval = pts_sens[i]
                              , rad_eval = 16
                              , box_size = 32
                              , pt_min = 25
                              , rad_max = rad_sens[j]
                              , k_glob = 3 
                              , k_loc = 3  
                              , type = 7)
    
    outs_iden[i,j] <- sum(sens_data2$outs)
    outs_sparse[i,j] <- sum(is.na(sens_data2$out_loc_lb)+1-1)
  }
}

# creating plot
setEPS()
postscript("outlier_sens.eps")


cols <- topo.colors(766)
par(mar =c(4.5,4.5,2,12)+0.1) 

data <- matrix(0,length(pts_sens)*length(rad_sens),4)

for(i in 1:length(rad_sens)){
  inds <- seq((i-1)*length(rad_sens)+1,i*length(rad_sens), by=1)
  
  data[c(inds),1] <- outs_iden[,i]
  data[c(inds),2] <- outs_sparse[,i]
  data[c(inds),3] <- c(10,25,50,100,200)
  data[c(inds),4] <- rad_sens[i]
}

dataplot <- data.frame(data)

colnames(dataplot) <- c("iden", "sparse", "pts", "rad")

dataplot$cols <- cols[dataplot$iden+1]
dataplot$cex <- sqrt(8919 - dataplot$sparse)/30

plot(dataplot$pts
     , log(dataplot$rad, base=2)
     , log='x'
     , cex = 1.02*sqrt(8919)/30
     , col = "black"
     , pch = 19
     , ylab = "Max Radius (km)"
     , xlab = "Points to Evaluate"
     , xaxt = "n"
     , yaxt = "n"
     , ylim = c(1.5,6.5)
     , xlim = c(7,300)
     )

points(dataplot$pts
     , log(dataplot$rad, base=2)
     , cex = 0.98*sqrt(8919)/30
     , col = "white"
     , pch = 19
)

points(dataplot$pts # x value
    , log(dataplot$rad, base=2)   # y value
    , cex = dataplot$cex
    , col = dataplot$cols
    , pch = 19
  )

axis(1, at=c(10,25,50,100,200), labels=c("10", "25", "50", "100", "200") )
axis(2, at=c(2,3,4,5,6), labels=c("4", "8", "16", "32", "64") )



par(xpd = TRUE)
# legend(x = 300
#        , y = 6 # location
#        , legend=c("% Tested", "(size)", "20", "40", "60", "80", "100", "% Outs of All Data", "(color)", "0", "2", "4", "6", "8") # legend entries
#        , pch = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 19, 19, 19, 19, 19)
#        , col=c(NA, NA, NA, NA, NA, NA, NA, NA, NA, cols[1], cols[180], cols[359], cols[538], cols[717] )
#        , ncol = 2
#        , pt.cex = c(0, 0, 0, 0, 0, 0, 0, 0, 0, sqrt(8919/5)/30, sqrt(8919*2/5)/30, sqrt(8919*3/5)/30, sqrt(8919*4/5)/30, sqrt(8919)/30)
# ) # colors

# legend(x = 300
#        , y = 60 # location
#        , legend=c("% Tested", "20", "40", "60", "80", "100") # legend entries
#        , pch = c( NA, 19, 19, 19, 19, 19)
#        , col=c( NA, "gray", "gray", "gray", "gray", "gray")
#        , ncol = 1
#        , pt.cex = c(0, sqrt(8919/5)/30, sqrt(8919*2/5)/30, sqrt(8919*3/5)/30, sqrt(8919*4/5)/30, sqrt(8919)/30)
# ) # colors


legend(x = 400
       , y = 6.5 # location
       , legend=c("# Outs/ # All (%)", "0", "1", "2", "3", "4",  "5", "6", "7", "8") # legend entries
       , pch = c(NA, 19, 19, 19, 19, 19, 19, 19, 19, 19)
       , col=c(NA, cols[1], cols[91], cols[179], cols[269], cols[358],  cols[447], cols[537], cols[626], cols[714]  )
       , ncol = 1
) # colors

text(x=400
     , y = 1+ c(2.5,2.1, 1.7, 1.1)
     , labels=c("Point size is ", "  proportional to %", "  of data tested.", "Black circle is 100%.")
     , adj = c(0,0)
     )


dev.off()
