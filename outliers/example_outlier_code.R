# code to generate dataset to test outlier algorithm

# packages/libraries/scripts----
library(rgdal)
library(data.table)
library(fields) # for colorbar
library(aqfig)
library(assertthat)

# set working directory, will need to change depending on user
setwd("/Users/calvinwhealton/GitHub/Geothermal_Codes")

# loading outlier functions from a script
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
colnames(conv_data) <- c("long", "lat")

# converting coordinates and transforming from WGS84 to UTM 17N
coordinates(conv_data) <- c("long", "lat")
proj4string(conv_data) <- CRS('+init=epsg:4326')
conv_data <- spTransform(conv_data, CRS('+init=epsg:26917'))

cornell_data$x_coord = conv_data@coords[,1]
cornell_data$y_coord = conv_data@coords[,2]

# calculating the harrison gradient----
cornell_data$harr_grad <- 1000*(cornell_data$bht_c - 16.512 + 0.018268*cornell_data$calc_depth - 
                             (2.3449*10^(-6))*cornell_data$calc_depth^2 - 9)/cornell_data$calc_depth

cornell_data$test <- cornell_data$harr_grad

# adding a column for 'censor' to the dataset (needed for the local points algorithm)----
cornell_data$censor <- cornell_data$test
pt_test$censor = pt_test$test

# creating datasets with NAs for error check----
cornell_dataNAx <- cornell_data
cornell_dataNAy <- cornell_data
cornell_dataNAtest <- cornell_data

cornell_dataNAx$x_coord[c(1,2,3)] <- NA
cornell_dataNAy$y_coord[c(1,2,3)] <- NA
cornell_dataNAtest$test[c(10,20,30)] <- NA

## testing the individual functions
# testing global function----
cornell_data1 <- outlier_glob(X=cornell_data, k_glob=3, type=7)

glob_outs <- sum(cornell_data1$out_glob_lo) + sum(cornell_data1$out_glob_hi)

# calculating manually
glob_quant <- as.numeric(quantile(cornell_data$harr_grad,c(0.25,0.5,0.75)))
glob_lb <- glob_quant[1] - 3*(glob_quant[2]-glob_quant[1])
glob_ub <- glob_quant[3] + 3*(glob_quant[3]-glob_quant[2])

glob_lo <- as.numeric(length(which(cornell_data$harr_grad < glob_lb)))
glob_hi <- as.numeric(length(which(cornell_data$harr_grad > glob_ub)))

if(are_equal(glob_outs, sum(glob_hi, glob_lo)) == FALSE){
  print('Test Failed. Check Global Outlier Function.')
}

# testing gridded local function----
grid_test2 <- outlier_loc_grid(X=grid_test,box_size=32000,pt_min=25,k_loc=3,type=7)

# testing radius-based local function-----
rad_test2 <- outlier_loc_rad(X=rad_test,rad_eval=16000,pt_min=25,k_loc=3,type=7)

# testing point-based local function-----
pt_test2 <- outlier_loc_pts(X=pt_test,pt_eval=25,rad_max=32000,k_loc=3,type=7,min_val=-Inf,max_val=Inf,rank=FALSE)

## testing whole algorithm
# testing NA errors in data ----
test_xNA <- outlier_iden(X=cornell_dataNAx
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16000
                         , box_size = 32000
                         , pt_min = 25
                         , rad_max = 16000
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)

test_yNA <- outlier_iden(X=cornell_dataNAy
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16000
                         , box_size = 32000
                         , pt_min = 25
                         , rad_max = 16000
                         , k_glob = 3 
                         , k_loc = 3  
                         , type = 7)

test_NAtest <- outlier_iden(X=cornell_dataNAtest
                         , algo = 1
                         , outcri = 1
                         , pt_eval = 25
                         , rad_eval = 16000
                         , box_size = 32000
                         , pt_min = 25
                         , rad_max = 16000
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
                         , rad_eval = 16000
                         , box_size = 32000
                         , pt_min = 25
                         , rad_max = 16000
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
                             , rad_eval = 16000
                             , box_size = 32000
                             , pt_min = 25
                             , rad_max = 16000
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
                             , rad_eval = 16000
                             , box_size = 32000
                             , pt_min = 25
                             , rad_max = 16000
                             , k_glob = 3 
                             , k_loc = 3
                             , type = 7)

# running global outlier function----
cd_glob <- cornell_data
cd_glob <- outlier_iden(X=cd_glob
                             , algo = 1
                             , outcri = 3
                             , pt_eval = 25
                             , rad_eval = 16000
                             , box_size = 32000
                             , pt_min = 25
                             , rad_max = 16000
                             , k_glob = 3 
                             , k_loc = 3  
                             , type = 7)

# running local outlier functions----
# gridded
cd_loc_grid <- cornell_data
cd_loc_grid <- outlier_iden(X=cd_loc_grid
                        , algo = 3
                        , outcri = 1
                        , pt_eval = 25
                        , rad_eval = 16000
                        , box_size = 32000
                        , pt_min = 25
                        , rad_max = 16000
                        , k_glob = 3 
                        , k_loc = 3  
                        , type = 7)

cd_loc_rad <- cornell_data
cd_loc_rad <- outlier_iden(X=cd_loc_rad
                            , algo = 2
                            , outcri = 1
                            , pt_eval = 25
                            , rad_eval = 16000
                            , box_size = 32000
                            , pt_min = 25
                            , rad_max = 16000
                            , k_glob = 3 
                            , k_loc = 3  
                            , type = 7)

cd_loc_pt <- cornell_data
cd_loc_pt <- outlier_iden(X=cd_loc_pt
                           , algo = 1
                           , outcri = 1
                           , pt_eval = 25
                           , rad_eval = 16000
                           , box_size = 32000
                           , pt_min = 25
                           , rad_max = 16000
                           , k_glob = 3 
                           , k_loc = 3  
                           , type = 7)

# running local and global outlier functions----
cd_loc_glob_grid <- cornell_data
cd_loc_glob_grid <- outlier_iden(X=cd_loc_glob_grid
                            , algo = 3
                            , outcri = 2
                            , pt_eval = 25
                            , rad_eval = 16000
                            , box_size = 32000
                            , pt_min = 25
                            , rad_max = 16000
                            , k_glob = 3 
                            , k_loc = 3  
                            , type = 7)

cd_loc_glob_rad <- cornell_data
cd_loc_glob_rad <- outlier_iden(X=cd_loc_glob_rad
                           , algo = 2
                           , outcri = 2
                           , pt_eval = 25
                           , rad_eval = 16000
                           , box_size = 32000
                           , pt_min = 25
                           , rad_max = 16000
                           , k_glob = 3 
                           , k_loc = 3  
                           , type = 7)

cd_loc_glob_pt <- cornell_data
cd_loc_glob_pt <- outlier_iden(X=cd_loc_glob_pt
                          , algo = 1
                          , outcri = 2
                          , pt_eval = 25
                          , rad_eval = 16000
                          , box_size = 32000
                          , pt_min = 25
                          , rad_max = 16000
                          , k_glob = 3 
                          , k_loc = 3  
                          , type = 7)

# running the wrapper function select_out_algo----
cd_wrap = cornell_data
# Making the names of the input data fields different from the defaults
colnames(cd_wrap[which(colnames(cd_wrap) == 'test')]) = 'NewName'
colnames(cd_wrap[which(colnames(cd_wrap) == 'x_coord')]) = 'XName'
colnames(cd_wrap[which(colnames(cd_wrap) == 'y_coord')]) = 'YName'
wrap_loc_pt <- select_out_algo(Data = cd_wrap, 
                  OutVarName = 'TESTED', InVarName = 'NewName', 
                  X_coordName = 'XName', Y_coordName = 'YName', 
                  CensorName = 'calc_depth', 
                  Threshold = 0, 
                  algo = 1, outcri = 1, 
                  pt_eval = 25, rad_eval = 16000, 
                  box_size = 32000, 
                  pt_min = 25, rad_max = 16000,
                  k_glob = 3, k_loc = 3, 
                  type = 7, 
                  min_val = 0, max_val = 7000, rank = FALSE)

assert_that(colnames(wrap_loc_pt$NotOutliers[which(colnames(cd_wrap) == 'test')]) == 'TESTED')
assert_that(colnames(wrap_loc_pt$NotOutliers[which(colnames(cd_wrap) == 'x_coord')]) == 'XName')
assert_that(colnames(wrap_loc_pt$NotOutliers[which(colnames(cd_wrap) == 'y_coord')]) == 'YName')
assert_that(colnames(wrap_loc_pt$NotOutliers[which(colnames(cd_wrap) == 'calc_depth')]) == 'calc_depth')
assert_that(nrow(wrap_loc_pt$NotOutliers) == nrow(cd_loc_pt[which(cd_loc_pt$outs == 0),]))
assert_that(nrow(wrap_loc_pt$Outliers) == nrow(cd_loc_pt[which(cd_loc_pt$outs == 1),]))

# sensitivity analysis for local points algorithm----
pts_sens <- c(10,25,50,100,200) # number of points for local neighborhood
rad_sens <- c(4000,8000,16000,32000,64000) # maximum size of radius

outs_iden <- matrix(0,length(pts_sens),length(rad_sens)) # matrix to hold number of outliers
outs_sparse <- matrix(0,length(pts_sens),length(rad_sens)) # number of points in sparse areas

#Calculate Outliers
for(i in 1:length(pts_sens)){
  for(j in 1:length(rad_sens)){
    sens_data <- cornell_data
    
    sens_data2 <- outlier_iden(X=sens_data
                               , algo = 1
                               , outcri = 1
                               , pt_eval = pts_sens[i]
                               , rad_eval = 16000
                               , box_size = 32000
                               , pt_min = 25
                               , rad_max = rad_sens[j]
                               , k_glob = 3 
                               , k_loc = 3  
                               , type = 7)
    
    outs_iden[i,j] <- sum(sens_data2$outs)
    outs_sparse[i,j] <- sum(sens_data2$out_loc_error)
  }
}

# creating plot
setEPS()
postscript(file = "outlier_sens.eps", title = "Sensitivity Outliers Local Points", width = 5, height = 5)

#Make color ramp
Pal = colorRampPalette(c('red', 'orange', 'yellow', 'green', 'blue', 'purple'))
cols <- Pal(max(outs_iden)+1)

#Set plotting margins
par(mar =c(3,3,0,9)+0.1)

data <- matrix(0,length(pts_sens)*length(rad_sens),4)

for(i in 1:length(rad_sens)){
  inds <- seq((i-1)*length(rad_sens)+1,i*length(rad_sens), by=1)
  
  data[c(inds),1] <- outs_iden[,i]
  data[c(inds),2] <- outs_sparse[,i]
  data[c(inds),3] <- pts_sens
  data[c(inds),4] <- rad_sens[i]/1000 #Convert to km for plotting position
}

dataplot <- data.frame(data)

colnames(dataplot) <- c("iden", "sparse", "pts", "rad")

#Changing the plotting location for pts = 10 to 12.5 for equal spacing in log base 2
dataplot$pts[dataplot$pts == 10] = 12.5

#Assigning color and size of points
dataplot$cols <- cols[dataplot$iden+1]
dataplot$cex <- sqrt(nrow(cornell_data) - dataplot$sparse)/30

plot(dataplot$pts
     , log(dataplot$rad, base=2)
     , log='x'
     , cex = 1.02*sqrt(nrow(cornell_data))/30
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
       , cex = 0.98*sqrt(nrow(cornell_data))/30
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
#Saving commented out legends for size of points legend entries
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

legend(x = 300
       , y = 6 # location
       , legend=c("# Outs/ # All (%)", seq(0,8,1)) # legend entries
       , pch = c(NA, rep(19,9))
       , col = c(NA, cols[round(max(outs_iden)/(max(outs_iden)/nrow(cornell_data))*seq(0,0.08,0.01),0) + 1])
       , ncol = 1
)
text(x=300
     , y = 0.75 + c(2.5,2.1, 1.7, 1.1)
     , labels=c("Point size is ", "  proportional to %", "  of data tested.", "Black circle is 100%.")
     , adj = c(0,0)
)

dev.off()
