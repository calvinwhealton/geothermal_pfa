# script to generate interpolation tables
# interpolation is of the play fairway metric variacne
# as a function of two distribution parameters for the risk factor
# Real and log space variances of the metrics are calculated.
# Log space is denoted by "_ls" in the variable name.

# Load library for writing out *.xlsx files.
library(xlsx)

# settting working directory
setwd('C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa')

# location of error interoplation tables
wd_error_interp <- paste(getwd(), '/Rasters/Rasters_in/error_interp_tabs', sep='') #Output from make_interp_tables.R

### for seismic worm angle to stress ----
mean_seSt <- seq(0,180,by=2) # range of means.
std_seSt <- seq(0,310,by=10) # range of standard deviations
nMC_seSt <- 100000 # number of monte carlo replicates

set.seed(10) # setting seed

# initializing matrix
seSt_var3 <- matrix(0,length(mean_seSt),length(std_seSt))
seSt_var5 <- matrix(0,length(mean_seSt),length(std_seSt))
seSt_var3_ls <- matrix(0,length(mean_seSt),length(std_seSt))
seSt_var5_ls <- matrix(0,length(mean_seSt),length(std_seSt))

# thresholds
seis_stress_min <- 0.001 # to avoid problems with numerically zero values
seis_stress_max <- 25
seis_stress_thresh3<- c(seis_stress_min,c(8,16),seis_stress_max)
seis_stress_thresh5<- c(seis_stress_min,c(5,10,15,20),seis_stress_max)

# set critical angles:
critical_ang1 = 65.2
critical_ang2 = 114.8

# loops to fill in matrix
for(i in 1:length(mean_seSt)){
  for(j in 1:length(std_seSt)){
    
    # creating random values
    # normal approximaton is wrong, but using von Mises was not working
    rand <- rnorm(nMC_seSt,mean_seSt[i],std_seSt[j])
    
    # Old method for RiskAng:
    # rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    # rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    
    # New Method:
    # Find the randomly generated values that are negative and add 360 degrees to them until they are all positive.
    while (length(rand[which(rand < 0)]) != 0){
      rand[which(rand < 0)] = rand[which(rand < 0)] + 360
    }
    # Now convert all values to [0,180]
    while (length(rand[which(rand > 180)]) != 0){
      rand[which(rand > 180)] = rand[which(rand > 180)] - 180
    }
    # Now take all angles and convert them to a risk angle
    a1 = abs(rand - critical_ang1)
    a2 = abs(rand - critical_ang2)
    # Bind them
    b = rbind(a1, a2)
    # Assign the minimum value (most risky) to those points
    rand = apply(b, 2, min)
    
    
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_seSt)
    pfm3[rand < seis_stress_thresh3[1]] <- 3
    pfm3[intersect(which(rand >= seis_stress_thresh3[1]),which(rand < seis_stress_thresh3[2]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh3[1]),which(rand < seis_stress_thresh3[2]))] - seis_stress_thresh3[1])/(seis_stress_thresh3[2]-seis_stress_thresh3[1])
    pfm3[intersect(which(rand >= seis_stress_thresh3[2]),which(rand < seis_stress_thresh3[3]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh3[2]),which(rand < seis_stress_thresh3[3]))] - seis_stress_thresh3[2])/(seis_stress_thresh3[3]-seis_stress_thresh3[2])
    pfm3[intersect(which(rand >= seis_stress_thresh3[3]),which(rand < seis_stress_thresh3[4]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh3[3]),which(rand < seis_stress_thresh3[4]))] - seis_stress_thresh3[3])/(seis_stress_thresh3[4]-seis_stress_thresh3[3])
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_seSt)
    pfm5[rand < seis_stress_thresh5[1]] <- 5
    pfm5[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_stress_thresh5[1]),which(rand < seis_stress_thresh5[2]))] - seis_stress_thresh5[1])/(seis_stress_thresh5[2]-seis_stress_thresh5[1])
    pfm5[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_stress_thresh5[2]),which(rand < seis_stress_thresh5[3]))] - seis_stress_thresh5[2])/(seis_stress_thresh5[3]-seis_stress_thresh5[2])
    pfm5[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_stress_thresh5[3]),which(rand < seis_stress_thresh5[4]))] - seis_stress_thresh5[3])/(seis_stress_thresh5[4]-seis_stress_thresh5[3])
    pfm5[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_stress_thresh5[4]),which(rand < seis_stress_thresh5[5]))] - seis_stress_thresh5[4])/(seis_stress_thresh5[5]-seis_stress_thresh5[4])
    pfm5[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_stress_thresh5[5]),which(rand < seis_stress_thresh5[6]))] - seis_stress_thresh5[5])/(seis_stress_thresh5[6]-seis_stress_thresh5[5])
    
    pfm5 <- 5- pfm5
    pfm3 <- 3 - pfm3
    
    # adding values to matrix
    seSt_var3[i,j] <- var(pfm3)
    seSt_var5[i,j] <- var(pfm5)
    
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    seSt_var3_ls[i,j] <- var(log(pfm3))
    seSt_var5_ls[i,j] <- var(log(pfm5))
    
  }
}
rm(a1, a2, b, i, j, rand, pfm3, pfm5)

### for seismic distance to earthquake ----
mean_seEq <- seq(0,26000,by=200) # range of means
std_seEq <- seq(0,2600,by=100) # range of standard deviations
nMC_seEq <- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrix
seEq_var3 <- matrix(0,length(mean_seEq),length(std_seEq))
seEq_var5 <- matrix(0,length(mean_seEq),length(std_seEq))
seEq_var3_ls <- matrix(0,length(mean_seEq),length(std_seEq))
seEq_var5_ls <- matrix(0,length(mean_seEq),length(std_seEq))

# thresholds
seis_eq_min <- 0.001 # to avoid problems with numerically zero values
seis_eq_max <- 25
seis_eq_thresh3<- 10^3*c(seis_eq_min,c(8,16),seis_eq_max)
seis_eq_thresh5<- 10^3*c(seis_eq_min,c(5,10,15,20),seis_eq_max)


# loops to fill-in matrix
for(i in 1:length(mean_seEq)){
  for(j in 1:length(std_seEq)){
    
    # creating random values
    # normal approximaton is wrong, but negative value not not reasonable
    rand <- abs(rnorm(nMC_seEq,mean_seEq[i],std_seEq[j]))
        
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_seEq)
    pfm3[rand < seis_eq_thresh3[1]] <-3
    pfm3[rand > seis_eq_thresh3[4]] <- 0
    pfm3[intersect(which(rand >= seis_eq_thresh3[1]),which(rand < seis_eq_thresh3[2]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh3[1]),which(rand < seis_eq_thresh3[2]))] - seis_eq_thresh3[1])/(seis_eq_thresh3[2]-seis_eq_thresh3[1])
    pfm3[intersect(which(rand >= seis_eq_thresh3[2]),which(rand < seis_eq_thresh3[3]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh3[2]),which(rand < seis_eq_thresh3[3]))] - seis_eq_thresh3[2])/(seis_eq_thresh3[3]-seis_eq_thresh3[2])
    pfm3[intersect(which(rand >= seis_eq_thresh3[3]),which(rand < seis_eq_thresh3[4]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh3[3]),which(rand < seis_eq_thresh3[4]))] - seis_eq_thresh3[3])/(seis_eq_thresh3[4]-seis_eq_thresh3[3])
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_seEq)
    pfm5[rand < seis_eq_thresh5[1]] <-5
    pfm5[rand > seis_eq_thresh5[6]] <-0
    pfm5[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] <- 5- (rand[intersect(which(rand >= seis_eq_thresh5[1]),which(rand < seis_eq_thresh5[2]))] - seis_eq_thresh5[1])/(seis_eq_thresh5[2]-seis_eq_thresh5[1])
    pfm5[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] <- 4- (rand[intersect(which(rand >= seis_eq_thresh5[2]),which(rand < seis_eq_thresh5[3]))] - seis_eq_thresh5[2])/(seis_eq_thresh5[3]-seis_eq_thresh5[2])
    pfm5[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] <- 3- (rand[intersect(which(rand >= seis_eq_thresh5[3]),which(rand < seis_eq_thresh5[4]))] - seis_eq_thresh5[3])/(seis_eq_thresh5[4]-seis_eq_thresh5[3])
    pfm5[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] <- 2- (rand[intersect(which(rand >= seis_eq_thresh5[4]),which(rand < seis_eq_thresh5[5]))] - seis_eq_thresh5[4])/(seis_eq_thresh5[5]-seis_eq_thresh5[4])
    pfm5[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] <- 1- (rand[intersect(which(rand >= seis_eq_thresh5[5]),which(rand < seis_eq_thresh5[6]))] - seis_eq_thresh5[5])/(seis_eq_thresh5[6]-seis_eq_thresh5[5])
    
    # Because the values are reversed, switch them here.
    pfm5 <- 5- pfm5
    pfm3 <- 3 - pfm3
    
    # adding values to matrix
    seEq_var3[i,j] <- var(pfm3)
    seEq_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    seEq_var3_ls[i,j] <- var(log(pfm3))
    seEq_var5_ls[i,j] <- var(log(pfm5))
    
  }
}
rm(i, j, rand, pfm3, pfm5)

### for reservoir ----
mean_re_rfc <- seq(-7,9.75,0.25) # range of means in log base e space
mean_re_RPIw <- seq(-14,3,0.25) # range of means in log base e space
mean_re_RPIg <- seq(-14,6.5,0.25) # range of means in log base e space
cv_re <- seq(0,0.4,by=0.05) # range of coefficient of variations
nMC_re <- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrix
re_rfc_var3 <- matrix(0,length(mean_re_rfc),length(cv_re))
re_rfc_var5 <- matrix(0,length(mean_re_rfc),length(cv_re))
re_rfc_var3_ls <- matrix(0,length(mean_re_rfc),length(cv_re))
re_rfc_var5_ls <- matrix(0,length(mean_re_rfc),length(cv_re))

re_RPIw_var3 <- matrix(0,length(mean_re_RPIw),length(cv_re))
re_RPIw_var5 <- matrix(0,length(mean_re_RPIw),length(cv_re))
re_RPIw_var3_ls <- matrix(0,length(mean_re_RPIw),length(cv_re))
re_RPIw_var5_ls <- matrix(0,length(mean_re_RPIw),length(cv_re))

re_RPIg_var3 <- matrix(0,length(mean_re_RPIg),length(cv_re))
re_RPIg_var5 <- matrix(0,length(mean_re_RPIg),length(cv_re))
re_RPIg_var3_ls <- matrix(0,length(mean_re_RPIg),length(cv_re))
re_RPIg_var5_ls <- matrix(0,length(mean_re_RPIg),length(cv_re))

# thresholds
res_rfc_min <- 0.001 # minimum for reservoir (both 3-color and 5-color)
res_rfc_max <- 10000 # maximum for reservoir (both 3-color and 5-color)
res_rfc_thresh3 <- c(res_rfc_min,c(10,100),res_rfc_max)
res_rfc_thresh5 <- c(res_rfc_min,c(1,10,100,1000),res_rfc_max)

res_RPIw_min <- 0.000001 # minimum for reservoir (both 3-color and 5-color)
res_RPIw_max <- 15    # maximum for reservoir (both 3-color and 5-color)
res_RPIw_thresh3 <- c(res_RPIw_min,c(0.1,1.0),res_RPIw_max)
res_RPIw_thresh5 <- c(res_RPIw_min,c(0.01,0.1,1,10),res_RPIw_max)

res_RPIg_min <- 0.000001 # minimum for reservoir (both 3-color and 5-color)
res_RPIg_max <- 15     # maximum for reservoir (both 3-color and 5-color)
res_RPIg_thresh3 <- c(res_RPIg_min,c(0.1,1.0),res_RPIg_max)
res_RPIg_thresh5 <- c(res_RPIg_min,c(0.01,0.1,1,10),res_RPIg_max)


# loops to fill-in matrix
for(i in 1:length(mean_re_rfc)){
  for(j in 1:length(cv_re)){
    
    # creating random values
    # using log-normal distribution
    sigma2 <- log(cv_re[j]^2 + 1)
    mu <- mean_re_rfc[i] - sigma2/2
    rand <- rnorm(nMC_re,mu,sqrt(sigma2))
    
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_re)
    pfm3[rand < log(res_rfc_thresh3[1])] <-3
    pfm3[rand > log(res_rfc_thresh3[4])] <- 0
    pfm3[intersect(which(rand >= log(res_rfc_thresh3[1])),which(rand < log(res_rfc_thresh3[2])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh3[1])),which(rand < log(res_rfc_thresh3[2])))] - log(res_rfc_thresh3[1]))/(log(res_rfc_thresh3[2])-log(res_rfc_thresh3[1]))
    pfm3[intersect(which(rand >= log(res_rfc_thresh3[2])),which(rand < log(res_rfc_thresh3[3])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh3[2])),which(rand < log(res_rfc_thresh3[3])))] - log(res_rfc_thresh3[2]))/(log(res_rfc_thresh3[3])-log(res_rfc_thresh3[2]))
    pfm3[intersect(which(rand >= log(res_rfc_thresh3[3])),which(rand < log(res_rfc_thresh3[4])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh3[3])),which(rand < log(res_rfc_thresh3[4])))] - log(res_rfc_thresh3[3]))/(log(res_rfc_thresh3[4])-log(res_rfc_thresh3[3]))
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_re)
    pfm5[rand < log(res_rfc_thresh5[1])] <-5
    pfm5[rand > log(res_rfc_thresh5[6])] <- 0
    pfm5[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_rfc_thresh5[1])),which(rand < log(res_rfc_thresh5[2])))] - log(res_rfc_thresh5[1]))/(log(res_rfc_thresh5[2])-log(res_rfc_thresh5[1]))
    pfm5[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_rfc_thresh5[2])),which(rand < log(res_rfc_thresh5[3])))] - log(res_rfc_thresh5[2]))/(log(res_rfc_thresh5[3])-log(res_rfc_thresh5[2]))
    pfm5[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_rfc_thresh5[3])),which(rand < log(res_rfc_thresh5[4])))] - log(res_rfc_thresh5[3]))/(log(res_rfc_thresh5[4])-log(res_rfc_thresh5[3]))
    pfm5[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_rfc_thresh5[4])),which(rand < log(res_rfc_thresh5[5])))] - log(res_rfc_thresh5[4]))/(log(res_rfc_thresh5[5])-log(res_rfc_thresh5[4]))
    pfm5[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_rfc_thresh5[5])),which(rand < log(res_rfc_thresh5[6])))] - log(res_rfc_thresh5[5]))/(log(res_rfc_thresh5[6])-log(res_rfc_thresh5[5]))
    
    # Because the values are reversed when they are calculated. Switch back here.
    pfm5 <- 5- pfm5
    pfm3 <- 3 - pfm3
    
    # adding values to matrix
    re_rfc_var3[i,j] <- var(pfm3)
    re_rfc_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    re_rfc_var3_ls[i,j] <- var(log(pfm3))
    re_rfc_var5_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5, sigma2, mu)

for(i in 1:length(mean_re_RPIw)){
  for(j in 1:length(cv_re)){
    
    # creating random values
    # using log-normal distribution
    sigma2 <- log(cv_re[j]^2 + 1)
    mu <- mean_re_RPIw[i] - sigma2/2
    rand <- rnorm(nMC_re,mu,sqrt(sigma2))
    
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_re)
    pfm3[rand < log(res_RPIw_thresh3[1])] <-3
    pfm3[rand > log(res_RPIw_thresh3[4])] <- 0
    pfm3[intersect(which(rand >= log(res_RPIw_thresh3[1])),which(rand < log(res_RPIw_thresh3[2])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh3[1])),which(rand < log(res_RPIw_thresh3[2])))] - log(res_RPIw_thresh3[1]))/(log(res_RPIw_thresh3[2])-log(res_RPIw_thresh3[1]))
    pfm3[intersect(which(rand >= log(res_RPIw_thresh3[2])),which(rand < log(res_RPIw_thresh3[3])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh3[2])),which(rand < log(res_RPIw_thresh3[3])))] - log(res_RPIw_thresh3[2]))/(log(res_RPIw_thresh3[3])-log(res_RPIw_thresh3[2]))
    pfm3[intersect(which(rand >= log(res_RPIw_thresh3[3])),which(rand < log(res_RPIw_thresh3[4])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh3[3])),which(rand < log(res_RPIw_thresh3[4])))] - log(res_RPIw_thresh3[3]))/(log(res_RPIw_thresh3[4])-log(res_RPIw_thresh3[3]))
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_re)
    pfm5[rand < log(res_RPIw_thresh5[1])] <-5
    pfm5[rand > log(res_RPIw_thresh5[6])] <- 0
    pfm5[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIw_thresh5[1])),which(rand < log(res_RPIw_thresh5[2])))] - log(res_RPIw_thresh5[1]))/(log(res_RPIw_thresh5[2])-log(res_RPIw_thresh5[1]))
    pfm5[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIw_thresh5[2])),which(rand < log(res_RPIw_thresh5[3])))] - log(res_RPIw_thresh5[2]))/(log(res_RPIw_thresh5[3])-log(res_RPIw_thresh5[2]))
    pfm5[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIw_thresh5[3])),which(rand < log(res_RPIw_thresh5[4])))] - log(res_RPIw_thresh5[3]))/(log(res_RPIw_thresh5[4])-log(res_RPIw_thresh5[3]))
    pfm5[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIw_thresh5[4])),which(rand < log(res_RPIw_thresh5[5])))] - log(res_RPIw_thresh5[4]))/(log(res_RPIw_thresh5[5])-log(res_RPIw_thresh5[4]))
    pfm5[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIw_thresh5[5])),which(rand < log(res_RPIw_thresh5[6])))] - log(res_RPIw_thresh5[5]))/(log(res_RPIw_thresh5[6])-log(res_RPIw_thresh5[5]))
    
    # Because the values are reversed when they are calculated. Switch back here.
    pfm5 <- 5- pfm5
    pfm3 <- 3 - pfm3
    
    # adding values to matrix
    re_RPIw_var3[i,j] <- var(pfm3)
    re_RPIw_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    re_RPIw_var3_ls[i,j] <- var(log(pfm3))
    re_RPIw_var5_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5, sigma2, mu)

for(i in 1:length(mean_re_RPIg)){
  for(j in 1:length(cv_re)){
    
    # creating random values
    # using log-normal distribution
    sigma2 <- log(cv_re[j]^2 + 1)
    mu <- mean_re_RPIg[i] - sigma2/2
    rand <- rnorm(nMC_re,mu,sqrt(sigma2))
    
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_re)
    pfm3[rand < log(res_RPIg_thresh3[1])] <-3
    pfm3[rand > log(res_RPIg_thresh3[4])] <- 0
    pfm3[intersect(which(rand >= log(res_RPIg_thresh3[1])),which(rand < log(res_RPIg_thresh3[2])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh3[1])),which(rand < log(res_RPIg_thresh3[2])))] - log(res_RPIg_thresh3[1]))/(log(res_RPIg_thresh3[2])-log(res_RPIg_thresh3[1]))
    pfm3[intersect(which(rand >= log(res_RPIg_thresh3[2])),which(rand < log(res_RPIg_thresh3[3])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh3[2])),which(rand < log(res_RPIg_thresh3[3])))] - log(res_RPIg_thresh3[2]))/(log(res_RPIg_thresh3[3])-log(res_RPIg_thresh3[2]))
    pfm3[intersect(which(rand >= log(res_RPIg_thresh3[3])),which(rand < log(res_RPIg_thresh3[4])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh3[3])),which(rand < log(res_RPIg_thresh3[4])))] - log(res_RPIg_thresh3[3]))/(log(res_RPIg_thresh3[4])-log(res_RPIg_thresh3[3]))
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_re)
    pfm5[rand < log(res_RPIg_thresh5[1])] <-5
    pfm5[rand > log(res_RPIg_thresh5[6])] <- 0
    pfm5[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] <- 5- (rand[intersect(which(rand >= log(res_RPIg_thresh5[1])),which(rand < log(res_RPIg_thresh5[2])))] - log(res_RPIg_thresh5[1]))/(log(res_RPIg_thresh5[2])-log(res_RPIg_thresh5[1]))
    pfm5[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] <- 4- (rand[intersect(which(rand >= log(res_RPIg_thresh5[2])),which(rand < log(res_RPIg_thresh5[3])))] - log(res_RPIg_thresh5[2]))/(log(res_RPIg_thresh5[3])-log(res_RPIg_thresh5[2]))
    pfm5[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] <- 3- (rand[intersect(which(rand >= log(res_RPIg_thresh5[3])),which(rand < log(res_RPIg_thresh5[4])))] - log(res_RPIg_thresh5[3]))/(log(res_RPIg_thresh5[4])-log(res_RPIg_thresh5[3]))
    pfm5[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] <- 2- (rand[intersect(which(rand >= log(res_RPIg_thresh5[4])),which(rand < log(res_RPIg_thresh5[5])))] - log(res_RPIg_thresh5[4]))/(log(res_RPIg_thresh5[5])-log(res_RPIg_thresh5[4]))
    pfm5[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] <- 1- (rand[intersect(which(rand >= log(res_RPIg_thresh5[5])),which(rand < log(res_RPIg_thresh5[6])))] - log(res_RPIg_thresh5[5]))/(log(res_RPIg_thresh5[6])-log(res_RPIg_thresh5[5]))
    
    # Because the values are reversed when they are calculated. Switch back here.
    pfm5 <- 5- pfm5
    pfm3 <- 3 - pfm3
    
    # adding values to matrix
    re_RPIg_var3[i,j] <- var(pfm3)
    re_RPIg_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    re_RPIg_var3_ls[i,j] <- var(log(pfm3))
    re_RPIg_var5_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5, sigma2, mu)

### for thermal depth to 80 degC----
mean_thd80 <- seq(1690,6490,by=200) # range of means
std_thd80 <- seq(40,1840,by=50) # range of standard deviations
nMC_thd80<- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrices
thd80_var3 <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var5 <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var3_ls <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var5_ls <- matrix(0,length(mean_thd80),length(std_thd80))

thd80_var3_old <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var5_old <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var3_old_ls <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var5_old_ls <- matrix(0,length(mean_thd80),length(std_thd80))

# thresholds - new values
therm_thresh3 <- rev(c(5000,3750,2350,1000)) #Look into changing these. 5000 is < 15 C/km and $14.5 mil/well. 4000 is about $10 mil/well, 17 C/km
therm_thresh5 <- rev(c(5000,4000,3000,2500,2000,1000)) #3000 is about $6.4 mil/well and 23 C/km. 2500 is about $4.8 mil/well and 28C/km.

# thresholds - old values
therm_thresh3_old <- rev(c(8750,3000,2000,500))
therm_thresh5_old <- rev(c(8750,4000,3000,2300,1500,500))

# looping over to fill-in matrices
for(i in 1:length(mean_thd80)){
  for(j in 1:length(std_thd80)){
    
    # generating random values
    rand <- rnorm(nMC_thd80,mean_thd80[i],std_thd80[j])
      
    # play fairway 3
    pfm3 <- rep(0,nMC_thd80)
    pfm3[rand < therm_thresh3[1]] <- 3
    pfm3[rand > therm_thresh3[4]] <- 0
    pfm3[intersect(which(rand >= therm_thresh3[1]),which(rand < therm_thresh3[2]))] <- 3 - (rand[intersect(which(rand >= therm_thresh3[1]),which(rand < therm_thresh3[2]))] - therm_thresh3[1])/(therm_thresh3[2]-therm_thresh3[1])
    pfm3[intersect(which(rand >= therm_thresh3[2]),which(rand < therm_thresh3[3]))] <- 2- (rand[intersect(which(rand >= therm_thresh3[2]),which(rand < therm_thresh3[3]))] - therm_thresh3[2])/(therm_thresh3[3]-therm_thresh3[2])
    pfm3[intersect(which(rand >= therm_thresh3[3]),which(rand < therm_thresh3[4]))] <- 1- (rand[intersect(which(rand >= therm_thresh3[3]),which(rand < therm_thresh3[4]))] - therm_thresh3[3])/(therm_thresh3[4]-therm_thresh3[3])
    
    # play fairway 5
    pfm5 <- rep(0,nMC_thd80)
    pfm5[rand < therm_thresh5[1]] <- 5
    pfm5[rand > therm_thresh5[6]] <- 0
    pfm5[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5[1]),which(rand < therm_thresh5[2]))] - therm_thresh5[1])/(therm_thresh5[2]-therm_thresh5[1])
    pfm5[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5[2]),which(rand < therm_thresh5[3]))] - therm_thresh5[2])/(therm_thresh5[3]-therm_thresh5[2])
    pfm5[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5[3]),which(rand < therm_thresh5[4]))] - therm_thresh5[3])/(therm_thresh5[4]-therm_thresh5[3])
    pfm5[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5[4]),which(rand < therm_thresh5[5]))] - therm_thresh5[4])/(therm_thresh5[5]-therm_thresh5[4])
    pfm5[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5[5]),which(rand < therm_thresh5[6]))] - therm_thresh5[5])/(therm_thresh5[6]-therm_thresh5[5])

    thd80_var3[i,j] <- var(pfm3)
    thd80_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    thd80_var3_ls[i,j] <- var(log(pfm3))
    thd80_var5_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5)

# old Thresholds
for(i in 1:length(mean_thd80)){
  for(j in 1:length(std_thd80)){
    
    # generating random values
    rand <- rnorm(nMC_thd80,mean_thd80[i],std_thd80[j])
    
    # play fairway 3
    pfm3 <- rep(0,nMC_thd80)
    pfm3[rand < therm_thresh3_old[1]] <- 3
    pfm3[rand > therm_thresh3_old[4]] <- 0
    pfm3[intersect(which(rand >= therm_thresh3_old[1]),which(rand < therm_thresh3_old[2]))] <- 3 - (rand[intersect(which(rand >= therm_thresh3_old[1]),which(rand < therm_thresh3_old[2]))] - therm_thresh3_old[1])/(therm_thresh3_old[2]-therm_thresh3_old[1])
    pfm3[intersect(which(rand >= therm_thresh3_old[2]),which(rand < therm_thresh3_old[3]))] <- 2- (rand[intersect(which(rand >= therm_thresh3_old[2]),which(rand < therm_thresh3_old[3]))] - therm_thresh3_old[2])/(therm_thresh3_old[3]-therm_thresh3_old[2])
    pfm3[intersect(which(rand >= therm_thresh3_old[3]),which(rand < therm_thresh3_old[4]))] <- 1- (rand[intersect(which(rand >= therm_thresh3_old[3]),which(rand < therm_thresh3_old[4]))] - therm_thresh3_old[3])/(therm_thresh3_old[4]-therm_thresh3_old[3])
    
    # play fairway 5
    pfm5 <- rep(0,nMC_thd80)
    pfm5[rand < therm_thresh5_old[1]] <- 5
    pfm5[rand > therm_thresh5_old[6]] <- 0
    pfm5[intersect(which(rand >= therm_thresh5_old[1]),which(rand < therm_thresh5_old[2]))] <- 5 - (rand[intersect(which(rand >= therm_thresh5_old[1]),which(rand < therm_thresh5_old[2]))] - therm_thresh5_old[1])/(therm_thresh5_old[2]-therm_thresh5_old[1])
    pfm5[intersect(which(rand >= therm_thresh5_old[2]),which(rand < therm_thresh5_old[3]))] <- 4- (rand[intersect(which(rand >= therm_thresh5_old[2]),which(rand < therm_thresh5_old[3]))] - therm_thresh5_old[2])/(therm_thresh5_old[3]-therm_thresh5_old[2])
    pfm5[intersect(which(rand >= therm_thresh5_old[3]),which(rand < therm_thresh5_old[4]))] <- 3- (rand[intersect(which(rand >= therm_thresh5_old[3]),which(rand < therm_thresh5_old[4]))] - therm_thresh5_old[3])/(therm_thresh5_old[4]-therm_thresh5_old[3])
    pfm5[intersect(which(rand >= therm_thresh5_old[4]),which(rand < therm_thresh5_old[5]))] <- 2- (rand[intersect(which(rand >= therm_thresh5_old[4]),which(rand < therm_thresh5_old[5]))] - therm_thresh5_old[4])/(therm_thresh5_old[5]-therm_thresh5_old[4])
    pfm5[intersect(which(rand >= therm_thresh5_old[5]),which(rand < therm_thresh5_old[6]))] <- 1- (rand[intersect(which(rand >= therm_thresh5_old[5]),which(rand < therm_thresh5_old[6]))] - therm_thresh5_old[5])/(therm_thresh5_old[6]-therm_thresh5_old[5])
    
    thd80_var3_old[i,j] <- var(pfm3)
    thd80_var5_old[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    thd80_var3_old_ls[i,j] <- var(log(pfm3))
    thd80_var5_old_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5)

### for utilization----
mean_util <- c(seq(5,65,by=2),900,1100) # range of means
std_util_pct <- seq(1,10,by=0.5) # range of standard deviations
nMC_util<- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrices
util_var3 <- matrix(0,length(mean_util),length(std_util_pct))
util_var5 <- matrix(0,length(mean_util),length(std_util_pct))
util_var3_ls <- matrix(0,length(mean_util),length(std_util_pct))
util_var5_ls <- matrix(0,length(mean_util),length(std_util_pct))

# thresholds
util_thresh3 <- c(5,13.5,16,25)
util_thresh5 <- c(5,12,13.5,16,20,25)


# looping over to fill-in matrices
for(i in 1:length(mean_util)){
  for(j in 1:length(std_util_pct)){
    
    # generating random values
    rand <- rnorm(nMC_util,mean_util[i],mean_util[i]*std_util_pct[j]/100)
    
    # play fairway 3
    pfm3 <- rep(0,nMC_util)
    pfm3[rand < util_thresh3[1]] <- 3
    pfm3[rand > util_thresh3[4]] <- 0
    pfm3[intersect(which(rand >= util_thresh3[1]),which(rand < util_thresh3[2]))] <- 3 - (rand[intersect(which(rand >= util_thresh3[1]),which(rand < util_thresh3[2]))] - util_thresh3[1])/(util_thresh3[2]-util_thresh3[1])
    pfm3[intersect(which(rand >= util_thresh3[2]),which(rand < util_thresh3[3]))] <- 2 - (rand[intersect(which(rand >= util_thresh3[2]),which(rand < util_thresh3[3]))] - util_thresh3[2])/(util_thresh3[3]-util_thresh3[2])
    pfm3[intersect(which(rand >= util_thresh3[3]),which(rand < util_thresh3[4]))] <- 1 - (rand[intersect(which(rand >= util_thresh3[3]),which(rand < util_thresh3[4]))] - util_thresh3[3])/(util_thresh3[4]-util_thresh3[3])
    
    # play fairway 5
    pfm5 <- rep(0,nMC_util)
    pfm5[rand < util_thresh5[1]] <- 5
    pfm5[rand > util_thresh5[6]] <- 0
    pfm5[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] <- 5 - (rand[intersect(which(rand >= util_thresh5[1]),which(rand < util_thresh5[2]))] - util_thresh5[1])/(util_thresh5[2]-util_thresh5[1])
    pfm5[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] <- 4 - (rand[intersect(which(rand >= util_thresh5[2]),which(rand < util_thresh5[3]))] - util_thresh5[2])/(util_thresh5[3]-util_thresh5[2])
    pfm5[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] <- 3 - (rand[intersect(which(rand >= util_thresh5[3]),which(rand < util_thresh5[4]))] - util_thresh5[3])/(util_thresh5[4]-util_thresh5[3])
    pfm5[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] <- 2 - (rand[intersect(which(rand >= util_thresh5[4]),which(rand < util_thresh5[5]))] - util_thresh5[4])/(util_thresh5[5]-util_thresh5[4])
    pfm5[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] <- 1 - (rand[intersect(which(rand >= util_thresh5[5]),which(rand < util_thresh5[6]))] - util_thresh5[5])/(util_thresh5[6]-util_thresh5[5])
    
    util_var3[i,j] <- var(pfm3)
    util_var5[i,j] <- var(pfm5)
    
    # Log Space
    
    # Set small values to a constant before computing the variance of the logs.
    pfm5[which(pfm5 < 0.2)] <- 0.2
    pfm3[which(pfm3 < 0.2)] <- 0.2
    
    util_var3_ls[i,j] <- var(log(pfm3))
    util_var5_ls[i,j] <- var(log(pfm5))
  }
}
rm(i, j, rand, pfm3, pfm5)

# write output tables:
setwd(wd_error_interp)

write.xlsx(util_var3
           ,'ut_slcoh_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(util_var5
           ,'ut_slcoh_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(util_var3_ls
           ,'ut_slcoh_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(util_var5_ls
           ,'ut_slcoh_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_rfc_var3
           ,'re_rfc_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_rfc_var5
           ,'re_rfc_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_rfc_var3_ls
           ,'re_rfc_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_rfc_var5_ls
           ,'re_rfc_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_RPIw_var3
           ,'re_RPIw_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_RPIw_var5
           ,'re_RPIw_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_RPIw_var3_ls
           ,'re_RPIw_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_RPIw_var5_ls
           ,'re_RPIw_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_RPIg_var3
           ,'re_RPIg_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_RPIg_var5
           ,'re_RPIg_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(re_RPIg_var3_ls
           ,'re_RPIg_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(re_RPIg_var5_ls
           ,'re_RPIg_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(seSt_var3
           ,'seSt_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(seSt_var5
           ,'seSt_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(seSt_var3_ls
           ,'seSt_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(seSt_var5_ls
           ,'seSt_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(seEq_var3
           ,'seEq_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(seEq_var5
           ,'seEq_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(seEq_var3_ls
           ,'seEq_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(seEq_var5_ls
           ,'seEq_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(thd80_var3
           ,'th_pfvar3.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(thd80_var5
           ,'th_pfvar5.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(thd80_var3_ls
           ,'th_pfvar3_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(thd80_var5_ls
           ,'th_pfvar5_ls.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(thd80_var3_old
           ,'th_pfvar3_old.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(thd80_var5_old
           ,'th_pfvar5_old.xlsx'
           ,col.names=F
           ,row.names=F)

write.xlsx(thd80_var3_old_ls
           ,'th_pfvar3_old_ls.xlsx'
           ,col.names=F
           ,row.names=F)
write.xlsx(thd80_var5_old_ls
           ,'th_pfvar5_old_ls.xlsx'
           ,col.names=F
           ,row.names=F)