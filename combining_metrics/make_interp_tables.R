# script to generate interpolation tables
# interpolation is of the play fairway metric variacne
# as a function of two distribution parameters for the risk factor

# settting working directory
setwd('/Users/calvinwhealton/GitHub/geothermal_pfa/error_interp_tabs')

### for seismic worm angle to stress ----
mean_seSt <- seq(0,72,by=3) # range of means
std_seSt <- seq(0,310,by=10) # range of standard deviations
nMC_seSt <- 100000 # number of monte carlo replicates

set.seed(10) # setting seed

# initializing matrix
seSt_var3 <- matrix(0,length(mean_seSt),length(std_seSt))
seSt_var5 <- matrix(0,length(mean_seSt),length(std_seSt))

# loops to fill-in matrix
for(i in 1:length(mean_seSt)){
  for(j in 1:length(std_seSt)){
    
    # creating random values
    # normal approximaton is wrong, but using von Mises was not working
    # abs value used because -1 is equivalent to +1
    rand <- abs(rnorm(nMC_seSt,mean_seSt[i],std_seSt[j]))
    
    # max distance to angle is 65.2 deg
    # process to convert everything to interval [0, 65.2]
    # is to take values that are on (65.2,2*65.2) and reverse them
    # and them put them on the scale of [0, 65.2]
    # the reason is that being close to 2*65.2 is being close to another failure (0)
    # the values greater than 2*65.2 have 2*65.2 subtracted and the new values have the previous
    # step performed on them
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    rand[intersect(which(rand > 65.2),which(rand < 2*65.2))] <- 2*65.2 - rand[intersect(which(rand > 65.2),which(rand < 2*65.2))]
    rand[which(rand > 2*65.2)] <- rand[which(rand > 2*65.2)] - 2*65.2
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_seSt)
    pfm3[rand < 0.01] <-3
    pfm3[intersect(which(rand >= 0.01),which(rand < 8))] <- 3- (rand[intersect(which(rand >= 0.01),which(rand < 8))] - 0.01)/(8-0.01)
    pfm3[intersect(which(rand >= 8),which(rand < 16))] <- 2- (rand[intersect(which(rand >= 8),which(rand < 16))] - 8)/8
    pfm3[intersect(which(rand >= 16),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 16),which(rand < 25))] - 16)/(25-16)
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_seSt)
    pfm5[rand < 0.01] <-5
    pfm5[intersect(which(rand >= 0.01),which(rand < 5))] <- 5- (rand[intersect(which(rand >= 0.01),which(rand < 5))] - 0.01)/(5-0.01)
    pfm5[intersect(which(rand >= 5),which(rand < 10))] <- 4- (rand[intersect(which(rand >= 5),which(rand < 10))] - 5)/5
    pfm5[intersect(which(rand >= 10),which(rand < 15))] <- 3- (rand[intersect(which(rand >= 10),which(rand < 15))] - 10)/5
    pfm5[intersect(which(rand >= 15),which(rand < 20))] <- 2- (rand[intersect(which(rand >= 15),which(rand < 20))] - 15)/5
    pfm5[intersect(which(rand >= 20),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 20),which(rand < 25))] - 20)/5
    
    # adding values to matrix
    seSt_var3[i,j] <- var(pfm3)
    seSt_var5[i,j] <- var(pfm5)
    
  }
}


### for seismic distance to earthquake ----
mean_seEq <- seq(0,26000,by=200) # range of means
std_seEq <- seq(0,2600,by=100) # range of standard deviations
nMC_seEq <- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrix
seEq_var3 <- matrix(0,length(mean_seEq),length(std_seEq))
seEq_var5 <- matrix(0,length(mean_seEq),length(std_seEq))

# loops to fill-in matrix
for(i in 1:length(mean_seEq)){
  for(j in 1:length(std_seEq)){
    
    # creating random values
    # normal approximaton is wrong, but negative value not not reasonable
    rand <- abs(rnorm(nMC_seEq,mean_seEq[i],std_seEq[j]))
        
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_seEq)
    pfm3[rand < 1] <-3
    pfm3[rand > 25000] <- 0
    pfm3[intersect(which(rand >= 1),which(rand < 8000))] <- 3- (rand[intersect(which(rand >= 1),which(rand < 8000))] - 1)/(8000-1)
    pfm3[intersect(which(rand >= 8000),which(rand < 16000))] <- 2- (rand[intersect(which(rand >= 8000),which(rand < 16000))] - 8000)/8000
    pfm3[intersect(which(rand >= 16000),which(rand < 25000))] <- 1- (rand[intersect(which(rand >= 16000),which(rand < 25000))] - 16000)/(25000-16000)
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_seSt)
    pfm5[rand < 1] <-5
    pfm5[rand > 25000] <-0
    pfm5[intersect(which(rand >= 1),which(rand < 5000))] <- 5- (rand[intersect(which(rand >= 1),which(rand < 5000))] - 1)/(5000-1)
    pfm5[intersect(which(rand >= 5000),which(rand < 10000))] <- 4- (rand[intersect(which(rand >= 5000),which(rand < 10000))] - 5000)/5000
    pfm5[intersect(which(rand >= 10000),which(rand < 15000))] <- 3- (rand[intersect(which(rand >= 10000),which(rand < 15000))] - 10000)/5000
    pfm5[intersect(which(rand >= 15000),which(rand < 20000))] <- 2- (rand[intersect(which(rand >= 15000),which(rand < 20000))] - 15000)/5000
    pfm5[intersect(which(rand >= 20000),which(rand < 25000))] <- 1- (rand[intersect(which(rand >= 20000),which(rand < 25000))] - 20000)/5000
    
    # adding values to matrix
    seEq_var3[i,j] <- var(pfm3)
    seEq_var5[i,j] <- var(pfm5)
    
  }
}


### for reservoir ----
mean_re <- seq(-9.5,6.25,0.25) # range of means in log-space
cv_re <- seq(0,2,by=0.1) # range of coefficient of variations
nMC_re <- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrix
re_var3 <- matrix(0,length(mean_re),length(cv_re))
re_var5 <- matrix(0,length(mean_re),length(cv_re))

# loops to fill-in matrix
for(i in 1:length(mean_re)){
  for(j in 1:length(cv_re)){
    
    # creating random values
    # using log-normal distribution
    sigma2 <- log(cv_re[j]^2 + 1)
    mu <- mean_re[i] - sigma2/2
    rand <- rnorm(nMC_re,mu,sqrt(sigma2))
    
    # calculating play fairway 3
    pfm3 <- rep(0,nMC_re)
    pfm3[rand < log(3*10^-5)] <-3
    pfm3[rand > log(301)] <- 0
    pfm3[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.1)))] <- 3- (rand[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.1)))] - log(3*10^-5))/(log(0.1)-log(3*10^-5))
    pfm3[intersect(which(rand >= log(0.1)),which(rand < log(1)))] <- 2- (rand[intersect(which(rand >= log(0.1)),which(rand < log(1)))] - log(0.1))/(log(1)-log(0.1))
    pfm3[intersect(which(rand >= log(1)),which(rand < log(301)))] <- 1- (rand[intersect(which(rand >= log(1)),which(rand < log(301)))] - log(1))/(log(301)-log(1))
    
    # calculating play fairway 5
    pfm5 <- rep(0,nMC_re)
    pfm5[rand < log(3*10^-5)] <-5
    pfm5[rand > log(301)] <- 0
    pfm5[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.01)))] <- 5- (rand[intersect(which(rand >= log(3*10^-5)),which(rand < log(0.01)))] - log(3*10^-5))/(log(0.01)-log(3*10^-5))
    pfm5[intersect(which(rand >= log(0.01)),which(rand < log(0.1)))] <- 4- (rand[intersect(which(rand >= log(0.01)),which(rand < log(0.1)))] - log(0.01))/(log(0.1)-log(0.01))
    pfm5[intersect(which(rand >= log(0.1)),which(rand < log(1)))] <- 3- (rand[intersect(which(rand >= log(0.1)),which(rand < log(1)))] - log(0.1))/(log(1)-log(0.1))
    pfm5[intersect(which(rand >= log(1)),which(rand < log(10)))] <- 2- (rand[intersect(which(rand >= log(1)),which(rand < log(10)))] - log(1))/(log(10)-log(1))
    pfm5[intersect(which(rand >= log(10)),which(rand < log(301)))] <- 1- (rand[intersect(which(rand >= log(10)),which(rand < log(301)))] - log(10))/(log(301)-log(10))
    
    # adding values to matrix
    re_var3[i,j] <- var(pfm3)
    re_var5[i,j] <- var(pfm5)
    
  }
}



### for thermal depth to 80 degC----
mean_thd80 <- seq(750,6350,by=200) # range of means
std_thd80 <- seq(40,1740,by=50) # range of standard deviations
nMC_thd80<- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrices
thd80_var3 <- matrix(0,length(mean_thd80),length(std_thd80))
thd80_var5 <- matrix(0,length(mean_thd80),length(std_thd80))

# looping over to fill-in matrices
for(i in 1:length(mean_thd80)){
  for(j in 1:length(std_thd80)){
    
    # generating random values
    rand <- rnorm(nMC_thd80,mean_thd80[i],std_thd80[j])
      
    # play fairway 3
    pfm3 <- rep(0,nMC_thd80)
    pfm3[rand < 500] <- 3
    pfm3[rand > 8750] <- 0
    pfm3[intersect(which(rand >= 500),which(rand < 2000))] <- 3 - (rand[intersect(which(rand >= 500),which(rand < 2000))] - 500)/(2000-500)
    pfm3[intersect(which(rand >= 2000),which(rand < 3000))] <- 2- (rand[intersect(which(rand >= 2000),which(rand < 3000))] - 2000)/(3000-2000)
    pfm3[intersect(which(rand >= 3000),which(rand < 8750))] <- 1- (rand[intersect(which(rand >= 3000),which(rand < 8750))] - 3000)/(8750-3000)
    
    # play fairway 5
    pfm5 <- rep(0,nMC_thd80)
    pfm5[rand < 500] <- 5
    pfm5[rand > 8750] <- 0
    pfm5[intersect(which(rand >= 500),which(rand < 1500))] <- 5 - (rand[intersect(which(rand >= 500),which(rand < 1500))] - 500)/(1500-500)
    pfm5[intersect(which(rand >= 1500),which(rand < 2300))] <- 4- (rand[intersect(which(rand >= 1500),which(rand < 2300))] - 1500)/(2300-1500)
    pfm5[intersect(which(rand >= 2300),which(rand < 3000))] <- 3- (rand[intersect(which(rand >= 2300),which(rand < 3000))] - 2300)/(3000-2300)
    pfm5[intersect(which(rand >= 3000),which(rand < 4000))] <- 2- (rand[intersect(which(rand >= 3000),which(rand < 4000))] - 3000)/(4000-3000)
    pfm5[intersect(which(rand >= 4000),which(rand < 8750))] <- 1- (rand[intersect(which(rand >= 4000),which(rand < 8750))] - 4000)/(8750-3000)

    thd80_var3[i,j] <- var(pfm3)
    thd80_var5[i,j] <- var(pfm5)
  }
}




### for utilization----
mean_util <- c(seq(5,65,by=2),900,1100) # range of means
std_util_pct <- seq(1,10,by=0.5) # range of standard deviations
nMC_util<- 100000 # number of monte carlo

set.seed(10) # setting seed

# initializing matrices
util_var3 <- matrix(0,length(mean_util),length(std_util_pct))
util_var5 <- matrix(0,length(mean_util),length(std_util_pct))

# looping over to fill-in matrices
for(i in 1:length(mean_util)){
  for(j in 1:length(std_util_pct)){
    
    # generating random values
    rand <- rnorm(nMC_util,mean_util[i],mean_util[i]*std_util_pct[j]/100)
    
    # play fairway 3
    pfm3 <- rep(0,nMC_util)
    pfm3[rand < 5] <- 3
    pfm3[rand > 25] <- 0
    pfm3[intersect(which(rand >= 5),which(rand < 13.5))] <- 3 - (rand[intersect(which(rand >= 5),which(rand < 13.5))] - 5)/(13.5-5)
    pfm3[intersect(which(rand >= 13.5),which(rand < 16))] <- 2- (rand[intersect(which(rand >= 13.5),which(rand < 16))] - 13.5)/(16-13.5)
    pfm3[intersect(which(rand >= 16),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 16),which(rand < 25))] - 16)/(25-16)
    
    # play fairway 5
    pfm5 <- rep(0,nMC_util)
    pfm5[rand < 5] <- 5
    pfm5[rand > 25] <- 0
    pfm5[intersect(which(rand >= 5),which(rand < 12))] <- 5 - (rand[intersect(which(rand >= 5),which(rand < 12))] - 5)/(12-5)
    pfm5[intersect(which(rand >= 12),which(rand < 13.5))] <- 4- (rand[intersect(which(rand >= 12),which(rand < 13.5))] - 12)/(13.5-12)
    pfm5[intersect(which(rand >= 13.5),which(rand < 16))] <- 3- (rand[intersect(which(rand >= 13.5),which(rand < 16))] - 13.5)/(16-13.5)
    pfm5[intersect(which(rand >= 16),which(rand < 20))] <- 2- (rand[intersect(which(rand >= 16),which(rand < 20))] - 16)/(20-16)
    pfm5[intersect(which(rand >= 20),which(rand < 25))] <- 1- (rand[intersect(which(rand >= 20),which(rand < 25))] - 20)/(25-20)
    
    util_var3[i,j] <- var(pfm3)
    util_var5[i,j] <- var(pfm5)
  }
}
setwd('/Users/calvinwhealton/GitHub/geothermal_pfa/error_interp_tabs')
write.xlsx(pfm3,'ut_slcoh3.xlsx')
write.xlsx(pfm5,'ut_slcoh5.xlsx')