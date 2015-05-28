# function to generate random data for conditional FE multinomial logit model
# inputs are:
#   n     = number of individuals
#   t     = number of time periods
#   coeff = coefficients for each non-base state [each row is a state]
#   msd   = mean and standard devaition of data
#   fe    = fixed-effect bounds (+/- fe), assumed uniform

# outputs are:
#   x = initial matrix of data for all individuals (with and without changes)
#   y = initial vector of states selected
#   xChange  = independent variables for individuals w/ a state change
#   yChange  = state variables for individual w/ a state change
#   idt   = IDs and times kept for analysis
#   nkept = number of individuals kept because they changed state
#   ndrop = number of individuals dropped because they did not change state
#   util = calculated utility for each state
#   error = error terms for each state

genRandData <- function(n,t,coeff,msd,fe,seed){
  
  # setting the seed
  set.seed(seed)
  
  # finding number of non-base states
  states <- nrow(coeff)
  
  # vector of state indices
  st_ind <- rep(0,states+1)
  
  # creating vector of state indices
  for(i in 1:states){
    st_ind[i+1] <- st_ind[i]*t + 1
  }
  
  # initializing matrix to hold ID-time indices
  idt <- matrix(0,n*t,2)
  
  # generatng ID-time indices
  for(i in 1:n){
    idt[(i-1)*t:i*t,1] <- i # ID column
    idt[(i-1)*t+1:i*t,2] <- 1:t
  }
  
  # number of independent variables
  vars <- ncol(coeff)
  
  # initializing matrix to hold independent variables
  x <- matrix(0,n*t,dim)
  
  # x values randomly generated, assumed normally distributed
  for(i in 1:dim){
    x[,i] <- rnorm(n,msd[1,i],msd[2,i])
  }
  
  # initializing matrix to hold fixed-effects
  u <- matrix(0,n,dim)
  
  # generating fixed-effect, uniform on interval [-fe, fe]
  for(i in 1:dim){
    u[,i] <- runif(n,-fe,fe)
  }
  
  # initializing matrix to hold error terms
  e <- matrix(0,n*t,dim)
  
  # generating Gumbel (Extreme Value Type I) error terms
  # for all states-time-individuals
  for(i in 1:dim){
    e[,i] <- -1*log(-1*log(runif(n*t)))
  }
  
  # initializing matrix to hold utility for each variable state
  util <- matrix(0,n*t,states)
  
  # calculating utility for each state-time-individual
  for(i in 1:states){
    util[,i] <- x %*% coeff[i,] + u[,i] + e[,i]
  }
  
  # initializing matrix to hold selected states
  y <- matrix(0,n*t,1)
  
  # assinging state based on utility
  for(i in 1:n*t){
    
    # base-state (y=0) is assigned 0 coefficients (0 relative utility)
    if(max(util[i,]) <= 0){
      y[i] <- 0
    }
    else{
      y[i] <- st_ind[which(util[i,] %in% max(util[i,]))+1]
    }  
  }
  
  # initializing vector to hold the indices kept
  ind_keep <- c()
  
  # cleaning the resulting data to eliminate non-changing individuals
  for(i in 1:n){
    if(unique(y[t*(i-1)+1:i*t]) != 1)
    ind_keep <- c(ind_keep,c(t*(i-1)+1:i*t))
  }
  
  # creating new variables for only those individuals that changed state
  xChange <- x[ind_keep,]
  yChange <- y[ind_keep]
  
  # calculations related to keeping/dropping variables
  num_kept <- length(ind_keep)/t
  num_drop <- n - num_kept
  
  # formattting calculations for output
  dataOutput <- list(  "x" = x
                     , "y" = y
                     , "xChange" = xChange
                     , "yChange" = yChange
                     , "num_kept" = num_kept
                     , "num_drop" = num_drop
                     , "util" = util
                     , "error" = e)
  
  # returning output of function
  return(dataOutput)
}