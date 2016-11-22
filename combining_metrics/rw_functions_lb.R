# function to solve Weibull distribution parameters
solveWeibull <- function(x # scale (lambda) and shape (k). 
){
  #lb is for the lower bound of the Weibull Dist.
  #Mean
  f1 <- (lb + x[1]*gamma(1 + 1/x[2])) - mom[1]
  #Variance
  f2 <- (x[1]^2)*(gamma(1+2/x[2]) -(gamma(1+1/x[2]))^2) - mom[2]
  
  return(c(f1,f2))
}

# function to solve Weibull distribution parameters
solveWeibull_k <- function(x # shape (k). 
){
  #lb is for the lower bound of the Weibull Dist.
  #Mean
  f1 <- (lb + u*gamma(1 + 1/x)) - mom[1]
  
  return(c(f1))
}

# function to solve Weibull distribution parameters
solveWeibull_u <- function(x # scale (lambda). 
){
  #lb is for the lower bound of the Weibull Dist.
  #Mean
  f1 <- (lb + x*gamma(1 + 1/k)) - mom[1]

  return(c(f1))
}

solveQuant <- function(x # value of the minimum
){
  
  #Summation of 4 variables 
  ssss <- 0
  for(j in ind_params){
    ssss <- ssss + ifelse(is.na(params[(j-1)*2+1]),0,((x - lb_mat[i,j])/params[(j-1)*2+1])^params[2*j])
  }
  
  #ps is defined outside of the function. This is the 100p percentile.
  f <- log(1-ps) + ssss
  return(f)
} 