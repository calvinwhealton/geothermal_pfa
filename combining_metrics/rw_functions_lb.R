# function to solve Weibull distribution parameters
solveWeibull <- function(x # scale (lambda) and shape (k). 
){
  if(right != Inf){
    #Use a truncated Weibull distribution
    #lb is for the lower bound of the Weibull Dist.
    #Mean
    f1 <- (lb + x[1]/(1 - exp(-((right - lb)/x[1])^x[2]))*pgamma(((right - lb)/x[1])^x[2], (1 + 1/x[2]))*gamma(1 + 1/x[2])) - mom[1]
    #Variance
    f2 <- ((x[1]^2)/(1 - exp(-((right - lb)/x[1])^x[2]))*pgamma(((right - lb)/x[1])^x[2], (1 + 2/x[2]))*gamma(1+2/x[2]) - (x[1]^2)/((1 - exp(-((right - lb)/x[1])^x[2]))^2)*(pgamma(((right - lb)/x[1])^x[2], (1 + 1/x[2]))*gamma(1+1/x[2]))^2) - mom[2]    
  }else{
    #Use the regular Weibull distribution
    #lb is for the lower bound of the Weibull Dist.
    #Mean
    f1 <- (lb + x[1]*gamma(1 + 1/x[2])) - mom[1]
    #Variance
    f2 <- ((x[1]^2)*(gamma(1+2/x[2]) - (gamma(1+1/x[2]))^2)) - mom[2]
  }
  
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
  if (right != Inf){
    #Mean for truncated Weibull
    f1 <- (lb + x[1]/(1 - exp(-((right - lb)/x[1])^k))*pgamma(((right - lb)/x[1])^k, (1 + 1/k))*gamma(1 + 1/k)) - mom[1]
  }else{
    #Mean
    f1 <- (lb + x*gamma(1 + 1/k)) - mom[1]
  }
  return(c(f1))
}

solveQuant <- function(x # value of the minimum
){
  if (right != Inf){
    #Use truncated Weibull
    #Summation of 4 variables 
    ssss <- 0
    for(j in ind_params){
      ssss <- ssss + ifelse(is.na(params[(j-1)*2+1]),0,log(1 - pweibull((x - lb_mat[i,j]), params[2*j], params[(j-1)*2+1])/pweibull((right - lb_mat[i,j]), params[2*j], params[(j-1)*2+1])))
    }
    
    #ps is defined outside of the function. This is the 100p percentile.
    f <- log(1 - ps) - ssss
    
  }else{
    #Summation of 4 variables 
    ssss <- 0
    for(j in ind_params){
      ssss <- ssss + ifelse(is.na(params[(j-1)*2+1]),0,((x - lb_mat[i,j])/params[(j-1)*2+1])^params[2*j])
    }
    
    #ps is defined outside of the function. This is the 100p percentile.
    f <- log(1-ps) + ssss
  }
  return(f)
}

# function to solve the generalized Beta distribution parameters
# assuming that the upper and lower bounds are known.
solveBeta <- function(x # alpha and beta 
){
  #lb and ub define the lower and upper bounds of the generalized Beta distribution.
  #Using the method of moments
  #Mean
  f1 <- ((ub - lb)*(x[1]/(x[1] + x[2])) + lb) - mom[1]
  #Variance
  f2 <- (x[1]*x[2]*(ub - lb)^2)/(((x[1] + x[2])^2)*(x[1] + x[2] + 1)) - mom[2]
  
  return(c(f1,f2))
}

solveBeta_beta <- function(x # alpha 
){
  #The lb and ub matrices define the lower and upper bounds of the generalized Beta distribution.
  #Using the method of moments
  #Mean
  f1 <- (ub*alpha - lb*x[1])/(alpha + x[1]) - mom[1]
  
  return(c(f1))
}

solveBeta_alpha <- function(x # beta 
){
  #The lb and ub matrices define the lower and upper bounds of the generalized Beta distribution.
  #Using the method of moments
  #Mean
  f1 <- (ub*x[1] - lb*beta)/(x[1] + beta) - mom[1]
  
  return(c(f1))
}

#Library for the 4 parameter beta CDF
library(ExtDist)

solveQuant_Beta <- function(x #value of the minimum
){
  #Multiplication of the CDFs for each variable 
  mmmm <- 1
  for(j in ind_params){
    mmmm <- mmmm * ifelse(is.na(params[(j-1)*2+1]), 1, (1 - pBeta_ab(x, shape1 = params[(j-1)*2+1], shape2 = params[j*2], a = lb_mat[i,j], b = ub_mat[i,j])))
  }
  
  #ps is defined outside of the function. This is the 100p percentile.
  f <- (1-ps) - mmmm
  return(f)
}