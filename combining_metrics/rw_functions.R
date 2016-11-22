# function to solve Weibull distribution parameters
solveWeibull <- function(x # scale (lambda) and shape (k)
){
  #Mean
  f1 <- x[1]*gamma(1 + 1/x[2]) - mom[1]
  #Variance
  f2 <- (x[1]^2)*(gamma(1+2/x[2]) -(gamma(1+1/x[2]))^2) - mom[2]
  return(c(f1,f2))
}

mom <- c(10,20)

library(rootSolve)

ss <- multiroot(solveWeibull,start=c(1,1))

solveQuant <- function(x # value of the minimum
){
  
  #Summation of 4 variables 
  ssss <- 0
  for(i in 1:(length(params)/2)){
    ssss <- ssss + ifelse(is.na(params[(i-1)*2+1]),0,(x/params[(i-1)*2+1])^params[2*i])
  }
  
  #ps is defined outside of the function. This is the 100ps percentile.
  f <- log(1-ps) + ssss
  return(f)
} 