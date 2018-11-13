# script to calculate bias in "mean" values with asymmetric boxplot outlier detection
# different values of algorithm k parameter used
# lognormal used

# values of k in boxplot detection algorithm
k_vals <- seq(0.5,5,0.1)

# matrices to store mean of censored and uncensored values
# 10,000 replicates for each value of k
means <- matrix(NA,nrow=10000,ncol=length(k_vals))
means2 <- matrix(NA,nrow=10000,ncol=length(k_vals))

# loop over k
for(i in 1:length(k_vals)){
  
  # loop over replicates
  for(j in 1:10000){
    # generating sample
    # lognormal, mu=0, sigma=1, => mode = 1 and mean = exp(0+1/2)
    x1 <- exp(rnorm(50))
    # lower and upper bounds
    lb <- quantile(x1,0.5)-k_vals[i]*(quantile(x1,0.5)-quantile(x1,0.25))
    ub <- quantile(x1,0.5)+k_vals[i]*(quantile(x1,0.75)-quantile(x1,0.5))
    
    # calculation, multiplying TRUE/FALSE gives 0 or 1
    means[j,i] <- mean(x1[which((x1 > lb)*(x1 < ub)==1)])
    means2[j,i] <- mean(x1)
  }
}

# making plots of results
plot(x=k_vals
     ,y=colMeans(means)
     ,xlab="k"
     ,ylab="mean"
     ,las=1
     ,ylim=c(0.9,1.7)
     ,main='Mean and 50% Distribution after Outlier Detection \nfor Various Outlier Detection k'
     ,pch=19
     )

# adding lines for central 50% of distribution of means of samples
for(i in 1:length(k_vals)){
  
  lines(y=quantile(means[,i],c(0.25,0.75))
        ,x=rep(k_vals[i],2)
        ,lty=1
        ,lwd=1
        ,col='gray')
  
}

# theoretical mean value for the samples
# from lognormal calculation
# will not be updated is parameters changed above!!!!!!!!!!
lines(y=rep(exp(0+1/2),length(k_vals))
      ,x=k_vals
      ,col='dodgerblue'
      ,lwd=2
      )
# legend
legend('bottomright'
       ,c('Mean of Censored Means','50% Interval of Censored Means','Theoretical Mean')
       ,lty=c(NA,1,1)
       ,pch=c(19,NA,NA)
       ,col=c('black','gray','dodgerblue')
       ,lwd=2)

