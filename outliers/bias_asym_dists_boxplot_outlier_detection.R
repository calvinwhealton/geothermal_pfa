# script to calculate bias in "mean" values with asymmetric boxplot outlier detection
# different values of algorithm k parameter used
# lognormal used

# values of k in boxplot detection algorithm
k_vals <- seq(0.2,5,0.1)

# log-space variance values
# coefficient of variation values fixed
cv_rs <- seq(0.1,1,0.1)
ls_var <- log(1+cv_rs^2)

# matrices to store mean of censored and uncensored values
# 10,000 replicates for each value of k
means_cen <- array(NA,dim=c(length(ls_var),length(k_vals),10000))
means_uncen <- array(NA,dim=c(length(ls_var),length(k_vals),10000))

set.seed(10)

# loop over logspace variance
for(h in 1:length(ls_var)){
  
  # loop over k
  for(i in 1:length(k_vals)){
    
    # loop over replicates
    for(j in 1:10000){
      # generating sample
      # lognormal, mu=0, sigma=1, => mode = 1 and mean = exp(0+1/2)
      x1 <- exp(rnorm(50,mean=0,sd=sqrt(ls_var[h])))
      
      # lower and upper bounds
      # lb=q
      lb <- quantile(x1,0.25)-k_vals[i]*(quantile(x1,0.5)-quantile(x1,0.25))
      ub <- quantile(x1,0.75)+k_vals[i]*(quantile(x1,0.75)-quantile(x1,0.5))
      
      # calculation, multiplying TRUE/FALSE gives 0 or 1
      means_cen[h,i,j] <- mean(x1[which((x1 > lb)*(x1 < ub)==1)])
      means_uncen[h,i,j] <- mean(x1)
    }
  }
}


# # making plots of results
# plot(x=k_vals
#      ,y=colMeans(means)
#      ,xlab="k"
#      ,ylab="mean"
#      ,las=1
#      ,ylim=c(0.9,1.7)
#      ,main='Mean and 50% Distribution after Outlier Detection \nfor Various Outlier Detection k'
#      ,pch=19
#      )
# 
# # adding lines for central 50% of distribution of means of samples
# for(i in 1:length(k_vals)){
#   
#   lines(y=quantile(means[,i],c(0.25,0.75))
#         ,x=rep(k_vals[i],2)
#         ,lty=1
#         ,lwd=1
#         ,col='gray')
#   
# }
# 
# # theoretical mean value for the samples
# # from lognormal calculation
# # will not be updated is parameters changed above!!!!!!!!!!
# lines(y=rep(exp(0+1/2),length(k_vals))
#       ,x=k_vals
#       ,col='dodgerblue'
#       ,lwd=2
#       )
# # legend
# legend('bottomright'
#        ,c('Mean of Censored Means','50% Interval of Censored Means','Theoretical Mean')
#        ,lty=c(NA,1,1)
#        ,pch=c(19,NA,NA)
#        ,col=c('black','gray','dodgerblue')
#        ,lwd=2)
# 




# plot of RMSE for the different coefficient of variation values
MSE_uncen <- matrix(NA,nrow=length(ls_var),ncol=length(k_vals))
MSE_cen <- matrix(NA,nrow=length(ls_var),ncol=length(k_vals))

BIAS_uncen <- matrix(NA,nrow=length(ls_var),ncol=length(k_vals))
BIAS_cen <- matrix(NA,nrow=length(ls_var),ncol=length(k_vals))

# loop over the log space variance
for(h in 1:length(ls_var)){
  
  # loop over the k parameter values
  for(i in 1:length(k_vals)){
    
    BIAS_uncen[h,i] <- exp(0+ls_var[h]/2) - mean(means_uncen[h,i,])
    BIAS_cen[h,i] <- exp(0+ls_var[h]/2) - mean(means_cen[h,i,])
    
    
    MSE_uncen[h,i] <- mean((exp(0+ls_var[h]/2) - means_uncen[h,i,])^2)
    MSE_cen[h,i] <- mean((exp(0+ls_var[h]/2) - means_cen[h,i,])^2)
  }
}

plot(y=NA,x=NA
     ,xlim=range(k_vals)
     ,ylim=c(min(MSE_cen),max(MSE_cen))
     ,las=1
     ,ylab='Mean Squared Error'
     ,xlab='Outlier Algorithm Parameter k'
     #,log='y'
     ,main='Mean Squared Error of Computed Mean \nwith and without Censoring'
     )

for(i in 1:length(ls_var)){
  
  lines(x=k_vals
        ,y=MSE_cen[i,]
        ,col=brewer.pal(10,'BrBG')[i]
        ,lty=1
        ,lwd=2
        )
  
  lines(x=k_vals
        ,y=MSE_uncen[i,]
        ,col=brewer.pal(10,'BrBG')[i]
        ,lty=2
        ,lwd=2
  )

  
}

legend('topright'
       ,c('CV = 0.1','CV = 0.2','CV = 0.3','CV = 0.4','CV = 0.5','CV = 0.6','CV = 0.7','CV = 0.8','CV = 0.9','CV = 1.0')
       ,lwd=2
       ,lty=1
       ,col=brewer.pal(10,'BrBG')
       ,ncol=2
       ,title='Color: LN Dist\'n, mu=0')

legend('right'
       ,c('Censored','Uncensored')
       ,lwd=2
       ,lty=c(1,2)
       ,col='gray'
       ,ncol=1
       ,title='Line Style')

plot(y=NA,x=NA
     ,xlim=range(k_vals)
     ,ylim=range(MSE_cen)
     ,las=1
     ,ylab='Squared Values'
     ,xlab='Outlier Algorithm Parameter k'
     #,log='y'
     ,main='MSE and Squared Bias for Censored Samples'
)


for(i in 1:length(ls_var)){
  
  lines(x=k_vals
        ,y=MSE_cen[i,]
        ,col=brewer.pal(10,'BrBG')[i]
        ,lty=1
        ,lwd=2
  )
  
  
}

for(i in 1:length(ls_var)){
  
  lines(x=k_vals
        ,y=BIAS_cen[i,]^2
        ,col=brewer.pal(10,'BrBG')[i]
        ,lty=3
        ,lwd=2
  )
  
}

legend('topright'
       ,c('CV = 0.1','CV = 0.2','CV = 0.3','CV = 0.4','CV = 0.5','CV = 0.6','CV = 0.7','CV = 0.8','CV = 0.9','CV = 1.0')
       ,lwd=2
       ,lty=1
       ,col=brewer.pal(10,'BrBG')
       ,ncol=2
       ,title='Color: LN Dist\'n, mu=0')

legend('right'
       ,c('MSE','Squared Bias')
       ,lwd=2
       ,lty=c(1,3)
       ,col='gray'
       ,ncol=1
       ,title='Line Style')
