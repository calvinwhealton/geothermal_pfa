freqMin_SpecialSeis <- function(x # x values to use in integration
                                ,params # distribution parameters, a vector of 8 values in order of [thermal reservoir seismic utilization]
                                ,lb # lower bound parameters, vector of length = length(params) / 2
                                ,ub # upper bound parameters, vector of length = length(params) / 2
                                ){
  
  # matrix to hold products of cdfs of the other variables
  cdf_prods <- matrix(1,length(x),length(lb))
  
  # matrix to hold pdf values
  pdf_vals <- matrix(0,length(x),length(lb))
  
  if(min(is.na(params)) == 1){
    # condition for when all parameters are NA, return NAs
    pmin_inds <- NA
    
  }else if(is.na(params[7] && (min(is.na(params[1:6])) == 0))){
    # condition of utilization parameters are NA but others defined
    # means utilization was 0 with no uncertainty
    # always will be minimum or tied for minimum
    # can set the ratios so that utilization is the only problem
    
    #FIXME: This is a hack. Does not allow for other factors to be 0.(max(is.na(params[1:6])) == 0)
    pmin_inds <- c(0,0,0,1)
    
  }else{
    # loop for calcuations of pdf and product cdfs
    for(i in 1:length(lb)){
      
      # pdf values
      if (i == 3 & lb[i] == 2.5){
        #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
        if(is.na(params[(i-1)*2+1])==F){
          pdf_vals[,i] <-dtnorm(x, params[(i-1)*2+1], sqrt(params[i*2]), lb[i], ub[i])
        }
      }else{
        #All Weibulls
        if(is.na(params[(i-1)*2+1])==F){
          pdf_vals[,i] <- dweibull(x, params[i*2], params[(i-1)*2+1])
        }
      }
      
      # product of other cdfs
      inds_calc <- seq(1,length(lb),1)[-i]
      # loop for products
      for(j in inds_calc){
        if (j == 3 & lb[j] == 2.5){
          
          #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <- cdf_prods[,i]*(1 - (ptnorm(x, params[(j-1)*2+1], sqrt(params[j*2]), lb[j], ub[j])))
          }
        }else{
          #Multiplication of Weibulls
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <-cdf_prods[,i]*(1 - (pweibull(x, params[2*j], params[(j-1)*2+1])))
          }
        }
      }
    }
    
    # product of pdf and other 1-cdf values
    pdf_cdfs <- pdf_vals*cdf_prods
    
    # initializing values of indices
    pmin_inds <- rep(NA,length(lb))
    
    # assuming constant x spacing and calculating the dx value
    dx <- (max(x)-min(x))/(length(x)-1)
    
    # using trapezoidal integration
    for(i in 1:length(lb)){
      pmin_inds[i] <- dx*(sum(pdf_cdfs[c(-1,-length(x)),i])+0.5*sum(pdf_cdfs[c(1,length(x)),i]))
    }
  }
  
  return(pmin_inds)
}

MeanMin_SpecialSeis <- function(x # x values to use in integration
                                ,params # distribution parameters, a vector of 8 values in order of [thermal reservoir seismic utilization]
                                ,lb # lower bound parameters, vector of length = length(params) / 2
                                ,ub # upper bound parameters, vector of length = length(params) / 2
){
  
  # matrix to hold products of cdfs of the other variables
  cdf_prods <- matrix(1,length(x),length(lb))
  
  # matrix to hold pdf values
  xpdf_vals <- matrix(0,length(x),length(lb))
  
  if(min(is.na(params)) == 1){
    # condition for when all parameters are NA, return NAs
    MeanMin_inds <- NA
    
  }else if(is.na(params[7] && (min(is.na(params[1:6])) == 0))){
    # condition of utilization parameters are NA but others defined
    # means utilization was 0 with no uncertainty
    # always will be minimum or tied for minimum
    # can set the ratios so that utilization is the only problem
    
    #FIXME: This is a hack. Says all factors are 0.
    MeanMin_inds <- c(0,0,0,0)
    
  }else{
    # loop for calcuations of x*pdf and product cdfs
    for(i in 1:length(lb)){
      
      # pdf values
      if (i == 3 & lb[i] == 2.5){
        #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
        if(is.na(params[(i-1)*2+1])==F){
          xpdf_vals[,i] <-x*dtnorm(x, params[(i-1)*2+1], sqrt(params[i*2]), lb[i], ub[i])
        }
      }else{
        #All Weibulls
        if(is.na(params[(i-1)*2+1])==F){
          xpdf_vals[,i] <- x*dweibull(x, params[i*2], params[(i-1)*2+1])
        }
      }
      
      # product of other cdfs
      inds_calc <- seq(1,length(lb),1)[-i]
      # loop for products
      for(j in inds_calc){
        if (j == 3 & lb[j] == 2.5){
          
          #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <- cdf_prods[,i]*(1 - (ptnorm(x, params[(j-1)*2+1], sqrt(params[j*2]), lb[j], ub[j])))
          }
        }else{
          #Multiplication of Weibulls
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <-cdf_prods[,i]*(1 - (pweibull(x, params[2*j], params[(j-1)*2+1])))
          }
        }
      }
    }
    
    # product of pdf and other 1-cdf values
    xpdf_cdfs <- xpdf_vals*cdf_prods
    
    # initializing values of indices
    MeanMin_inds <- rep(NA,length(lb))
    
    # assuming constant x spacing and calculating the dx value
    dx <- (max(x)-min(x))/(length(x)-1)
    
    # using trapezoidal integration
    for(i in 1:length(lb)){
      MeanMin_inds[i] <- dx*(sum(xpdf_cdfs[c(-1,-length(x)),i])+0.5*sum(xpdf_cdfs[c(1,length(x)),i]))
    }
  }
  
  return(MeanMin_inds)
}

VarMin_SpecialSeis <- function(x # x values to use in integration
                                ,params # distribution parameters, a vector of 8 values in order of [thermal reservoir seismic utilization]
                                ,lb # lower bound parameters, vector of length = length(params) / 2
                                ,ub # upper bound parameters, vector of length = length(params) / 2
                                ,SiteMean #approximate mean of the minimum for sites 
){
  
  # matrix to hold products of cdfs of the other variables
  cdf_prods <- matrix(1,length(x),length(lb))
  
  # matrix to hold pdf values
  xminusmu_pdf_vals <- matrix(0,length(x),length(lb))
  
  if(min(is.na(params)) == 1){
    # condition for when all parameters are NA, return NAs
    VarMin_inds <- NA
    
  }else if(is.na(params[7] && (min(is.na(params[1:6])) == 0))){
    # condition of utilization parameters are NA but others defined
    # means utilization was 0 with no uncertainty
    # always will be minimum or tied for minimum
    # can set the ratios so that utilization is the only problem
    
    #FIXME: This is a hack. Says all factors are 0 variance.
    VarMin_inds <- c(0,0,0,0)
    
  }else{
    # loop for calcuations of ((x - mu)^2)*pdf and product cdfs
    for(i in 1:length(lb)){
      
      # pdf values
      if (i == 3 & lb[i] == 2.5){
        #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
        if(is.na(params[(i-1)*2+1])==F){
          xminusmu_pdf_vals[,i] <- ((x - SiteMean)^2)*dtnorm(x, params[(i-1)*2+1], sqrt(params[i*2]), lb[i], ub[i])
        }
      }else{
        #All Weibulls
        if(is.na(params[(i-1)*2+1])==F){
          xminusmu_pdf_vals[,i] <- ((x - SiteMean)^2)*dweibull(x, params[i*2], params[(i-1)*2+1])
        }
      }
      
      # product of other cdfs
      inds_calc <- seq(1,length(lb),1)[-i]
      # loop for products
      for(j in inds_calc){
        if (j == 3 & lb[j] == 2.5){
          
          #Seismic variable with a lower bound of 2.5 uses the doubly truncated normal.
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <- cdf_prods[,i]*(1 - (ptnorm(x, params[(j-1)*2+1], sqrt(params[j*2]), lb[j], ub[j])))
          }
        }else{
          #Multiplication of Weibulls
          if(is.na(params[(j-1)*2+1])==F){
            cdf_prods[,i] <-cdf_prods[,i]*(1 - (pweibull(x, params[2*j], params[(j-1)*2+1])))
          }
        }
      }
    }
    
    # product of pdf and other 1-cdf values
    xminusmu_pdf_cdfs <- xminusmu_pdf_vals*cdf_prods
    
    # initializing values of indices
    VarMin_inds <- rep(NA,length(lb))
    
    # assuming constant x spacing and calculating the dx value
    dx <- (max(x)-min(x))/(length(x)-1)
    
    # using trapezoidal integration
    for(i in 1:length(lb)){
      VarMin_inds[i] <- dx*(sum(xminusmu_pdf_cdfs[c(-1,-length(x)),i])+0.5*sum(xminusmu_pdf_cdfs[c(1,length(x)),i]))
    }
  }
  
  return(VarMin_inds)
}

################################
# running calculation
# calcd_freq_min <- matrix(0,nrow=nrow(dataParams),ncol=ncol(dataParams)/2)
# 
# for(i in 1:nrow(dataParams)){
#   calcd_freq_min[i,] <- freqMin_SpecialSeis(x=seq(0.01,5,0.01)
#                                             ,params=dataParams[i,]
#                                             ,lb=lb_mat[i,]
#                                             ,ub=ub_mat[i,])
#   
# }