# code to generate Type I Error Rates for the Asymmetric Boxplot algorithm

# packages/libraries----
library(foreach)
library(doParallel)
library(doRNG)
library(parallel)

# set working directory, will need to change depending on user----
setwd("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers")

#Define the outlier identification algorithm----
outID2 = function(data #dataframe to test for outliers
                  ){
  if (PI != FALSE){
    #Perfect information. Use population quantiles
    if (PI == 'Norm'){
      Blow = qnorm(0.25) - k*(qnorm(0.5) - qnorm(0.25))
      Bup = qnorm(0.75) + k*(qnorm(0.75) - qnorm(0.5))
      
      #Calculate the percentage of the outliers in the data
      outPct = pnorm(Blow) + (1 - pnorm(Bup))
    }else if (PI == 'Beta1'){
      Blow = qbeta(0.25,1,1) - k*(qbeta(0.5,1,1) - qbeta(0.25,1,1))
      Bup = qbeta(0.75,1,1) + k*(qbeta(0.75,1,1) - qbeta(0.5,1,1))
      
      #Calculate the percentage of the outliers in the data
      outPct = pbeta(Blow,1,1) + (1 - pbeta(Bup,1,1))
    }else if (PI == 'Beta2'){
      Blow = qbeta(0.25,2,2) - k*(qbeta(0.5,2,2) - qbeta(0.25,2,2))
      Bup = qbeta(0.75,2,2) + k*(qbeta(0.75,2,2) - qbeta(0.5,2,2))
      
      #Calculate the percentage of the outliers in the data
      outPct = pbeta(Blow,2,2) + (1 - pbeta(Bup,2,2))
    }else if (PI == 'T2'){
      Blow = qt(0.25,2) - k*(qt(0.5,2) - qt(0.25,2))
      Bup = qt(0.75,2) + k*(qt(0.75,2) - qt(0.5,2))
      
      #Calculate the percentage of the outliers in the data
      outPct = pt(Blow,2) + (1 - pt(Bup,2))
    }else if (PI == 'T4'){
      Blow = qt(0.25,4) - k*(qt(0.5,4) - qt(0.25,4))
      Bup = qt(0.75,4) + k*(qt(0.75,4) - qt(0.5,4))
      
      #Calculate the percentage of the outliers in the data
      outPct = pt(Blow,4) + (1 - pt(Bup,4))
    }else if (PI == 'T6'){
      Blow = qt(0.25,6) - k*(qt(0.5,6) - qt(0.25,6))
      Bup = qt(0.75,6) + k*(qt(0.75,6) - qt(0.5,6))
      
      #Calculate the percentage of the outliers in the data
      outPct = pt(Blow,6) + (1 - pt(Bup,6))
    }else if (PI == 'T8'){
      Blow = qt(0.25,8) - k*(qt(0.5,8) - qt(0.25,8))
      Bup = qt(0.75,8) + k*(qt(0.75,8) - qt(0.5,8))
      
      #Calculate the percentage of the outliers in the data
      outPct = pt(Blow,8) + (1 - pt(Bup,8))
    }else if (PI == 'T10'){
      Blow = qt(0.25,10) - k*(qt(0.5,10) - qt(0.25,10))
      Bup = qt(0.75,10) + k*(qt(0.75,10) - qt(0.5,10))
      
      #Calculate the percentage of the outliers in the data
      outPct = pt(Blow,10) + (1 - pt(Bup,10))
    }
  }else{
    #Get the quantiles of the data
    quants = quantile(data)
    
    #Get the lower and upper bounds based on the data
    Blow = quants[2] - k*(quants[3] - quants[2])
    Bup = quants[4] + k*(quants[4] - quants[3])  
    
    #Assign data columns for the outliers
    #outID = vector('numeric', length(data)) #all zeros
    #outID[which(data > Bup | data < Blow)] = 1 #Outlier
    
    #Calculate the percentage of the outliers in the data
    outPct = length(which(data > Bup | data < Blow))/length(data)
  }
  
  return(outPct)
}

#Set the parameters of the data to be tested----
#Create a vector of k values to be tested
kvec = seq(1,5,0.5)

#Specify the number of replicates
rep = 100000

#Specify the sample size
N = 25

#Specify the number of distributions being tested
dists = 8

#Create a dataframe to store the Type I Errors
Results = as.data.frame(matrix(0, nrow = length(kvec), ncol = dists))
colnames(Results) = c('Normal', 'Beta1', 'Beta2', 'T2', 'T4', 'T6', 'T8', 'T10')
rownames(Results) = kvec
ResultsLarge = Results 

#Set the random seed
set.seed(7523)

#Set the perfect information parameter to FALSE
PI = FALSE

#Calculate the Type I Error rates for several distributions - Serial----
for (i in 1:length(kvec)){
  #Generate the random values for the distributions
  normDist = as.data.frame(matrix(replicate(rep, rnorm(N, 0, 1)), nrow = N, ncol = rep))
  betaDist1 = as.data.frame(matrix(replicate(rep, rbeta(N, 1, 1)), nrow = N, ncol = rep))
  betaDist2 = as.data.frame(matrix(replicate(rep, rbeta(N, 2, 2)), nrow = N, ncol = rep))
  tDist2 = as.data.frame(matrix(replicate(rep, rt(N, 2)), nrow = N, ncol = rep))
  tDist4 = as.data.frame(matrix(replicate(rep, rt(N, 4)), nrow = N, ncol = rep))
  tDist6 = as.data.frame(matrix(replicate(rep, rt(N, 6)), nrow = N, ncol = rep))
  tDist8 = as.data.frame(matrix(replicate(rep, rt(N, 8)), nrow = N, ncol = rep))
  tDist10 = as.data.frame(matrix(replicate(rep, rt(N, 10)), nrow = N, ncol = rep))
  
  k = kvec[i] #used in outID2
  
  #Gather the outliers
  normOut = apply(normDist, 2, outID2)
  beta1Out = apply(betaDist1, 2, outID2)
  beta2Out = apply(betaDist2, 2, outID2)
  t2Out = apply(tDist2, 2, outID2)
  t4Out = apply(tDist4, 2, outID2)
  t6Out = apply(tDist6, 2, outID2)
  t8Out = apply(tDist8, 2, outID2)
  t10Out = apply(tDist10, 2, outID2)
  
  ErrorMat = rbind(normOut,beta1Out, beta2Out, t2Out, t4Out,t6Out,t8Out,t10Out)
  
  #Summarize the Type 1 Error rates for this k value
  Results[i,] = t(apply(ErrorMat, 1, FUN = mean))*100
}

#Calculate the Type I Error rates for several distributions - Parallel----
#Detect number of computer cores:
cores = detectCores() - 1 

for (i in 1:length(kvec)){
  #Generate the random values for the distributions
  normDist = as.data.frame(matrix(replicate(rep, rnorm(N, 0, 1)), nrow = N, ncol = rep))
  betaDist1 = as.data.frame(matrix(replicate(rep, rbeta(N, 1, 1)), nrow = N, ncol = rep))
  betaDist2 = as.data.frame(matrix(replicate(rep, rbeta(N, 2, 2)), nrow = N, ncol = rep))
  tDist2 = as.data.frame(matrix(replicate(rep, rt(N, 2)), nrow = N, ncol = rep))
  tDist4 = as.data.frame(matrix(replicate(rep, rt(N, 4)), nrow = N, ncol = rep))
  tDist6 = as.data.frame(matrix(replicate(rep, rt(N, 6)), nrow = N, ncol = rep))
  tDist8 = as.data.frame(matrix(replicate(rep, rt(N, 8)), nrow = N, ncol = rep))
  tDist10 = as.data.frame(matrix(replicate(rep, rt(N, 10)), nrow = N, ncol = rep))
  
  k = kvec[i] #used in outID2
  cl<-makeCluster(cores) #Needs to be in the loop for the seed to make reproducible results
  clusterExport(cl = cl, varlist = c("PI", "k", "outID2"))
  
  #Gather the outliers
  normOut = parCapply(cl, normDist, outID2)
  beta1Out = parCapply(cl, betaDist1, outID2)
  beta2Out = parCapply(cl, betaDist2, outID2)
  t2Out = parCapply(cl, tDist2, outID2)
  t4Out = parCapply(cl, tDist4, outID2)
  t6Out = parCapply(cl, tDist6, outID2)
  t8Out = parCapply(cl, tDist8, outID2)
  t10Out = parCapply(cl, tDist10, outID2)
  stopCluster(cl)
  
  #Summarize the Type 1 Error rates for this k value
  Results[i,] = cbind(mean(normOut),mean(beta1Out), mean(beta2Out), mean(t2Out), mean(t4Out), mean(t6Out), mean(t8Out), mean(t10Out))*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2,k)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out)

#Perfect Information----
for (i in 1:length(kvec)){
  #Generate the random values for the distributions
  normDist = 1
  betaDist1 = 2
  betaDist2 = 3
  tDist2 = 4
  tDist4 = 5
  tDist6 = 6
  tDist8 = 7
  tDist10 = 8
  
  k = kvec[i] #used in outID2
  
  #Gather the outliers
  PI = 'Norm'
  normOut = outID2(normDist)
  PI = 'Beta1'
  beta1Out = outID2(betaDist1)
  PI = 'Beta2'
  beta2Out = outID2(betaDist2)
  PI = 'T2'
  t2Out = outID2(tDist2)
  PI = 'T4'
  t4Out = outID2(tDist4)
  PI = 'T6'
  t6Out = outID2(tDist6)
  PI = 'T8'
  t8Out = outID2(tDist8)
  PI = 'T10'
  t10Out = outID2(tDist10)
  
  #Summarize the Type 1 Error rates for this k value
  ResultsLarge[i,] = cbind(normOut,beta1Out, beta2Out, t2Out, t4Out, t6Out, t8Out, t10Out)*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2,k,i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out)

#Export the results
write.csv(Results, file =  "TypeIError_25_100k.csv", row.names = TRUE)
write.csv(ResultsLarge, file =  "TypeIError_PI.csv", row.names = TRUE)

#Parallel - Slow ----
#Function for testing the outlier identification
outID = function(data, #dataframe to test for outliers
                 vec,  #column name (string) containing the data to be tested.
                 k     #value defining the outlier bound locations.
){
  #Get the quantiles of the data
  quants = quantile(data[,vec])
  
  #Get the lower and upper bounds based on the data
  Blow = quants[2] - k*(quants[3] - quants[2])
  Bup = quants[4] + k*(quants[4] - quants[3])
  
  #Assign data columns for the outliers
  data$outID = 0 #Not an outlier
  data[which(data[,vec] > Bup | data[,vec] < Blow), 'outID'] = 1 #Outlier
  
  return(data)
}
#Detect number of computer cores:
cores = detectCores() - 1 
cl<-makeCluster(cores)
registerDoParallel(cl)

#Calculate the Type I Error rates for several distributions
#Note that the random seed is set using the doRNG package.
Results = t(foreach(i = 1:length(kvec), .combine = data.frame) %do% {
  ErrorMat = foreach(j=1:rep, .combine = data.frame, .options.RNG = 7532) %dorng% {  
    #Generate the random values for the distributions
    normDist = as.data.frame(rnorm(N, 0, 1))
    betaDist1 = as.data.frame(rbeta(N, 1, 1))
    betaDist2 = as.data.frame(rbeta(N, 2, 2))
    tDist2 = as.data.frame(rt(N, 2))
    tDist4 = as.data.frame(rt(N, 4))
    tDist6 = as.data.frame(rt(N, 6))
    tDist8 = as.data.frame(rt(N, 8))
    tDist10 = as.data.frame(rt(N, 10))
    colnames(normDist) = colnames(betaDist1) = colnames(betaDist2) = colnames(tDist2) = colnames(tDist4)= colnames(tDist6)= colnames(tDist8)= colnames(tDist10) = 'test'
    
    #Gather the outliers
    normOut = outID(normDist, 'test', kvec[i])
    beta1Out = outID(betaDist1, 'test', kvec[i])
    beta2Out = outID(betaDist2, 'test', kvec[i])
    t2Out = outID(tDist2, 'test', kvec[i])
    t4Out = outID(tDist4, 'test', kvec[i])
    t6Out = outID(tDist6, 'test', kvec[i])
    t8Out = outID(tDist8, 'test', kvec[i])
    t10Out = outID(tDist10, 'test', kvec[i])
    
    #Calculate the percentage of outliers for each distribution and store in the ErrorMat
    data.frame(c(length(which(normOut$outID == 1))/N*100, 
                 length(which(beta1Out$outID == 1))/N*100, 
                 length(which(beta2Out$outID == 1))/N*100, 
                 length(which(t2Out$outID == 1))/N*100, 
                 length(which(t4Out$outID == 1))/N*100, 
                 length(which(t6Out$outID == 1))/N*100, 
                 length(which(t8Out$outID == 1))/N*100, 
                 length(which(t10Out$outID == 1))/N*100))
  }
  
  #Summarize the Type 1 Error rates for this k value
  apply(ErrorMat, 1, FUN = mean)
})
colnames(Results) = c('Normal', 'Beta1', 'Beta2', 'T2', 'T4', 'T6', 'T8', 'T10')
rownames(Results) = kvec
stopCluster(cl)