# code to generate Type I Error Rates for the Asymmetric Boxplot algorithm

# packages/libraries----
library(foreach)
library(doParallel)
library(doRNG)
library(parallel)
library(Hmisc)

# set working directory, will need to change depending on user----
setwd("C:\\Users\\Jared\\Documents\\Cornell\\Research\\Publications\\DOE Grant\\CombiningRiskFactorCode\\geothermal_pfa\\outliers")
#setwd("C:\\Users\\jds485\\Documents\\TypeIErr")

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
    }else if (PI == 'Gamma11'){
      Blow = qgamma(0.25,1,1) - k*(qgamma(0.5,1,1) - qgamma(0.25,1,1))
      Bup = qgamma(0.75,1,1) + k*(qgamma(0.75,1,1) - qgamma(0.5,1,1))
      
      #Calculate the percentage of the outliers in the data
      outPct = pgamma(Blow,1,1) + (1 - pgamma(Bup,1,1))
    }else if (PI == 'Gamma22'){
      Blow = qgamma(0.25,2,2) - k*(qgamma(0.5,2,2) - qgamma(0.25,2,2))
      Bup = qgamma(0.75,2,2) + k*(qgamma(0.75,2,2) - qgamma(0.5,2,2))
      
      #Calculate the percentage of the outliers in the data
      outPct = pgamma(Blow,2,2) + (1 - pgamma(Bup,2,2))
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

#Specify the number of distributions being tested
dists = 10

#Create a dataframe to store the Type I Errors
Results_n25 = as.data.frame(matrix(0, nrow = length(kvec), ncol = dists))
colnames(Results_n25) = c('Normal', 'Beta1', 'Beta2', 'T2', 'T4', 'T6', 'T8', 'T10', 'Gamma11', 'Gamma22')
rownames(Results_n25) = kvec
#Standard Deviation of Errors
SDs = Results_n25
#Other Sample Sizes
ResultsLarge = Results_n200 = Results_n100 = Results_n50 = Results_n25

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
  Results_n25[i,] = t(apply(ErrorMat, 1, FUN = mean))*100
}

#Calculate the Type I Error rates for several distributions - Parallel----
#Detect number of computer cores:
cores = detectCores() - 1 

# N = 25 sample size----
N = 25
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
  gDist11 = as.data.frame(matrix(replicate(rep, rgamma(N, 1, 1)), nrow = N, ncol = rep))
  gDist22 = as.data.frame(matrix(replicate(rep, rgamma(N, 2, 2)), nrow = N, ncol = rep))
  
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
  g11Out = parCapply(cl, gDist11, outID2)
  g22Out = parCapply(cl, gDist22, outID2)
  stopCluster(cl)
  
  #Summarize the Type 1 Error rates for this k value
  Results_n25[i,] = cbind(mean(normOut),mean(beta1Out), mean(beta2Out), mean(t2Out), mean(t4Out), mean(t6Out), mean(t8Out), mean(t10Out), mean(g11Out), mean(g22Out))*100
  #Summarize the Standard Deviation of this k value
  SDs[i,] = cbind(sd(normOut),sd(beta1Out), sd(beta2Out), sd(t2Out), sd(t4Out), sd(t6Out), sd(t8Out), sd(t10Out), sd(g11Out), sd(g22Out))*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2, gDist11, gDist22, k, i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out, g11Out, g22Out)

# N = 50 sample size----
N = 50
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
  gDist11 = as.data.frame(matrix(replicate(rep, rgamma(N, 1, 1)), nrow = N, ncol = rep))
  gDist22 = as.data.frame(matrix(replicate(rep, rgamma(N, 2, 2)), nrow = N, ncol = rep))
  
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
  g11Out = parCapply(cl, gDist11, outID2)
  g22Out = parCapply(cl, gDist22, outID2)
  stopCluster(cl)
  
  #Summarize the Type 1 Error rates for this k value
  Results_n50[i,] = cbind(mean(normOut),mean(beta1Out), mean(beta2Out), mean(t2Out), mean(t4Out), mean(t6Out), mean(t8Out), mean(t10Out), mean(g11Out), mean(g22Out))*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2, gDist11, gDist22, k, i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out, g11Out, g22Out)

# N = 100 sample size----
N = 100
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
  gDist11 = as.data.frame(matrix(replicate(rep, rgamma(N, 1, 1)), nrow = N, ncol = rep))
  gDist22 = as.data.frame(matrix(replicate(rep, rgamma(N, 2, 2)), nrow = N, ncol = rep))
  
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
  g11Out = parCapply(cl, gDist11, outID2)
  g22Out = parCapply(cl, gDist22, outID2)
  stopCluster(cl)
  
  #Summarize the Type 1 Error rates for this k value
  Results_n100[i,] = cbind(mean(normOut),mean(beta1Out), mean(beta2Out), mean(t2Out), mean(t4Out), mean(t6Out), mean(t8Out), mean(t10Out), mean(g11Out), mean(g22Out))*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2, gDist11, gDist22, k, i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out, g11Out, g22Out)

# N = 200 sample size----
N = 200
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
  gDist11 = as.data.frame(matrix(replicate(rep, rgamma(N, 1, 1)), nrow = N, ncol = rep))
  gDist22 = as.data.frame(matrix(replicate(rep, rgamma(N, 2, 2)), nrow = N, ncol = rep))
  
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
  g11Out = parCapply(cl, gDist11, outID2)
  g22Out = parCapply(cl, gDist22, outID2)
  stopCluster(cl)
  
  #Summarize the Type 1 Error rates for this k value
  Results_n200[i,] = cbind(mean(normOut),mean(beta1Out), mean(beta2Out), mean(t2Out), mean(t4Out), mean(t6Out), mean(t8Out), mean(t10Out), mean(g11Out), mean(g22Out))*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2, gDist11, gDist22, k, i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out, g11Out, g22Out)

#Perfect Information (Large Sample Size)----
for (i in 1:length(kvec)){
  #Enter a dummy value to run the function
  normDist = 1
  betaDist1 = 2
  betaDist2 = 3
  tDist2 = 4
  tDist4 = 5
  tDist6 = 6
  tDist8 = 7
  tDist10 = 8
  gDist11 = 9
  gDist22 = 10
  
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
  PI = 'Gamma11'
  g11Out = outID2(gDist11)
  PI = 'Gamma22'
  g22Out = outID2(gDist22)
  
  #Summarize the Type 1 Error rates for this k value
  ResultsLarge[i,] = cbind(normOut,beta1Out, beta2Out, t2Out, t4Out, t6Out, t8Out, t10Out, g11Out, g22Out)*100
}
rm(betaDist1, betaDist2, normDist, tDist10, tDist8, tDist6, tDist4, tDist2, gDist11, gDist22, k, i)
rm(beta1Out, beta2Out, normOut, t10Out, t8Out, t6Out, t4Out, t2Out, g11Out, g22Out)

#Export the results----
write.csv(Results_n25, file =  "TypeIError_25_100k.csv", row.names = TRUE)
write.csv(Results_n50, file =  "TypeIError_50_100k.csv", row.names = TRUE)
write.csv(Results_n100, file =  "TypeIError_100_100k.csv", row.names = TRUE)
write.csv(Results_n200, file =  "TypeIError_200_100k.csv", row.names = TRUE)
write.csv(ResultsLarge, file =  "TypeIError_PI.csv", row.names = TRUE)

#Save the results----
#save.image("~/TypeIErr/TypeIError_AllDists.RData")

#Plot the results----
#Colors for distributions
cols = c('black', 'red', 'purple', 'springgreen', 'white', 'green', 'white', 'darkgreen', 'blue', 'skyblue')
#Line types for sample size
ltys = seq(1,5,1)
#Plot characters for sample size
pchs = c(15, 18, 17, 16)

#T and Normal on same plot, Beta on another
sets = rbind(1,2)
layout(sets)
png('TypeIError.png', res = 600, width = 7, height = 9, units = 'in')
layout(sets)
par(mar = c(4,4.5,0.7,0.7), xaxs = 'i', yaxs = 'i')
#Normal and T
matplot(x = as.numeric(rownames(Results_n25[,c(1,4,6,8)])), y = Results_n25[,c(1,4,6,8)],
        type = 'o', pch = pchs[1], lty = ltys[3], col = cols[c(1,4,6,8)],
        xlim = c(1,5), ylim = c(0,30), 
        xlab = 'k', ylab = 'Type I Error (%)', 
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 0.7)
#rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
# matplot(x = as.numeric(rownames(Results_n25[,c(1,4,6,8)])), y = Results_n25[,c(1,4,6,8)],
#         type = 'o', pch = pchs[1], lty = ltys[3], col = cols[c(1,4,6,8)],
#         xlim = c(1,5), ylim = c(0,30), 
#         xlab = 'k', ylab = 'Type I Error (%)', 
#         cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, add = TRUE, cex = 0.7)
minor.tick(ny=5,nx=2)
matplot(x = as.numeric(rownames(Results_n50[,c(1,4,6,8)])), y = Results_n50[,c(1,4,6,8)],
        type = 'o', pch = pchs[2], lty = ltys[4], col = cols[c(1,4,6,8)],
        xlim = c(1,5), ylim = c(0,30), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE)
matplot(x = as.numeric(rownames(Results_n100[,c(1,4,6,8)])), y = Results_n100[,c(1,4,6,8)],
        type = 'o', pch = pchs[3], lty = ltys[2], col = cols[c(1,4,6,8)],
        xlim = c(1,5), ylim = c(0,30), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
# matplot(x = as.numeric(rownames(Results_n200[,c(1,4,6,8)])), y = Results_n200[,c(1,4,6,8)],
#         type = 'o', pch = 16, lty = ltys[5], col = cols[c(1,4,6,8)],
#         xlim = c(1,5), ylim = c(0,30), 
#         xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
matplot(x = as.numeric(rownames(ResultsLarge[,c(1,4,6,8)])), y = ResultsLarge[,c(1,4,6,8)],
        type = 'o', pch = pchs[4], lty = ltys[1], col = cols[c(1,4,6,8)],
        xlim = c(1,5), ylim = c(0,30), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
#legend('topright', legend = c('Normal', 't(2)', 't(4)','t(6)','t(8)','t(10)'), col = cols[c(1,4,6,8)], lty = 1)
legend('right', legend = c('t(2)','t(6)','t(10)', 'Normal'), col = cols[c(4,6,8,1)], lty = 1, cex = 1.2)
legend('topright', title = 'Number of Points', legend = c('25', '50','100','Perfect'), col = 'black', lty = c(3, 4, 2, 1), pch = pchs, horiz = TRUE, cex = 1.2)

#Beta
matplot(x = as.numeric(rownames(Results_n25[,c(2,3,9,10)])), y = Results_n25[,c(2,3,9,10)],
        type = 'o', pch = pchs[1], lty = ltys[3], col = cols[c(2,3,9,10)],
        xlim = c(1,5), ylim = c(0,20), 
        xlab = 'k', ylab = 'Type I Error (%)', 
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 0.7)
# rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "black")
# matplot(x = as.numeric(rownames(Results_n25[,c(2,3,9,10)])), y = Results_n25[,c(2,3,9,10)],
#         type = 'o', pch = pchs[1], lty = ltys[3], col = cols[c(2,3,9,10)],
#         xlim = c(1,5), ylim = c(0,20), 
#         xlab = 'k', ylab = 'Type I Error (%)', 
#         cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, add = TRUE, cex = 0.7)
minor.tick(ny=5,nx=2)
matplot(x = as.numeric(rownames(Results_n50[,c(2,3,9,10)])), y = Results_n50[,c(2,3,9,10)],
        type = 'o', pch = pchs[2], lty = ltys[4], col = cols[c(2,3,9,10)],
        xlim = c(1,5), ylim = c(0,20), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
matplot(x = as.numeric(rownames(Results_n100[,c(2,3,9,10)])), y = Results_n100[,c(2,3,9,10)],
        type = 'o', pch = pchs[3], lty = ltys[2], col = cols[c(2,3,9,10)],
        xlim = c(1,5), ylim = c(0,20), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
# matplot(x = as.numeric(rownames(Results_n200[,c(2,3,9,10)])), y = Results_n200[,c(2,3,9,10)],
#         type = 'o', pch = 16, lty = ltys[5], col = cols[c(2,3,9,10)],
#         xlim = c(1,5), ylim = c(0,20), 
#         xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
matplot(x = as.numeric(rownames(ResultsLarge[,c(2,3,9,10)])), y = ResultsLarge[,c(2,3,9,10)],
        type = 'o', pch = pchs[4], lty = ltys[1], col = cols[c(2,3,9,10)],
        xlim = c(1,5), ylim = c(0,20), 
        xlab = '', ylab = '', axes = FALSE, add = TRUE, cex = 0.7)
legend('right', legend = c('Beta(1,1)', 'Beta(2,2)', 'Gamma(1)', 'Gamma(2)'), col = cols[c(2,3,9,10)], lty = 1, cex = 1.2)
legend('topright', title = 'Number of Points', legend = c('25', '50','100','Perfect'), col = 'black', lty = c(3, 4, 2, 1), pch = pchs, horiz = TRUE, cex = 1.2)
#text(x = 1.1, y = 18, 'B', cex = 1.5)
dev.off()

#Parallel - Slow, use the implementation above. Saving for reference only. ----
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