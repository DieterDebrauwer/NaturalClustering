# R script accompanying Section 7.2 (Example 2: nested Archimedean copula) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we consider Expected shortfall for Setting 3 from Table 6.

library("HAC")
library(copula)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)

#Determine copula parameters based on Kendall's tau values
par1=iTau(gumbelCopula(dim=2),0.8)
par2=iTau(gumbelCopula(dim=2),0.5)
par3=iTau(gumbelCopula(dim=2),0)
#Setup HAC (nested Archimdean copula)
Truetree = list(list("X1", "X2", par1), list("X3", "X4", par2), par3)
Truemodel = hac(type = 1, tree = Truetree)





###### Generate grid of weights for single-step optimization
generate_combinations <- function(current, remaining, max_sum) {
  if (length(remaining) == 0) {
    last_value <- max_sum - sum(current)
    if (last_value > 0) {
      return(list(c(current, last_value)))
    } else {
      return(list())
    }
  }
  
  result <- list()
  w <- seq(0.01, 0.99, 0.01)
  for (value in w) {
    new_sum <- sum(current) + value
    if (new_sum < max_sum) {
      new_current <- c(current, value)
      result <- c(result, generate_combinations(new_current, remaining[-1], max_sum))
    }
  }
  
  return(result)
}

# Generate combinations
combinations <- generate_combinations(numeric(0), rep(1, 3), 1)
combinations <- do.call(rbind, combinations)
combinations <- as.matrix(combinations)
Grid=combinations
dim(Grid)


# Pre-generate 500 copula samples (each of size 400) 
# so that all optimization methods use the same input data.
# This ensures comparability across simulations.
set.seed(1)
copulaSamples <- vector("list", 500)
for (j in 1:500){
  copulaSamples[[j]] <- rHAC(400, Truemodel)
}



## Function: apply_ES
# Purpose: Estimate Expected Shortfall at level tau for a whole grid of portfolio weights.
# Input:
#   - tau: risk level (here 0.95)
#   - weight_matrix: candidate portfolio weights (each row is a weight vector)
#   - GiantSampleMatrix: large simulated sample (which we will generate by sampling from estimated copula (see further))
# Output: vector of ES estimates (one per row of weight_matrix).
apply_ES <- function(tau, weight_matrix, GiantSampleMatrix) {
  results <- apply(weight_matrix, 1, function(weight) {
    Zcol <- rowSums(t(t(GiantSampleMatrix) * weight))
    
    # Compute empirical CDF
    FecdfEvaluated <- rank(Zcol) / length(Zcol)
    
    # Define d_tau, partialLoss, lambda, and findroot
    d_tau=function(x,tau){1/(1-tau)*as.numeric(x>=tau)*as.numeric(x<=1)*as.numeric(x>=0)}
    
    
    partialLoss <- function(z, t) {
      -2 * (z - t)
    }
    
    DFZ <- matrix(d_tau(FecdfEvaluated, tau), nrow = 1)
    
    lambda <- function(t, tau) {
      return(1 / length(DFZ) * DFZ %*% matrix(partialLoss(Zcol, t), ncol = 1))
    }
    
    findroot <- function(tau) {
      uniroot(function(t) {
        lambda(t, tau)
      }, lower = min(Zcol), upper = max(Zcol), extendInt = "yes")$root
    }
    
    findroot(tau)
  })
  return(results)
}













###################################################################################
###################################################################################
#################   Single-step method                      #######################
###################################################################################
###################################################################################

nCores <- detectCores() - 2
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)

weightEstimateMatrix <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  # Code inside the for-loop
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))

  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))

  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)

  U=cbind(u1,u2,u3,u4)
  colnames(U)=NULL

  estimaterHAC=estimate.copula(X=U,type=1,hac=Truemodel,method=1)
  CopulaSample <- rHAC(30000, hac(type=1,tree=estimaterHAC$tree))


  GiantSampleMatrix <- cbind(qexp(CopulaSample[, 1], rate = rate1),
                             qnorm(CopulaSample[, 2], mean = muhat1, sd = sd1),
                             qnorm(CopulaSample[, 3], mean = muhat2, sd = sd2),
                             qnorm(CopulaSample[, 4], mean = muhat3, sd = sd3))

  weight_matrix <- Grid
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  weightEstimate <- Grid[which.min(riskPerWeight), ]
  c(weightEstimate, min(riskPerWeight))
}

stopCluster(cl_outer)
toc()


#save(weightEstimateMatrix,file="simulationStudyHACextremilesExtraSampleSize50.RData")


###############################################
###############################################
##Optimization 1
###############################################
###############################################
#Generate Grid
Grid2d=round(cbind(seq(0.01,0.99,0.01),1-seq(0.01,0.99,0.01)),digits=2)

nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt1 <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))
  
  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))
  
  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)
  
  
  
  
  ### cluster 1
  fitCopula1 <- fitCopula(gumbelCopula(dim = 2), cbind(u1, u2), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula1)))
  GiantSampleMatrix <- cbind(qexp(CopulaSample[, 1], rate = rate1),
                             qnorm(CopulaSample[, 2], mean = muhat1, sd = sd1))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster1 <- Grid2d[which.min(riskPerWeight), ]
  
  
  ### cluster 2
  fitCopula2 <- fitCopula(gumbelCopula(dim = 2), cbind(u3, u4), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula2)))
  GiantSampleMatrix <- cbind(qnorm(CopulaSample[, 1], mean = muhat2, sd = sd2),
                             qnorm(CopulaSample[, 2], mean = muhat3, sd = sd3))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster2 <- Grid2d[which.min(riskPerWeight), ]
  
  #between clusters
  variable1=rowSums(t(cluster1*t(sample[,1:2])))
  variable2=rowSums(t(cluster2*t(sample[,3:4])))
  u5=pobs(variable1)
  u6=pobs(variable2)
  betaCop  <- empCopula(cbind(u5,u6),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  Opt1Between <- Grid2d[which.min(riskPerWeight), ]
  ExtrOpt1=min(riskPerWeight)
  
  c(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2]),ExtrOpt1,cluster1,cluster2)
}


stopCluster(cl_outer)


###############################################
###############################################
##Optimization 2
###############################################
###############################################
nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt2 <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))
  
  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))
  
  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)
  
  
  
  
  ### cluster 1
  fitCopula1 <- fitCopula(gumbelCopula(dim = 2), cbind(u1, u2), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula1)))
  GiantSampleMatrix <- cbind(qexp(CopulaSample[, 1], rate = rate1),
                             qnorm(CopulaSample[, 2], mean = muhat1, sd = sd1))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster1 <- Grid2d[which.min(riskPerWeight), ]
  
  
  ### cluster 2
  fitCopula2 <- fitCopula(gumbelCopula(dim = 2), cbind(u3, u4), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula2)))
  GiantSampleMatrix <- cbind(qnorm(CopulaSample[, 1], mean = muhat2, sd = sd2),
                             qnorm(CopulaSample[, 2], mean = muhat3, sd = sd3))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster2 <- Grid2d[which.min(riskPerWeight), ]
  
  #Across clusters
  variable1=rowSums(t(c(1/2,1/2)*t(sample[,1:2])))
  variable2=rowSums(t(c(1/2,1/2)*t(sample[,3:4])))
  u5=pobs(variable1)
  u6=pobs(variable2)
  betaCop  <- empCopula(cbind(u5,u6),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  Opt1Between <- Grid2d[which.min(riskPerWeight), ]
  
  variable1=rowSums(t(cluster1*t(sample[,1:2])))
  variable2=rowSums(t(cluster2*t(sample[,3:4])))
  u5=pobs(variable1)
  u6=pobs(variable2)
  betaCop  <- empCopula(cbind(u5,u6),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2))
  weight_matrix <- matrix(data=c(Opt1Between),nrow=1,ncol=2)
  ESOpt2 <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  c(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2]),ESOpt2)
}



###############################################
###############################################
##Optimization 3 
###############################################
###############################################
#Estimate expected shortfall parametrically for exponential data
ESMarginsExp=function(data){
  tau=0.95
  d_tau=function(x,tau){1/(1-tau)*as.numeric(x>=tau)*as.numeric(x<=1)*as.numeric(x>=0)}
  partialLoss <- function(z, t) {
    -2 * (z - t)
  }
  cdfEvaluated=pexp(data,rate=1/mean(data))
  
  lambda=function(t,tau){
    cdfEvaluated[which(cdfEvaluated==1)]=0
    DFZ=matrix(d_tau(cdfEvaluated,tau),nrow=1)
    return(1/length(DFZ)*DFZ%*%matrix(partialLoss(data,t),ncol=1))}
  findroot=function(tau){
    return(uniroot(function(t){lambda(t,tau)},lower=min(data),upper=max(data),extendInt = "yes")$root)
  }
  EstimatedExtr=findroot(tau)
  return(EstimatedExtr)
}


#Estimate expected shortfall parametrically for exponential data
ESMarginsNorm=function(data){
  tau=0.95
  d_tau=function(x,tau){1/(1-tau)*as.numeric(x>=tau)*as.numeric(x<=1)*as.numeric(x>=0)}
  partialLoss <- function(z, t) {
    -2 * (z - t)
  }
  cdfEvaluated=pnorm(data,mean=mean(data),sd=sd(data))
  
  lambda=function(t,tau){
    cdfEvaluated[which(cdfEvaluated==1)]=0
    DFZ=matrix(d_tau(cdfEvaluated,tau),nrow=1)
    return(1/length(DFZ)*DFZ%*%matrix(partialLoss(data,t),ncol=1))}
  findroot=function(tau){
    return(uniroot(function(t){lambda(t,tau)},lower=min(data),upper=max(data),extendInt = "yes")$root)
  }
  EstimatedExtr=findroot(tau)
  return(EstimatedExtr)
}



nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt3Representative <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))
  
  #representative cluster 1
  toMinimize1=c(ESMarginsExp(sample[,1]),ESMarginsNorm(sample[,2]))
  for (i in 1:2){
    sum=0
    for (k in 3:4){
      sum=sum+abs(cor(sample[,i],sample[,k],method="kendall"))
    }
    toMinimize1[i]=toMinimize1[i]*sum
  }
  clusterRepresentative1=which.min(toMinimize1)
  
  
  #representative cluster 2
  toMinimize2=c(ESMarginsNorm(sample[,3]),ESMarginsNorm(sample[,4]))
  for (i in 3:4){
    sum=0
    for (k in 1:2){
      sum=sum+abs(cor(sample[,i],sample[,k],method="kendall"))
    }
    toMinimize2[i-2]=toMinimize2[i-2]*sum
  }
  clusterRepresentative2=which.min(toMinimize2)+2
  
  
  
  
  
  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))
  
  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)
  
  
  
  if (clusterRepresentative1==1){
    U1=u1
  } else {    
    U1=u2
  }
  
  if (clusterRepresentative2==3){
    U2=u3
  } else {    
    U2=u4
  }
  
  
  betaCop  <- empCopula(cbind(U1,U2),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=matrix(NA,nrow=30000,ncol=2)
  if (clusterRepresentative1==1){
    GiantSampleMatrix[,1]=qexp(CopulaSample[,1],rate=rate1)
  } else {  
    GiantSampleMatrix[,1]=qnorm(CopulaSample[,1],mean = muhat1, sd = sd1)
  }
  
  if (clusterRepresentative2==3){
    GiantSampleMatrix[,2]=qnorm(CopulaSample[,2],mean = muhat2, sd = sd2)
  } else {  
    GiantSampleMatrix[,2]=qnorm(CopulaSample[,2],mean = muhat3, sd = sd3)
  }
  
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  Opt3Representative <- Grid2d[which.min(riskPerWeight), ]
  ExtrOpt3Representative=min(riskPerWeight)
  Opt3RepresentativeWeight=c(0,0,0,0)
  Opt3RepresentativeWeight[clusterRepresentative1]=Opt3Representative[1]
  Opt3RepresentativeWeight[clusterRepresentative2]=Opt3Representative[2]
  
  c(Opt3RepresentativeWeight,ExtrOpt3Representative)}






###############################################
###############################################
##Optimization 4
###############################################
###############################################

nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt4 <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))
  
  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))
  
  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)
  
  
  
  
  ### cluster 1
  fitCopula1 <- fitCopula(gumbelCopula(dim = 2), cbind(u1, u2), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula1)))
  GiantSampleMatrix <- cbind(qexp(CopulaSample[, 1], rate = rate1),
                             qnorm(CopulaSample[, 2], mean = muhat1, sd = sd1))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster1 <- Grid2d[which.min(riskPerWeight), ]
  
  
  ### cluster 2
  fitCopula2 <- fitCopula(gumbelCopula(dim = 2), cbind(u3, u4), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula2)))
  GiantSampleMatrix <- cbind(qnorm(CopulaSample[, 1], mean = muhat2, sd = sd2),
                             qnorm(CopulaSample[, 2], mean = muhat3, sd = sd3))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster2 <- Grid2d[which.min(riskPerWeight), ]
  
  #between clusters
  variable1=rowSums(t(c(1/2,1/2)*t(sample[,1:2])))
  variable2=rowSums(t(c(1/2,1/2)*t(sample[,3:4])))
  u5=pobs(variable1)
  u6=pobs(variable2)
  betaCop  <- empCopula(cbind(u5,u6),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  Opt4Between <- Grid2d[which.min(riskPerWeight), ]
  ExtrOpt4=min(riskPerWeight)
  
  c(c(c(1/2,1/2)*Opt4Between[1],c(1/2,1/2)*Opt4Between[2]),ExtrOpt4)}







###############################################
###############################################
##Optimization 5
###############################################
###############################################


nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt5 <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 2/3),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.5))
  
  rate1 <- 1 / mean(sample[, 1])
  muhat1 <- mean(sample[, 2])
  muhat2 <- mean(sample[, 3])
  muhat3 <- mean(sample[, 4])
  n <- length(sample[, 1])
  sd1 <- sqrt((n - 1) / n * var(sample[, 2]))
  sd2 <- sqrt((n - 1) / n * var(sample[, 3]))
  sd3 <- sqrt((n - 1) / n * var(sample[, 4]))
  
  u1 <- pexp(sample[, 1], rate = rate1)
  u2 <- pnorm(sample[, 2], mean = muhat1, sd = sd1)
  u3 <- pnorm(sample[, 3], mean = muhat2, sd = sd2)
  u4 <- pnorm(sample[, 4], mean = muhat3, sd = sd3)
  
  
  
  
  ### cluster 1
  fitCopula1 <- fitCopula(gumbelCopula(dim = 2), cbind(u1, u2), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula1)))
  GiantSampleMatrix <- cbind(qexp(CopulaSample[, 1], rate = rate1),
                             qnorm(CopulaSample[, 2], mean = muhat1, sd = sd1))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster1 <- Grid2d[which.min(riskPerWeight), ]
  
  
  ### cluster 2
  fitCopula2 <- fitCopula(gumbelCopula(dim = 2), cbind(u3, u4), method = "mpl", estimate.variance = F)
  CopulaSample=rCopula(30000,gumbelCopula(dim=2,coef(fitCopula2)))
  GiantSampleMatrix <- cbind(qnorm(CopulaSample[, 1], mean = muhat2, sd = sd2),
                             qnorm(CopulaSample[, 2], mean = muhat3, sd = sd3))
  weight_matrix <- Grid2d
  riskPerWeight <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  cluster2 <- Grid2d[which.min(riskPerWeight), ]
  
  
  variable1=rowSums(t(cluster1*t(sample[,1:2])))
  variable2=rowSums(t(cluster2*t(sample[,3:4])))
  u5=pobs(variable1)
  u6=pobs(variable2)
  betaCop  <- empCopula(cbind(u5,u6),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2))
  weight_matrix <- matrix(data=c(1/2,1/2),nrow=1,ncol=2)
  ExtrOpt5 <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  c(c(cluster1*1/2,cluster2*1/2),ExtrOpt5)}






# save(weightEstimateMatrixOpt1,weightEstimateMatrixOpt2,weightEstimateMatrixOpt3Representative,
#      weightEstimateMatrixOpt4,weightEstimateMatrixOpt5,
#      file="simulationHACsetting3ESSampleSize50.RData")








######### Comparing different methods in terms of MSE, Bias, ....
# load("simulationHACsetting3ESSampleSize400.RData")
# load("simulationHACsetting3ESSingleStepSampleSize400.RData")

trueExtrsingle=1.998166
trueWeightSingle=c(0.01,0.40,0.56,0.03)

trueExtrOpt1=1.987322
trueWeightOpt1=c(0.0041,0.4059,0.5841,0.0059)

trueExtrOpt2=2.034554
trueWeightOpt2=c(0.003,0.297,0.693,0.007)

trueExtrOpt3Representative=1.979807
trueWeightOpt3Representative=c(0,0.41,0.59,0)

trueExtrOpt4=2.380255
trueWeightOpt4=c(0.15,0.15,0.35,0.35)

trueExtrOpt5=2.022879
trueWeightOpt5=c(0.005,0.495,0.495,0.005) 






#Mean Euclidean distance
SingleStepWeights=weightEstimateMatrix[,1:4]
Opt1Weights=weightEstimateMatrixOpt1[,1:4]
Opt2Weights=weightEstimateMatrixOpt2[,1:4]
Opt3RepresentativeWeights=weightEstimateMatrixOpt3Representative[,1:4]
Opt4Weights=weightEstimateMatrixOpt4[,1:4]
Opt5Weights=weightEstimateMatrixOpt5[,1:4]

SingleStepExtr=weightEstimateMatrix[,5]
Opt1Extr=weightEstimateMatrixOpt1[,5]
Opt2Extr=weightEstimateMatrixOpt2[,5]
Opt3RepresentativeExtr=weightEstimateMatrixOpt3Representative[,5]
Opt4Extr=weightEstimateMatrixOpt4[,5]
Opt5Extr=weightEstimateMatrixOpt5[,5]

round(mean(sqrt(rowSums((t(t(SingleStepWeights)-trueWeightSingle))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt1Weights)-trueWeightOpt1))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt2Weights)-trueWeightOpt2))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt3RepresentativeWeights)-trueWeightOpt3Representative))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt4Weights)-trueWeightOpt4))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt5Weights)-trueWeightOpt5))^2))),digits=3)

#risk bias
round(mean((SingleStepExtr-rep(trueExtrsingle,500))),digits=3)
round(mean((Opt1Extr-rep(trueExtrOpt1,500))),digits=3)
round(mean((Opt2Extr-rep(trueExtrOpt2,500))),digits=3)
round(mean((Opt3RepresentativeExtr-rep(trueExtrOpt3Representative,500))),digits=3)
round(mean((Opt4Extr-rep(trueExtrOpt4,500))),digits=3)
round(mean((Opt5Extr-rep(trueExtrOpt5,500))),digits=3)

#risk variance
round(var(SingleStepExtr),digits=3)
round(var(Opt1Extr),digits=3)
round(var(Opt2Extr),digits=3)
round(var(Opt3RepresentativeExtr),digits=3)
round(var(Opt4Extr),digits=3)
round(var(Opt5Extr),digits=3)

#MSE risk
round(mean((SingleStepExtr-rep(trueExtrsingle,500))^2),digits=3)
round(mean((Opt1Extr-rep(trueExtrOpt1,500))^2),digits=3)
round(mean((Opt2Extr-rep(trueExtrOpt2,500))^2),digits=3)
round(mean((Opt3RepresentativeExtr-rep(trueExtrOpt3Representative,500))^2),digits=3)
round(mean((Opt4Extr-rep(trueExtrOpt4,500))^2),digits=3)
round(mean((Opt5Extr-rep(trueExtrOpt5,500))^2),digits=3)

