# R script accompanying Section 7.2 (Example 2: nested Archimedean copula) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we consider Expected shortfall for Setting 1 from Table 6.


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
    FecdfEvaluated <- rank(Zcol) / (length(Zcol)+1)
    
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

nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)

weightEstimateMatrix <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  # Code inside the for-loop
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
  
  
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
#save(weightEstimateMatrix, file="simulationStudyHACES14SampleSize50.Rdata")




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
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
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
  ESOpt1=min(riskPerWeight)
  
  c(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2]),ESOpt1)
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
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
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
  EstimatedES=findroot(tau)
  return(EstimatedES)
}


#Estimate expected shortfall parametrically for normal data
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
  EstimatedES=findroot(tau)
  return(EstimatedES)
}



Grid2d=round(cbind(seq(0.01,0.99,0.01),1-seq(0.01,0.99,0.01)),digits=2)


nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt3Representative <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
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
  ESOpt3Representative=min(riskPerWeight)
  Opt3RepresentativeWeight=c(0,0,0,0)
  Opt3RepresentativeWeight[clusterRepresentative1]=Opt3Representative[1]
  Opt3RepresentativeWeight[clusterRepresentative2]=Opt3Representative[2]
  
  c(Opt3RepresentativeWeight,ESOpt3Representative)}






###############################################
###############################################
##Optimization 4
###############################################
###############################################
#Generate grid
Grid2d=round(cbind(seq(0.01,0.99,0.01),1-seq(0.01,0.99,0.01)),digits=2)


nCores <- detectCores() - 3
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)
weightEstimateMatrixOpt4 <- foreach(j = 1:500, .combine = rbind, .packages = c("HAC", "copula", "stats")) %dopar% {
  CopulaSample <- copulaSamples[[j]]
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
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
  ESOpt4=min(riskPerWeight)
  
  c(c(c(1/2,1/2)*Opt4Between[1],c(1/2,1/2)*Opt4Between[2]),ESOpt4)}







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
  sample <- cbind(qexp(CopulaSample[, 1], rate = 1.4),
                  qnorm(CopulaSample[, 2], sd = 1.5),
                  qnorm(CopulaSample[, 3], sd = 1.25),
                  qnorm(CopulaSample[, 4], sd = 1.25))
  
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
  ESOpt5 <- apply_ES(0.95, weight_matrix, GiantSampleMatrix)
  c(c(cluster1*1/2,cluster2*1/2),ESOpt5)}






#save(weightEstimateMatrixOpt1,weightEstimateMatrixOpt2,weightEstimateMatrixOpt3Representative,weightEstimateMatrixOpt4,weightEstimateMatrixOpt5,file="simulationStudyHACES14SampleSize400Opt1To5.Rdata")


######### Comparing different methods in terms of MSE, Bias, ....
# load("simulationStudyHACES14SampleSize50.Rdata")
# load("simulationStudyHACES14SampleSize50Opt1To5.Rdata")

load("simulationStudyHACES14SampleSize400.Rdata")
load("simulationStudyHACES14SampleSize400Opt1To5.Rdata")

trueESsingle=1.903419
trueWeightSingle=c(0.35,0.11,0.27,0.27)

trueESOpt1=1.906987
trueWeightOpt1=c(0.4752,0.0048,0.2600,0.2600)

trueESOpt2=1.913347
trueWeightOpt2=c(0.4356,0.0044,0.2800,0.2800)

trueESOpt3Representative=1.951604
trueWeightOpt3Representative=c(0.5,0,0.5,0) 

trueESOpt4=1.906898
trueWeightOpt4=c(0.22,0.22,0.28,0.28)

trueESOpt5=1.908894
trueWeightOpt5=c(0.495,0.005,0.250,0.250) 






#Mean Euclidean distance
SingleStepWeights=weightEstimateMatrix[,1:4]
Opt1Weights=weightEstimateMatrixOpt1[,1:4]
Opt2Weights=weightEstimateMatrixOpt2[,1:4]
Opt3RepresentativeWeights=weightEstimateMatrixOpt3Representative[,1:4]
Opt4Weights=weightEstimateMatrixOpt4[,1:4]
Opt5Weights=weightEstimateMatrixOpt5[,1:4]

SingleStepES=weightEstimateMatrix[,5]
Opt1ES=weightEstimateMatrixOpt1[,5]
Opt2ES=weightEstimateMatrixOpt2[,5]
Opt3RepresentativeES=weightEstimateMatrixOpt3Representative[,5]
Opt4ES=weightEstimateMatrixOpt4[,5]
Opt5ES=weightEstimateMatrixOpt5[,5]

round(mean(sqrt(rowSums((t(t(SingleStepWeights)-trueWeightSingle))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt1Weights)-trueWeightOpt1))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt2Weights)-trueWeightOpt2))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt3RepresentativeWeights)-trueWeightOpt3Representative))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt4Weights)-trueWeightOpt4))^2))),digits=3)
round(mean(sqrt(rowSums((t(t(Opt5Weights)-trueWeightOpt5))^2))),digits=3)

#risk bias
round(mean((SingleStepES-rep(trueESsingle,500))),digits=3)
round(mean((Opt1ES-rep(trueESOpt1,500))),digits=3)
round(mean((Opt2ES-rep(trueESOpt2,500))),digits=3)
round(mean((Opt3RepresentativeES-rep(trueESOpt3Representative,500))),digits=3)
round(mean((Opt4ES-rep(trueESOpt4,500))),digits=3)
round(mean((Opt5ES-rep(trueESOpt5,500))),digits=3)

#risk variance
round(var(SingleStepES),digits=3)
round(var(Opt1ES),digits=3)
round(var(Opt2ES),digits=3)
round(var(Opt3RepresentativeES),digits=3)
round(var(Opt4ES),digits=3)
round(var(Opt5ES),digits=3)

#MSE risk
round(mean((SingleStepES-rep(trueESsingle,500))^2),digits=3)
round(mean((Opt1ES-rep(trueESOpt1,500))^2),digits=3)
round(mean((Opt2ES-rep(trueESOpt2,500))^2),digits=3)
round(mean((Opt3RepresentativeES-rep(trueESOpt3Representative,500))^2),digits=3)
round(mean((Opt4ES-rep(trueESOpt4,500))^2),digits=3)
round(mean((Opt5ES-rep(trueESOpt5,500))^2),digits=3)





