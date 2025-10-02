# R script accompanying Section 7.2 (Example 2: nested Archimedean copula) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we consider Expected shortfall for misspecification 1 from the supplementary material





library(NMOF)
library(MASS)
library(corpcor)
library(copula)
library(parallel)
library(foreach)
library(doParallel)
###################################################################################
########################## setting 1 with misspecification     #################### 
###################################################################################
### Use wrong clusters: (X1,X5,X7)    (X2,X4,X6)   (X3,X8)
between=0.1
group1=0.8
group2=0.7
group3=0.6
correlations=c(group1,group1,c(0.2,0.2,0.2,0.2,0.2),
               group1,rep(between,5),
               rep(between,5),
               group2,group2,rep(between,2),
               group2,rep(between,2),
               rep(between,2),
               group3)

# Number of variables (size of the matrix)
d <- 8
cor_matrix <- diag(1, d)

# Fill the lower triangular part
cor_matrix[lower.tri(cor_matrix)] <- correlations

# Mirror the lower triangle to the upper triangle
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
cor_matrixBig=cor_matrix
# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs


###############################################
###############################################
## Single-step
###############################################
###############################################
set.seed(1)
SingleStepWeights=matrix(NA,nrow=500,ncol=8)
SingleStepES=rep(NA,500)

for (j in 1:500){
  sample=mvrnorm(n=50,mu=rep(0,8),Sigma=cov_matrix)
  EstimatedCov=cov(sample)
  minvariance=minvar(EstimatedCov)
  singlestep=c(minvariance[1:8])
  varianceSingle=t(singlestep)%*%EstimatedCov%*%singlestep
  ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)
  SingleStepWeights[j,]=singlestep
  SingleStepES[j]=ESSingle
}

round(colMeans(SingleStepWeights),digits=2)
mean(SingleStepES)



###############################################
###############################################
##Optimization 1: with nonparametric stage 2
###############################################
###############################################


###### Generate weigths
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
combinations <- generate_combinations(numeric(0), rep(1, 2), 1)
combinations <- do.call(rbind, combinations)
combinations <- as.matrix(combinations)
Grid=combinations
dim(Grid)





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










set.seed(1)
Opt1Weights=matrix(NA,nrow=500,ncol=8)
Opt1ES=rep(NA,500)

nCores <- detectCores() - 2
cl_outer <- makeCluster(nCores)
registerDoParallel(cl_outer)

for (j in 1:500){
  sample=mvrnorm(n=400,mu=rep(0,8),Sigma=cov_matrix)
  
  ### cluster 1
  cov_matrix1 <- cov(sample[,c(1,5,7)])
  minvariance=minvar(cov_matrix1)
  cluster1=c(minvariance[1:3])
  
  ### cluster 2
  cov_matrix2 <- cov(sample[,c(2,4,6)])
  minvariance=minvar(cov_matrix2)
  cluster2=c(minvariance[1:3])
  
  ### cluster 3
  cov_matrix3 <- cov(sample[,c(3,8)])
  minvariance=minvar(cov_matrix3)
  cluster3=c(minvariance[1:2])
  
  #between clusters
  variable1=rowSums(t(cluster1*t(sample[,c(1,5,7)])))
  variable2=rowSums(t(cluster2*t(sample[,c(2,4,6)])))
  variable3=rowSums(t(cluster3*t(sample[,c(3,8)])))
  
  
  
  
  
  u5=pobs(variable1)
  u6=pobs(variable2)
  u7=pobs(variable3)
  betaCop  <- empCopula(cbind(u5,u6,u7),smoothing="beta")   
  CopulaSample  <-  rCopula(30000, copula = betaCop)
  GiantSampleMatrix=toEmpMargins(CopulaSample,x=cbind(variable1,variable2,variable3))
  weight_matrix <- Grid
  
  
  riskPerWeight <- foreach(j = 1:dim(Grid)[1], .combine = rbind, .packages = c("copula", "stats")) %dopar% {
    c(apply_ES(0.95, matrix(Grid[j,],nrow=1), GiantSampleMatrix))}
  
  Opt1Between <- Grid[which.min(riskPerWeight), ]
  Opt1Weights[j,]=c(cluster1[1]*Opt1Between[1],cluster2[1]*Opt1Between[2],
                    cluster3[1]*Opt1Between[3],cluster2[2]*Opt1Between[2],
                    cluster1[2]*Opt1Between[1],cluster2[3]*Opt1Between[2],
                    cluster1[3]*Opt1Between[1],cluster3[2]*Opt1Between[3])
  ESOpt1=min(riskPerWeight)
  Opt1ES[j]=ESOpt1
}

round(colMeans(Opt1Weights),digits=2)
mean(Opt1ES)

stopCluster(cl_outer)



###############################################
###############################################
##Optimization 2
###############################################
###############################################
set.seed(1)
Opt2Weights=matrix(NA,nrow=500,ncol=8)
Opt2ES=rep(NA,500)

for (j in 1:500){
  sample=mvrnorm(n=400,mu=rep(0,8),Sigma=cov_matrix)
  cov_matrixTotal=cov(sample)
  ### cluster 1
  cov_matrix1 <- cov(sample[,c(1,5,7)])
  minvariance=minvar(cov_matrix1)
  cluster1=c(minvariance[1:3])
  
  ### cluster 2
  cov_matrix2 <- cov(sample[,c(2,4,6)])
  minvariance=minvar(cov_matrix2)
  cluster2=c(minvariance[1:3])
  
  ### cluster 3
  cov_matrix3 <- cov(sample[,c(3,8)])
  minvariance=minvar(cov_matrix3)
  cluster3=c(minvariance[1:2])
  
  cluster1Eq=rep(1/3,3)
  cluster2Eq=rep(1/3,3)
  cluster3Eq=rep(1/2,2)
  
  
  #between clusters
  variable1=rowSums(t(cluster1Eq*t(sample[,c(1,5,7)])))
  variable2=rowSums(t(cluster2Eq*t(sample[,c(2,4,6)])))
  variable3=rowSums(t(cluster3Eq*t(sample[,c(3,8)])))
  
  cov_matrixBetweenClusters=cov(cbind(variable1,variable2,variable3))
  minvariance=minvar(cov_matrixBetweenClusters)
  Opt2Between=c(minvariance[1:3])
  Opt2Weights[j,]=c(cluster1[1]*Opt2Between[1],cluster2[1]*Opt2Between[2],
                    cluster3[1]*Opt2Between[3],cluster2[2]*Opt2Between[2],
                    cluster1[2]*Opt2Between[1],cluster2[3]*Opt2Between[2],
                    cluster1[3]*Opt2Between[1],cluster3[2]*Opt2Between[3])
  varianceOpt2=t(Opt2Weights[j,])%*%cov_matrixTotal%*%Opt2Weights[j,]
  Opt2ES[j]=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)
}

round(colMeans(Opt2Weights),digits=2)
mean(Opt2ES)







###############################################
###############################################
##Optimization 3
###############################################
###############################################


set.seed(1)
Opt3RepresentativesWeights=matrix(0,nrow=500,ncol=8)
Opt3RepresentativesES=rep(NA,500)
Opt3RepresentativesESCluster1=rep(NA,500)
Opt3RepresentativesESCluster2=rep(NA,500)
Opt3RepresentativesESCluster3=rep(NA,500)

for (k in 1:500){
  sample=mvrnorm(n=400,mu=rep(0,8),Sigma=cov_matrix)
  cov_matrixTotal=cov(sample)
  cor_matrixKendall=cor(sample,method="kendall")
  std_devs=apply(sample,2,sd)
  cor_matrixEstimate=cor(sample)
  
  #representative cluster 1
  toMinimize1=std_devs[c(1,5,7)]
  cluster1Numbers=c(1,5,7)
  for (i in 1:3){
    m=cluster1Numbers[i]
    sum=0
    for (j in c(2,3,4,6,8)){
      sum=sum+abs(cor_matrixKendall[m,j])
    }
    toMinimize1[i]=toMinimize1[i]*sum
  }
  clusterRepresentative1=cluster1Numbers[which.min(toMinimize1)]
  Opt3RepresentativesESCluster1[k]=clusterRepresentative1
  
  #representative cluster 2
  toMinimize2=std_devs[c(2,4,6)]
  cluster2Numbers=c(2,4,6)
  for (i in 1:3){
    m=cluster2Numbers[i]
    sum=0
    for (j in c(1,3,5,7,8)){
      sum=sum+abs(cor_matrixKendall[m,j])
    }
    toMinimize2[i]=toMinimize2[i]*sum
  }
  clusterRepresentative2=cluster2Numbers[which.min(toMinimize2)]
  Opt3RepresentativesESCluster2[k]=clusterRepresentative2
  
  #representative cluster 3
  toMinimize3=std_devs[c(3,8)]
  cluster3Numbers=c(3,8)
  for (i in 1:2){
    m=cluster3Numbers[i]
    sum=0
    for (j in c(1,2,4,5,6,7)){
      sum=sum+abs(cor_matrixKendall[m,j])
    }
    toMinimize3[i]=toMinimize3[i]*sum
  }
  clusterRepresentative3=cluster3Numbers[which.min(toMinimize3)]
  Opt3RepresentativesESCluster3[k]=clusterRepresentative3
  
  correlationsRepresentatives=c(cor_matrixEstimate[clusterRepresentative1,clusterRepresentative2],
                                cor_matrixEstimate[clusterRepresentative1,clusterRepresentative3],
                                cor_matrixEstimate[clusterRepresentative2,clusterRepresentative3])
  cor_matrixRepresentatives <- diag(1, 3)
  cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
  cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
  std_devsRepresentatives=std_devs[c(clusterRepresentative1,clusterRepresentative2,clusterRepresentative3)]
  diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
  cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives
  minvarianceRepresentatives=minvar(cov_matrixRepresentatives)
  
  
  WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
  varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
  ESRepresentativesEstimate=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)
  
  Opt3RepresentativesWeights[k,clusterRepresentative1]=WeigthsRepresentatives[1]
  Opt3RepresentativesWeights[k,clusterRepresentative2]=WeigthsRepresentatives[2]
  Opt3RepresentativesWeights[k,clusterRepresentative3]=WeigthsRepresentatives[3]
  Opt3RepresentativesES[k]=ESRepresentativesEstimate
}




###############################################
###############################################
##Optimization 4
###############################################
###############################################

set.seed(1)
Opt4Weights=matrix(NA,nrow=500,ncol=8)
Opt4ES=rep(NA,500)

for (j in 1:500){
  sample=mvrnorm(n=400,mu=rep(0,8),Sigma=cov_matrix)
  cov_matrixTotal=cov(sample)
  
  cluster1Eq=rep(1/3,3)
  cluster2Eq=rep(1/3,3)
  cluster3Eq=rep(1/2,2)
  
  
  variable1=rowSums(t(cluster1Eq*t(sample[,c(1,5,7)])))
  variable2=rowSums(t(cluster2Eq*t(sample[,c(2,4,6)])))
  variable3=rowSums(t(cluster3Eq*t(sample[,c(3,8)])))
  
  cov_matrixBetweenClusters=cov(cbind(variable1,variable2,variable3))
  minvariance=minvar(cov_matrixBetweenClusters)
  Opt2Between=c(minvariance[1:3])
  Opt4Weights[j,]=c(cluster1Eq[1]*Opt2Between[1],cluster2Eq[1]*Opt2Between[2],
                    cluster3Eq[1]*Opt2Between[3],cluster2Eq[2]*Opt2Between[2],
                    cluster1Eq[2]*Opt2Between[1],cluster2Eq[3]*Opt2Between[2],
                    cluster1Eq[3]*Opt2Between[1],cluster3Eq[2]*Opt2Between[3])
  varianceOpt4=t(Opt4Weights[j,])%*%cov_matrixTotal%*%Opt4Weights[j,]
  Opt4ES[j]=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)
}







###############################################
###############################################
##Optimization 5
###############################################
###############################################
set.seed(1)
Opt5Weights=matrix(NA,nrow=500,ncol=8)
Opt5ES=rep(NA,500)

for (j in 1:500){
  sample=mvrnorm(n=400,mu=rep(0,8),Sigma=cov_matrix)
  cov_matrixTotal=cov(sample)
  ### cluster 1
  cov_matrix1 <- cov(sample[,c(1,5,7)])
  minvariance=minvar(cov_matrix1)
  cluster1=c(minvariance[1:3])
  
  ### cluster 2
  cov_matrix2 <- cov(sample[,c(2,4,6)])
  minvariance=minvar(cov_matrix2)
  cluster2=c(minvariance[1:3])
  
  ### cluster 3
  cov_matrix3 <- cov(sample[,c(3,8)])
  minvariance=minvar(cov_matrix3)
  cluster3=c(minvariance[1:2])
  
  
  weightsOrderedByClusterOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
  Opt5Weights[j,]=c(weightsOrderedByClusterOpt5[1],weightsOrderedByClusterOpt5[4],weightsOrderedByClusterOpt5[7],
                    weightsOrderedByClusterOpt5[5],weightsOrderedByClusterOpt5[2],weightsOrderedByClusterOpt5[6],
                    weightsOrderedByClusterOpt5[3],weightsOrderedByClusterOpt5[8])
  
  varianceOpt5=t(Opt5Weights[j,])%*%cov_matrixTotal%*%Opt5Weights[j,]
  Opt5ES[j]=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)
  
}



