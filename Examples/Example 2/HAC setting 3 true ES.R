# R script accompanying Section 5.2 (Example 2: nested Archimdean copula) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we consider Expected shortfall for Setting 3 from Table 6.

library("HAC")
library(copula)
library(tictoc)
library("foreach")

set.seed(1)
#Determine copula parameters based on Kendall's tau values
par1=iTau(gumbelCopula(dim=2),0.8)
par2=iTau(gumbelCopula(dim=2),0.5)
par3=iTau(gumbelCopula(dim=2),0)
#Setup HAC (nested Archimdean copula)
tree = list(list("X1", "X2", par1), list("X3", "X4", par2), par3)
model = hac(type = 1, tree = tree)
model
plot(model, cex = 0.8, circles = 0.35)

#Draw large sample from copula and convert using inverse probability integral transform
CopulaSample <- rHAC(1000000,model)
sample=cbind(qexp(CopulaSample[,1],rate=2/3),qnorm(CopulaSample[,2],sd=1.5),qnorm(CopulaSample[,3],sd=1.25),qnorm(CopulaSample[,4],sd=1.5))








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





#Computes (Estimated based on very large sample) Expected shortfall.
#INPUT: (tau= Expected shortfall level, weight= portfolio weight)
#OUTPUT: Expected shortfall
ES=function(tau,weight){
  Zcol=weight[1]*sample[,1]+weight[2]*sample[,2]+weight[3]*sample[,3]+weight[4]*sample[,4]
  FecdfEvaluated=rank(Zcol)/(length(Zcol)+1)
  d_tau=function(x,tau){1/(1-tau)*as.numeric(x>=tau)*as.numeric(x<=1)*as.numeric(x>=0)}
  partialLoss=function(z,t){-2*(z-t)}
  DFZ=matrix(d_tau(FecdfEvaluated,tau),nrow=1)
  lambda=function(t,tau){
    return(1/length(DFZ)*DFZ%*%matrix(partialLoss(Zcol,t),ncol=1))}
  findroot=function(tau){
    return(uniroot(function(t){lambda(t,tau)},lower=min(Zcol),upper=max(Zcol),extendInt = "yes")$root)
  }
  return(findroot(tau))
}



#################################################################################
####################          Single step      ##################################
#################################################################################

#Iterate over grid
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
riskPerWeigth=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,Grid[i,])
parallel::stopCluster(cl)
toc()
#Takes more than 12 hours to run due to large grid size

#Optimal weight and minimal risk
Grid[which.min(riskPerWeigth),]#0.01 0.40 0.56 0.03
min(riskPerWeigth)#1.998166









#################################################################################
####################      optimization    1    ##################################
#################################################################################
###### Generate weights for two-step optimizations
w=seq(0.01, 0.99, 0.01)
Grid=cbind(w,1-w)

##########Optimize Cluster 1
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
Cluster1Risks=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(Grid[i,],0,0))
parallel::stopCluster(cl)
toc()
Optimal1Cluster1=Grid[which.min(Cluster1Risks),]
Optimal1Cluster1

##########Optimize Cluster 2
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
Cluster2Risks=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(0,0,Grid[i,]))
parallel::stopCluster(cl)
toc()
Optimal1Cluster2=Grid[which.min(Cluster2Risks),]
Optimal1Cluster2

###Optimize Across clusters
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
betweenClustersRisks=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(Grid[i,1]*Optimal1Cluster1,Grid[i,2]*Optimal1Cluster2))
parallel::stopCluster(cl)
toc()

Grid[which.min(betweenClustersRisks),]
##Total result
weight=c(Grid[which.min(betweenClustersRisks),1]*Optimal1Cluster1,Grid[which.min(betweenClustersRisks),2]*Optimal1Cluster2)
#Weight: 0.0041 0.4059 0.5841 0.0059 
#ES: 1.987322
risk=ES(0.95,c(Grid[which.min(betweenClustersRisks),1]*Optimal1Cluster1,Grid[which.min(betweenClustersRisks),2]*Optimal1Cluster2))



#################################################################################
####################      optimization    2    ##################################
#################################################################################
###Optimize Across clusters with uniform weights within clusters
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
betweenClustersRisksOpt2=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(Grid[i,1]*c(1/2,1/2),Grid[i,2]*c(1/2,1/2)))
parallel::stopCluster(cl)
toc()
Grid[which.min(betweenClustersRisksOpt2),]

##Total result
weight2=c(Grid[which.min(betweenClustersRisksOpt2),1]*Optimal1Cluster1,Grid[which.min(betweenClustersRisksOpt2),2]*Optimal1Cluster2)
#
risk2=ES(0.95,c(Grid[which.min(betweenClustersRisksOpt2),1]*Optimal1Cluster1,Grid[which.min(betweenClustersRisksOpt2),2]*Optimal1Cluster2))
#
#Weight: 0.003 0.297 0.693 0.007
#Risk: 2.034554



#################################################################################
####################            optimization    3         #######################
#################################################################################
#Variables 2 and 3 are cluster representatives


tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
riskOpt3Representative=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(0,Grid[i,1],Grid[i,2],0))
parallel::stopCluster(cl)
toc()
Grid[which.min(riskOpt3Representative),]
min(riskOpt3Representative)
#Weight: 0 0.41 0.59 0 
#Risk: 1.979807



#################################################################################
####################      optimization    4    ##################################
#################################################################################
###Across cluster weight using uniform within cluster weights
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
betweenClustersRisksOpt4=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(Grid[i,1]*c(1/2,1/2),Grid[i,2]*c(1/2,1/2)))
parallel::stopCluster(cl)
toc()
Grid[which.min(betweenClustersRisksOpt4),]

##Total result
c(Grid[which.min(betweenClustersRisksOpt4),1]*c(1/2,1/2),Grid[which.min(betweenClustersRisksOpt4),2]*c(1/2,1/2))
#Weight: 0.15 0.15 0.35 0.35
ES(0.95,c(Grid[which.min(betweenClustersRisksOpt4),1]*c(1/2,1/2),Grid[which.min(betweenClustersRisksOpt4),2]*c(1/2,1/2)))
#Risk:2.380255




#################################################################################
####################      optimization    5    ##################################
#################################################################################
########## Cluster 1 optimization 
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
Cluster1Risks=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(Grid[i,],0,0))
parallel::stopCluster(cl)
toc()
Optimal1Cluster1=Grid[which.min(Cluster1Risks),]
Optimal1Cluster1 #
##########Cluster 2 optimization 
tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
Cluster2Risks=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% ES(0.95,c(0,0,Grid[i,]))
parallel::stopCluster(cl)
toc()
Optimal1Cluster2=Grid[which.min(Cluster2Risks),]
Optimal1Cluster2 #


##Total result combine using uniform across cluster weights
c(1/2*Optimal1Cluster1,1/2*Optimal1Cluster2)
#Weight: 0.005 0.495 0.495 0.005 
#ES: 2.022879
ES(0.95,c(1/2*Optimal1Cluster1,1/2*Optimal1Cluster2))



#Summary
#sample size 1 million and set.seed(1)
#Single step:  
#OPT 1:  0.0041 0.4059 0.5841 0.0059  with risk 1.987322
#OPT 2:  0.003 0.297 0.693 0.007 with risk 2.034554
#OPT 3:  0.0041 0.4059 0.5841 0.0059 with risk 1.987322
#OPT 4:  0.15 0.15 0.35 0.35 with risk 2.380255
#OPT 5:  0.005 0.495 0.495 0.005 with risk 2.022879


