# R script accompanying Section 9.2 (Stock portfolio application) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we implement the two step optimizations.

library(quantmod)
library(copula)
library(stats)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
library(latex2exp)

start_date <- "2024-01-01"
end_date <- "2025-01-01"


## Fetch daily adjusted closing prices from Yahoo Finance and compute daily returns

#################################################################################
##########################     Technology         ###############################
#################################################################################

#Apple
apple_df <- getSymbols('AAPL', src='yahoo', auto.assign=FALSE)
apple <- exp(as.numeric(diff(log(apple_df$AAPL.Close[start_date <= index(apple_df) & index(apple_df) < end_date])))[-1])

#Microsoft 
microsoft_df <- getSymbols('MSFT', src='yahoo', auto.assign=FALSE)
microsoft <- exp(as.numeric(diff(log(microsoft_df$MSFT.Close[start_date <= index(microsoft_df) & index(microsoft_df) < end_date])))[-1])

#Google
google_df <- getSymbols('GOOGL', src='yahoo', auto.assign=FALSE)
google <- exp(as.numeric(diff(log(google_df$GOOGL.Close[start_date <= index(google_df) & index(google_df) < end_date])))[-1])

#Adobe
adobe_df <- getSymbols('ADBE', src='yahoo', auto.assign=FALSE)
adobe <- exp(as.numeric(diff(log(adobe_df$ADBE.Close[start_date <= index(adobe_df) & index(adobe_df) < end_date])))[-1])



#################################################################################
##########################     Healthcare         ###############################
#################################################################################


#Johnson & Johnson 
johnson_df <- getSymbols('JNJ', src='yahoo', auto.assign=FALSE)
johnson <- exp(as.numeric(diff(log(johnson_df$JNJ.Close[start_date <= index(johnson_df) & index(johnson_df) < end_date])))[-1])

#Pfizer
pfizer_df <- getSymbols('PFE', src='yahoo', auto.assign=FALSE)
pfizer <- exp(as.numeric(diff(log(pfizer_df$PFE.Close[start_date <= index(pfizer_df) & index(pfizer_df) < end_date])))[-1])

#Merck 
merck_df <- getSymbols('MRK', src='yahoo', auto.assign=FALSE)
merck <- exp(as.numeric(diff(log(merck_df$MRK.Close[start_date <= index(merck_df) & index(merck_df) < end_date])))[-1])

#abbvie 
abbvie_df <- getSymbols('ABBV', src='yahoo', auto.assign=FALSE)
abbvie <- exp(as.numeric(diff(log(abbvie_df$ABBV.Close[start_date <= index(abbvie_df) & index(abbvie_df) < end_date])))[-1])




#################################################################################
##################               Finance               ##########################
#################################################################################

#JPMorgan Chase
jpmorgan_df <- getSymbols('JPM', src='yahoo', auto.assign=FALSE)
jpmorgan<- exp(as.numeric(diff(log(jpmorgan_df$JPM.Close[start_date <= index(jpmorgan_df) & index(jpmorgan_df) < end_date])))[-1])

#Goldman Sachs
goldman_df <- getSymbols('GS', src='yahoo', auto.assign=FALSE)
goldman <- exp(as.numeric(diff(log(goldman_df$GS.Close[start_date <= index(goldman_df) & index(goldman_df) < end_date])))[-1])

#Bank of America
bank_df <- getSymbols('BAC', src='yahoo', auto.assign=FALSE)
bank <- exp(as.numeric(diff(log(bank_df$BAC.Close[start_date <= index(bank_df) & index(bank_df) < end_date])))[-1])

#Morgan Stanley
stanley_df <- getSymbols('MS', src='yahoo', auto.assign=FALSE)
stanley <- exp(as.numeric(diff(log(stanley_df$MS.Close[start_date <= index(stanley_df) & index(stanley_df) < end_date])))[-1])


#To work with \widetilde{R}_t as explained in paper.
df=data.frame(apple,microsoft,google,adobe,johnson,pfizer,merck,abbvie,jpmorgan,goldman,bank,stanley)
df=df-1



# Implements the full nonparametric estimator (NP) from Debrauwer and Gijbels (2024) in the case of Extremiles.
# tau =  level (0.95 here).
# Input: portfolio-weighted returns. Output: estimated extremile risk.
ExtrNPNP=function(tau,data){
  library(copula)
  set.seed(1)
  X=data
  u=pobs(X)
  ecop.orig  <- empCopula(u,smoothing="beta")   #Empirical Beta copula
  U.orig  <-  rCopula(100000, copula = ecop.orig)
  GiantSampleMatrix=toEmpMargins(U.orig,x=X)
  Zcol=rowSums(GiantSampleMatrix)
  FecdfEvaluated=rank(Zcol)/(length(Zcol)+1)
  d_tau=function(x,tau){((log(0.5)/log(1-tau)*(1-x)^(log(0.5)/log(1-tau)-1))*as.numeric(tau>0)*as.numeric(tau<=0.5)+log(0.5)/log(tau)*x^(log(0.5)/log(tau)-1)*as.numeric(tau>0.5)*as.numeric(tau<=1))*as.numeric(x<1)}
  partialLoss=function(z,t,delta){-2*(z-t)}
  lambda=function(t,tau,delta){
    FecdfEvaluated[which(FecdfEvaluated==1)]=0
    DFZ=matrix(d_tau(FecdfEvaluated,tau),nrow=1)
    return(1/length(DFZ)*DFZ%*%matrix(partialLoss(Zcol,t,delta),ncol=1))}
  findroot=function(tau,delta){
    return(uniroot(function(t){lambda(t,tau,delta)},lower=min(Zcol),upper=max(Zcol),extendInt = "yes")$root)
  }
  
  EstimatedExtr=findroot(tau,5)#This loss function does not have a delta so we fill arbitrary value
  return(EstimatedExtr)
  
}


##Helper which only requires the weights
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(-df)%*%diag(x)))}



##We implement every Optimization method separately


###############################################
###############################################
##Optimization 1
###############################################
###############################################


###### Generate weights for within cluster optimization
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


##Per cluster we determine the optimal weights (stage 1 of Optimization 1)


############################### cluster 1
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(-df[,1:4])%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
cluster1=Grid[which.min(NPNPrisk),]
cluster1#0.55 0.40 0.02 0.03
min(NPNPrisk)#0.02026971




############################### cluster 2
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(-df[,5:8])%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
cluster2=Grid[which.min(NPNPrisk),]
#cluster 2=0.59 0.15 0.21 0.05
cluster2
min(NPNPrisk)#0.0142816



############################### cluster 3
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(-df[,9:12])%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
cluster3=Grid[which.min(NPNPrisk),]
cluster3# 0.18 0.11 0.56 0.15
min(NPNPrisk)#0.02048599




### Generate grid of weights for across cluster optimization 
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



##Across cluster optimization (stage 2 of Optimization 1)

variable1=rowSums(t(cluster1*t(df[,1:4])))
variable2=rowSums(t(cluster2*t(df[,5:8])))
variable3=rowSums(t(cluster3*t(df[,9:12])))

risk=function(x){
  return(ExtrNPNP(0.95,data=cbind(-variable1,-variable2,-variable3)%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
Opt1Between=Grid[which.min(NPNPrisk),]#0.28 0.60 0.12
ExtrOpt1=min(NPNPrisk)#0.01192556
Opt1Between
ExtrOpt1



#Final weight:
#c(0.28*c(0.55,0.40,0.02,0.03), 0.60*c(0.59,0.15,0.21,0.05) ,0.12*c(0.18, 0.11, 0.56, 0.15))

#Risk: 0.01192556















###############################################
###############################################
##Optimization 2
###############################################
###############################################

#Stage 1 of Optimization 2: uniform within cluster weights
variable1=rowSums(t(rep(1/4,4)*t(df[,1:4])))
variable2=rowSums(t(rep(1/4,4)*t(df[,5:8])))
variable3=rowSums(t(rep(1/4,4)*t(df[,9:12])))


##Determine across cluster weights
risk=function(x){
  return(ExtrNPNP(0.95,data=cbind(-variable1,-variable2,-variable3)%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
Opt2Between=Grid[which.min(NPNPrisk),]#0.21 0.59 0.20
ExtrOpt2=min(NPNPrisk)#0.0130436


##Determine final risk, using again the optimized within cluster weights 
variable1=rowSums(t(cluster1*t(df[,1:4])))
variable2=rowSums(t(cluster2*t(df[,5:8])))
variable3=rowSums(t(cluster3*t(df[,9:12])))


risk=function(x){
  return(ExtrNPNP(0.95,data=cbind(-variable1,-variable2,-variable3)%*%diag(x)))}
ExtrOpt2FinalRisk=risk(c(Opt2Between))#0.01206293

FinalWeightOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
# 0.1155 0.0840 0.0042 0.0063 0.3481 0.0885 0.1239 0.0295 0.0360 0.0220 0.1120 0.0300



###############################################
###############################################
##Optimization 3
###############################################
###############################################

#Estimate the 95% Extremile of the margins
ExtrMargins=function(data){
  tau=0.95
  d_tau <- function(x, tau) {
    (log(0.5) / log(tau) * x^(log(0.5) / log(tau) - 1) * as.numeric(tau > 0.5) * as.numeric(tau <= 1)) * as.numeric(x <= 1)
  }
  
  partialLoss <- function(z, t) {
    -2 * (z - t)
  }
  cdfEvaluated=rank(data)/(length(data)+1)
  
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


##Find the representative for each cluster.


#representative cluster 1
toMinimize1=c(ExtrMargins(-df[,1]),ExtrMargins(-df[,2]),ExtrMargins(-df[,3]),ExtrMargins(-df[,4]))
for (i in 1:4){
  sum=0
  for (k in 5:12){
    sum=sum+abs(cor(df[,i],df[,k],method="kendall"))
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)
#1




#representative cluster 2
toMinimize2=c(ExtrMargins(-df[,5]),ExtrMargins(-df[,6]),ExtrMargins(-df[,7]),ExtrMargins(-df[,8]))
for (i in 5:8){
  sum=0
  for (k in c(1:4,9:12)){
    sum=sum+abs(cor(df[,i],df[,k],method="kendall"))
  }
  toMinimize2[i-4]=toMinimize2[i-4]*sum
}
clusterRepresentative2=which.min(toMinimize2)+4
#7




#representative cluster 3
toMinimize3=c(ExtrMargins(-df[,9]),ExtrMargins(-df[,10]),ExtrMargins(-df[,11]),ExtrMargins(-df[,12]))
for (i in 9:12){
  sum=0
  for (k in c(1:8)){
    sum=sum+abs(cor(df[,i],df[,k],method="kendall"))
  }
  toMinimize3[i-8]=toMinimize3[i-8]*sum
}
clusterRepresentative3=which.min(toMinimize3)+8
#11


##Optimized portfolio using the representatives
risk=function(x){
  return(ExtrNPNP(0.95,data=cbind(-df[,clusterRepresentative1],-df[,clusterRepresentative2],-df[,clusterRepresentative3])%*%diag(x)))}

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
weightsOpt3=Grid[which.min(NPNPrisk),]#0.35 0.30 0.35
min(NPNPrisk)#0.0146372
#Total portfolio weight: 0.35  0  0  0  0  0  0.30  0  0  0  0.35  0 





###############################################
###############################################
##Optimization 4
###############################################
###############################################
#From the first stage of Optimization 2 we know the optimal across cluster weights and corresponding risk
#when using uniform within cluster weights.
c(0.21*rep(1/4,4),0.59*rep(1/4,4),0.20*rep(1/4,4))
#Final weight: 0.0525 0.0525 0.0525 0.0525 0.1475 0.1475 0.1475 0.1475 0.0500 0.0500 0.0500 0.0500
#with risk 0.0130436


###############################################
###############################################
##Optimization 5
###############################################
###############################################
#From stage 1 of optimization 1 we know the optimal within cluster weights
cluster1=c(0.55,0.40,0.02,0.03) 
cluster2=c(0.59,0.15,0.21,0.05)
cluster3=c(0.18,0.11,0.56,0.15)



variable1=rowSums(t(cluster1*t(df[,1:4])))
variable2=rowSums(t(cluster2*t(df[,5:8])))
variable3=rowSums(t(cluster3*t(df[,9:12])))


##determine across cluster weights
risk=function(x){
  return(ExtrNPNP(0.95,data=cbind(-variable1,-variable2,-variable3)%*%diag(x)))}

ExtrOpt5=risk(c(1/3,1/3,1/3))#0.013001

#weight:
round(c(1/3*cluster1 ,1/3*cluster2, 1/3*cluster3),digits=3)
#0.183333333 0.133333333 0.006666667 0.010000000 0.196666667 0.050000000 0.070000000
#0.016666667 0.060000000 0.036666667 0.186666667 0.050000000









#Summary:
#Opt1:
#weight: 0.1540 0.1120 0.0056 0.0084 0.3540 0.0900 0.1260 0.0300 0.0216 0.0132 0.0672 0.0180 
#risk: 0.01192556  

#Opt2:
#weight: 0.1155 0.0840 0.0042 0.0063 0.3481 0.0885 0.1239 0.0295 0.0360 0.0220 0.1120 0.0300
#risk: 0.01206293   

#Opt3: 
#weight: 0.35  0  0  0  0  0  0.30  0  0  0  0.35  0 
#risk: 0.0146372

#Opt4: 
#weight: 0.0525 0.0525 0.0525 0.0525 0.1475 0.1475 0.1475 0.1475 0.0500 0.0500 0.0500 0.0500
#risk: 0.0130436

#Opt5:
#weight: 0.183 0.133 0.007 0.010 0.197 0.050 0.070 0.017 0.060 0.037 0.187 0.050 (rounded)
#risk: 0.013001










#Data study 

#Box tests
Box.test(apple,lag=3,type="Ljung-Box") 
Box.test(microsoft,lag=3,type="Ljung-Box") 
Box.test(google,lag=3,type="Ljung-Box") 
Box.test(adobe,lag=3,type="Ljung-Box") 

Box.test(johnson,lag=3,type="Ljung-Box")
Box.test(pfizer,lag=3,type="Ljung-Box")
Box.test(merck,lag=3,type="Ljung-Box")
Box.test(abbvie,lag=3,type="Ljung-Box")

Box.test(jpmorgan,lag=3,type="Ljung-Box")
Box.test(goldman,lag=3,type="Ljung-Box")
Box.test(bank,lag=3,type="Ljung-Box")
Box.test(stanley,lag=3,type="Ljung-Box")




#Estimated risks log-returns
ExtrMargins=function(data){
  tau=0.95
  d_tau <- function(x, tau) {
    (log(0.5) / log(tau) * x^(log(0.5) / log(tau) - 1) * as.numeric(tau > 0.5) * as.numeric(tau <= 1)) * as.numeric(x <= 1)
  }
  
  partialLoss <- function(z, t) {
    -2 * (z - t)
  }
  cdfEvaluated=rank(data)/(length(data))
  
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


round(apply(-df,2,ExtrMargins),digits=4)

#Mean and variance
round(apply(df,2,mean),digits=4)
round(apply(df,2,var),digits=5)



# dependence 


#Correlations matrix and plots
df=data.frame(apple,microsoft,google,adobe,johnson,pfizer,merck,abbvie,jpmorgan,goldman,bank,stanley)
pairs(1-df,pch=16)
cor_matrix=cor(df,method="kendall")
colnames(cor_matrix)=c("Apple","Microsoft","Google","Adobe","Johnson & Johnson","Pfizer","Merck","AbbVie","JPMorgan Chase","Goldman Sachs","Bank of America", "Morgan Stanley")
rownames(cor_matrix)=c("Apple","Microsoft","Google","Adobe","Johnson & Johnson","Pfizer","Merck","AbbVie","JPMorgan Chase","Goldman Sachs","Bank of America", "Morgan Stanley")
library(corrplot)
corrplot(cor_matrix, method="color", tl.cex=0.8, addCoef.col="black")


















########################################################################## 
# Backtest: apply portfolio weights (Opt1–Opt5 and single-step benchmark)
# to out-of-sample period Jan–Jun 2025.
# Returns are plotted as 'losses' (positive values) vs 'gains' (negative).
########################################################################## 
start_date <- "2025-01-01"
end_date <- "2025-06-13"



#################################################################################
##########################     Technology         ###############################
#################################################################################

#Apple
apple_df <- getSymbols('AAPL', src='yahoo', auto.assign=FALSE)
apple <- exp(as.numeric(diff(log(apple_df$AAPL.Close[start_date <= index(apple_df) & index(apple_df) <= end_date])))[-1])

apple_close <- Cl(apple_df)[paste0(start_date, "/", end_date)]
gross_returns <- apple_close[-1] / lag(apple_close, k = 1)[-1]

#Microsoft 
microsoft_df <- getSymbols('MSFT', src='yahoo', auto.assign=FALSE)
microsoft <- exp(as.numeric(diff(log(microsoft_df$MSFT.Close[start_date <= index(microsoft_df) & index(microsoft_df) <= end_date])))[-1])

#Google
google_df <- getSymbols('GOOGL', src='yahoo', auto.assign=FALSE)
google <- exp(as.numeric(diff(log(google_df$GOOGL.Close[start_date <= index(google_df) & index(google_df) <= end_date])))[-1])

#Adobe
adobe_df <- getSymbols('ADBE', src='yahoo', auto.assign=FALSE)
adobe <- exp(as.numeric(diff(log(adobe_df$ADBE.Close[start_date <= index(adobe_df) & index(adobe_df) <= end_date])))[-1])



#################################################################################
##########################     Healthcare         ###############################
#################################################################################


#Johnson & Johnson 
johnson_df <- getSymbols('JNJ', src='yahoo', auto.assign=FALSE)
johnson <- exp(as.numeric(diff(log(johnson_df$JNJ.Close[start_date <= index(johnson_df) & index(johnson_df) <= end_date])))[-1])

#Pfizer
pfizer_df <- getSymbols('PFE', src='yahoo', auto.assign=FALSE)
pfizer <- exp(as.numeric(diff(log(pfizer_df$PFE.Close[start_date <= index(pfizer_df) & index(pfizer_df) <= end_date])))[-1])

#Merck 
merck_df <- getSymbols('MRK', src='yahoo', auto.assign=FALSE)
merck <- exp(as.numeric(diff(log(merck_df$MRK.Close[start_date <= index(merck_df) & index(merck_df) <= end_date])))[-1])

#abbvie 
abbvie_df <- getSymbols('ABBV', src='yahoo', auto.assign=FALSE)
abbvie <- exp(as.numeric(diff(log(abbvie_df$ABBV.Close[start_date <= index(abbvie_df) & index(abbvie_df) <= end_date])))[-1])




#################################################################################
##################               Finance               ##########################
#################################################################################

#JPMorgan Chase
jpmorgan_df <- getSymbols('JPM', src='yahoo', auto.assign=FALSE)
jpmorgan<- exp(as.numeric(diff(log(jpmorgan_df$JPM.Close[start_date <= index(jpmorgan_df) & index(jpmorgan_df) <= end_date])))[-1])

#Goldman Sachs
goldman_df <- getSymbols('GS', src='yahoo', auto.assign=FALSE)
goldman <- exp(as.numeric(diff(log(goldman_df$GS.Close[start_date <= index(goldman_df) & index(goldman_df) <= end_date])))[-1])

#Bank of America
bank_df <- getSymbols('BAC', src='yahoo', auto.assign=FALSE)
bank <- exp(as.numeric(diff(log(bank_df$BAC.Close[start_date <= index(bank_df) & index(bank_df) <= end_date])))[-1])

#Morgan Stanley
stanley_df <- getSymbols('MS', src='yahoo', auto.assign=FALSE)
stanley <- exp(as.numeric(diff(log(stanley_df$MS.Close[start_date <= index(stanley_df) & index(stanley_df) <= end_date])))[-1])





df=data.frame(apple,microsoft,google,adobe,johnson,pfizer,merck,abbvie,jpmorgan,goldman,bank,stanley)

##Fill in appropriate weights
weightSingle=c(0.14,0.06,0.03,0.04,0.42,0,0.13,0.05,0.06,0,0,0.07)
weightOpt1=c(0.1540,0.1120,0.0056,0.0084,0.3540,0.0900,0.1260,0.0300,0.0216,0.0132,0.0672,0.0180)
weightOpt2=c(0.1155,0.0840,0.0042,0.0063,0.3481,0.0885,0.1239,0.0295,0.0360,0.0220,0.1120,0.0300)
weightOpt3=c(0.35,0,0,0,0,0,0.30,0,0,0,0.35,0)
weightOpt4=c(0.0525,0.0525,0.0525,0.0525,0.1475,0.1475,0.1475,0.1475,0.0500,0.0500,0.0500,0.0500)

cluster1=c(0.55,0.40,0.02,0.03) 
cluster2=c(0.59,0.15,0.21,0.05)
cluster3=c(0.18,0.11,0.56,0.15)
weightOpt5=c(1/3*cluster1 ,1/3*cluster2, 1/3*cluster3)






returnsSingle=1-as.matrix(df)%*%matrix(weightSingle,nrow=12,ncol=1)
returns1=1-as.matrix(df)%*%matrix(weightOpt1,nrow=12,ncol=1)
returns2=1-as.matrix(df)%*%matrix(weightOpt2,nrow=12,ncol=1)
returns3=1-as.matrix(df)%*%matrix(weightOpt3,nrow=12,ncol=1)
returns4=1-as.matrix(df)%*%matrix(weightOpt4,nrow=12,ncol=1)
returns5=1-as.matrix(df)%*%matrix(weightOpt5,nrow=12,ncol=1)


##Construction of Figure with performance of different portfolios in backtest. Postive values are losses.
par(mar = c(5.1, 4.3, 4.1, 2.1))
plot(index(gross_returns),returnsSingle,type="l",ylim=c(-0.15,0.15),lwd=3,xaxt="n",xlab="",ylab=TeX("$\\tilde{R}_t$"))
pretty_dates <- pretty(index(gross_returns))
axis(1, at = pretty_dates, labels = format(pretty_dates, "%b %d %Y"), cex.axis = 0.8)
ylims <- par("usr")[3:4]
xpos <- par("usr")[1] -10
# "Gains" (for negative values)
text(x = xpos, y = mean(c(ylims[1], 0)), labels = "Gains", srt = 90, adj = 0.5, xpd = TRUE)

# "Losses" (for positive values)
text(x = xpos,  y = mean(c(0, ylims[2])), labels = "Losses", srt = 90, adj = 0.5, xpd = TRUE)
lines(index(gross_returns),returns1,col="red",lwd=3)
lines(index(gross_returns),returns2,col="blue",lwd=3)
lines(index(gross_returns),returns3,col="forestgreen",lwd=4,lty=3)
lines(index(gross_returns),returns4,col="orange",lwd=3)
lines(index(gross_returns),returns5,col="purple",lwd=3,lty=6)
legend("topleft",col=c("black","red","blue","forestgreen","orange","purple"),
       legend=c("Single-step","Optimization 1","Optimization 2","Optimization 3","Optimization 4","Optimization 5"),
       lty=c(1,1,1,3,1,6),lwd=c(3,3,3,4,3,3),bty="n")
abline(h=0)
par(mar = c(5.1, 4.1, 4.1, 2.1))#default



##Construction of Figure with performance of individual stocks. Stocks from the same cluster have the same color.
par(mar = c(5.1, 4.3, 4.1, 2.1))
plot(index(gross_returns), 1-coredata(gross_returns), type = "l",
     ylim = c(-0.15, 0.15), lwd = 3, ylab = TeX("$\\tilde{R}_t$"), xlab = "", col = "black", xaxt = "n")
pretty_dates <- pretty(index(gross_returns))
axis(1, at = pretty_dates, labels = format(pretty_dates, "%b %d %Y"), cex.axis = 0.8)
ylims <- par("usr")[3:4]
xpos <- par("usr")[1] -10
# "Gains" (for negative values)
text(x = xpos, y = mean(c(ylims[1], 0)), labels = "Gains", srt = 90, adj = 0.5, xpd = TRUE)

# "Losses" (for positive values)
text(x = xpos,  y = mean(c(0, ylims[2])), labels = "Losses", srt = 90, adj = 0.5, xpd = TRUE)


lines(index(gross_returns),1-df[,2], col = "black",lwd=3)
lines(index(gross_returns),1-df[,3], col = "black",lwd=3)
lines(index(gross_returns),1-df[,4], col = "black",lwd=3)

lines(index(gross_returns),1-df[,5], col = "#FF8C00",lwd=3)
lines(index(gross_returns),1-df[,6], col = "#FF8C00",lwd=3)
lines(index(gross_returns),1-df[,7], col = "#FF8C00",lwd=3)
lines(index(gross_returns),1-df[,8], col = "#FF8C00",lwd=3)

lines(index(gross_returns),1-df[,9], col = "#006400",lwd=3,lty=5)
lines(index(gross_returns),1-df[,10], col = "#006400",lwd=3,lty=5)
lines(index(gross_returns),1-df[,11], col = "#006400",lwd=3,lty=5)
lines(index(gross_returns),1-df[,12], col = "#006400",lwd=3,lty=5)
legend("topleft",col=c("black","#FF8C00","#006400"),lwd=c(3,3,3),lty=c(1,1,2),legend=c("Technology","Healthcare","Finance"),bty="n")
abline(h=0)
par(mar = c(5.1, 4.1, 4.1, 2.1))#default













