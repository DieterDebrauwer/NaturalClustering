# R script accompanying Section 9.2 (Stock portfolio application) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we implement the Single step optimization
#This script takes some time to run due to the grid sizes.

library(quantmod)
library(copula)
library(stats)
library(tictoc)
library(parallel)
library(foreach)
library(doParallel)


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


#################################################################################
#################################################################################

###### Generate grid of weigths
generate_combinations <- function(current, remaining, max_sum) {
  if (length(remaining) == 0) {
    last_value <- max_sum - sum(current)
    if (last_value %in% seq(0.05, 0.5, 0.05)) {
      return(list(c(current, last_value)))
    } else {
      return(list())
    }
  }
  
  result <- list()
  w <- seq(0.05, 0.5, 0.05)
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
combinations <- generate_combinations(numeric(0), rep(1, 11), 1)
combinations <- do.call(rbind, combinations)
combinations <- as.matrix(combinations)
Grid=combinations
dim(Grid)#1056




###############################################
###############################################
#single-step
###############################################
###############################################
##To work with \widetilde{R}_t as explained in paper.

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
  ecop.orig  <- empCopula(u,smoothing="beta")   
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
  
  EstimatedExtr=findroot(tau,5)
  return(EstimatedExtr)
  
}




##Helper which only requires the weights
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(-df)%*%diag(x)))}


#As explained in the paper we follow a procedure where we first consider a coarse grid and 
#use this as a starting point to make furter refinements. 

#First coarse grid
tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(Grid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(Grid[i,])
parallel::stopCluster(cl)
toc()
SingleStepCoarse=as.numeric(Grid[which.min(NPNPrisk),])

#Inital optimal grid point
as.numeric(Grid[which.min(NPNPrisk),])# 0.05 0.05 0.05 0.05 0.30 0.05 0.10 0.05 0.05 0.05 0.05 0.15


#Make a small grid around first optimal point
stepsize=0.02
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepCoarse[1:4],SingleStepCoarse[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]


tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()

#Second optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
# 0.07 0.07 0.03 0.03 0.38 0.03 0.12 0.05 0.03 0.03 0.03 0.13



#Refine again
stepsize=0.02
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepFiner[1:4],SingleStepFiner[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]
dim(finerGrid)




tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()
#Third optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
# 0.09 0.05 0.05 0.03 0.42 0.01 0.14 0.05 0.03 0.01 0.01 0.11


#Refine grid further
stepsize=0.02
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepFiner[1:4],SingleStepFiner[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]
dim(finerGrid)


tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()
#Fourth optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
# 0.11 0.07 0.03 0.03 0.40 0.01 0.14 0.05 0.05 0.01 0.01 0.09




#Refine further
stepsize=0.01
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepFiner[1:4],SingleStepFiner[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]
dim(finerGrid)


tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()
#Fifth optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
# 0.12 0.06 0.04 0.04 0.42 0.00 0.14 0.05 0.05 0.00 0.00 0.08


#Refine 
stepsize=0.01
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepFiner[1:4],SingleStepFiner[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]
dim(finerGrid)

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()
#Sixth optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
#0.13 0.07 0.03 0.04 0.42 0.00 0.13 0.05 0.06 0.00 0.00 0.07



#Final refinement
stepsize=0.01
stepnum=1
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)
FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,
                      basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(SingleStepFiner[1:4],SingleStepFiner[6:12]),
                        each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:4],1-rowSums(finerGrid),finerGrid[,5:11])
finerGridPos=apply(finerGrid,1,function(row) all(row>0))
finerGrid=finerGrid[finerGridPos,]
dim(finerGrid)

tic()
nCores=parallel::detectCores()-2
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()
#Last optimal grid point
SingleStepFiner=as.numeric(finerGrid[which.min(NPNPrisk),])
#0.14 0.06 0.03 0.04 0.42 0.00 0.13 0.05 0.06 0.00 0.00 0.07





#           Overview of grid points          and                  corresponding risk 
#0.05 0.05 0.05 0.05 0.30 0.05 0.10 0.05 0.05 0.05 0.05 0.15        0.0126286
#0.07 0.07 0.03 0.03 0.38 0.03 0.12 0.05 0.03 0.03 0.03 0.13        0.01202582
#0.09 0.05 0.05 0.03 0.42 0.01 0.14 0.05 0.03 0.01 0.01 0.11        0.01177677
#0.11 0.07 0.03 0.03 0.40 0.01 0.14 0.05 0.05 0.01 0.01 0.09        0.01169021
#0.12 0.06 0.04 0.04 0.42 0.00 0.14 0.05 0.05 0.00 0.00 0.08        0.01163352
#0.13 0.07 0.03 0.04 0.42 0.00 0.13 0.05 0.06 0.00 0.00 0.07        0.01161869
#0.14 0.06 0.03 0.04 0.42 0.00 0.13 0.05 0.06 0.00 0.00 0.07        0.01160853



