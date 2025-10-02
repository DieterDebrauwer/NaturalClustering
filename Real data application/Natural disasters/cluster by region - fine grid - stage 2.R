# R script accompanying Section 9.1 (US disaster data 1980–2024) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#We consider clustering by region type.

#In this script we go over a finer grid and optimize across each cluster (so stage 2)
#Data is retrieved from https://www.ncei.noaa.gov/access/billions/

library(tictoc)
library(copula)
library(doParallel)

###############################################
###############################################
##Optimization 1
###############################################
###############################################
#Central
drought=c(138.7,0,0,44.8,0,0,15,0,332.7,35.8,0,102.2,0,5.7,0,8.4,0,0,9.8,32.1,7,0,67.9,23.3,0,26.2,0.7,26.3,39.4,0,0,66.4,312.6,24.8,0,0,0,0,3.8,0,0,0,16.1,17.8,0)
flooding=c(0,0,0,0,0,40.1,0,0,0,0,0,0,0,445.3,0,rep(0,13),75.2,0,64.6,31.3,0,29.2,0,0.5,25.1,22.6,0,104.4,0,0,30.6,6.3,0)
freeze=c(0,0,0,18.4,rep(0,23),29.4, rep(0,9),0.5,rep(0,7))
severeStorm=c(0,0,31.8,0,2.3,19.7,0,0,0,0,0,13.6,5.1,6.8,6.6,4.6,0,21.1,23.9,7.1,0,92.6,67.6,119,16.3,0,148.7,0,44.1,67.8,14.8,344.1,218.1,46.8,60.4,60.7,36.2,75.8,27.6,93.8,237,165.6,62.6,312.1,189.7)
tropicalCyclone=c(rep(0,15),3.3,4.2,rep(0,6),3.4,18.5,4.3,0,0,127.9,0,0,1.3,13.2,0,0,0,0,6.3,0,0,5.1,5.1,0,0,0)
wildfire=c(rep(0,36),36.3,rep(0,8))
winterStorm=c(0,0,13.1,0,0,13.3,0,0,0,6.4,0,0,0.7,25.8,23.5,0,9.6,0,0,14.3,1.2,rep(0,9),0.9,14.8,0,0,15.8,13.2,0,0,1.6,0,0,26.6,46.5,0,13.7)

C=0.37*drought+0.42*flooding+0.21*severeStorm



#East North Central
drought=c(206.3,0,0,0,0,0,0,0,397,43.4,0,35.4,0,0,0,2.8,0,0,0,0,18.8,0,0,98.9,0,13.6,27.4,49,53,0,0,19.5,258.5,165.9,0,0,0,0,0,0,23.4,59,24.7,67.2,0)
flooding=c(rep(0,6),53.8,rep(0,6),805.5,0,0,0,85.7,rep(0,10),450,0,0,44.8,0,0,44.1,0,0,0,0,242.6,0,0,0,22.6,0)
freeze=c(0,0,0,10.7,rep(0,23),0.3,rep(0,17))
severeStorm=c(0,0,7.3,rep(0,8),8.9,0,32.1,6.3,0,0,0,206,0,0,45.8,0,25.6,17.7,0,93.4,0,168,0.6,32.7,159,16.3,69,42.6,13,2.2,226.1,59.4,7.1,539.4,125.3,351.3,255.3,71.2)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,23),2.7,0,0,0,1.6,rep(0,17))
winterStorm=c(0,0,14.2,0,0,14.3,rep(0,13),2.5,rep(0,14),4.3,5.1,rep(0,6),15.9,0,4.2)

ENC=0.55*drought+0.18*flooding+0.27*severeStorm



#North east
drought=c(0,0,0,0,0,0,0,0,26.4,0,0,3.9,0,2.8,0,2.4,0,0,0,11.7,0,0,29.3,0,0,0,0,9.1,6.1,0,0,21.2,18.4,0,0,0,2,0,0,0,0,0,0,0,0)
flooding=c(rep(0,5),5.2,rep(0,10),23.1,rep(0,9),36.5,0,0,0,43.1,0,0,0,4.7,0.7,0,rep(0,6),36.1,0)
freeze=c(0,0,0,6.8,rep(0,41))
severeStorm=c(0,0,2,0,18.7,16.4,0,0,0,0,0,0.8,0,0,0,0,0,9.2,18,0,0,0.8,12.7,7.1,1.1,0,0,41.6,7.8,4.9,12.7,8.9,11,1.4,30.1,10.1,8.7,9.2,27,31.1,21.9,10,8.6,32.2,28.7)
tropicalCyclone=c(rep(0,5),42.6,rep(0,5),60.6,rep(0,4),21,0,0,46.6,0,14.5,0,55.3,65.5,0,0,0,4.2,0,0,276.6,1392.4,0,0,0,0,0,1.5,0,75.6,405.6,0,0,0)
wildfire=c(rep(0,45))
winterStorm=c(0,0,9.6,0,0,9.8,0,0,0,2.9,0,0,95.4,62.1,34.4,0,78,0,46.1,31.6,rep(0,10),15,36.9,0,0,12.3,42.3,0,0,46.7,0,0,1.6,37.6,28.4,3.5)

NE=0.02*tropicalCyclone+0.98*winterStorm


#North west (NW)
drought=c(0,0,0,0,0,0,0,0,179.4,0,0,48.8,0,0,0,0,0,0,0,0,31.8,0,68.6,83.2,0,37.5,0,21.7,71,94.8,0,0,14.6,60.4,67.8,78.9,0,0,11.3,0,6.5,81.8,67.9,62.4,0)
flooding=c(0,0,0,14.2,rep(0,12),207.3,116.3,rep(0,27))
freeze=c(0,0,0,32.8,rep(0,41))
severeStorm=c(rep(0,6),9.3,rep(0,21),8.3,rep(0,16))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),42.5,0,0,0,44.8,rep(0,5),64.5,0,28.7,26.9,0,0,39.4,44.8,13.4,0,0,12.2,47.6,0,0,39.3,22.8,101.2,55.7,0,236.4,133.5,34.9,0,0)
winterStorm=c(rep(0,41),45.4,14.1,0,123)

NW=0.37*drought+0.40*flooding+0.23*wildfire


#South (S)
drought=c(276.6,0,0,101.4,0,0,70.7,0,315.3,72.7,0,35.1,0,0,0,25.1,88.7,0,99.2,18.1,110.4,0,86.2,25.2,0,4.1,105.9,6.6,58.9,81.1,0,193.1,210.1,104.1,61.6,0,2.2,0,68.1,0,55.3,0,254.9,172.9,0)
flooding=c(95.8,0,0,160.5,0,rep(0,5),77.1,0,0,50.7,62.1,0,0,0,50.3,0,rep(0,8),6.7,0,1,65.6,0,0,0,71.9,438,21.5,0,152.4,0,33.3,0,0,0)
freeze=c(0,0,0,37.6,21.6,rep(0,23),rep(0,9),0.3,rep(0,7))
severeStorm=c(0,38.4,45.9,0,2,0,0,0,0,27.4,0,6.6,144,5.2,43.8,357.5,0,7.4,11.6,83.9,34,8.4,6.7,141.6,1.5,26.2,23.3,21.2,80.9,110.6,100.1,235.2,166.4,174.7,74.1,113.4,267.6,146.2,67.6,146.9,195,176.8,91,409.1,307.3)
tropicalCyclone=c(77.7,0,0,317.6,0,180,0,0,0,47.1,0,0,68,0,0,3,0,0,99.8,0,0,375.9,87.3,0,70.1,5538.3,0,0,1132.5,0,0,1.7,80.7,0,0,0,0,3434.6,0,131.2,738.7,1276.7,0,0,0)
wildfire=c(rep(0,14),7.8,rep(0,8),1.1,0,0,13.7,0,16.5,6.4,0,38.1,1.8,0,0,0.8,rep(0,9))
winterStorm=c(0,0,13.8,0,0,15.2,0,0,0,31.9,0,0,0,23.6,117.8,0,10.2,0,0,12.6,20.8,rep(0,10),19.1,0,0,3,rep(0,6),522.9,28.8,0,17.4)

S=0.36*drought+0.30*flooding+0.09*severeStorm+0.01*tropicalCyclone+0.24*winterStorm


#South East (SE)
drought=c(250.1,0,0,130,0,0,63.7,0,39.7,0,0,0,0,55.9,0,11.9,0,0,61.2,41,39.6,0,38,0,0,0,17.8,29.4,27.1,4.7,0,51.2,9.7,0,0,0,2.8,0,0,0,0,0,0,1.7,0)
flooding=c(0,0,0,4.8,0,55.1,0,rep(0,13),33.8,0,0,0,0,0,2.5,0,0,23.7,0.9,0,0,0,0,45.2,0.5,0.4,0,0,0,0,0,23,0)
freeze=c(0,60.6,0,97.1,0,97.7,0,0,0,132.4,rep(0,17),15.7,rep(0,9),20.9,rep(0,7))
severeStorm=c(0,8.6,1.2,0,18.3,0,0,0,0,13.8,0,11.3,51.7,0,0,0,0,3.5,11.8,8,0,0,13.8,22.3,1.1,7.1,10.4,7.6,31.1,34.1,3.7,212.7,19.7,24.3,28.5,13.2,19.6,29.7,44.4,10.5,63.6,26.8,30.6,85,50)
tropicalCyclone=c(rep(0,5),76.5,0,0,0,473.9,0,0.6,1401.9,0,49.3,252.7,191.8,0,94.1,199.6,0,0.2,12.3,116,1579.8,872.1,0,0,5,0,0,86.3,16.3,0,0,0,219.2,855.1,1000.4,31.7,239.8,34.2,1901.4,55.9,0)
wildfire=c(rep(0,22),14.7,0,0,0,0.9,5.9,1.8,0,0,1.2,0,0,0,0,2.8,rep(0,8))
winterStorm=c(0,0,10.8,0,0,16.9,0,0,0,6.5,0,0,0.7,156.9,38,0,16.9,0,0,12.8,10.1,rep(0,9),3.8,0,0,0,19.3,8.9,0,0,17.5,0,0,1.6,31,0,1.6)

SE=0.16*drought+0.29*freeze+0.14*severeStorm+0.01*tropicalCyclone+0.40*winterStorm

#South West (SW)
drought=c(0,0,0,0,0,0,0,0,0,39.5,0,0,0,0,0,4.3,44.4,0,0,0,52.8,0,68.2,74.9,0,0,56.2,6.3,36,12.6,0,54.2,45.2,62.7,13.8,1.2,0,0,7.4,0,41.8,57,97,17.9,0)
flooding=c(0,0,0,384,0,rep(0,28),119.7,0,2.6,0,0,0,0,0,4.9,0,0,0)
freeze=c(0,0,0,13.1,rep(0,41))
severeStorm=c(0,0,13.5,0,139.8,0,11.2,0,0,0,193.8,0,0,0,4.8,0,rep(0,5),1.5,0,1.7,79.4,0,0,0,28.5,90.8,380,36,115.4,0,129,45.1,95.2,223,323.8,92.7,0,23.7,12.1,263.1,90.3)
tropicalCyclone=c(rep(0,28),19.2,rep(0,16))
wildfire=c(rep(0,10),11.7,0,0,0,41.4,rep(0,5),47.3,0,50.3,9.5,0,0,17.7,10.8,14.5,15.3,0,33.4,62.3,0,0,0.8,4.6,2.2,69.1,0,95.5,204.9,63.2,0,0)
winterStorm=c(rep(0,41),7.4,7.7,0,7.2)

SW=0.69*drought+0.12*flooding+0.10*severeStorm+0.09*wildfire


#West (W)
drought=c(0,0,0,0,0,0,0,0,0,3.4,0,0,0,0,0,0,0,0,0,0,7.8,0,5.7,0.5,0,0,0,4.8,2.7,5.6,0,0,5.9,7.5,35.5,112.7,99.5,0,4.2,0,9.2,35.2,67.2,0,0)
flooding=c(0,0,0,91.9,0,rep(0,10),156.4,0,136,0,rep(0,18),45.7,0,0,0,28.9,0,112.1,0)
freeze=c(rep(0,10),271.2,rep(0,7),139,rep(0,8),55.9,rep(0,17))
severeStorm=c(rep(0,6),45.1,rep(0,11),32.9,0,rep(0,8),17.8,rep(0,8),10.5,rep(0,7))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),35.1,240.3,0,91.3,11.3,rep(0,7),4.9,147.5,0,0,14.7,73.1,13.9,10.9,0,2.3,6.7,0,0,69.5,16.5,486.9,647.2,118.4,343.5,136.5,35.2,0,0)
winterStorm=c(rep(0,41),1.1,0,0,0)

W=0.49*drought+0.32*flooding+0.18*freeze+0.01*wildfire


#West North Central (WNC)
drought=c(3340.8,0,0,0,0,0,0,0,4092.4,605.5,0,0,0,0,0,0,0,0,0,0,267.1,0,855.5,498.2,0,58,595.2,19.4,307.5,0,0,0,1604.9,448.9,0,47.2,0,623.2,0,0,193.5,932.8,807.6,478.3,0)
flooding=c(0,0,0,19.2,rep(0,9),1646.3,0,0,1.1,1227.5,0,rep(0,9),29.7,0,0,345.4,rep(0,7),1237.5,0,0,0,0,0)
freeze=c(0,0,0,31,rep(0,23),3.3,rep(0,17))
severeStorm=c(0,0,10.5,0,19,0,13.8,0,0,0,0,0,0,50.9,29.2,0,0,0,13.7,8.9,0,194.4,0,88.3,66.1,0,5.2,0,57.5,34.8,0,183.1,42.6,268.6,536.9,45.9,328.5,360.6,120.4,187.7,118.3,111.4,741.8,156.5,212.5)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),31.6,rep(0,9),142,0,0,65.4,0,0,79.7,31.9,13.9,0,0,7.6,80.7,0,0,18.6,10.1,204.1,96.1,0,45.7,119.8,10.9,0,0)
winterStorm=c(0,0,9.3,rep(0,39),7.8,0,0)

WNC=0.01*drought+0.04*flooding+0.01*severeStorm+0.94*wildfire



##Jitter 
set.seed(1)
df=cbind(C,ENC,NE,NW,S,SE,SW,W,WNC)
zero_counts <- apply(df, 2, function(x) sum(x == 0))

# 1. Identify the positions of zeros
zero_positions <- which(df == 0, arr.ind = TRUE)
# 2. Generate random numbers between 0 and 0.1
random_numbers <- runif(nrow(zero_positions), min = 0, max = 0.1)
# 3. Replace zeros with the generated random numbers
df[zero_positions] <- random_numbers



stepsize=0.02
stepnum=2
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)

FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(0.05,0.05,0.05,0.05,0.05,0.25,0.05,0.05), each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:7],1-rowSums(finerGrid),finerGrid[,8])



###Remove all negative weights
# Define the function to check rows
is_valid_row <- function(row, tolerance = 1e-6) {
  all(row > tolerance)
}

# Apply the function to each row
valid_rows <- apply(finerGrid, 1, is_valid_row)

# Subset the dataframe to keep only valid rows
df_filtered <- finerGrid[valid_rows, ]

finerGrid=df_filtered
finerGrid=round(finerGrid,digits=2)
dim(finerGrid)#390625







##ExtrNPNP: estimate extremile at level tau using nonparametric copula and nonparametric margins
ExtrNPNP=function(tau,data){
  library(copula)
  set.seed(1)
  X=data
  u=pobs(X)
  ecop.orig  <- empCopula(u,smoothing="beta")   
  U.orig  <-  rCopula(100000, copula = ecop.orig)
  
  GiantSampleMatrix=toEmpMargins(U.orig,x=X)
  Zcol=rowSums(GiantSampleMatrix)
  Fecdf=ecdf(Zcol)
  FecdfEvaluated=Fecdf(Zcol)
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



#helper function 
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(df)%*%diag(x)))}



tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()


finerGrid[which.min(NPNPrisk),]
min(NPNPrisk)



###############################################
###############################################
##Optimization 2
###############################################
###############################################
#Central
drought=c(138.7,0,0,44.8,0,0,15,0,332.7,35.8,0,102.2,0,5.7,0,8.4,0,0,9.8,32.1,7,0,67.9,23.3,0,26.2,0.7,26.3,39.4,0,0,66.4,312.6,24.8,0,0,0,0,3.8,0,0,0,16.1,17.8,0)
flooding=c(0,0,0,0,0,40.1,0,0,0,0,0,0,0,445.3,0,rep(0,13),75.2,0,64.6,31.3,0,29.2,0,0.5,25.1,22.6,0,104.4,0,0,30.6,6.3,0)
freeze=c(0,0,0,18.4,rep(0,23),29.4, rep(0,9),0.5,rep(0,7))
severeStorm=c(0,0,31.8,0,2.3,19.7,0,0,0,0,0,13.6,5.1,6.8,6.6,4.6,0,21.1,23.9,7.1,0,92.6,67.6,119,16.3,0,148.7,0,44.1,67.8,14.8,344.1,218.1,46.8,60.4,60.7,36.2,75.8,27.6,93.8,237,165.6,62.6,312.1,189.7)
tropicalCyclone=c(rep(0,15),3.3,4.2,rep(0,6),3.4,18.5,4.3,0,0,127.9,0,0,1.3,13.2,0,0,0,0,6.3,0,0,5.1,5.1,0,0,0)
wildfire=c(rep(0,36),36.3,rep(0,8))
winterStorm=c(0,0,13.1,0,0,13.3,0,0,0,6.4,0,0,0.7,25.8,23.5,0,9.6,0,0,14.3,1.2,rep(0,9),0.9,14.8,0,0,15.8,13.2,0,0,1.6,0,0,26.6,46.5,0,13.7)

C=1/3*(drought+flooding+severeStorm)



#East North Central
drought=c(206.3,0,0,0,0,0,0,0,397,43.4,0,35.4,0,0,0,2.8,0,0,0,0,18.8,0,0,98.9,0,13.6,27.4,49,53,0,0,19.5,258.5,165.9,0,0,0,0,0,0,23.4,59,24.7,67.2,0)
flooding=c(rep(0,6),53.8,rep(0,6),805.5,0,0,0,85.7,rep(0,10),450,0,0,44.8,0,0,44.1,0,0,0,0,242.6,0,0,0,22.6,0)
freeze=c(0,0,0,10.7,rep(0,23),0.3,rep(0,17))
severeStorm=c(0,0,7.3,rep(0,8),8.9,0,32.1,6.3,0,0,0,206,0,0,45.8,0,25.6,17.7,0,93.4,0,168,0.6,32.7,159,16.3,69,42.6,13,2.2,226.1,59.4,7.1,539.4,125.3,351.3,255.3,71.2)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,23),2.7,0,0,0,1.6,rep(0,17))
winterStorm=c(0,0,14.2,0,0,14.3,rep(0,13),2.5,rep(0,14),4.3,5.1,rep(0,6),15.9,0,4.2)

ENC=1/3*(drought+flooding+severeStorm)



#North east
drought=c(0,0,0,0,0,0,0,0,26.4,0,0,3.9,0,2.8,0,2.4,0,0,0,11.7,0,0,29.3,0,0,0,0,9.1,6.1,0,0,21.2,18.4,0,0,0,2,0,0,0,0,0,0,0,0)
flooding=c(rep(0,5),5.2,rep(0,10),23.1,rep(0,9),36.5,0,0,0,43.1,0,0,0,4.7,0.7,0,rep(0,6),36.1,0)
freeze=c(0,0,0,6.8,rep(0,41))
severeStorm=c(0,0,2,0,18.7,16.4,0,0,0,0,0,0.8,0,0,0,0,0,9.2,18,0,0,0.8,12.7,7.1,1.1,0,0,41.6,7.8,4.9,12.7,8.9,11,1.4,30.1,10.1,8.7,9.2,27,31.1,21.9,10,8.6,32.2,28.7)
tropicalCyclone=c(rep(0,5),42.6,rep(0,5),60.6,rep(0,4),21,0,0,46.6,0,14.5,0,55.3,65.5,0,0,0,4.2,0,0,276.6,1392.4,0,0,0,0,0,1.5,0,75.6,405.6,0,0,0)
wildfire=c(rep(0,45))
winterStorm=c(0,0,9.6,0,0,9.8,0,0,0,2.9,0,0,95.4,62.1,34.4,0,78,0,46.1,31.6,rep(0,10),15,36.9,0,0,12.3,42.3,0,0,46.7,0,0,1.6,37.6,28.4,3.5)

NE=1/2*(tropicalCyclone+winterStorm)


#North west (NW)
drought=c(0,0,0,0,0,0,0,0,179.4,0,0,48.8,0,0,0,0,0,0,0,0,31.8,0,68.6,83.2,0,37.5,0,21.7,71,94.8,0,0,14.6,60.4,67.8,78.9,0,0,11.3,0,6.5,81.8,67.9,62.4,0)
flooding=c(0,0,0,14.2,rep(0,12),207.3,116.3,rep(0,27))
freeze=c(0,0,0,32.8,rep(0,41))
severeStorm=c(rep(0,6),9.3,rep(0,21),8.3,rep(0,16))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),42.5,0,0,0,44.8,rep(0,5),64.5,0,28.7,26.9,0,0,39.4,44.8,13.4,0,0,12.2,47.6,0,0,39.3,22.8,101.2,55.7,0,236.4,133.5,34.9,0,0)
winterStorm=c(rep(0,41),45.4,14.1,0,123)

NW=1/3*(drought+flooding+wildfire)


#South (S)
drought=c(276.6,0,0,101.4,0,0,70.7,0,315.3,72.7,0,35.1,0,0,0,25.1,88.7,0,99.2,18.1,110.4,0,86.2,25.2,0,4.1,105.9,6.6,58.9,81.1,0,193.1,210.1,104.1,61.6,0,2.2,0,68.1,0,55.3,0,254.9,172.9,0)
flooding=c(95.8,0,0,160.5,0,rep(0,5),77.1,0,0,50.7,62.1,0,0,0,50.3,0,rep(0,8),6.7,0,1,65.6,0,0,0,71.9,438,21.5,0,152.4,0,33.3,0,0,0)
freeze=c(0,0,0,37.6,21.6,rep(0,23),rep(0,9),0.3,rep(0,7))
severeStorm=c(0,38.4,45.9,0,2,0,0,0,0,27.4,0,6.6,144,5.2,43.8,357.5,0,7.4,11.6,83.9,34,8.4,6.7,141.6,1.5,26.2,23.3,21.2,80.9,110.6,100.1,235.2,166.4,174.7,74.1,113.4,267.6,146.2,67.6,146.9,195,176.8,91,409.1,307.3)
tropicalCyclone=c(77.7,0,0,317.6,0,180,0,0,0,47.1,0,0,68,0,0,3,0,0,99.8,0,0,375.9,87.3,0,70.1,5538.3,0,0,1132.5,0,0,1.7,80.7,0,0,0,0,3434.6,0,131.2,738.7,1276.7,0,0,0)
wildfire=c(rep(0,14),7.8,rep(0,8),1.1,0,0,13.7,0,16.5,6.4,0,38.1,1.8,0,0,0.8,rep(0,9))
winterStorm=c(0,0,13.8,0,0,15.2,0,0,0,31.9,0,0,0,23.6,117.8,0,10.2,0,0,12.6,20.8,rep(0,10),19.1,0,0,3,rep(0,6),522.9,28.8,0,17.4)

S=1/5*(drought+flooding+severeStorm+tropicalCyclone+winterStorm)


#South East (SE)
drought=c(250.1,0,0,130,0,0,63.7,0,39.7,0,0,0,0,55.9,0,11.9,0,0,61.2,41,39.6,0,38,0,0,0,17.8,29.4,27.1,4.7,0,51.2,9.7,0,0,0,2.8,0,0,0,0,0,0,1.7,0)
flooding=c(0,0,0,4.8,0,55.1,0,rep(0,13),33.8,0,0,0,0,0,2.5,0,0,23.7,0.9,0,0,0,0,45.2,0.5,0.4,0,0,0,0,0,23,0)
freeze=c(0,60.6,0,97.1,0,97.7,0,0,0,132.4,rep(0,17),15.7,rep(0,9),20.9,rep(0,7))
severeStorm=c(0,8.6,1.2,0,18.3,0,0,0,0,13.8,0,11.3,51.7,0,0,0,0,3.5,11.8,8,0,0,13.8,22.3,1.1,7.1,10.4,7.6,31.1,34.1,3.7,212.7,19.7,24.3,28.5,13.2,19.6,29.7,44.4,10.5,63.6,26.8,30.6,85,50)
tropicalCyclone=c(rep(0,5),76.5,0,0,0,473.9,0,0.6,1401.9,0,49.3,252.7,191.8,0,94.1,199.6,0,0.2,12.3,116,1579.8,872.1,0,0,5,0,0,86.3,16.3,0,0,0,219.2,855.1,1000.4,31.7,239.8,34.2,1901.4,55.9,0)
wildfire=c(rep(0,22),14.7,0,0,0,0.9,5.9,1.8,0,0,1.2,0,0,0,0,2.8,rep(0,8))
winterStorm=c(0,0,10.8,0,0,16.9,0,0,0,6.5,0,0,0.7,156.9,38,0,16.9,0,0,12.8,10.1,rep(0,9),3.8,0,0,0,19.3,8.9,0,0,17.5,0,0,1.6,31,0,1.6)

SE=1/5*(drought+freeze+severeStorm+tropicalCyclone+winterStorm)


#South West (SW)
drought=c(0,0,0,0,0,0,0,0,0,39.5,0,0,0,0,0,4.3,44.4,0,0,0,52.8,0,68.2,74.9,0,0,56.2,6.3,36,12.6,0,54.2,45.2,62.7,13.8,1.2,0,0,7.4,0,41.8,57,97,17.9,0)
flooding=c(0,0,0,384,0,rep(0,28),119.7,0,2.6,0,0,0,0,0,4.9,0,0,0)
freeze=c(0,0,0,13.1,rep(0,41))
severeStorm=c(0,0,13.5,0,139.8,0,11.2,0,0,0,193.8,0,0,0,4.8,0,rep(0,5),1.5,0,1.7,79.4,0,0,0,28.5,90.8,380,36,115.4,0,129,45.1,95.2,223,323.8,92.7,0,23.7,12.1,263.1,90.3)
tropicalCyclone=c(rep(0,28),19.2,rep(0,16))
wildfire=c(rep(0,10),11.7,0,0,0,41.4,rep(0,5),47.3,0,50.3,9.5,0,0,17.7,10.8,14.5,15.3,0,33.4,62.3,0,0,0.8,4.6,2.2,69.1,0,95.5,204.9,63.2,0,0)
winterStorm=c(rep(0,41),7.4,7.7,0,7.2)

SW=1/4*(drought+flooding+severeStorm+wildfire)


#West (W)
drought=c(0,0,0,0,0,0,0,0,0,3.4,0,0,0,0,0,0,0,0,0,0,7.8,0,5.7,0.5,0,0,0,4.8,2.7,5.6,0,0,5.9,7.5,35.5,112.7,99.5,0,4.2,0,9.2,35.2,67.2,0,0)
flooding=c(0,0,0,91.9,0,rep(0,10),156.4,0,136,0,rep(0,18),45.7,0,0,0,28.9,0,112.1,0)
freeze=c(rep(0,10),271.2,rep(0,7),139,rep(0,8),55.9,rep(0,17))
severeStorm=c(rep(0,6),45.1,rep(0,11),32.9,0,rep(0,8),17.8,rep(0,8),10.5,rep(0,7))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),35.1,240.3,0,91.3,11.3,rep(0,7),4.9,147.5,0,0,14.7,73.1,13.9,10.9,0,2.3,6.7,0,0,69.5,16.5,486.9,647.2,118.4,343.5,136.5,35.2,0,0)
winterStorm=c(rep(0,41),1.1,0,0,0)

W=1/4*(drought+flooding+freeze+wildfire)


#West North Central (WNC)
drought=c(3340.8,0,0,0,0,0,0,0,4092.4,605.5,0,0,0,0,0,0,0,0,0,0,267.1,0,855.5,498.2,0,58,595.2,19.4,307.5,0,0,0,1604.9,448.9,0,47.2,0,623.2,0,0,193.5,932.8,807.6,478.3,0)
flooding=c(0,0,0,19.2,rep(0,9),1646.3,0,0,1.1,1227.5,0,rep(0,9),29.7,0,0,345.4,rep(0,7),1237.5,0,0,0,0,0)
freeze=c(0,0,0,31,rep(0,23),3.3,rep(0,17))
severeStorm=c(0,0,10.5,0,19,0,13.8,0,0,0,0,0,0,50.9,29.2,0,0,0,13.7,8.9,0,194.4,0,88.3,66.1,0,5.2,0,57.5,34.8,0,183.1,42.6,268.6,536.9,45.9,328.5,360.6,120.4,187.7,118.3,111.4,741.8,156.5,212.5)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),31.6,rep(0,9),142,0,0,65.4,0,0,79.7,31.9,13.9,0,0,7.6,80.7,0,0,18.6,10.1,204.1,96.1,0,45.7,119.8,10.9,0,0)
winterStorm=c(0,0,9.3,rep(0,39),7.8,0,0)

WNC=1/4*(drought+flooding+severeStorm+wildfire)

##Jitter
set.seed(1)
df=cbind(C,ENC,NE,NW,S,SE,SW,W,WNC)
zero_counts <- apply(df, 2, function(x) sum(x == 0))

# 1. Identify the positions of zeros
zero_positions <- which(df == 0, arr.ind = TRUE)
# 2. Generate random numbers between 0 and 0.1
random_numbers <- runif(nrow(zero_positions), min = 0, max = 0.1)
# 3. Replace zeros with the generated random numbers
df[zero_positions] <- random_numbers





stepsize=0.02
stepnum=2
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)

FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(0.05,0.05,0.05,0.05,0.05,0.15,0.05,0.05), each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:3],1-rowSums(finerGrid),finerGrid[,4:8])



###Remove all negative weights
# Define the function to check rows
is_valid_row <- function(row, tolerance = 1e-6) {
  all(row > tolerance)
}

# Apply the function to each row
valid_rows <- apply(finerGrid, 1, is_valid_row)

# Subset the dataframe to keep only valid rows
df_filtered <- finerGrid[valid_rows, ]

finerGrid=df_filtered
finerGrid=round(finerGrid,digits=2)
dim(finerGrid)#390625




tic()
nCores=parallel::detectCores()-1
cl <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(cl)
NPNPrisk=foreach(i=1:dim(finerGrid)[1], .combine='c',.export=ls(envir=globalenv())) %dopar% risk(finerGrid[i,])
parallel::stopCluster(cl)
toc()


finerGrid[which.min(NPNPrisk),]#0.03 0.01 0.01 0.7 0.01 0.01 0.19 0.03 0.01
min(NPNPrisk)#61.82329
#Corresponding risk: (fill in weight from here into risk function from opt 1!!!!!!!)
risk(c(finerGrid[which.min(NPNPrisk),]))#













###############################################
###############################################
##Optimization 3
###############################################
###############################################

ExtrMargins=function(data){
  tau=0.95
  d_tau <- function(x, tau) {
    (log(0.5) / log(tau) * x^(log(0.5) / log(tau) - 1) * as.numeric(tau > 0.5) * as.numeric(tau <= 1)) * as.numeric(x <= 1)
  }
  
  partialLoss <- function(z, t) {
    -2 * (z - t)
  }
  cdfEvaluated=rank(data)/length(data)
  
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


#Central
drought=c(138.7,0,0,44.8,0,0,15,0,332.7,35.8,0,102.2,0,5.7,0,8.4,0,0,9.8,32.1,7,0,67.9,23.3,0,26.2,0.7,26.3,39.4,0,0,66.4,312.6,24.8,0,0,0,0,3.8,0,0,0,16.1,17.8,0)
flooding=c(0,0,0,0,0,40.1,0,0,0,0,0,0,0,445.3,0,rep(0,13),75.2,0,64.6,31.3,0,29.2,0,0.5,25.1,22.6,0,104.4,0,0,30.6,6.3,0)
freeze=c(0,0,0,18.4,rep(0,23),29.4, rep(0,9),0.5,rep(0,7))
severeStorm=c(0,0,31.8,0,2.3,19.7,0,0,0,0,0,13.6,5.1,6.8,6.6,4.6,0,21.1,23.9,7.1,0,92.6,67.6,119,16.3,0,148.7,0,44.1,67.8,14.8,344.1,218.1,46.8,60.4,60.7,36.2,75.8,27.6,93.8,237,165.6,62.6,312.1,189.7)
tropicalCyclone=c(rep(0,15),3.3,4.2,rep(0,6),3.4,18.5,4.3,0,0,127.9,0,0,1.3,13.2,0,0,0,0,6.3,0,0,5.1,5.1,0,0,0)
wildfire=c(rep(0,36),36.3,rep(0,8))
winterStorm=c(0,0,13.1,0,0,13.3,0,0,0,6.4,0,0,0.7,25.8,23.5,0,9.6,0,0,14.3,1.2,rep(0,9),0.9,14.8,0,0,15.8,13.2,0,0,1.6,0,0,26.6,46.5,0,13.7)

C=cbind(drought,flooding,severeStorm)



#East North Central
drought=c(206.3,0,0,0,0,0,0,0,397,43.4,0,35.4,0,0,0,2.8,0,0,0,0,18.8,0,0,98.9,0,13.6,27.4,49,53,0,0,19.5,258.5,165.9,0,0,0,0,0,0,23.4,59,24.7,67.2,0)
flooding=c(rep(0,6),53.8,rep(0,6),805.5,0,0,0,85.7,rep(0,10),450,0,0,44.8,0,0,44.1,0,0,0,0,242.6,0,0,0,22.6,0)
freeze=c(0,0,0,10.7,rep(0,23),0.3,rep(0,17))
severeStorm=c(0,0,7.3,rep(0,8),8.9,0,32.1,6.3,0,0,0,206,0,0,45.8,0,25.6,17.7,0,93.4,0,168,0.6,32.7,159,16.3,69,42.6,13,2.2,226.1,59.4,7.1,539.4,125.3,351.3,255.3,71.2)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,23),2.7,0,0,0,1.6,rep(0,17))
winterStorm=c(0,0,14.2,0,0,14.3,rep(0,13),2.5,rep(0,14),4.3,5.1,rep(0,6),15.9,0,4.2)

ENC=cbind(drought,flooding,severeStorm)



#North east
drought=c(0,0,0,0,0,0,0,0,26.4,0,0,3.9,0,2.8,0,2.4,0,0,0,11.7,0,0,29.3,0,0,0,0,9.1,6.1,0,0,21.2,18.4,0,0,0,2,0,0,0,0,0,0,0,0)
flooding=c(rep(0,5),5.2,rep(0,10),23.1,rep(0,9),36.5,0,0,0,43.1,0,0,0,4.7,0.7,0,rep(0,6),36.1,0)
freeze=c(0,0,0,6.8,rep(0,41))
severeStorm=c(0,0,2,0,18.7,16.4,0,0,0,0,0,0.8,0,0,0,0,0,9.2,18,0,0,0.8,12.7,7.1,1.1,0,0,41.6,7.8,4.9,12.7,8.9,11,1.4,30.1,10.1,8.7,9.2,27,31.1,21.9,10,8.6,32.2,28.7)
tropicalCyclone=c(rep(0,5),42.6,rep(0,5),60.6,rep(0,4),21,0,0,46.6,0,14.5,0,55.3,65.5,0,0,0,4.2,0,0,276.6,1392.4,0,0,0,0,0,1.5,0,75.6,405.6,0,0,0)
wildfire=c(rep(0,45))
winterStorm=c(0,0,9.6,0,0,9.8,0,0,0,2.9,0,0,95.4,62.1,34.4,0,78,0,46.1,31.6,rep(0,10),15,36.9,0,0,12.3,42.3,0,0,46.7,0,0,1.6,37.6,28.4,3.5)

NE=cbind(tropicalCyclone,winterStorm)


#North west (NW)
drought=c(0,0,0,0,0,0,0,0,179.4,0,0,48.8,0,0,0,0,0,0,0,0,31.8,0,68.6,83.2,0,37.5,0,21.7,71,94.8,0,0,14.6,60.4,67.8,78.9,0,0,11.3,0,6.5,81.8,67.9,62.4,0)
flooding=c(0,0,0,14.2,rep(0,12),207.3,116.3,rep(0,27))
freeze=c(0,0,0,32.8,rep(0,41))
severeStorm=c(rep(0,6),9.3,rep(0,21),8.3,rep(0,16))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),42.5,0,0,0,44.8,rep(0,5),64.5,0,28.7,26.9,0,0,39.4,44.8,13.4,0,0,12.2,47.6,0,0,39.3,22.8,101.2,55.7,0,236.4,133.5,34.9,0,0)
winterStorm=c(rep(0,41),45.4,14.1,0,123)

NW=cbind(drought,flooding,wildfire)


#South (S)
drought=c(276.6,0,0,101.4,0,0,70.7,0,315.3,72.7,0,35.1,0,0,0,25.1,88.7,0,99.2,18.1,110.4,0,86.2,25.2,0,4.1,105.9,6.6,58.9,81.1,0,193.1,210.1,104.1,61.6,0,2.2,0,68.1,0,55.3,0,254.9,172.9,0)
flooding=c(95.8,0,0,160.5,0,rep(0,5),77.1,0,0,50.7,62.1,0,0,0,50.3,0,rep(0,8),6.7,0,1,65.6,0,0,0,71.9,438,21.5,0,152.4,0,33.3,0,0,0)
freeze=c(0,0,0,37.6,21.6,rep(0,23),rep(0,9),0.3,rep(0,7))
severeStorm=c(0,38.4,45.9,0,2,0,0,0,0,27.4,0,6.6,144,5.2,43.8,357.5,0,7.4,11.6,83.9,34,8.4,6.7,141.6,1.5,26.2,23.3,21.2,80.9,110.6,100.1,235.2,166.4,174.7,74.1,113.4,267.6,146.2,67.6,146.9,195,176.8,91,409.1,307.3)
tropicalCyclone=c(77.7,0,0,317.6,0,180,0,0,0,47.1,0,0,68,0,0,3,0,0,99.8,0,0,375.9,87.3,0,70.1,5538.3,0,0,1132.5,0,0,1.7,80.7,0,0,0,0,3434.6,0,131.2,738.7,1276.7,0,0,0)
wildfire=c(rep(0,14),7.8,rep(0,8),1.1,0,0,13.7,0,16.5,6.4,0,38.1,1.8,0,0,0.8,rep(0,9))
winterStorm=c(0,0,13.8,0,0,15.2,0,0,0,31.9,0,0,0,23.6,117.8,0,10.2,0,0,12.6,20.8,rep(0,10),19.1,0,0,3,rep(0,6),522.9,28.8,0,17.4)

S=cbind(drought,flooding,severeStorm,tropicalCyclone,winterStorm)


#South East (SE)
drought=c(250.1,0,0,130,0,0,63.7,0,39.7,0,0,0,0,55.9,0,11.9,0,0,61.2,41,39.6,0,38,0,0,0,17.8,29.4,27.1,4.7,0,51.2,9.7,0,0,0,2.8,0,0,0,0,0,0,1.7,0)
flooding=c(0,0,0,4.8,0,55.1,0,rep(0,13),33.8,0,0,0,0,0,2.5,0,0,23.7,0.9,0,0,0,0,45.2,0.5,0.4,0,0,0,0,0,23,0)
freeze=c(0,60.6,0,97.1,0,97.7,0,0,0,132.4,rep(0,17),15.7,rep(0,9),20.9,rep(0,7))
severeStorm=c(0,8.6,1.2,0,18.3,0,0,0,0,13.8,0,11.3,51.7,0,0,0,0,3.5,11.8,8,0,0,13.8,22.3,1.1,7.1,10.4,7.6,31.1,34.1,3.7,212.7,19.7,24.3,28.5,13.2,19.6,29.7,44.4,10.5,63.6,26.8,30.6,85,50)
tropicalCyclone=c(rep(0,5),76.5,0,0,0,473.9,0,0.6,1401.9,0,49.3,252.7,191.8,0,94.1,199.6,0,0.2,12.3,116,1579.8,872.1,0,0,5,0,0,86.3,16.3,0,0,0,219.2,855.1,1000.4,31.7,239.8,34.2,1901.4,55.9,0)
wildfire=c(rep(0,22),14.7,0,0,0,0.9,5.9,1.8,0,0,1.2,0,0,0,0,2.8,rep(0,8))
winterStorm=c(0,0,10.8,0,0,16.9,0,0,0,6.5,0,0,0.7,156.9,38,0,16.9,0,0,12.8,10.1,rep(0,9),3.8,0,0,0,19.3,8.9,0,0,17.5,0,0,1.6,31,0,1.6)

SE=cbind(drought,freeze,severeStorm,tropicalCyclone,winterStorm)

#South West (SW)
drought=c(0,0,0,0,0,0,0,0,0,39.5,0,0,0,0,0,4.3,44.4,0,0,0,52.8,0,68.2,74.9,0,0,56.2,6.3,36,12.6,0,54.2,45.2,62.7,13.8,1.2,0,0,7.4,0,41.8,57,97,17.9,0)
flooding=c(0,0,0,384,0,rep(0,28),119.7,0,2.6,0,0,0,0,0,4.9,0,0,0)
freeze=c(0,0,0,13.1,rep(0,41))
severeStorm=c(0,0,13.5,0,139.8,0,11.2,0,0,0,193.8,0,0,0,4.8,0,rep(0,5),1.5,0,1.7,79.4,0,0,0,28.5,90.8,380,36,115.4,0,129,45.1,95.2,223,323.8,92.7,0,23.7,12.1,263.1,90.3)
tropicalCyclone=c(rep(0,28),19.2,rep(0,16))
wildfire=c(rep(0,10),11.7,0,0,0,41.4,rep(0,5),47.3,0,50.3,9.5,0,0,17.7,10.8,14.5,15.3,0,33.4,62.3,0,0,0.8,4.6,2.2,69.1,0,95.5,204.9,63.2,0,0)
winterStorm=c(rep(0,41),7.4,7.7,0,7.2)

SW=cbind(drought,flooding,severeStorm,wildfire)


#West (W)
drought=c(0,0,0,0,0,0,0,0,0,3.4,0,0,0,0,0,0,0,0,0,0,7.8,0,5.7,0.5,0,0,0,4.8,2.7,5.6,0,0,5.9,7.5,35.5,112.7,99.5,0,4.2,0,9.2,35.2,67.2,0,0)
flooding=c(0,0,0,91.9,0,rep(0,10),156.4,0,136,0,rep(0,18),45.7,0,0,0,28.9,0,112.1,0)
freeze=c(rep(0,10),271.2,rep(0,7),139,rep(0,8),55.9,rep(0,17))
severeStorm=c(rep(0,6),45.1,rep(0,11),32.9,0,rep(0,8),17.8,rep(0,8),10.5,rep(0,7))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),35.1,240.3,0,91.3,11.3,rep(0,7),4.9,147.5,0,0,14.7,73.1,13.9,10.9,0,2.3,6.7,0,0,69.5,16.5,486.9,647.2,118.4,343.5,136.5,35.2,0,0)
winterStorm=c(rep(0,41),1.1,0,0,0)

W=cbind(drought,flooding,freeze,wildfire)


#West North Central (WNC)
drought=c(3340.8,0,0,0,0,0,0,0,4092.4,605.5,0,0,0,0,0,0,0,0,0,0,267.1,0,855.5,498.2,0,58,595.2,19.4,307.5,0,0,0,1604.9,448.9,0,47.2,0,623.2,0,0,193.5,932.8,807.6,478.3,0)
flooding=c(0,0,0,19.2,rep(0,9),1646.3,0,0,1.1,1227.5,0,rep(0,9),29.7,0,0,345.4,rep(0,7),1237.5,0,0,0,0,0)
freeze=c(0,0,0,31,rep(0,23),3.3,rep(0,17))
severeStorm=c(0,0,10.5,0,19,0,13.8,0,0,0,0,0,0,50.9,29.2,0,0,0,13.7,8.9,0,194.4,0,88.3,66.1,0,5.2,0,57.5,34.8,0,183.1,42.6,268.6,536.9,45.9,328.5,360.6,120.4,187.7,118.3,111.4,741.8,156.5,212.5)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),31.6,rep(0,9),142,0,0,65.4,0,0,79.7,31.9,13.9,0,0,7.6,80.7,0,0,18.6,10.1,204.1,96.1,0,45.7,119.8,10.9,0,0)
winterStorm=c(0,0,9.3,rep(0,39),7.8,0,0)

WNC=cbind(drought,flooding,severeStorm,wildfire)




set.seed(1)
AllColumns=cbind(C,ENC,NE,NW,S,SE,SW,W,WNC)
zero_counts <- apply(AllColumns, 2, function(x) sum(x == 0))

# 1. Identify the positions of zeros
zero_positions <- which(AllColumns == 0, arr.ind = TRUE)
# 2. Generate random numbers between 0 and 0.1
random_numbers <- runif(nrow(zero_positions), min = 0, max = 0.1)
# 3. Replace zeros with the generated random numbers
AllColumns[zero_positions] <- random_numbers



#representative cluster C
toMinimize1=apply(AllColumns[,1:3],2,ExtrMargins)
for (i in 1:3){
  sum=0
  for (k in 4:33){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster ENC
toMinimize2=apply(AllColumns[,4:6],2,ExtrMargins)
for (i in 4:6){
  sum=0
  for (k in c(1:3,7:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3



#representative cluster NE
toMinimize3=apply(AllColumns[,7:8],2,ExtrMargins)
for (i in 7:8){
  sum=0
  for (k in c(1:6,9:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6




#representative cluster NW
toMinimize4=apply(AllColumns[,9:11],2,ExtrMargins)
for (i in 9:11){
  sum=0
  for (k in c(1:8,12:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize4[i-8]=toMinimize4[i-8]*sum
}
clusterRepresentative4=which.min(toMinimize4)+8





#representative cluster S
toMinimize5=apply(AllColumns[,12:16],2,ExtrMargins)
for (i in 12:16){
  sum=0
  for (k in c(1:11,17:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize5[i-11]=toMinimize5[i-11]*sum
}
clusterRepresentative5=which.min(toMinimize5)+11



#representative cluster SE
toMinimize6=apply(AllColumns[,17:21],2,ExtrMargins)
for (i in 17:21){
  sum=0
  for (k in c(1:16,22:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize6[i-16]=toMinimize6[i-16]*sum
}
clusterRepresentative6=which.min(toMinimize6)+16



#representative cluster SW
toMinimize7=apply(AllColumns[,22:25],2,ExtrMargins)
for (i in 22:25){
  sum=0
  for (k in c(1:21,26:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize7[i-21]=toMinimize7[i-21]*sum
}
clusterRepresentative7=which.min(toMinimize7)+21



#representative cluster W
toMinimize8=apply(AllColumns[,26:29],2,ExtrMargins)
for (i in 26:29){
  sum=0
  for (k in c(1:25,30:33)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize8[i-25]=toMinimize8[i-25]*sum
}
clusterRepresentative8=which.min(toMinimize8)+25



#representative cluster WNC
toMinimize9=apply(AllColumns[,30:33],2,ExtrMargins)
for (i in 30:33){
  sum=0
  for (k in c(1:29)){
    sum=sum+abs(cor(AllColumns[,i],AllColumns[,k],method="kendall"))
  }
  toMinimize9[i-29]=toMinimize9[i-29]*sum
}
clusterRepresentative9=which.min(toMinimize9)+29




representatives=AllColumns[,c(clusterRepresentative1,clusterRepresentative2,clusterRepresentative3,
                              clusterRepresentative4,clusterRepresentative5,clusterRepresentative6,
                              clusterRepresentative7,clusterRepresentative8,clusterRepresentative9)]


stepsize=0.02
stepnum=2
basicSeq=seq(-stepsize*stepnum,stepsize*stepnum,stepsize)

FullCombi=expand.grid(basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq,basicSeq)
finerGrid=FullCombi+rep(c(0.05,0.05,0.15,0.05,0.10,0.05,0.15,0.15), each = nrow(FullCombi))
finerGrid=cbind(finerGrid[,1:2],1-rowSums(finerGrid),finerGrid[,c(3:8)])
dim(finerGrid)



apply_Extr <- function(tau, weight_matrix, GiantSampleMatrix) {
  results <- apply(weight_matrix, 1, function(weight) {
    Zcol <- rowSums(t(t(GiantSampleMatrix) * weight))
    
    # Compute empirical CDF
    FecdfEvaluated <- rank(Zcol) / length(Zcol)
    
    # Define d_tau, partialLoss, lambda, and findroot
    d_tau <- function(x, tau) {
      (log(0.5) / log(tau) * x^(log(0.5) / log(tau) - 1) * as.numeric(tau > 0.5) * as.numeric(tau <= 1)) * as.numeric(x <= 1)
    }
    
    partialLoss <- function(z, t) {
      -2 * (z - t)
    }
    
    DFZ <- matrix(d_tau(FecdfEvaluated, tau), nrow = 1)
    
    lambda <- function(t, tau) {
      FecdfEvaluated[FecdfEvaluated == 1] <- 0
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



betaCop  <- empCopula(apply(representatives,2,pobs),smoothing="beta")   
CopulaSample  <-  rCopula(30000, copula = betaCop)
GiantSampleMatrix=toEmpMargins(CopulaSample,x=representatives)
weight_matrix <- finerGrid
riskPerWeight <- apply_Extr(0.95, weight_matrix, GiantSampleMatrix)
Opt3Representative <- finerGrid[which.min(riskPerWeight), ]#0.01 0.01 0.33 0.13 0.01 0.14 0.07 0.13 0.17
ExtrOpt3Representative=min(riskPerWeight)#43.20952
Opt3RepresentativeWeight=rep(0,33)
Opt3RepresentativeWeight[clusterRepresentative1]=as.numeric(Opt3Representative[1])
Opt3RepresentativeWeight[clusterRepresentative2]=as.numeric(Opt3Representative[2])
Opt3RepresentativeWeight[clusterRepresentative3]=as.numeric(Opt3Representative[3])
Opt3RepresentativeWeight[clusterRepresentative4]=as.numeric(Opt3Representative[4])
Opt3RepresentativeWeight[clusterRepresentative5]=as.numeric(Opt3Representative[5])
Opt3RepresentativeWeight[clusterRepresentative6]=as.numeric(Opt3Representative[6])
Opt3RepresentativeWeight[clusterRepresentative7]=as.numeric(Opt3Representative[7])
Opt3RepresentativeWeight[clusterRepresentative8]=as.numeric(Opt3Representative[8])
Opt3RepresentativeWeight[clusterRepresentative9]=as.numeric(Opt3Representative[9])

###############################################
###############################################
##Optimization 4
###############################################
###############################################
finerGrid[which.min(NPNPrisk),]#0.03 0.01 0.01 0.7 0.01 0.01 0.19 0.03 0.01
min(NPNPrisk)#61.82329
#This is obtained from optimization 2







###############################################
###############################################
##Optimization 5
###############################################
###############################################
#Central
drought=c(138.7,0,0,44.8,0,0,15,0,332.7,35.8,0,102.2,0,5.7,0,8.4,0,0,9.8,32.1,7,0,67.9,23.3,0,26.2,0.7,26.3,39.4,0,0,66.4,312.6,24.8,0,0,0,0,3.8,0,0,0,16.1,17.8,0)
flooding=c(0,0,0,0,0,40.1,0,0,0,0,0,0,0,445.3,0,rep(0,13),75.2,0,64.6,31.3,0,29.2,0,0.5,25.1,22.6,0,104.4,0,0,30.6,6.3,0)
freeze=c(0,0,0,18.4,rep(0,23),29.4, rep(0,9),0.5,rep(0,7))
severeStorm=c(0,0,31.8,0,2.3,19.7,0,0,0,0,0,13.6,5.1,6.8,6.6,4.6,0,21.1,23.9,7.1,0,92.6,67.6,119,16.3,0,148.7,0,44.1,67.8,14.8,344.1,218.1,46.8,60.4,60.7,36.2,75.8,27.6,93.8,237,165.6,62.6,312.1,189.7)
tropicalCyclone=c(rep(0,15),3.3,4.2,rep(0,6),3.4,18.5,4.3,0,0,127.9,0,0,1.3,13.2,0,0,0,0,6.3,0,0,5.1,5.1,0,0,0)
wildfire=c(rep(0,36),36.3,rep(0,8))
winterStorm=c(0,0,13.1,0,0,13.3,0,0,0,6.4,0,0,0.7,25.8,23.5,0,9.6,0,0,14.3,1.2,rep(0,9),0.9,14.8,0,0,15.8,13.2,0,0,1.6,0,0,26.6,46.5,0,13.7)

C=0.37*drought+0.42*flooding+0.21*severeStorm



#East North Central
drought=c(206.3,0,0,0,0,0,0,0,397,43.4,0,35.4,0,0,0,2.8,0,0,0,0,18.8,0,0,98.9,0,13.6,27.4,49,53,0,0,19.5,258.5,165.9,0,0,0,0,0,0,23.4,59,24.7,67.2,0)
flooding=c(rep(0,6),53.8,rep(0,6),805.5,0,0,0,85.7,rep(0,10),450,0,0,44.8,0,0,44.1,0,0,0,0,242.6,0,0,0,22.6,0)
freeze=c(0,0,0,10.7,rep(0,23),0.3,rep(0,17))
severeStorm=c(0,0,7.3,rep(0,8),8.9,0,32.1,6.3,0,0,0,206,0,0,45.8,0,25.6,17.7,0,93.4,0,168,0.6,32.7,159,16.3,69,42.6,13,2.2,226.1,59.4,7.1,539.4,125.3,351.3,255.3,71.2)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,23),2.7,0,0,0,1.6,rep(0,17))
winterStorm=c(0,0,14.2,0,0,14.3,rep(0,13),2.5,rep(0,14),4.3,5.1,rep(0,6),15.9,0,4.2)

ENC=0.55*drought+0.18*flooding+0.27*severeStorm



#North east
drought=c(0,0,0,0,0,0,0,0,26.4,0,0,3.9,0,2.8,0,2.4,0,0,0,11.7,0,0,29.3,0,0,0,0,9.1,6.1,0,0,21.2,18.4,0,0,0,2,0,0,0,0,0,0,0,0)
flooding=c(rep(0,5),5.2,rep(0,10),23.1,rep(0,9),36.5,0,0,0,43.1,0,0,0,4.7,0.7,0,rep(0,6),36.1,0)
freeze=c(0,0,0,6.8,rep(0,41))
severeStorm=c(0,0,2,0,18.7,16.4,0,0,0,0,0,0.8,0,0,0,0,0,9.2,18,0,0,0.8,12.7,7.1,1.1,0,0,41.6,7.8,4.9,12.7,8.9,11,1.4,30.1,10.1,8.7,9.2,27,31.1,21.9,10,8.6,32.2,28.7)
tropicalCyclone=c(rep(0,5),42.6,rep(0,5),60.6,rep(0,4),21,0,0,46.6,0,14.5,0,55.3,65.5,0,0,0,4.2,0,0,276.6,1392.4,0,0,0,0,0,1.5,0,75.6,405.6,0,0,0)
wildfire=c(rep(0,45))
winterStorm=c(0,0,9.6,0,0,9.8,0,0,0,2.9,0,0,95.4,62.1,34.4,0,78,0,46.1,31.6,rep(0,10),15,36.9,0,0,12.3,42.3,0,0,46.7,0,0,1.6,37.6,28.4,3.5)

NE=0.02*tropicalCyclone+0.98*winterStorm


#North west (NW)
drought=c(0,0,0,0,0,0,0,0,179.4,0,0,48.8,0,0,0,0,0,0,0,0,31.8,0,68.6,83.2,0,37.5,0,21.7,71,94.8,0,0,14.6,60.4,67.8,78.9,0,0,11.3,0,6.5,81.8,67.9,62.4,0)
flooding=c(0,0,0,14.2,rep(0,12),207.3,116.3,rep(0,27))
freeze=c(0,0,0,32.8,rep(0,41))
severeStorm=c(rep(0,6),9.3,rep(0,21),8.3,rep(0,16))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),42.5,0,0,0,44.8,rep(0,5),64.5,0,28.7,26.9,0,0,39.4,44.8,13.4,0,0,12.2,47.6,0,0,39.3,22.8,101.2,55.7,0,236.4,133.5,34.9,0,0)
winterStorm=c(rep(0,41),45.4,14.1,0,123)

NW=0.37*drought+0.40*flooding+0.23*wildfire


#South (S)
drought=c(276.6,0,0,101.4,0,0,70.7,0,315.3,72.7,0,35.1,0,0,0,25.1,88.7,0,99.2,18.1,110.4,0,86.2,25.2,0,4.1,105.9,6.6,58.9,81.1,0,193.1,210.1,104.1,61.6,0,2.2,0,68.1,0,55.3,0,254.9,172.9,0)
flooding=c(95.8,0,0,160.5,0,rep(0,5),77.1,0,0,50.7,62.1,0,0,0,50.3,0,rep(0,8),6.7,0,1,65.6,0,0,0,71.9,438,21.5,0,152.4,0,33.3,0,0,0)
freeze=c(0,0,0,37.6,21.6,rep(0,23),rep(0,9),0.3,rep(0,7))
severeStorm=c(0,38.4,45.9,0,2,0,0,0,0,27.4,0,6.6,144,5.2,43.8,357.5,0,7.4,11.6,83.9,34,8.4,6.7,141.6,1.5,26.2,23.3,21.2,80.9,110.6,100.1,235.2,166.4,174.7,74.1,113.4,267.6,146.2,67.6,146.9,195,176.8,91,409.1,307.3)
tropicalCyclone=c(77.7,0,0,317.6,0,180,0,0,0,47.1,0,0,68,0,0,3,0,0,99.8,0,0,375.9,87.3,0,70.1,5538.3,0,0,1132.5,0,0,1.7,80.7,0,0,0,0,3434.6,0,131.2,738.7,1276.7,0,0,0)
wildfire=c(rep(0,14),7.8,rep(0,8),1.1,0,0,13.7,0,16.5,6.4,0,38.1,1.8,0,0,0.8,rep(0,9))
winterStorm=c(0,0,13.8,0,0,15.2,0,0,0,31.9,0,0,0,23.6,117.8,0,10.2,0,0,12.6,20.8,rep(0,10),19.1,0,0,3,rep(0,6),522.9,28.8,0,17.4)

S=0.36*drought+0.30*flooding+0.09*severeStorm+0.01*tropicalCyclone+0.24*winterStorm


#South East (SE)
drought=c(250.1,0,0,130,0,0,63.7,0,39.7,0,0,0,0,55.9,0,11.9,0,0,61.2,41,39.6,0,38,0,0,0,17.8,29.4,27.1,4.7,0,51.2,9.7,0,0,0,2.8,0,0,0,0,0,0,1.7,0)
flooding=c(0,0,0,4.8,0,55.1,0,rep(0,13),33.8,0,0,0,0,0,2.5,0,0,23.7,0.9,0,0,0,0,45.2,0.5,0.4,0,0,0,0,0,23,0)
freeze=c(0,60.6,0,97.1,0,97.7,0,0,0,132.4,rep(0,17),15.7,rep(0,9),20.9,rep(0,7))
severeStorm=c(0,8.6,1.2,0,18.3,0,0,0,0,13.8,0,11.3,51.7,0,0,0,0,3.5,11.8,8,0,0,13.8,22.3,1.1,7.1,10.4,7.6,31.1,34.1,3.7,212.7,19.7,24.3,28.5,13.2,19.6,29.7,44.4,10.5,63.6,26.8,30.6,85,50)
tropicalCyclone=c(rep(0,5),76.5,0,0,0,473.9,0,0.6,1401.9,0,49.3,252.7,191.8,0,94.1,199.6,0,0.2,12.3,116,1579.8,872.1,0,0,5,0,0,86.3,16.3,0,0,0,219.2,855.1,1000.4,31.7,239.8,34.2,1901.4,55.9,0)
wildfire=c(rep(0,22),14.7,0,0,0,0.9,5.9,1.8,0,0,1.2,0,0,0,0,2.8,rep(0,8))
winterStorm=c(0,0,10.8,0,0,16.9,0,0,0,6.5,0,0,0.7,156.9,38,0,16.9,0,0,12.8,10.1,rep(0,9),3.8,0,0,0,19.3,8.9,0,0,17.5,0,0,1.6,31,0,1.6)

SE=0.16*drought+0.29*freeze+0.14*severeStorm+0.01*tropicalCyclone+0.40*winterStorm

#South West (SW)
drought=c(0,0,0,0,0,0,0,0,0,39.5,0,0,0,0,0,4.3,44.4,0,0,0,52.8,0,68.2,74.9,0,0,56.2,6.3,36,12.6,0,54.2,45.2,62.7,13.8,1.2,0,0,7.4,0,41.8,57,97,17.9,0)
flooding=c(0,0,0,384,0,rep(0,28),119.7,0,2.6,0,0,0,0,0,4.9,0,0,0)
freeze=c(0,0,0,13.1,rep(0,41))
severeStorm=c(0,0,13.5,0,139.8,0,11.2,0,0,0,193.8,0,0,0,4.8,0,rep(0,5),1.5,0,1.7,79.4,0,0,0,28.5,90.8,380,36,115.4,0,129,45.1,95.2,223,323.8,92.7,0,23.7,12.1,263.1,90.3)
tropicalCyclone=c(rep(0,28),19.2,rep(0,16))
wildfire=c(rep(0,10),11.7,0,0,0,41.4,rep(0,5),47.3,0,50.3,9.5,0,0,17.7,10.8,14.5,15.3,0,33.4,62.3,0,0,0.8,4.6,2.2,69.1,0,95.5,204.9,63.2,0,0)
winterStorm=c(rep(0,41),7.4,7.7,0,7.2)

SW=0.69*drought+0.12*flooding+0.10*severeStorm+0.09*wildfire


#West (W)
drought=c(0,0,0,0,0,0,0,0,0,3.4,0,0,0,0,0,0,0,0,0,0,7.8,0,5.7,0.5,0,0,0,4.8,2.7,5.6,0,0,5.9,7.5,35.5,112.7,99.5,0,4.2,0,9.2,35.2,67.2,0,0)
flooding=c(0,0,0,91.9,0,rep(0,10),156.4,0,136,0,rep(0,18),45.7,0,0,0,28.9,0,112.1,0)
freeze=c(rep(0,10),271.2,rep(0,7),139,rep(0,8),55.9,rep(0,17))
severeStorm=c(rep(0,6),45.1,rep(0,11),32.9,0,rep(0,8),17.8,rep(0,8),10.5,rep(0,7))
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),35.1,240.3,0,91.3,11.3,rep(0,7),4.9,147.5,0,0,14.7,73.1,13.9,10.9,0,2.3,6.7,0,0,69.5,16.5,486.9,647.2,118.4,343.5,136.5,35.2,0,0)
winterStorm=c(rep(0,41),1.1,0,0,0)

W=0.49*drought+0.32*flooding+0.18*freeze+0.01*wildfire


#West North Central (WNC)
drought=c(3340.8,0,0,0,0,0,0,0,4092.4,605.5,0,0,0,0,0,0,0,0,0,0,267.1,0,855.5,498.2,0,58,595.2,19.4,307.5,0,0,0,1604.9,448.9,0,47.2,0,623.2,0,0,193.5,932.8,807.6,478.3,0)
flooding=c(0,0,0,19.2,rep(0,9),1646.3,0,0,1.1,1227.5,0,rep(0,9),29.7,0,0,345.4,rep(0,7),1237.5,0,0,0,0,0)
freeze=c(0,0,0,31,rep(0,23),3.3,rep(0,17))
severeStorm=c(0,0,10.5,0,19,0,13.8,0,0,0,0,0,0,50.9,29.2,0,0,0,13.7,8.9,0,194.4,0,88.3,66.1,0,5.2,0,57.5,34.8,0,183.1,42.6,268.6,536.9,45.9,328.5,360.6,120.4,187.7,118.3,111.4,741.8,156.5,212.5)
tropicalCyclone=c(rep(0,45))
wildfire=c(rep(0,10),31.6,rep(0,9),142,0,0,65.4,0,0,79.7,31.9,13.9,0,0,7.6,80.7,0,0,18.6,10.1,204.1,96.1,0,45.7,119.8,10.9,0,0)
winterStorm=c(0,0,9.3,rep(0,39),7.8,0,0)

WNC=0.01*drought+0.04*flooding+0.01*severeStorm+0.94*wildfire




set.seed(1)
df=cbind(C,ENC,NE,NW,S,SE,SW,W,WNC)
zero_counts <- apply(df, 2, function(x) sum(x == 0))

# 1. Identify the positions of zeros
zero_positions <- which(df == 0, arr.ind = TRUE)
# 2. Generate random numbers between 0 and 0.1
random_numbers <- runif(nrow(zero_positions), min = 0, max = 0.1)
# 3. Replace zeros with the generated random numbers
df[zero_positions] <- random_numbers




##Extremile: NP+NP 
ExtrNPNP=function(tau,data){
  library(copula)
  set.seed(1)
  X=data
  u=pobs(X)
  ecop.orig  <- empCopula(u,smoothing="beta")   
  U.orig  <-  rCopula(100000, copula = ecop.orig)
  
  GiantSampleMatrix=toEmpMargins(U.orig,x=X)
  Zcol=rowSums(GiantSampleMatrix)
  Fecdf=ecdf(Zcol)
  FecdfEvaluated=Fecdf(Zcol)
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




#NP+NP 
risk=function(x){
  return(ExtrNPNP(0.95,data=as.matrix(df)%*%diag(x)))}

risk(rep(1/9,9))#58.58659





