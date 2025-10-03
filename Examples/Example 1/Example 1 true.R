# R script accompanying Section 5.1 (Example 1: multivariate normal distribution) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we compute the optimal portfolios for the five settings of Table 1, 
#and the misspecification models from the supplementary material


library(NMOF)


###################################################################################
############### Setting 1: high in clusters, low between  
###################################################################################
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 1##
##################

### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 3##
##################
#representative cluster 1
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

toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6





correlationsRepresentatives=c(0.1,0.1,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 4##
##################
#Between: 0.3473537 0.4789003 0.1737460 #via optimization 2

weightsOpt4=c(0.3473537*c(1/3,1/3,1/3),0.4789003*c(1/3,1/3,1/3),0.1737460*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 5##
##################
weightsOpt5=c(1/3*c(0.5,0.0,0.5), 1/3*c(1/3,1/3,1/3), 1/3*c(0.5,0.5)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)


######
#Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)

#double step: Optimization 4
round(weightsOpt4,digits=2)
round(ESOpt4,digits=2)

#double step: Optimization 5
round(weightsOpt5,digits=2)
round(ESOpt5,digits=2)











###################################################################################
############### Setting 2: low in clusters, low between clusters
###################################################################################
between=0.1
group1=0.3
group2=0.3
group3=0.2
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)







##################
##Optimization 3##
##################
between=0.1
group1=0.3
group2=0.3
group3=0.2
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

#representative cluster 1
toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6





correlationsRepresentatives=c(0.1,0.1,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 4##
##################
#Between:  #via optimization 2

weightsOpt4=c(Opt2Between[1]*c(1/3,1/3,1/3),Opt2Between[2]*c(1/3,1/3,1/3),Opt2Between[3]*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 5##
##################

weightsOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)



######
##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)

#double step: Optimization 4
round(weightsOpt4,digits = 2)
round(ESOpt4,digits=3)

#double step: Optimization 5
round(weightsOpt5,digits = 2)
round(ESOpt5,digits=3)











###################################################################################
############### Setting 3: high in clusters, high between clusters
###################################################################################
between=0.1
group1=0.8
group2=0.7
group3=0.6
correlations=c(group1,group1,c(0.5,0.5,0.4,0.4,0.3),
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)



##################
##Optimization 3##
##################
between=0.1
group1=0.8
group2=0.7
group3=0.6
correlations=c(group1,group1,c(0.5,0.5,0.4,0.4,0.3),
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

#representative cluster 1
toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6





correlationsRepresentatives=c(0.1,0.1,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)












##################
##Optimization 4##
##################
weightsOpt4=c(Opt2Between[1]*c(1/3,1/3,1/3),Opt2Between[2]*c(1/3,1/3,1/3),Opt2Between[3]*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 5##
##################
weightsOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)



######
##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)

#double step: Optimization 4
round(weightsOpt4,digits = 2)
round(ESOpt4,digits=3)

#double step: Optimization 5
round(weightsOpt5,digits = 2)
round(ESOpt5,digits=3)









###################################################################################
############### Setting 4: low in clusters, high in between clusters
###################################################################################
between=0.1
group1=0.3
group2=0.3
group3=0.2
correlations=c(group1,group1,c(0.5,0.5,0.4,0.4,0.3),
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 3##
##################
between=0.1
group1=0.3
group2=0.3
group3=0.2
correlations=c(group1,group1,c(0.5,0.5,0.4,0.4,0.3),
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

#representative cluster 1
toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6





correlationsRepresentatives=c(0.1,0.1,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 4##
##################
weightsOpt4=c(Opt2Between[1]*c(1/3,1/3,1/3),Opt2Between[2]*c(1/3,1/3,1/3),Opt2Between[3]*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 5##
##################
weightsOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)


######
##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)

#double step: Optimization 4
round(weightsOpt4,digits = 2)
round(ESOpt4,digits=3)


#double step: Optimization 5
round(weightsOpt5,digits = 2)
round(ESOpt5,digits=3)





###################################################################################
############### Setting 5
###################################################################################
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,5,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 3##
##################
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

#representative cluster 1
toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  print(sum)
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6





correlationsRepresentatives=c(0.1,0.1,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)



##################
##Optimization 4##
##################
weightsOpt4=c(Opt2Between[1]*c(1/3,1/3,1/3),Opt2Between[2]*c(1/3,1/3,1/3),Opt2Between[3]*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

##################
##Optimization 5##
##################
weightsOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)


######
##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)

#double step: Optimization 4
round(weightsOpt4,digits = 2)
round(ESOpt4,digits=3)

#double step: Optimization 5
round(weightsOpt5,digits = 2)
round(ESOpt5,digits=3)




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

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[upper.tri(cor_matrix)] <- c(cor_matrixBig[1,5],cor_matrixBig[1,7],cor_matrixBig[5,7])
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1],std_devs[5],std_devs[7])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[upper.tri(cor_matrix)] <- c(cor_matrixBig[2,4],cor_matrixBig[2,6],cor_matrixBig[4,6])
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[2],std_devs[4],std_devs[6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[upper.tri(cor_matrix)] <- c(cor_matrixBig[3,8])
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[3],std_devs[8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[5]^2+cluster1[3]^2*std_devs[7]^2+
  2*cluster1[1]*cluster1[2]*cov_matrix[1,5]+
  2*cluster1[1]*cluster1[3]*cov_matrix[1,7]+
  2*cluster1[2]*cluster1[3]*cov_matrix[5,7]

varS2=cluster2[1]^2*std_devs[2]^2+cluster2[2]^2*std_devs[4]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*cov_matrix[2,4]+
  2*cluster2[1]*cluster2[3]*cov_matrix[2,6]+
  2*cluster2[2]*cluster2[3]*cov_matrix[4,6]


varS3=cluster3[1]^2*std_devs[3]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*cov_matrix[3,8]


#Between cluster 1 and cluster 2 Sums
combinations=expand.grid(c(1,5,7),c(2,4,6))
covarianceS1S2=0
weights=c(cluster1[1],cluster2[1],0,cluster2[2],cluster1[2],cluster2[3],cluster1[3])
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



#Between cluster 1 and cluster 3 Sums
combinations=expand.grid(c(1,5,7),c(3,8))
weights=c(cluster1[1],0,cluster3[1],0,cluster1[2],0,cluster1[3],cluster3[2])

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}


#Between cluster 2 and cluster 3 Sums
combinations=expand.grid(c(2,4,6),c(3,8))
covarianceS2S3=0
weights=c(0,cluster2[1],cluster3[1],cluster2[2],0,cluster2[3],0,cluster3[2])
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)
weightsOrderedByCluster=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])
weightsOpt1=c(weightsOrderedByCluster[1],weightsOrderedByCluster[4],weightsOrderedByCluster[7],
              weightsOrderedByCluster[5],weightsOrderedByCluster[2],weightsOrderedByCluster[6],
              weightsOrderedByCluster[3],weightsOrderedByCluster[8])




##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)
#between clusters
varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[5]^2+cluster1Eq[3]^2*std_devs[7]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*cov_matrix[1,5]+
  2*cluster1Eq[1]*cluster1Eq[3]*cov_matrix[1,7]+
  2*cluster1Eq[2]*cluster1Eq[3]*cov_matrix[5,7]

varS2=cluster2Eq[1]^2*std_devs[2]^2+cluster2Eq[2]^2*std_devs[4]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*cov_matrix[2,4]+
  2*cluster2Eq[1]*cluster2Eq[3]*cov_matrix[2,6]+
  2*cluster2Eq[2]*cluster2Eq[3]*cov_matrix[4,6]


varS3=cluster3Eq[1]^2*std_devs[3]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*cov_matrix[3,8]


#Between cluster 1 and cluster 2 Sums
combinations=expand.grid(c(1,5,7),c(2,4,6))
covarianceS1S2=0
weights=c(cluster1Eq[1],cluster2Eq[1],0,cluster2Eq[2],cluster1Eq[2],cluster2Eq[3],cluster1Eq[3])
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



#Between cluster 1 and cluster 3 Sums
combinations=expand.grid(c(1,5,7),c(3,8))
weights=c(cluster1Eq[1],0,cluster3Eq[1],0,cluster1Eq[2],0,cluster1Eq[3],cluster3Eq[2])

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



#Between cluster 2 and cluster 3 Sums
combinations=expand.grid(c(2,4,6),c(3,8))
covarianceS2S3=0
weights=c(0,cluster2Eq[1],cluster3Eq[1],cluster2Eq[2],0,cluster2Eq[3],0,cluster3Eq[2])
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOrderedByCluster=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
weightsOpt2=c(weightsOrderedByCluster[1],weightsOrderedByCluster[4],weightsOrderedByCluster[7],
              weightsOrderedByCluster[5],weightsOrderedByCluster[2],weightsOrderedByCluster[6],
              weightsOrderedByCluster[3],weightsOrderedByCluster[8])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)





##################
##Optimization 3##
##################
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

#representative cluster 1
toMinimize1=std_devs[c(1,5,7)]
cluster1Numbers=c(1,5,7)
for (i in 1:3){
  m=cluster1Numbers[i]
  sum=0
  for (j in c(2,3,4,6,8)){
    sum=sum+abs(asin(cor_matrix[m,j])*2/pi)
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=cluster1Numbers[which.min(toMinimize1)]


#representative cluster 2
toMinimize2=std_devs[c(2,4,6)]
cluster2Numbers=c(2,4,6)
for (i in 1:3){
  m=cluster2Numbers[i]
  sum=0
  for (j in c(1,3,5,7,8)){
    sum=sum+abs(asin(cor_matrix[m,j])*2/pi)
  }
  toMinimize2[i]=toMinimize2[i]*sum
}
clusterRepresentative2=cluster2Numbers[which.min(toMinimize2)]



#representative cluster 3
toMinimize3=std_devs[c(3,8)]
cluster3Numbers=c(3,8)
for (i in 1:2){
  m=cluster3Numbers[i]
  sum=0
  for (j in c(1,2,4,5,6,7)){
    sum=sum+abs(asin(cor_matrix[m,j])*2/pi)
  }
  toMinimize3[i]=toMinimize3[i]*sum
}
clusterRepresentative3=cluster3Numbers[which.min(toMinimize3)]


#representatives: 7,4,8


correlationsRepresentatives=c(0.1,0.6,0.1)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(3,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)
WeightsRepresentativesOpt3=c(0,0,0,WeigthsRepresentatives[2],0,0,WeigthsRepresentatives[1],WeigthsRepresentatives[3])



##################
##Optimization 4##
##################
weightsOrderedByCluster=c(c(1/3,1/3,1/3)*Opt2Between[1],c(1/3,1/3,1/3)*Opt2Between[2],c(1/2,1/2)*Opt2Between[3])
weightsOpt4=c(weightsOrderedByCluster[1],weightsOrderedByCluster[4],weightsOrderedByCluster[7],
              weightsOrderedByCluster[5],weightsOrderedByCluster[2],weightsOrderedByCluster[6],
              weightsOrderedByCluster[3],weightsOrderedByCluster[8])
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)

round(weightsOpt4,digits=2)

##################
##Optimization 5##
##################
weightsOrderedByClusterOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
weightsOpt5=c(weightsOrderedByClusterOpt5[1],weightsOrderedByClusterOpt5[4],weightsOrderedByClusterOpt5[7],
              weightsOrderedByClusterOpt5[5],weightsOrderedByClusterOpt5[2],weightsOrderedByClusterOpt5[6],
              weightsOrderedByClusterOpt5[3],weightsOrderedByClusterOpt5[8])

varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)


######
##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
round(weightsOpt1,digits=2)
round(ESOpt1,digits=3)

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(weightsOpt3,digits = 2)
round(ESOpt3,digits=3)







###################################################################################
##################################### misspecification 2  #########################
###################################################################################
between=0.5
group1=0.5
group2=0.5
group3=0.5
correlations=c(group1,group1,rep(between,5),
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

# Display the correlation matrix
print(cor_matrix)
eigen(cor_matrix)$values

# Compute standard deviations from variances
std_devs=c(2,2.5,2,2,2,2,3,3)
# Compute the diagonal matrix of standard deviations
diag_std_devs <- diag(std_devs)

# Compute the covariance matrix
cov_matrix <- diag_std_devs %*% cor_matrix %*% diag_std_devs

minvariance=minvar(cov_matrix)
singlestep=c(minvariance[1:8])
varianceSingle=t(singlestep)%*%cov_matrix%*%singlestep
ESSingle=sqrt(varianceSingle)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 1##
##################
### cluster 1
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group1,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs1 <- c(std_devs[1:3])
diag_std_devs <- diag(std_devs1)
cov_matrix1 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix1)
cluster1=c(minvariance[1:3])


### cluster 2
d <- 3
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group2,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[4:6])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster2=c(minvariance[1:3])



### cluster 3
d <- 2
cor_matrix <- diag(1, d)
cor_matrix[lower.tri(cor_matrix)] <- rep(group3,3)
cor_matrix <- cor_matrix + t(cor_matrix) - diag(1, d)
std_devs2 <- c(std_devs[7:8])
diag_std_devs <- diag(std_devs2)
cov_matrix2 <- diag_std_devs %*% cor_matrix %*% diag_std_devs
minvariance=minvar(cov_matrix2)
cluster3=c(minvariance[1:2])


#between clusters
varS1=cluster1[1]^2*std_devs[1]^2+cluster1[2]^2*std_devs[2]^2+cluster1[3]^2*std_devs[3]^2+
  2*cluster1[1]*cluster1[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1[1]*cluster1[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1[2]*cluster1[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2[1]^2*std_devs[4]^2+cluster2[2]^2*std_devs[5]^2+cluster2[3]^2*std_devs[6]^2+
  2*cluster2[1]*cluster2[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2[1]*cluster2[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2[2]*cluster2[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3[1]^2*std_devs[7]^2+cluster3[2]^2*std_devs[8]^2+2*cluster3[1]*cluster3[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1,cluster2)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1,cluster3)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1,rep(0,3),cluster3)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2,cluster3)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2,cluster3)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt1Between=c(minvariance[1:3])
varianceOpt1=t(Opt1Between)%*%cov_matrixBetweenClusters%*%Opt1Between
ESOpt1=sqrt(varianceOpt1)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 2##
##################
#Equal weights in clusters
cluster1Eq=rep(1/3,3)
cluster2Eq=rep(1/3,3)
cluster3Eq=rep(1/2,2)

varS1=cluster1Eq[1]^2*std_devs[1]^2+cluster1Eq[2]^2*std_devs[2]^2+cluster1Eq[3]^2*std_devs[3]^2+
  2*cluster1Eq[1]*cluster1Eq[2]*correlations[1]*std_devs[1]*std_devs[2]+
  2*cluster1Eq[1]*cluster1Eq[3]*correlations[2]*std_devs[1]*std_devs[3]+
  2*cluster1Eq[2]*cluster1Eq[3]*correlations[8]*std_devs[2]*std_devs[3]

varS2=cluster2Eq[1]^2*std_devs[4]^2+cluster2Eq[2]^2*std_devs[5]^2+cluster2Eq[3]^2*std_devs[6]^2+
  2*cluster2Eq[1]*cluster2Eq[2]*correlations[19]*std_devs[4]*std_devs[5]+
  2*cluster2Eq[1]*cluster2Eq[3]*correlations[20]*std_devs[4]*std_devs[6]+
  2*cluster2Eq[2]*cluster2Eq[3]*correlations[23]*std_devs[5]*std_devs[6]


varS3=cluster3Eq[1]^2*std_devs[7]^2+cluster3Eq[2]^2*std_devs[8]^2+2*cluster3Eq[1]*cluster3Eq[2]*correlations[length(correlations)]*std_devs[7]*std_devs[8]


weights=c(cluster1Eq,cluster2Eq)
combinations=expand.grid(c(1:3),c(4:6))
covarianceS1S2=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S2=covarianceS1S2+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster1Eq,cluster3Eq)
combinations=expand.grid(c(1:3),c(7:8))
weights=c(cluster1Eq,rep(0,3),cluster3Eq)

covarianceS1S3=0
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS1S3=covarianceS1S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}



weights=c(cluster2Eq,cluster3Eq)
combinations=expand.grid(c(4:6),c(7:8))
covarianceS2S3=0
weights=c(rep(0,3),cluster2Eq,cluster3Eq)
for (i in 1:dim(combinations)[1]){
  row=as.numeric(combinations[i,])
  covarianceS2S3=covarianceS2S3+weights[row[1]]*weights[row[2]]*cov_matrix[row[1],row[2]]
}

cov_matrixBetweenClusters=matrix(c(varS1,covarianceS1S2,covarianceS1S3,
                                   covarianceS1S2,varS2,covarianceS2S3,
                                   covarianceS1S3,covarianceS2S3,varS3),nrow = 3, byrow = TRUE)
minvariance=minvar(cov_matrixBetweenClusters)
Opt2Between=c(minvariance[1:3])
weightsOpt2=c(cluster1*Opt2Between[1],cluster2*Opt2Between[2],cluster3*Opt2Between[3])
varianceOpt2=t(weightsOpt2)%*%cov_matrix%*%weightsOpt2
ESOpt2=sqrt(varianceOpt2)*dnorm(qnorm(0.95))/(1-0.95)














##################
##Optimization 3##
#################
between=0.5
group1=0.5
group2=0.5
group3=0.5
correlations=c(group1,group1,rep(between,5),
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


#representative cluster 1

toMinimize1=std_devs[1:3]
for (i in 1:3){
  sum=0
  for (j in 4:8){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize1[i]=toMinimize1[i]*sum
}
clusterRepresentative1=which.min(toMinimize1)


#representative cluster 2
toMinimize2=std_devs[4:6]
for (i in 4:6){
  sum=0
  for (j in c(1,2,3,7,8)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize2[i-3]=toMinimize2[i-3]*sum
}
clusterRepresentative2=which.min(toMinimize2)+3


#representative cluster 3
toMinimize3=std_devs[7:8]
for (i in 7:8){
  sum=0
  for (j in c(1,2,3,4,5,6)){
    sum=sum+asin(cor_matrix[i,j])*2/pi
  }
  toMinimize3[i-6]=toMinimize3[i-6]*sum
}
clusterRepresentative3=which.min(toMinimize3)+6




#representatives 1,4,7
correlationsRepresentatives=c(0.5,0.5,0.5)
cor_matrixRepresentatives <- diag(1, 3)
cor_matrixRepresentatives[lower.tri(cor_matrixRepresentatives)] <- correlationsRepresentatives
cor_matrixRepresentatives <- cor_matrixRepresentatives + t(cor_matrixRepresentatives) - diag(1, 3)
std_devsRepresentatives=c(2,2,3)
diag_std_devsRepresentatives <- diag(std_devsRepresentatives)
cov_matrixRepresentatives <- diag_std_devsRepresentatives %*% cor_matrixRepresentatives %*% diag_std_devsRepresentatives

minvarianceRepresentatives=minvar(cov_matrixRepresentatives)


WeigthsRepresentatives=c(minvarianceRepresentatives[1:3])
varianceRepresentatives=t(WeigthsRepresentatives)%*%cov_matrixRepresentatives%*%WeigthsRepresentatives
ESRepresentatives=sqrt(varianceRepresentatives)*dnorm(qnorm(0.95))/(1-0.95)




##################
##Optimization 4##
##################
weightsOpt4=c(Opt2Between[1]*c(1/3,1/3,1/3),Opt2Between[2]*c(1/3,1/3,1/3),Opt2Between[3]*c(1/2,1/2))
varianceOpt4=t(weightsOpt4)%*%cov_matrix%*%weightsOpt4
ESOpt4=sqrt(varianceOpt4)*dnorm(qnorm(0.95))/(1-0.95)


##################
##Optimization 5##
##################

weightsOpt5=c(1/3*c(cluster1), 1/3*c(cluster2), 1/3*c(cluster3)) #Via optimization 1
varianceOpt5=t(weightsOpt5)%*%cov_matrix%*%weightsOpt5
ESOpt5=sqrt(varianceOpt5)*dnorm(qnorm(0.95))/(1-0.95)




##Summary
####single step
round(singlestep,digits=2)
round(ESSingle,digits = 3)


#double step: Optimization 1
round(c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3]),digits = 2)
round(ESOpt1,digits=3)
weightsOpt1=c(cluster1*Opt1Between[1],cluster2*Opt1Between[2],cluster3*Opt1Between[3])

#double step: Optimization 2
round(weightsOpt2,digits=2)
round(ESOpt2,digits=3)

#double step: Optimization 3
round(WeigthsRepresentatives,digits = 2)
round(ESRepresentatives,digits=3)

