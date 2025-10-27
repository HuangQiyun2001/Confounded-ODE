load('ncov_result1.RData')
library(nCov2019)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)
set.seed(6)
selected<-setdiff(c(1:33),c(10,17))
Astar_0<-Ahat_0[selected,selected]
Astar=Ahat[selected,selected]

data<-data_select[selected,]
p=31
Int_X<-matrix(nrow=p, ncol=n)
for (i in 1:n){
  for (j in 1:p){
    data_0<-data[j,1:(i)]
    data_0[1]<-(data_0[1]+data_0[i])/2
    Int_X[j,i]<-sum(data_0[1:(i-1)])/n
  }
}


D_01<-matrix(0,nrow=p,ncol=p)
D_02<-matrix(0,nrow=p,ncol=p)
IntX_centered<-matrix(nrow=p, ncol=n)
for (i in 1:n){
  for (j in 1:p){
    IntX_centered[j,i]<-Int_X[j,i]-mean(Int_X[j,])
  }
}
data_centered<-matrix(nrow=p, ncol=n)
for (i in 1:n){
  for (j in 1:p){
    data_centered[j,i]<-data[j,i]-mean(data[j,])
  }
}

X_0<-rep(0,p)
Lambda_list<-0.001*1.5^c(1:15)
result<-Dhat_cv(data_centered, IntX_centered,D_01,D_02,X_0, 1e-5, 50,p,n,Lambda_list)
X0<-result[,1]
Dhat<-result[,2:(p+1)]
data_cha<-data[,2:n]-data[,1:(n-1)]
sigmahat<-c()
for (i in 1:p){
  sigmahat<-c(sigmahat,sum(data_cha[i,]^2)/(n-1))
}
Res<-matrix(nrow=p,ncol=n)
for (i in 1:n){
  Res[,i]<-data_centered[,i]-X0-Dhat %*% IntX_centered[,i]
}

Res_scaled=diag((sigmahat)^(-0.5))%*%Res
pca=prcomp(t(Res_scaled))
lambda_ratio<-(pca$sdev[1:ceiling(sqrt(p))]-pca$sdev[2:(ceiling(sqrt(p))+1)])/(pca$sdev[2:(ceiling(sqrt(p))+1)]-pca$sdev[3:(ceiling(sqrt(p))+2)])
dhat<-order(lambda_ratio,decreasing = TRUE)[1]
dhat
u=diag((sigmahat)^(0.5))%*%(pca$rotation[,1:dhat])
proj<-diag(p)-u %*% ginv(t(u) %*% u) %*% t(u)
proj_y= proj %*% data

save.image('ncov_data2.RData')


