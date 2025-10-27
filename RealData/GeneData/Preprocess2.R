load('Generesult1.RData')
library(GEOquery)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)
library(impute)
set.seed(6)
confounder<-c(3)
remain<-setdiff(c(1:p),confounder)
cfd2<-confounder
cnt=p
cnt=p
for (j in 1:p){
  for (i in 1:j){
    cnt=cnt+1
    if (i %in% confounder){
      cfd2<-c(cfd2,cnt)
    }
    else if (j %in% confounder){
      cfd2<-c(cfd2,cnt)
    }
  }
}
remain2<-setdiff(c(1:q),cfd2)
data=data[remain,]
Int_X<-Int_X[remain2,]
IntX_centered<-IntX_centered[remain2,]
data_centered<-data_centered[remain,]
Astar<-Ahat[remain,remain2]
Astar0<-Ahat_0[remain,remain2]
p<-12
q<-p*(p+3)/2
X_0<-rep(0,p)
D_01<-matrix(0,nrow=p,ncol=q)
D_02<-matrix(0,nrow=p,ncol=q)
result<-Dhat_cv(data_centered, IntX_centered,D_01,D_02,X_0, 1e-5, 50,p,q,n,Lambda_list)
X0<-result[,1]
Dhat<-result[,2:(q+1)]
data_cha<-data[,2:n]-data[,1:(n-1)]
sigmahat<-c()
for (i in 1:p){
  sigmahat<-c(sigmahat,sum(data_cha[i,]^2)/(n-1))
}
Res<-matrix(nrow=p,ncol=n)
for (i in 1:n){
  Res[,i]<-data_centered[,i]-X0-Dhat %*% IntX_centered[,i]
}
pca=prcomp(t(Res))

Res_scaled=diag((sigmahat)^(-0.5))%*%Res
pca=prcomp(t(Res_scaled))
lambda_ratio<-(pca$sdev[1:ceiling(sqrt(p))]-pca$sdev[2:(ceiling(sqrt(p))+1)])/(pca$sdev[2:(ceiling(sqrt(p))+1)]-pca$sdev[3:(ceiling(sqrt(p))+2)])
dhat<-order(lambda_ratio,decreasing = TRUE)[1]
dhat
u=diag((sigmahat)^(0.5))%*%(pca$rotation[,1:dhat])
proj<-diag(p)-u %*% ginv(t(u) %*% u) %*% t(u)
proj_y= proj %*% data
save.image('gene_data1.RData')
