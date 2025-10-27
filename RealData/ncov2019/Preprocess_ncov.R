library(nCov2019)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)
set.seed(2)
res <- query()
y<-res$historical
y_country<-y$table
country<-levels(factor(y_country$country))
country_selected<-c("Iceland","Finland","Sweden","Norway","Denmark","Poland",
                    "Czechia","Slovakia","Hungary","Germany",
                    "Austria","Switzerland","Liechtenstein","Netherlands","Belgium","Luxembourg","France","Monaco"
                    ,"Romania","Bulgaria","Serbia","Macedonia","Albania"
                    ,"Greece","Slovenia","Croatia","Bosnia","Italy"
                   ,"San Marino","Malta","Spain","Portugal","Andorra"
                    )

data_ncov<-matrix(nrow=33,ncol=1143)
date<-levels(factor(y_country$date))
rownames(data_ncov)<-country_selected
colnames(data_ncov)<-date
for (i in 1:229743){
  u<-y_country[i,1]
  v<-as.character(y_country[i,2])
  w<-y_country[i,3]
  if (u %in% country_selected){
  data_ncov[u,v]<-w
  }
}

save.image(file='datancov.RData')
load('datancov.RData')
data_select<-data_ncov[,1:100]
p<-33
n<-100
Stot<-sd(data_select)
for (i in 1:p){
  data_select[i,]<-data_select[i,]/Stot
}
data<-data_select
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

positive<-function(A,p){
  for (i in 1:p){
    for (j in 1:p){
      A[i,j]<-max(A[i,j],0)
    }
  }
  return (A)
}
norm1<-function(D){
  return(sum(abs(D)))
}
norm12<- function(D,p){
  s=0
  for (i in 1:p){
    s=s+sqrt(sum(D[i,])^2)
  }
  return(s)
}
norm_nuclear<-function(D,p){
  svd<-svd(D)
  v<-svd$d
  return (sum(v))
}
Loss_D<-function(D1,D2,S,R,lambda){
  L=0
  L=L+(norm((S-(D1+D2) %*% R),"f"))^2
  L=L+lambda[1]*norm1(D1)
  return(L)
}

Loss_D_new<-function(D1,D2,S,R,lambda,p){
  L=0
  L=L+(norm((S-(D1+D2) %*% R),"f"))^2/(2*ncol(R))
  L=L+lambda[1]*norm1(D1)
  for (i in 1:p){
    L=L+lambda[2]*sqrt(sum((D2[i,])^2))
  }
  return(L)
}

sgn<-function(D,p){
  for (i in 1:p){
    for (j in 1:p){
      if (D[i,j]>0){
        D[i,j]=1
      }
      else if (D[i,j]==0){
        D[i,j]=0
      }
      else if (D[i,j]<0){
        D[i,j]=-1
      }
    }
  }
  return (D)
}


D_estimation_new<-function(data, Int_X, D_01,D_02,X_0, thres, iter,lambda1,lambda2,p,n){
  D1_previous=D_01
  D2_previous=D_02
  D1_current=D_01
  D2_current=D_02
  X0_current=X_0
  X0_matrix<-matrix(nrow=p,ncol=n)
  for (l in 1:n){
    X0_matrix[,l]<-X0_current
  }
  loss_list=c()
  Loss_current<-Loss_D_new(D1_current,D2_current,data,Int_X,c(lambda1,lambda2),p)
  print(Loss_current)
  Loss_list<-c(Loss_current)
  for (k in 1:iter){
    Dhat<-matrix(nrow=p,ncol=p)
    X0<-numeric(p)
    for (i in 1:p){
      lm<-glmnet(t(Int_X),(data[i,]-D2_current[i,]%*%Int_X),alpha=1,lambda=lambda1,standardize = FALSE)
      Dhat[i,]<-as.numeric(lm$beta)
      X0[i]<-lm$a0
    }
    X0_current<-X0
    X0_matrix<-matrix(nrow=p,ncol=n)
    for (l in 1:n){
      X0_matrix[,l]<-X0_current
    }
    D1_current<-Dhat
    for (i in 1:p){
      model<-gglasso(t(Int_X),(data[i,]-D1_current[i,]%*%Int_X-rep(X0_current[i],n))
                     ,group=rep(1,p),lambda = lambda2,intercept = FALSE)
      Dhat[i,]<-as.numeric(model$beta)
    }
    D2_current<-Dhat
    Loss_previous<-Loss_current
    Loss_current<-Loss_D_new(D1_current,D2_current,data,Int_X,c(lambda1,lambda2),p)
    print(k)
    print(Loss_current)
    Loss_list<-c(Loss_list, Loss_current)
    if(abs(Loss_current-Loss_previous)<thres){
      break
    }
  }
  D_current<-D1_current+D2_current
  return(cbind(X0_current, D_current))
}
FS_confouuder<-function(Res, p, s){
  Cov_matrix<-matrix(nrow=p,ncol=p)
  for (i in 1:p){
    for (j in 1:p){
      Cov_matrix[i,j]=cov(Res[i,],Res[j,])
    }
  }
  for (i in 1:p){
    Cov_matrix[i,i]=0
  }
  max_stat<-c()
  for (i in 1:p){
    ord_i<-order(abs(Cov_matrix[i,]),decreasing = TRUE)[1:1]
    max_stat<-c(max_stat, sum(abs(Cov_matrix[i,])[ord_i]))
  }
  ord<-order(max_stat,decreasing = TRUE)
  return (ord[1:s])
}
Dhat_cv<-function(data, Int_X, D_01,D_02,X_0, thres, iter,p,n,Lambda_list){
  Loss_list<-c()
  groups <- cut(1:n, breaks = 5, labels = FALSE)
  for (lambda in Lambda_list){
    lambda1=lambda
    if (p<100){
      lambda2=lambda1*2/3
    }
    if (p>=100){
      lambda2=lambda1*3/2
    }
    Loss<-0
    for (iter in 1:5){
      index<-which(groups==iter)
      remain0<-which(groups!=iter)
      n_cv=length(remain0)
      result<-D_estimation_new(data[,remain0], Int_X[,remain0],D_01,D_02,X_0, thres, iter,lambda1,lambda2,p,n_cv)
      X0<-result[,1]
      Dhat<-result[,2:(p+1)]
      for (i in index){
        Loss=Loss+sum((data[,i]-Dhat%*%Int_X[,i]-X0)^2)
      }
    }
    Loss_list<-c(Loss_list,Loss)
  }
  lambda_star<-Lambda_list[order(Loss_list)[1]]
  print(lambda_star)
  lambda1=lambda_star
  if (p<100){
    lambda2=lambda1*2/3
  }
  if (p>=100){
    lambda2=lambda1*3/2
  }
  result<-D_estimation_new(data,Int_X,D_01,D_02,X_0,thres, iter,lambda1,lambda2,p,n)
  return(result)
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






Ahat_estimation<-function(Int_X, data, Lambda_list,n,p){
  Loss_list<-c()
  groups <- cut(1:n, breaks = 5, labels = FALSE)
  groups=sample(groups,size=n)
  for (lambda in Lambda_list){
    Loss<-0
    for (iter in 1:5){
      index<-which(groups==iter)
      remain0<-which(groups!=iter)
      X_0<-rep(0,p)
      Ahat_0<-matrix(nrow=p,ncol=p)
      for (i in 1:p){
        if (sd(data[i,remain0])==0){
          Ahat_0[i,]<-rep(0,p)
          X_0[i]<-0
        }
        if (sd(data[i,remain0])!=0){
          lm<-glmnet(t(Int_X[,remain0]),data[i,remain0],alpha=1,lambda=lambda,standardize = FALSE)
          Ahat_0[i,]<-as.numeric(lm$beta)
          X_0[i]<-as.numeric(lm$a0)
        }
      }
      
      for (i in index){
        Loss=Loss+sum((data[,i]-Ahat_0%*%Int_X[,i]-X_0)^2)
      }
    }
    Loss_list<-c(Loss_list,Loss)
  }
  lambda_star<-Lambda_list[order(Loss_list)[1]]
  print(lambda_star)
  Ahat_0<-matrix(nrow=p,ncol=p)
  for (i in 1:p){
    if (sd(data[i,remain0])==0){
      Ahat_0[i,]<-rep(0,p)
      X_0[i]<-0
    }
    if (sd(data[i,remain0])!=0){
      lm<-glmnet(t(Int_X),data[i,],alpha=1,lambda=lambda_star,standardize = FALSE)
      Ahat_0[i,]<-as.numeric(lm$beta)
    }
  }
  return(Ahat_0)
}
save.image('ncov_data1.RData')










