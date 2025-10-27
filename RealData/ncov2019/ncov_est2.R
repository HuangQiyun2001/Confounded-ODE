load('ncov_data2.RData')
library(nCov2019)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)
set.seed(6)
Lambda_list<-0.001*1.5^c(1:15)
Ahat_0<-Ahat_estimation(Int_X, data, Lambda_list,n,p)
Ahat<-Ahat_estimation(Int_X, proj_y, Lambda_list,n,p)
loss_base1=((norm((Ahat_0-Astar_0),'f')/p)^2)
loss_proposed=((norm((Ahat-Astar_0),'f')/p)^2)
loss_base2=((norm((Dhat-Astar_0),'f')/p)^2)
loss_base1
loss_base2
loss_proposed

loss_base1=((norm((Ahat_0-Astar),'f')/p)^2)
loss_proposed=((norm((Ahat-Astar),'f')/p)^2)
loss_base2=((norm((Dhat-Astar),'f')/p)^2)
loss_base1
loss_base2
loss_proposed

rownames(Ahat_0)<-rownames(data)
colnames(Ahat_0)<-rownames(data)
rownames(Ahat)<-rownames(data)
colnames(Ahat)<-rownames(data)


save.image(file='ncov_result2.RData')
