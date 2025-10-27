load('gene_data1.RData')
library(GEOquery)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)
library(impute)
set.seed(6)
Lambda_list<-0.001*1.5^c(1:15)
Ahat_0<-Ahat_estimation(Int_X, data, Lambda_list,n,p,q)
Ahat<-Ahat_estimation(Int_X, proj_y, Lambda_list,n,p,q)
norm(Ahat_0-Astar0,'f')^2/(p*q)
norm(Ahat-Astar0,'f')^2/(p*q)
norm(Dhat-Astar0,'f')^2/(p*q)
norm(Ahat_0-Astar,'f')^2/(p*q)
norm(Ahat-Astar,'f')^2/(p*q)
norm(Dhat-Astar,'f')^2/(p*q)

gene_select<-gene_select[remain]
rug_matrix<-matrix(0,nrow=p,ncol=p)
for (u in 1:p){
  for (i in 1:p){
    if (Ahat[u,i]!=0){
      rug_matrix[u,i]<-1
    }
  }
  cnt=p
  for (j in 1:p){
    for (i in 1:(j)){
      cnt=cnt+1
      if (Ahat[u,cnt]!=0){
        rug_matrix[u,i]<-1
        rug_matrix[u,j]<-1
      }
    }
  }
}
rownames(rug_matrix)<-gene_select
colnames(rug_matrix)<-gene_select
library(igraph)
vi<-rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if (rug_matrix[i,j]!=0){
      if (i==j){
        rug_matrix[i,i]=0
      }
      else {
        vi[i]=vi[i]+1
        vi[j]=vi[j]+1
      }
    }
  }
}
gene_selected<-(vi>=1)
adj_matrix<-t(rug_matrix)[gene_selected,gene_selected]

net <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "directed",   
  weighted = TRUE     
)


summary(net)
plot(net,
     edge.arrow.size = 0.5,  
     vertex.label.cex = 1.2,  
     vertex.size = 15,        
     vertex.color = "lightblue",
     edge.width = abs(E(net)$weight)  
)

rug_matrix0<-matrix(0,nrow=p,ncol=p)
for (u in 1:p){
  for (i in 1:p){
    if (Ahat_0[u,i]!=0){
      rug_matrix0[u,i]<-1
    }
  }
  cnt=p
  for (j in 1:p){
    for (i in 1:(j)){
      cnt=cnt+1
      if (Ahat_0[u,cnt]!=0){
        rug_matrix0[u,i]<-1
        rug_matrix0[u,j]<-1
      }
    }
  }
}
rownames(rug_matrix0)<-gene_select
colnames(rug_matrix0)<-gene_select
library(igraph)
vi<-rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if (rug_matrix0[i,j]!=0){
      if (i==j){
        rug_matrix[i,i]=0
      }
      vi[i]=vi[i]+1
      vi[j]=vi[j]+1
    }
  }
}
gene_selected<-(vi>=1)
adj_matrix<-t(rug_matrix0)[gene_selected,gene_selected]

net <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "directed",   
  weighted = TRUE     
)


summary(net)
plot(net,
     edge.arrow.size = 0.5,  
     vertex.label.cex = 1.2,  
     vertex.size = 15,        
     vertex.color = "lightblue",
     edge.width = abs(E(net)$weight)  
)

save.image('Generesult2.RData')

