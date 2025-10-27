load('ncov_data1.RData')
library(nCov2019)
library(Matrix)
library(cubature)
library(MASS)
library(glmnet)
library(gglasso)
library(sparsepca)

set.seed(2)
Lambda_list<-0.001*1.5^c(1:15)
Ahat<-Ahat_estimation(Int_X, proj_y, Lambda_list,n,p)
Ahat_0<-Ahat_estimation(Int_X, data, Lambda_list,n,p)


rownames(Ahat_0)<-rownames(data)
colnames(Ahat_0)<-rownames(data)
rownames(Ahat)<-rownames(data)
colnames(Ahat)<-rownames(data)

library(igraph)
rug_matrix<-Ahat
vi<-rep(0,p)
for (i in 1:p){
  for (j in 1:p){
    if (rug_matrix[i,j]!=0){
      if (i==j){
        rug_matrix[i,j]=0
      }
      else{
        vi[i]=vi[i]+1
        vi[j]=vi[j]+1
      }
    }
  }
}
gene_selected<-(vi>=1)
adj_matrix<-t(rug_matrix)[gene_selected,gene_selected]
p_select<-sum(gene_selected==TRUE)

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
edge_colors <- ifelse(E(net)$weight > 0, "purple", "purple")

plot(net,
     edge.arrow.size = 0.5,
     vertex.label.color = "black",
     vertex.frame.color = "gray",
     vertex.color = "lightblue",
     edge.color = edge_colors,
     vertex.size =20,
     edge.width = (log(abs(E(net)$weight)+1)*2+0.5),
     # 布局方式：圆形
)
save.image(file='ncov_result1.RData')

