library(parallel)


GetResult <- function(seed_r,n_star,p_star,noisesd){
  set.seed(seed_r)
  library(Matrix)
  library(cubature)
  library(MASS)
  library(glmnet)
  library(gglasso)
  library(deSolve)
  T=1
  n=n_star
  times <- seq(0, T, by = T/(n))[1:n]
  p=p_star
  q=2*p
  d=3
  confounder<-c(1,3,5)
  X0<-runif(p,4,6)
  X0[confounder]<-runif(d,10,15)
  l<-runif(1,-0.5,0.5)
  data<-matrix(nrow=p,ncol=n)
  noise<-matrix(nrow=p,ncol=n)
  A<-matrix(nrow=p,ncol=p)
  for (i in 1:(p/2)){
    A[(2*i),(2*i-1)]<-(-2*pi*i)
    A[(2*i-1),(2*i)]<-(2*pi*i)
    A[(2*i),(2*i)]<-0
    A[(2*i-1),(2*i-1)]<-0
  }
  ode_sim <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      dPrey <- alpha * Predator + beta * (cos(Prey/5)-1)
      dPredator <- delta * Prey + gamma * (cos(Predator/5)-1)
      return(list(c(dPrey, dPredator)))
    })
  }
  alpha1<-2*pi*c(1:(p/2))
  alpha2<-c(1:(p/2))+runif((p/2),-1,1)
  alpha3<--2*pi*c(1:(p/2))
  alpha4<-c(1:(p/2))+runif((p/2),-1,1)
  A_star<-matrix(0,nrow=p,ncol=q)
  for (i in 1:(p/2)){
    A_star[(2*i-1),(2*i)]<-alpha1[i]
    A_star[(2*i),(2*i-1)]<-alpha3[i]
  }
  for (i in 1:(p/2)){
    A_star[(2*i-1),(p+(2*i-1))]<-alpha2[i]
    A_star[(2*i),(p+(2*i))]<-alpha4[i]
  }
  sigma=noisesd
  for (i in 1:(p/2)){
    parameters <- c(alpha=alpha1[i],beta=alpha2[i],
                    delta=alpha3[i],gamma=alpha4[i])
    
    initial_state <- c(Prey=X0[(2*i-1)],Predator=X0[(2*i)])
    solution <- ode(
      y = initial_state,
      times = times,
      func = ode_sim,
      parms = parameters
    )
    data[(2*i-1),]<-solution[1:(n),2]
    data[(2*i),]<-solution[1:(n),3]
    noise[(2*i-1),]<-rnorm(n,0,sigma)
    noise[(2*i),]<-rnorm(n,0,sigma)
  }
  data_clean<-data
  data<-data_clean+noise

  loess_cv <- function(data, spans, k = 5) {
    set.seed(123)  
    n <- nrow(data)
    folds <- sample(rep(1:k, length.out = n))  
    cv_errors <- matrix(NA, nrow = k, ncol = length(spans))
    colnames(cv_errors) <- spans
    
    for (i in 1:k) {
      train_data <- data[folds != i, ]
      test_data <- data[folds == i, ]
      
      for (j in seq_along(spans)) {
        fit <- loess(y ~ x, data = train_data, span = spans[j], degree = 1)  
        pred <- predict(fit, newdata = test_data)
        cv_errors[i, j] <- mean((test_data$y - pred)^2, na.rm = TRUE)  
      }
    }
    
    avg_mse <- colMeans(cv_errors)
    best_span <- spans[which.min(avg_mse)]
    
    return(list(best_span = best_span, cv_errors = avg_mse))
  }
  
  data_smooth<-matrix(nrow=p,ncol=n)
  span<-n^(-0.25)*1.2^c(-10:0)
  for (i in 1:p){
    smooth_result<-loess_cv(data.frame(x=times,y=data[i,]),spans = span)
    bandwith<-smooth_result$best_span/(n^(1/5))
    smooth_result<-loess(y~x,data.frame(x=times,y=data[i,]),span=bandwith)
    data_smooth[i,]<-smooth_result$fitted
  }
  
  remain<-setdiff(c(1:p),confounder)
  remain2<-c(remain,remain+rep(p,p-d))
  p=p-d
  data=data[remain,]
  data_clean=data_clean[remain,]
  A_star0<-A_star
  A_star=A_star[remain,remain2]
  q<-2*p
  Covariate<-matrix(nrow=q,ncol=n)
  rownames(Covariate)<-c(1:q)
  for (i in 1:p){
    Covariate[i,]<-data[i,]
    Covariate[(i+p),]<-data[i,]^2
  }
  sigma2_hat<-c()
  for (i in 1:p){
    sumi<-sum((data[i,(2:n)]-data[i,(1:(n-1))])^2)
    sigma2_hat<-c(sigma2_hat, (sumi/(n-1)))
  }
  Int_X2<-matrix(nrow=q, ncol=n)
  for (i in 1:n){
    for (j in 1:p){
      data_0<-data[j,1:(i)]
      data_0[1]<-(data_0[1]+data_0[i])/2
      Int_X2[j,i]<-sum(data_0[1:(i-1)])*T/n
    }
    for (j in 1:p){
      data_0<-cos(data_smooth[j,1:(i)]/5)-rep(1,i)
      data_0[1]<-(data_0[1]+data_0[i])/2
      Int_X2[(j+p),i]<-sum(data_0[1:(i-1)])*T/n
    }
  }
  Int_X<-Int_X2
  IntX_centered<-matrix(nrow=q, ncol=n)
  for (i in 1:n){
    for (j in 1:q){
      IntX_centered[j,i]<-Int_X[j,i]-mean(Int_X[j,])
    }
  }
  data_centered<-matrix(nrow=p, ncol=n)
  for (i in 1:n){
    for (j in 1:p){
      data_centered[j,i]<-data[j,i]-mean(data[j,])
    }
  }
  
  
  
  Ahat_estimation<-function(Int_X, data, Lambda_list,n,p,q){
    Loss_list<-c()
    groups <- cut(1:n, breaks = 5, labels = FALSE)
    groups=sample(groups,size=n)
    for (lambda in Lambda_list){
      Loss<-0
      for (iter in 1:5){
        index<-which(groups==iter)
        remain0<-which(groups!=iter)
        X_0<-rep(0,p)
        Ahat_0<-matrix(0,nrow=p,ncol=q)
        for (i in 1:p){
          if (sd(data[i,remain0])!=0){
            lm<-glmnet(t(Int_X[,remain0]),data[i,remain0],alpha=1,lambda=lambda)
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
    Ahat_0<-matrix(0,nrow=p,ncol=q)
    for (i in 1:p){
      if (sd(data[i,remain0])!=0){
        lm<-glmnet(t(Int_X),data[i,],alpha=1,lambda=lambda_star)
        Ahat_0[i,]<-as.numeric(lm$beta)
      }
    }
    return(Ahat_0)
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
  
  D_estimation_new<-function(data, Int_X, D_01,D_02,X_0, thres, iter,lambda1,lambda2,p,q,n){
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
      Dhat<-matrix(nrow=p,ncol=q)
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
                       ,group=rep(1,q),lambda = lambda2,intercept = FALSE)
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
  Dhat_cv<-function(data, Int_X, D_01,D_02,X_0, thres, iter,p,q,n,Lambda_list){
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
        result<-D_estimation_new(data[,remain0], Int_X[,remain0],D_01,D_02,X_0, thres, iter,lambda1,lambda2,p,q,n_cv)
        X0<-result[,1]
        Dhat<-result[,2:(q+1)]
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
      lambda2=lambda*2/3
    }
    if (p>=100){
      lambda2=lambda1*3/2
    }
    result<-D_estimation_new(data,Int_X,D_01,D_02,X_0,thres, iter,lambda1,lambda2,p,q,n)
    return(result)
  }
  
  X_0<-rep(0,p)
  Lambda_list<-0.001*1.5^c(1:15)
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
  lambda_ratio<-lambda_ratio<-(pca$sdev[1:ceiling(sqrt(p))]-pca$sdev[2:(ceiling(sqrt(p))+1)])/(pca$sdev[2:(ceiling(sqrt(p))+1)]-pca$sdev[3:(ceiling(sqrt(p))+2)])
  dhat<-order(lambda_ratio,decreasing = TRUE)[1]
  dhat
  u=pca$rotation[,1:d]
  u2=pca$rotation[,1:dhat]
  proj<-diag(p)-u %*% ginv(t(u) %*% u) %*% t(u)
  proj_y= proj %*% data
  proj2<-diag(p)-u2 %*% ginv(t(u2) %*% u2) %*% t(u2)
  proj_y2= proj2 %*% data
  
  s_fs<-min(20,p)
  fs_result<-FS_confouuder(Res_scaled,p,s_fs)
  for (i in setdiff(c(1:p),fs_result)){
    Res_scaled[i,]<-rep(0,n)
    Res[i,]<-rep(0,n)
  }
  pca=prcomp(t(Res_scaled))
  lambda_ratio<-lambda_ratio<-(pca$sdev[1:ceiling(sqrt(p))]-pca$sdev[2:(ceiling(sqrt(p))+1)])/(pca$sdev[2:(ceiling(sqrt(p))+1)]-pca$sdev[3:(ceiling(sqrt(p))+2)])
  if (p>=100){
    dhat<-order(lambda_ratio,decreasing = TRUE)[1]
  }
  u_fs=pca$rotation[,1:d]
  proj_fs<-diag(p)-u_fs %*% ginv(t(u_fs) %*% u_fs) %*% t(u_fs)
  proj_y_fs= proj_fs %*% data
  u_fs=pca$rotation[,1:dhat]
  proj_fs2<-diag(p)-u_fs %*% ginv(t(u_fs) %*% u_fs) %*% t(u_fs)
  proj_y_fs2= proj_fs2 %*% data
  Ahat_estimation<-function(Int_X, data, Lambda_list,n,p,q){
    Loss_list<-c()
    groups <- cut(1:n, breaks = 5, labels = FALSE)
    groups=sample(groups,size=n)
    for (lambda in Lambda_list){
      Loss<-0
      for (iter in 1:5){
        index<-which(groups==iter)
        remain0<-which(groups!=iter)
        X_0<-rep(0,p)
        Ahat_0<-matrix(nrow=p,ncol=q)
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
    Ahat_0<-matrix(nrow=p,ncol=q)
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
  Lambda_list<-0.001*1.5^c(1:15)
  Ahat_0<-Ahat_estimation(Int_X, data, Lambda_list,n,p,q)
  print((norm((Ahat_0-A_star),'f')/norm(A_star,'f'))^2)
  print((norm((Dhat-A_star),'f')/norm(A_star,'f'))^2)
  if (p<100){
    Ahat<-Ahat_estimation(Int_X, proj_y, Lambda_list,n,p,q)
    Ahat2<-Ahat_estimation(Int_X, proj_y2, Lambda_list,n,p,q)
    print((norm((Ahat-A_star),'f')/norm(A_star,'f'))^2)
    print((norm((Ahat2-A_star),'f')/norm(A_star,'f'))^2)
  }
  if (p>=100){
    Ahat<-Ahat_estimation(Int_X, proj_y_fs, Lambda_list,n,p,q)
    Ahat2<-Ahat_estimation(Int_X, proj_y_fs2, Lambda_list,n,p,q)
    print((norm((Ahat-A_star),'f')/norm(A_star,'f'))^2)
    print((norm((Ahat2-A_star),'f')/norm(A_star,'f'))^2)
  }
  
  
  fs_acc<-0
  for (i in 1:min(6,p)){
    if (i %in% fs_result){
      fs_acc=fs_acc+1/min(6,p)
    }
  }
  u_star<-A_star0[remain,setdiff(c(1:(p+d)),remain)]
  
  
  proj_true<-diag(p)-u_star %*% ginv(t(u_star) %*% u_star) %*% t(u_star)
  
  error_B<-norm((proj-proj_true),'2')/norm(proj,'2')
  error_B_fs<-norm((proj_fs-proj_true),'2')/norm(proj,'2')
  proj_true_y=proj_true%*% data
  
  Ahat_o2<-Ahat_estimation(Int_X, proj_true_y, Lambda_list,n,p,q)
  print((norm((Ahat_o2-A_star),'f')/norm(A_star,'f'))^2)
  proj_A<-proj_true%*%A_star
  
  rs1<-((norm((Ahat_0-A_star),'f')/norm(A_star,'f'))^2)
  rs2<-((norm((Dhat-A_star),'f')/norm(A_star,'f'))^2)
  rs3<-((norm((Ahat-A_star),'f')/norm(A_star,'f'))^2)
  rs4<-((norm((Ahat2-A_star),'f')/norm(A_star,'f'))^2)
  rs5<-((norm((proj_A-A_star),'f')/norm(A_star,'f'))^2)
  rs6<-((norm((Ahat_o2-A_star),'f')/norm(A_star,'f'))^2)
  rs7<-dhat
  rs8<-fs_acc
  
  res <- list(error_0=rs1,error_D=rs2,error_our=rs3,error_our2=rs4,error_oracle=rs6,
              error_PA=rs5,dhat=rs7,fs_acc=rs8)
}

nrun <- 100  ##重复次数

for (noisesd in c(0.6)) {
  for (n_star in c(200,400)) {
    for (p_star in c(40,60,80,120,200)) {
      cat("n",":",n_star,"\n","p",':',p_star,"noise",noisesd,"\n")
      
      ncore <- 105  ##调用的kernel的个数
      
      # stopCluster(cl)
      res <- mclapply(
        X = 1:nrun,
        FUN = function(seed_r) GetResult(seed_r,n_star,p_star,noisesd),
        mc.cores = ncore,
        mc.preschedule = TRUE,  # 计算密集型通常更快；不均衡任务可改 FALSE
        mc.set.seed = TRUE      # 每个子进程独立 RNG，保证可重复
      )
      
      noisesdp <- round(noisesd, digits = 4)
      save(res,file = paste0("0910simuresults2","n",n_star,"p",p_star,"noise",noisesd))
    }
  }
}

results_final<-matrix(nrow=5,ncol=8)
colnames(results_final)<-c("error_0","error_D","error_our","error_our2","error_oracle","error_PA","dhat","fs_acc"
)
rownames(results_final)<-c(40,60,80,120,200)
ii<-1
n_star=400
noisesd=0.6
for (p_star in c(40,60,80,120,200)){
  file_name=paste0("0910simuresults2","n",n_star,"p",p_star,"noise",noisesd)
  load(file_name)
  results<-res
  error_0<-c()
  error_D<-c()
  error_our<-c()
  error_our2<-c()
  error_oracle<-c()
  error_PA<-c()
  dhat<-c()
  fs_acc<-c()
  for (i in 1:100){
    error_0<-c(error_0, (results[[i]])[1])
    error_D<-c(error_D, (results[[i]])[2])
    error_our<-c(error_our, (results[[i]])[3])
    error_our2<-c(error_our2, (results[[i]])[4])
    error_oracle<-c(error_oracle, (results[[i]])[5])
    error_PA<-c(error_PA, (results[[i]])[6])
    dhat<-c(dhat, (results[[i]])[7])
    fs_acc<-c(fs_acc, (results[[i]])[8])
  }
  aa<-round(mean(as.numeric(error_0)),3)
  bb<-round(sd(as.numeric(error_0))/10,3)
  uu1<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(error_D)),3)
  bb<-round(sd(as.numeric(error_D))/10,3)
  uu2<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(error_our)),3)
  bb<-round(sd(as.numeric(error_our))/10,3)
  uu3<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(error_our2)),3)
  bb<-round(sd(as.numeric(error_our2))/10,3)
  uu4<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(error_oracle)),3)
  bb<-round(sd(as.numeric(error_oracle))/10,3)
  uu5<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(error_PA)),3)
  bb<-round(sd(as.numeric(error_PA))/10,3)
  uu6<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(dhat)),3)
  bb<-round(sd(as.numeric(dhat))/10,3)
  uu7<-paste(aa,'(',bb,')')
  aa<-round(mean(as.numeric(fs_acc)),3)
  bb<-round(sd(as.numeric(fs_acc))/10,3)
  uu8<-paste(aa,'(',bb,')')
  results_final[ii,1]<-uu1
  results_final[ii,2]<-uu2
  results_final[ii,3]<-uu3
  results_final[ii,4]<-uu4
  results_final[ii,5]<-uu5
  results_final[ii,6]<-uu6
  results_final[ii,7]<-uu7
  results_final[ii,8]<-uu8
  ii=ii+1
}
show(results_final)
