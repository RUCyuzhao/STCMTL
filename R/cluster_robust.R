#' @title SCMTL algorithm(robust version)
#'
#' @param X \code{matrix} independent variable.
#' @param Y \code{matrix} dependent variable.
#' @param group \code{matrix} specify the structure of the task.
#' @param K cluster number.
#' @param ext.prop the proportion of extreme neighbors for each cluster. Without any prior information, set ext.prop=5/min(group[-1]-group[-length(group)]).
#' @param pure.prop  (optional) the expected proportion of pure tasks. By default pure.prop=0.5.
#' @param nlambda (optional) a user supplied lambda sequence. By default lambda=2^(-15,3,length=300).
#' @param alpha (optional) a 0-1 number to control the step size. By default alpha=1.
#' @param max_iter max iterations.
#' @param thresh convergence threshold for the algorithm.
#' @param show \code{TRUE/FALSE} whether show every step's result.
#' @param method2 the type of final procedure's algorithm. The value must be one of {lasso,scad}
#' @param nfold number of folds - default is 5.
#' @param cut threshold of outlier.
#'
#' @return The estimated cluster center matrix U, membership matrix V, coefficient matrix W and the outlier tasks set.
#' @export
#'
#' @examples
cluster_robust <- function(X,Y,group,K,ext.prop=NULL,nlambda=2^(seq(-15,3,length=100)),pure.prop = 0.5,alpha = 1,max_iter = 50,thresh=1e-2,show=F,method2 = "lasso",nfold = 5,cut = 3){
  T <- length(group)-1
  p <- NCOL(X)
  N <- NROW(X)
  Z <- matrix(0,nrow = N,ncol = p*T)
  for(i in 1:T){
    Z[group[i]:(group[i+1]-1),((i-1)*p+1):(i*p)] <- X[group[i]:(group[i+1]-1),]
  }
  y_multi <-matrix(rep(Y,T),nrow = N,byrow = F)
  I <- matrix(0,nrow = T,ncol = N)
  for(i in 1:T){
    I[i,group[i]:(group[i+1]-1)] <- 1/(group[i+1]-group[i])
  }

  # initialize
  beta0 <- matrix(nrow = p,ncol = T)
  beta_real <- matrix(nrow = p,ncol = T)
  lambda <- c()
  for(i in 1:T){
    lambda[i] <- cv.glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],nfolds = nfold,type.measure = "mse",standardize=FALSE,lambda = nlambda)$lambda.min/10
    fit <- glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],lambda = lambda[i],standardize = FALSE)
    beta0[,i]<- as.array(coef(fit))[1:p]
  }
  beta_real <- beta0

  #得到V，U的初值
  delta <- I%*%abs(y_multi-X%*%beta0)
  soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
  V <- t(soup.out$memberships[[1]])
  U_chu <- matrix(0,nrow = T,ncol = K)
  for(i in 1:T){
    U_chu[i,] <- lsei(A=t(V),B=t(delta)[,i],type=2)$X
  }
  re <- NMF_robust(delta,1000,1e-5,U_chu = U_chu,V_chu = V,K = 5,alpha = 1)
  ind_out <- which.max(re$L)
  delta_cluster <- delta[,-ind_out]
  if(sum(!is.na(ind_out))>0){
    soup.out <- SOUP(t(delta_cluster),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
    V <- t(soup.out$memberships[[1]])
    ind <- apply(V,2,function(x){return(sum(x^2)==0)})
    if(sum(ind,na.rm = T)>0 | sum(is.na(ind))>0){
      temp <- sample(1:K,1)
      V[temp,which(ind)] <- 1
    }
  }

  U <- matrix(0,nrow = p,ncol = K)
  for(i in 1:p){
    U[i,] <- lsei(A=t(V),B=t(beta0[,-ind_out])[,i],type = 2)$X
  }

  V_cal <- matrix(0,nrow = K ,ncol = T)
  V_cal[,-ind_out] <- V
  jilu <- c()
  choice <- permutations(K)
  rmse <- c(sum((Y-Z%*%as.vector(beta0))^2))
  rmse_cha <- c(rmse)
  V_cha <- c(sum(svd(V)$d^2))
  iter <- 0
  panding <- 1
  beta_temp <- matrix(nrow = p,ncol = T)
  l1 <- 0
  while((l1 != length(ind_out)|(rmse_cha[iter+1] >=thresh & alpha > 0 )) & iter<= max_iter){
    t1 <- Sys.time()
    l1 <- length(ind_out)
    V_cal1 <- V_cal
    V1 <- V
    U1 <- U
    beta_temp[,-ind_out] <- U1%*%V1
    beta_temp[,ind_out] <- beta0[,ind_out]
    for(i in setdiff(1:T,ind_out)){
      X_temp <- X[group[i]:(group[i+1]-1),]
      Y_temp <- Y[group[i]:(group[i+1]-1)]
      temp <- min((group[i+1]-group[i]),p/2)
      if(iter == 0)
        temp <- 5
      temp <- coordinate_decent(X_temp,Y_temp,max_iter = temp,lambda=lambda[i],beta_temp = beta_temp[,i] ,thresh = 1e-2)*alpha+(1-alpha)*beta_temp[,i]
      beta_temp[,i] <- temp[1:p]
    }
    delta <- I%*%abs(y_multi-X%*%beta_temp)
    soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
    V_chu <- t(soup.out$memberships[[1]])
    U_chu <- matrix(0,nrow = T,ncol = K)
    for(i in 1:T){
      U_chu[i,] <- lsei(A=t(V_chu),B=t(delta)[,i],type=2)$X
    }
    delta_mu <- apply(delta,1,mean)/K
    U_chu <- matrix(rep(delta_mu,K),ncol = K,byrow = FALSE)
    re <- NMF_robust(delta,1000,1e-5,U_chu = U_chu,V_chu = V_chu+1e-6,K = 5)
    # ind_out <- find_change(re$L,0.9,3)
    re$L[ind_out] <- 0
    if(max(re$L)>mean(re$L[-c(ind_out,which.max(re$L))])+cut*sd(re$L[-c(ind_out,which.max(re$L))])){
      ind_out <- c(ind_out,which.max(re$L))
    }
    sample_out <- c()
    p_out <- c()
    for(i in ind_out){
      sample_out <- c(sample_out,group[i]:(group[i+1]-1))
      p_out <- c(p_out,((i-1)*p+1):(i*p))
    }
    delta_cluster <- delta[,-ind_out]
    if(sum(!is.na(ind_out))>0){
      soup.out <- SOUP(t(delta_cluster),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
      V <- t(soup.out$memberships[[1]])
      ind <- apply(V,2,function(x){return(sum(x^2)==0)})
      if(sum(ind,na.rm = T)>0 | sum(is.na(ind))>0){
        temp <- sample(1:K,1)
        V[temp,which(ind)] <- 1
      }
    }
    V_cal <- matrix(0,nrow = K ,ncol = T)
    V_cal[,-ind_out] <- V
    V_cha[iter+2] <- sum(svd(V_cal[,-ind_out]-V_cal1[,-ind_out])$d^2)
    ind <- choice[1,]
    for(j in 1:NROW(choice)){
      if(V_cha[iter+2] > sum(svd(V_cal1[choice[j,],-ind_out]-V_cal[,-ind_out])$d^2) ){
        V_cha[iter+2] <- sum(svd(V_cal1[choice[j,],-ind_out]-V_cal[,-ind_out])$d^2)
        ind <- choice[j,]
      }
    }
    V_cha[iter+2] <- sum((V_cal1[ind,]-V_cal)^2)
    U <- U[,ind]
    t2 <- Sys.time()
    Y_bottom <- Z[,-p_out]%*%as.vector(U%*%V)
    for(i in 1:K){
      X_temp <- matrix(nrow = N,ncol = p)
      for(j in setdiff(1:T,ind_out)){
        X_temp[group[j]:(group[j+1]-1),] <- V_cal[i,j] * X[group[j]:(group[j+1]-1),]
      }
      X_temp <- X_temp[-sample_out,]
      Y_temp <- Y[-sample_out]-Y_bottom[-sample_out]+X_temp%*%U[,i]
      # fit_U <- glmnet(X_temp[,-1],Y_temp,lambda = mean(lambda)/rate2,family = "gaussian")
      # U[,i] <- as.array(coef(fit_U))[1:(p)]*alpha+(1-alpha)*U[,i]
      U[,i] <- lassoshooting::lassoshooting(X_temp,Y_temp,lambda = mean(lambda))$coefficients*alpha+(1-alpha)*U[,i]
    }
    beta_out <- matrix(0,nrow = p,ncol = T)
    beta_out[,ind_out] <- beta_real[,ind_out]
    rmse[iter+2] <- sum((Y-Z%*%as.vector(U%*%V_cal+beta_out))^2)
    rmse_cha[iter+2] <- abs(rmse[iter+2]-rmse[iter+1])/rmse[iter+1]
    iter <- iter+1
    if(show == TRUE){
      print(paste0("第",iter,"次循环 "," 收敛：",rmse_cha[iter+1]," 估计准确度：",rmse[iter+1]))
    }
  }
  # Y_bottom <- Z[,-p_out]%*%as.vector(U%*%V)
  for(i in 1:K){
    X_temp <- matrix(nrow = N,ncol = p)
    for(j in setdiff(1:T,ind_out)){
      X_temp[group[j]:(group[j+1]-1),] <- V_cal[i,j] * X[group[j]:(group[j+1]-1),]
    }
    X_temp <- X_temp[-sample_out,]
    Y_temp <- Y[-sample_out]-Y_bottom[-sample_out]+X_temp%*%U[,i]
    if(i==1){
      for(j in ind_out){
        lambda_out <- cv.glmnet(X[group[j]:(group[j+1]-1),-1],Y[group[j]:(group[j+1]-1)],nfolds = nfold,family="gaussian",standardize = FALSE)$lambda.min
        fit <- glmnet(X[group[j]:(group[j+1]-1),-1],Y[group[j]:(group[j+1]-1)],lambda = lambda_out,standardize = FALSE)
        beta_temp[,j] <- as.array(coef(fit))[1:(p)]
      }
    }
    if(method2=="lasso"){
      lambda_U <- cv.glmnet(X_temp[,-1],Y_temp,nfolds = nfold,family="gaussian",standardize = FALSE,lambda = nlambda)$lambda.min
      fit <- glmnet(X_temp[,-1],Y_temp,lambda = lambda_U,standardize = FALSE)
      U[,i]<- as.array(coef(fit))[1:(p)]
    }
    if(method2=="scad"){
      cvfit <- cv.ncvreg(X_temp[,-1],Y_temp,method="SCAD",nfolds = nfold,lambda = rev(nlambda),nlambda=length(nlambda))
      fit <- cvfit$fit
      U[,i] <- fit$beta[,cvfit$min]
    }
  }
  beta_out <- matrix(0,nrow = p,ncol = T)
  beta_out[,ind_out] <- beta_temp[,ind_out]
  rmse[iter+2] <- sum((Y-Z%*%as.vector(U%*%V_cal+beta_out))^2)
  l <- list(U,V_cal,sqrt(rmse),beta0,U%*%V_cal+beta_out,ind_out)
  names(l) <- c("U","V","RSS","beta0","beta","outlier_task")
  return(l)
}
