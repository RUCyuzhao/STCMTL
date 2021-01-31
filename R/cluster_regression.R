#' @title SCMTL Algorithm(linear regression)
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
#'
#' @return The estimated cluster center matrix U, membership matrix V and coefficient matrix W.
#' @import SOUP
#' @import glmnet
#' @import MASS
#' @import limSolve
#' @import ncvreg
#' @import e1071
#' @import lassoshooting
#' @export
#'
#' @examples  x_train <- as.matrix(school$train$x)
#' y_train <- school$train$y[,1]
#' group <- as.numeric(school$group[1,])
#' re <- cluster_regression(x_train,y_train,group=group,K=2,ext.prop = 5/(length(group)-1)
#' ,max_iter = 10,thresh = 1e-3,method2 = "scad",show = TRUE)

cluster_regression <- function(X,Y,group,K,ext.prop=NULL,nlambda=2^(seq(-15,3,length=100)),pure.prop = 0.5,alpha = 1,max_iter = 50,thresh=1e-2,show=F,method2 = "lasso",nfold = 5){
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
    I[i,group[i]:(group[i+1]-1)] <- 1
  }
  t1 <- Sys.time()
  beta0 <- matrix(nrow = p,ncol = T)
  lambda <- c()
  for(i in 1:T){
    lambda[i] <- cv.glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],nfolds = nfold,type.measure = "mse",standardize=FALSE,lambda = nlambda)$lambda.min/10
    fit <- glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],lambda = lambda[i],standardize = FALSE)
    beta0[,i]<- as.array(coef(fit))[1:p]
  }
  t_cv <- Sys.time()-t1
  delta <- I%*%abs(y_multi-X%*%beta0)
  soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
  V <- t(soup.out$memberships[[1]])
  ind <- apply(V,2,function(x){return(sum(x^2)==0)})
  if(sum(ind)>0){
    temp <- sample(1:K,1)
    V[temp,which(ind)] <- 1
  }
  U <- matrix(0,nrow = p,ncol = K)
  for(i in 1:p){
    U[i,] <- lsei(A=t(V),B=t(beta0)[,i],type = 2)$X
  }
  jilu <- c()
  choice <- permutations(K)
  rmse <- c(sum((Y-Z%*%as.vector(beta0))^2))
  V_cha <- c(sum(svd(V)$d^2))
  U_cha <- c(sum(svd(U)$d^2))
  W_cha <- c(sum(svd(U%*%V)$d^2))
  rmse_cha <- c(rmse)
  iter <- 0
  panding <- 1
  step1_time <- c()
  step2_time <- c()
  step3_time <- c()
  while(rmse_cha[iter+1] >=thresh & alpha > 0 & iter<= max_iter){
    t0 <- Sys.time()
    V1 <- V
    U1 <- U
    beta_temp <- U1%*%V1
    for(i in 1:T){
      X_temp <- X[group[i]:(group[i+1]-1),]
      Y_temp <- Y[group[i]:(group[i+1]-1)]
      temp <- min((group[i+1]-group[i]),p/2)
      if(iter == 0)
        temp <- 5
      temp <- coordinate_decent(X_temp,Y_temp,max_iter = temp,lambda=lambda[i],beta_temp = beta_temp[,i] ,thresh = 1e-2)*alpha+(1-alpha)*beta_temp[,i]
      beta_temp[,i] <- temp[1:p]
    }
    step1_time[iter+1] <- Sys.time()-t0
    t1 <- Sys.time()
    delta <- I%*%abs(y_multi-X%*%beta_temp)
    soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
    V <- t(soup.out$memberships[[1]])
    ind <- apply(V,2,function(x){return(sum(x^2)==0)})
    if(sum(ind)>0){
      temp <- sample(1:K,1)
      V[temp,which(ind)] <- 1
    }
    step2_time[iter+1] <- Sys.time()-t1
    t1 <- Sys.time()
    V_cha[iter+2] <- sum(svd(V1-V)$d^2)
    ind <- choice[1,]
    for(j in 1:NROW(choice)){
      if(V_cha[iter+2] > sum(svd(V1[choice[j,],]-V)$d^2) ){
        V_cha[iter+2] <- sum(svd(V1[choice[j,],]-V)$d^2)
        ind <- choice[j,]
      }
    }
    V_cha[iter+2] <- sqrt(V_cha[iter+2]/sum(svd(V1)$d^2))
    U <- U[,ind]
    U1 <- U
    Y_bottom <- Z%*%as.vector(U%*%V)
    for(i in 1:K){
      X_temp <- matrix(nrow = N,ncol = p)
      for(j in 1:T){
        X_temp[group[j]:(group[j+1]-1),] <- V[i,j] * X[group[j]:(group[j+1]-1),]
      }
      Y_temp <- Y-Y_bottom+X_temp%*%U[,i]

      U[,i] <- lassoshooting::lassoshooting(X_temp,Y_temp,lambda = mean(lambda))$coefficients*alpha+(1-alpha)*U[,i]

    }
    U_cha[iter+2] <- sqrt(sum(svd(U1-U)$d^2)/sum(svd(U1)$d^2))
    W_cha[iter+2] <- sqrt(sum(svd(beta_temp-U%*%V)$d^2)/sum(svd(beta_temp)$d^2))
    rmse[iter+2] <- sum((Y-Z%*%as.vector(U%*%V))^2)
    rmse_cha[iter+2] <- abs(rmse[iter+2]-rmse[iter+1])/rmse[iter+1]
    step3_time[iter+1] <- Sys.time()-t1
    iter <- iter+1

    if(show == TRUE){
      print(paste0("第",iter,"次循环 "," 目标函数收敛：",rmse_cha[iter+1]," 纠错步运行时间：",step1_time[iter]," Soup运行时间:",step2_time[iter]," U更新时间：",step3_time[iter]," 总时间：",Sys.time()-t0))
    }
  }
  t1 <- Sys.time()
  Y_bottom <- Z%*%as.vector(U%*%V)
  for(i in 1:K){
    X_temp <- matrix(nrow = N,ncol = p)
    for(j in 1:T){
      X_temp[group[j]:(group[j+1]-1),] <- V[i,j] * X[group[j]:(group[j+1]-1),]
    }
    Y_temp <- Y-Y_bottom+X_temp%*%U[,i]

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
  t_cv<- t_cv+Sys.time()-t1
  rmse <- c(rmse,sum((Y-Z%*%as.vector(U%*%V))^2))
  l <- list(U,V,V_cha,U_cha,W_cha,sqrt(rmse),beta0,step1_time,step2_time,step3_time,t_cv)
  names(l) <- c("U","V","V_cha","U_cha","W_cha","RSS","beta0","step1_time","step2_time","step3_time","CV_time")
  return(l)
}
