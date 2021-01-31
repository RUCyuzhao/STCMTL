#' @title SCMTL Algorithm(logit regression)
#' @description
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
#' @param nfold number of folds - default is 5.
#' @param thresh_inner convergence threshold of CDS.
#'
#' @return The estimated cluster center matrix U, membership matrix V and coefficient matrix W.
#' @export
#'
#' @examples
#' x_train <- as.matrix(USPS$train$x)
#' y_train <- USPS$train$y[,1]
#' group <- USPS$group
#' re <- cluster_classify(x_train,y_train,group=group,K=2,ext.prop = 5/(length(group)-1),max_iter = 10,thresh = 1e-3,show = TRUE,thresh_inner = 1e-6)

cluster_classify <- function(X,Y,group,K,ext.prop=NULL,pure.prop = 0.5,alpha = 1,max_iter = 50,thresh=1e-2,show=F,CV=TRUE,nfolds = 5,thresh_inner = 1e-3){
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
    lambda[i] <- cv.glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],nfolds = 5,type.measure = "auc",standardize=FALSE,family = "binomial",lambda = seq(5e-5,5e-1,length = 200))$lambda.min/10
    fit <- glmnet(X[group[i]:(group[i+1]-1),-1],Y[group[i]:(group[i+1]-1)],lambda = lambda[i],family = "binomial")
    beta0[,i]<- as.array(coef(fit))[1:p]
  }
  delta <- I%*%sigmoid(y_multi*X%*%beta0)
  soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
  V <- t(soup.out$memberships[[1]])
  ind <- apply(V,2,function(x){return(sum(x^2)==0)})
  if(sum(ind,na.rm = T)>0 | sum(is.na(ind))>0){
    temp <- sample(1:K,1)
    V[temp,which(ind)] <- 1
  }
  U <- matrix(0,nrow = p,ncol = K)
  for(i in 1:p){
    U[i,] <- lsei(A=t(V),B=t(beta0)[,i],type = 2)$X
  }
  t_cv <- Sys.time()-t1
  jilu <- c()
  choice <- permutations(K)
  le <- mean(sigmoid(Y*Z%*%as.vector(U%*%V)))
  V_cha <- c(sum(svd(V)$d^2))
  U_cha <- c(sum(svd(U)$d^2))
  L <- c(sum(log(1+exp(-Y*Z%*%as.vector(U%*%V)))))
  L_cha <- c(L)
  iter <- 0
  panding <- 1
  t_v <- c()
  t_u <- c()
  while(L_cha[iter+1] >=thresh & alpha > 0 & iter<= max_iter){
    t1 <- Sys.time()
    V1 <- V
    U1 <- U
    beta_temp <- U1%*%V1
    for(i in 1:T){
      X_temp <- X[group[i]:(group[i+1]-1),]
      Y_temp <- Y[group[i]:(group[i+1]-1)]
      temp <- 2
      if(iter == 0)
        temp <- 1
      temp <- coordinate_descent_logit(X_temp,Y_temp,max_iter = temp,lambda=lambda[i],beta_temp = beta_temp[,i] ,thresh = 1e-5,sigma = 0.02,cont = 0.5)*alpha+(1-alpha)*beta_temp[,i]
      beta_temp[,i] <- temp[1:p]
    }
    delta <- I%*%sigmoid(y_multi*X%*%beta_temp)

    soup.out <- SOUP(t(delta),Ks=K,type = "count",ext.prop = ext.prop,pure.prop = pure.prop)
    V <- t(soup.out$memberships[[1]])
    ind <- apply(V,2,function(x){return(sum(x^2)==0)})
    if(sum(ind,na.rm = T)>0 | sum(is.na(ind))>0){
      temp <- sample(1:K,1)
      V[temp,which(ind)] <- 1
    }
    t_v[iter+1] <- Sys.time()-t1
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

    for(i in 1:K){
      U[,i] <- CDS(X= X,Y_temp = Y,max_iter = 1000,thresh = thresh_inner,lambda = mean(lambda)/20,U = U,V=V,l = i,group = group,cont = 0.5,sigma = 0.02)
    }
    U_cha[iter+2] <- sqrt(sum(svd(U1-U)$d^2)/sum(svd(U1)$d^2))
    t_u[iter+1] <- Sys.time()-t1
    le[iter+2] <- mean(sign(sigmoid(Z%*%as.vector(U%*%V))-0.5)!=Y)
    L[iter+2] <- sum(log(1+exp(-Y*Z%*%as.vector(U%*%V))))
    L_cha[iter+2] <- abs(L[iter+2]-L[iter+1])/L[iter+1]
    iter <- iter+1
    if(show == TRUE){
      print(paste0("第",iter,"次循环 ","收敛：",L_cha[iter+1]," ER值：",le[iter+1]," V更新时间:",t_v[iter]," U更新时间：",t_u[iter]))
    }
  }
  l <- list(U,V,U_cha,V_cha,le,L,beta0,t_cv,t_u,t_v)
  names(l) <- c("U","V","U_cha","V_cha","ER","ob_function","beta0","time_cv","time_u","time_v")
  return(l)
}
