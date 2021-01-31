#' @title Updating matrix U(logit)
#' @description Using coordinate descent to update one of the vector of matrix U.
#' @param X_temp \code{matrix} dependent variable
#' @param Y_temp \code{matrix} independent variable
#' @param max_iter \code{double} max iterations
#' @param beta_temp \code{vector} the initialization of the coefficient vector.
#' @param thresh \code{double} convergence threshold for the algorithm.
#' @param lambda \code{double} penalty parameter.
#' @param cont \code{double} step for backforward line search.
#' @param sigma \code{double} threshold for backforward line search.
#' @param U \code{matrix} initialization of matrix U.
#' @param V \code{matrix} initialization of matrix V.
#' @param l which column of matrix U to update.
#' @param group the task structure.
#'
#' @return the updated lth column of matrix U.
#' @export
#'
#' @examples
CDS <- function(X,Y_temp,max_iter,beta_temp=NULL,thresh=1e-5,lambda = 0.001,U,V,l,group,cont,sigma){
  n <- length(Y_temp)
  p <- NCOL(X)
  if(is.null(beta_temp))
    beta_temp <- U[,l]
  k <- c()
  X_temp <- X
  for(j in 1:(length(group)-1)){
    X_temp[group[j]:(group[j+1]-1),] <- V[l,j] * X[group[j]:(group[j+1]-1),]
    if(NROW(V)>2){
      k[group[j]:(group[j+1]-1)] <- exp(Y_temp[group[j]:(group[j+1]-1)]*X_temp[group[j]:(group[j+1]-1),]%*%(U[,-l]%*%V[-l,j]))
    }else{
      k[group[j]:(group[j+1]-1)] <- exp(Y_temp[group[j]:(group[j+1]-1)]*X_temp[group[j]:(group[j+1]-1),]%*%(U[,-l]*V[-l,j]))
    }

  }
  j <- 0
  start_1 <- c()
  cha <- 1
  while(j <= max_iter & cha >thresh){
    start_1<- beta_temp
    ind_zero <- which(beta_temp==0)
    ind_sam <- sample(setdiff(1:p,ind_zero),p-length(ind_zero),replace = FALSE)
    for(ind in ind_sam){
      start_2 <- beta_temp
      del1 <- 1/n/lambda*sum((sigmoid_new(k,Y_temp*(X_temp%*%beta_temp))-1)*Y_temp*X_temp[,ind])
      del2 <- 1/n/lambda*t(X_temp[,ind])%*%((sigmoid_new(k,Y_temp*(X_temp%*%beta_temp)))*(1-sigmoid_new(k,Y_temp*(X_temp%*%beta_temp)))*X_temp[,ind])
      if(del1+1 <= del2 * beta_temp[ind]){
        start_2[ind] <- beta_temp[ind]-(del1+1)/del2
      }else if(del1-1 >= del2*beta_temp[ind]){
        start_2[ind] <- beta_temp[ind]-(del1-1)/del2
      }else{
        start_2[ind] <- 0
      }
      if(start_2[ind]!=0 ){
        temp <-  mean(log(1+exp(-Y_temp*X_temp%*%beta_temp)))+lambda*sum(abs(beta_temp))
        temp_beta <- start_2[ind]-beta_temp[ind]
        alpha <- 1
        beta_temp[ind] <- start_2[ind]-(1-alpha)*temp_beta
        cha1 <- mean(log(1+exp(-Y_temp*X_temp%*%start_2)))+lambda*sum(abs(start_2))- temp
        cha2 <- sigma*lambda*(del1*temp+abs(start_2[ind])-abs(beta_temp[ind]))
        m <- 1
        if(cha1>0 & cha2 <0){
          # print("出错了")
          alpha <- 0
          beta_temp[ind] <- start_2[ind]-(1-alpha)*temp_beta
        }else{
          while(cha1 > cha2 & alpha > 1e-5){
            alpha <- cont^m
            beta_temp[ind] <- start_2[ind]-(1-alpha)*temp_beta
            cha1 <- mean(log(1+exp(-Y_temp*X_temp%*%beta_temp)))+lambda*sum(abs(beta_temp))- temp
            m <- m+1
            cha2 <- alpha*cha2
          }
        }
        #print(paste0("m：",m))
      }else{
        beta_temp <- start_2
      }

    }
    cha <- sqrt(mean((beta_temp-start_1)^2))
    # print(paste0("收敛：",cha))
    j <- j+1
  }
  return(beta_temp)
}
