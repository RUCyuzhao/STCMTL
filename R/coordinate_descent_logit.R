#' @title CD algorithm(Logit)
#' @description Using coordinate decent algorithm to solve l1 regularized logistic regression.
#' @param X_temp \code{matrix} dependent variable
#' @param Y_temp \code{matrix} independent variable
#' @param max_iter \code{double} max iterations
#' @param beta_temp \code{vector} the initialization of the coefficient vector.
#' @param thresh \code{double} convergence threshold for the algorithm.
#' @param lambda \code{double} penalty parameter.
#' @param cont \code{double} step for backforward line search.
#' @param sigma \code{double} threshold for backforward line search.
#'
#' @return estimated coefficient vector.
#' @export
#'
#' @examples
coordinate_descent_logit <- function(X_temp,Y_temp,max_iter,beta_temp=NULL,thresh=1e-5,lambda = 0.001,cont,sigma){
  n <- length(Y_temp)
  p <- NCOL(X_temp)
  if(is.null(beta_temp)){
    beta_temp <- rep(1,p)
  }
  j <- 0
  cha <- 1
  while(j <= max_iter & cha >thresh){
    start_1 <- beta_temp
    for(ind in 1:p){
      start_2 <- beta_temp
      del1 <- 1/n/lambda*sum((sigmoid(Y_temp*(X_temp%*%beta_temp))-1)*Y_temp*X_temp[,ind])
      del2 <- 1/n/lambda*t(X_temp[,ind])%*%(sigmoid(Y_temp*(X_temp%*%beta_temp))*(1-sigmoid(Y_temp*(X_temp%*%beta_temp)))*X_temp[,ind])
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
        cha1 <- mean(log(1+exp(-Y_temp*X_temp%*%start_2)))+lambda*sum(abs(start_2))- temp
        cha2 <- sigma*lambda*(del1*temp+abs(start_2[ind])-abs(beta_temp[ind]))
        m <- 1
        if(cha1 <= cha2){
          beta_temp[ind] <- start_2[ind]
        }else{
          while(cha1 > cha2){
            alpha <- cont^m
            beta_temp[ind] <- start_2[ind]-(1-alpha)*temp_beta
            cha1 <- mean(log(1+exp(-Y_temp*X_temp%*%beta_temp)))+lambda*sum(abs(beta_temp))- temp
            m <- m+1
            cha2 <- alpha*cha2
          }
        }
        # print(paste0("m：",m))
      }else{
        beta_temp <- start_2
      }
    }
    cha <- sqrt(mean((beta_temp-start_1)^2))
    j <- j+1
  }
  return(beta_temp)
}
