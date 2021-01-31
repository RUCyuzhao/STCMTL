#' @title CD algorithm(Regression)
#' @description  Using coordinate decent algorithm to solve lasso problem.
#' @param X_temp \code{matrix} dependent variable.
#' @param Y_temp \code{matrix} independent variable.
#' @param max_iter \code{double} max iterations.
#' @param beta_temp \code{vector} the initialization of the coefficient vector.
#' @param thresh \code{double} convergence threshold for the algorithm.
#' @param lambda \code{double} penalty parameter.
#'
#' @return \code{vector} the estimated coefficient
#' @export
#'
#' @examples
coordinate_decent <- function(X_temp,Y_temp,max_iter,beta_temp=NULL,thresh=1e-5,lambda = 0.001){
  n <- length(Y_temp)
  p <- NCOL(X_temp)
  if(is.null(beta_temp)){
    beta_temp <- rep(1,p)
  }
  j <- 0
  start_1 <- c()
  start_2 <- c()
  cha <- 1
  while(j <= max_iter & cha >thresh){
    ind <- sample(which(beta_temp!=0),1)
    Z_ind <- sum(X_temp[,ind]^2)
    rho <- t(X_temp[,ind])%*%(Y_temp - X_temp %*% beta_temp+beta_temp[ind]*X_temp[,ind])
    start_1[j+1] <- beta_temp[ind]
    if(rho < -lambda/2*n){
      beta_temp[ind] <- (rho+lambda/2*n)/Z_ind
    }else if(rho > lambda/2*n){
      beta_temp[ind] <- (rho-lambda/2*n)/Z_ind
    }else{
      beta_temp[ind] <- 0
    }
    start_2[j+1] <- beta_temp[ind]
    if(j > 2*p){
      cha <- mean(abs(start_1[(j-10):(j+1)]-start_2[(j-10):(j+1)]))
    }
    j <- j+1
  }
  return(beta_temp)
}
