#' @title NMF_robust
#' @description A robust formulation of NMF using L21-norm.
#' @param delta \code{matrix} the input matrix needs to be decomposed.
#' @param max_iter max iterations.
#' @param thresh threshold of convergence
#' @param V_chu intialization of matrix V.
#' @param U_chu initializaion of matrix U.
#' @param K cluster number.
#' @param alpha damping term.
#'
#' @return the difference between estimatio and real data.
#' @export
#'
#' @examples
NMF_robust <- function(delta, max_iter,thresh,V_chu = NA,U_chu = NA,K,alpha = 0.6){
  if(sum(c(is.na(V_chu), is.na(U_chu)))!=0){
    T <- NCOL(delta)
    p_T <- NROW(delta)
    avg <- sqrt(mean(delta)/K)
    V_chu <- abs(avg*matrix(rnorm(T*K),nrow = K, ncol = T))
    U_chu <- abs(avg*matrix(rnorm(p_T*K),nrow = p_T, ncol = K))
  }
  iter <- 0
  cha1 <- 1e-4
  delta_hat <- U_chu %*% V_chu
  L <- sqrt(apply(delta-delta_hat,2,function(x){return(sum(x^2))}))
  cha2 <- sum(L)
  Ud <- U_chu
  V <- V_chu
  while(iter<=max_iter & abs(cha1-cha2)/cha1 > thresh){
    cha1 <- cha2
    D <- 1/L
    temp1 <- delta * D
    temp1 <- apply(temp1,2,function(x){
      if(sum(x==0)>0){
        x[which(x==0)] <- runif(length(which(x==0)),0,1e-6)
      }
      return(x)
    })
    temp2 <- delta_hat *D
    temp2 <- apply(temp2,2,function(x){
      if(sum(x==0)>0){
        x[which(x==0)] <- runif(length(which(x==0)),0,1e-6)
      }
      return(x)
    })
    Ud1 <- Ud
    Ud <- Ud1*(1-alpha) + Ud * (temp1%*%t(V))/(temp2%*%t(V))*(alpha)
    V <- V*(1-alpha) + V * (t(Ud1)%*%temp1)/(t(Ud1)%*%temp2)*(alpha)
    delta_hat <- Ud%*%V
    L <- sqrt(apply(delta-delta_hat,2,function(x){return(sum(x^2))}))
    cha2 <- sum(L)
    iter <- iter+1
  }
  return(list("U" = Ud,"V" = V,"L" = L))
}
