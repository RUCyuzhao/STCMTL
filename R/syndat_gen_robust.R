#' @title Synthetic dataset generator(sparse)
#' @description Generate the synthetic dataset(sparse scenario) in the article.
#' @param n each task's number.
#' @param K cluster number.
#' @param p variable number.
#' @param p0 useless variable number(set to zero).
#' @param pure_num the pure tasks number in each task group.
#' @param seed random number seed.
#' @param sam_num how many different samples do you want to generate.
#' @param out_num outlier task number.
#'
#' @return A large list contains train sample, test sample,coefficient vector and task structure group.
#' @export
#'
#' @examples
syndat_gen_robust <- function(n,K,p,p0,T,pure_num,seed=1,sam_num,out_num){
  set.seed(seed)
  U0 <- matrix(0,nrow = p,ncol=K)
  inter <- floor((p-p0)/(K+1))
  for(i in 1:K){
    U0[((i-1)*inter+1):((i+1)*inter),i] <- runif(2*inter,0.1,0.5)*sample(c(-1,1),inter,replace = T)
  }
  V0 <- matrix(0,nrow = K,ncol = T)
  for(i in 1:K){
    V0[i,((i-1)*pure_num+1):(i*pure_num)] <- 1
  }
  if(K*pure_num < (T-out_num)){
    for(i in (K*pure_num+1):(T-out_num)){
      ind <- sample(1:K,1)
      V0[ind,i] <- runif(1,0.5,1)
      V0[sample(setdiff(1:K,ind),1),i] <- runif(1,0,0.5)
    }
  }

  V0[,-c((T-out_num+1):T)] <- scale(V0[,-c((T-out_num+1):T)],center=F,scale = rowSums(t(V0[,-c((T-out_num+1):T)])))
  x <- list()
  y <- list()
  beta_chu <- U0%*%V0
  for(i in 1:out_num){
    beta_chu[sample(1:p,inter*2),T-out_num+i] <- runif(inter*2,0.5,1)*sample(c(-1,1),inter*2,replace = T)
  }
  beta_chu <- rbind(rep(0,T),beta_chu)
  group <- c(1,cumsum(rep(n*2,T))+1)
  xtrain_list <- list()
  ytrain_list <- list()
  xtest_list <- list()
  ytest_list <- list()
  for(j in 1:sam_num){
    x <- list()
    y <- list()
    set.seed(j)
    for(i in 1:T){
      x[[i]] <- cbind(rep(1,n*2),matrix(rnorm(n*2*p),nrow = n*2,ncol = p))
      y[[i]] <- x[[i]] %*%beta_chu[,i]+rnorm(n*2,sd=0.5)
    }
    X <- do.call(rbind,x)
    Y <- as.array(do.call(rbind,y))
    set.seed(j)
    ind1 <- c()
    ind2 <- c()
    ind3 <- c()
    for(i in 1:(length(group)-1)){
      sam1 <- sample(group[i]:(group[i+1]-1),n)
      ind1 <- c(ind1,sam1)
      sam2 <- sample(setdiff(c(group[i]:(group[i+1]-1)),sam1),n)
      ind2 <- c(ind2,sam2)
    }
    xtrain_list[[j]] <- X[ind1,]
    ytrain_list[[j]] <- Y[ind1,]
    xtest_list[[j]] <- X[ind2,]
    ytest_list[[j]] <- Y[ind2,]
  }
  return(list(train=list(xtrain_list,ytrain_list),test=list(xtest_list,ytest_list),group = c(1,cumsum(rep(n,T))+1),beta = beta_chu))
}
