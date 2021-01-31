#' @title Outlier task detection
#' @description Find the existence of outlier task through simple threshold.
#' @param L the F norm difference between W and UV, according to NMF_robust algorithm.
#' @param percent the minimum percent of classifiable tasks.
#' @param cut threshold of outlier.
#'
#' @return detected outlier tasks.
#' @export
#'
#' @examples
find_change <- function(L,percent = 0.8,cut = 4){
  paixu <- rank(L)
  ind <- min(paixu[which(L > quantile(L,percent))])
  while(ind < length(L) & mean(L[which(paixu<=ind)])+cut*sd(L[which(paixu<=ind)])>mean(L[which(paixu>ind)])){
    ind <- ind+1
  }
  return(which(paixu>=ind))
}
