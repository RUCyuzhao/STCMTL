#' @title New Sigmoid Function
#' @description A more flexible sigmoid function: \eqn{k/(k+exp(-x))}.
#' @param k a single number
#' @param x can be vecotr or number
#'
#' @return function value
#' @export
#'
sigmoid_new <- function(k,x){
  return(k/(k+exp(-x)))
}
