#' Function that computes the covariance matrix of the Multinomial model 
#' 
#' @param p a probability function
#' @returns A covariance matrix
#' 
#' @export

SigmaMultinomial <- function(p){
  return(diag(p) - p %*% t(p))
}