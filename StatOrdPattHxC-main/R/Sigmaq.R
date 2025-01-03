#' Function that computes the covariance matrix for a time series
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param ent type of entropy: S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)
#' @param beta parameter for the Tsallis and Rényi entropies
#'
#' @name Sigmaq
#'
#' @keywords internal
#'

utils::globalVariables("PDmatrix")

Sigmaq <- function(TS, emb, ent, beta){
  # Compute Sigma matrix
  S <- Sigma(TS, emb)

  # Delete null rows and null columns
  SR <- S[rowSums(S[,]) != 0,]
  SC <- SR[,colSums(SR[,]) != 0]

  # Compute vector of probabilities
  qvec <- OPprob(TS, emb)
  qvecPos <- qvec[qvec>0]

  # Compute matrix of partial derivatives
  J <- PDmatrix(q = qvecPos, ent, beta)

  return(J %*% SC %*% t(J))
}



