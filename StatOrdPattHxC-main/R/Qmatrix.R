#' Function that computes the Ordinal Patterns transitions matrices
#'
#' Let \eqn{\bm \pi=(\pi_1,\pi_2,\dots,\pi_n)} be the sequence of ordinal patterns
#' from the time series \eqn{\bm x=(x_1,x_2,\dots,x_{n-m+1})} computed over windows
#' of embedding dimension \eqn{m}. For each lag \eqn{\ell =1,2, \dots , m-1}, the ordinal
#' pattern transition matrix of order \eqn{\ell} is
#' \eqn{\textbf{Q}^{(\ell)} \in \mathbb{R}^{m! \times m!}}, with elements
#' \deqn{q_{ij}^{(\ell)} = \Pr\big(\pi_t=\pi^{(i)} \wedge \pi_{t+\ell}=\pi^{(j)}) = q_i \Pr(\pi_{t+\ell} = \pi^{(j)} \mid \pi_t=\pi^{(i)}\big),}
#' for \eqn{i,j=1,2,\dots ,m!}.
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param lag the lag
#' @returns an \eqn{m!\times m!} matrix of probabilities
#'
#' @export
#' 

Qmatrix <- function(TS, emb, lag){
  
  k <- factorial(emb)
  
  OP <- OPseq(TS, emb, lag)
  
  Qfreq <- matrix(0, nrow = k, ncol = k)
  
  for(i in 1:k){
    for(j in 1:k){
      s <- which(OP == i)
      r <- which(OP[s+1] == j)
      Qfreq[i,j] <- length(r)
    }
  }
  
  return(Qfreq / (length(OP) - 1))
}


