#' Statistical Complexity of a probability function
#'
#' Given a probability function \eqn{\bm p}, its statistical complexity
#' is the product of its entropy \eqn{H(\bm p)} and its distance to an
#' equilibrium distribution \eqn{\bm u}.
#'
#' The normalized Jensen-Shannon distance to the uniform distribution
#' \eqn{\bm u=(1/k, 1/k,\dots,1/k)}, where \eqn{k\geq 2} is the size of \eqn{p},
#' is \eqn{Q=Q'/\max{Q'}}, where
#' \deqn{Q'(\bm{p}, \bm{u}) = \sum_{\ell=1}^{k} \Big[\big(p(\ell)-\frac{1}{k}\big) \ln\big(k p(\ell)\big)\Big],}
#' and
#' \deqn{\max Q' = -2\Big[\frac{k!+1}{k!} \ln(k!+1) - 2 \ln(2k!) + \ln k!\Big].}
#'
#' @export
#'
#' @param p a probability function
#' @usage StatComplexity(p)
#' @returns a value in \eqn{[0,1]}
#' 
#' @examples
#' 
#' # Time series
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' y <- mov.av(x, order=50) # smoothed with moving averages
#' # Statistical Complexity of patterns of dimension 4
#' StatComplexity(OPprob(x, emb=4))
#' StatComplexity(OPprob(y, emb=4))
#'

StatComplexity <- function(p){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){
    
    k <- length(p)
    prestricted <- p[p>0] # Discard entries with zeroes
    
    Q0 <- -2 / (
      (k+1)/k * log(k+1) -
        2 * log(2*k) +
        log(k)
    )
    
    Div <- -sum(((p+1/k)/2) * log((p+1/k)/2)) + 
      sum((prestricted*log(prestricted)))/2 - log(k)/2 
    
    HS <- HShannon(p)
    
    return(Q0*Div*HS)
    
  } else {
    print("ERROR: Not a valid probability function")
  }
}