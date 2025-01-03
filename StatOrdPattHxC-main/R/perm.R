#' Function that finds all possible ordinal patterns of the
#' elements \eqn{0, 1, ... , m-1}
#'
#' @param v a vector \eqn{0, 1, \dots , m-1}
#'
#' @keywords internal

perm <- function(v) {
  n <- length(v)
  if (n == 1) v
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(v[-i])))
    X
  }
}

