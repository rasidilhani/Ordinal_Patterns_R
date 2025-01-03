#' Function that assigns the label i to an
#'
#' @param pat an OP, for example, pat = 0,1,2 if \eqn{x_1 < x_2 < x_3}
#' @keywords internal

pi_i <- function(pat){
  a <- length(pat) - 1
  op <- as.data.frame(perm(0:a))
  return(row.match(pat, op))
}


