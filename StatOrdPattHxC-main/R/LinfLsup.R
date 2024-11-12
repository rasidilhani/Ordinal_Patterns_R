#' Boundaries of the \eqn{H\times C} closed manifold
#'
#' A data frame with the lower and upper boundaries of the \eqn{H\times C} closed manifold
#' for \eqn{D=3,4,5,6}.
#'
#' @format ## `LinfLsup`
#' A data frame with 184936 rows and 4 variables:
#' \describe{
#'  \item{H}{entropy (value)}
#'  \item{C}{statistical complexity (value)}
#'  \item{Dimension}{embedding dimension (factor, "3", "4", "5", "6")}
#'  \item{Side}{Lower or Upper side of the feasible region (factor, "Lower" and "Upper")}
#' }
#' @source <https://dx.doi.org/10.1016/j.physa.2005.11.053>
#' 
#' Check this package's vignette for an example.
#' LinfLsup
"LinfLsup"
