#' Function that computes the p-value of the hypothesis test for two time series
#'
#' @param TSx a time series
#' @param TSy another time series
#' @param emb embedding dimension
#' @param ent type of entropy: S (Shannon, default), T (Tsallis, beta is required), R (Rényi, beta is required), or F (Fisher)
#' @param beta parameter for the Tsallis and Rényi entropy (default 0.9)
#'
#' @importFrom stats pnorm
#'
#' @export
#'

pval <- function(TSx, TSy, emb, ent="S", beta=0.9){

  # Compute the length of the OP sequences

  n_x <- length(TSx) - emb + 1
  n_y <- length(TSy) - emb + 1

  # Compute the OP probability vectors

  q_x <- OPprob(TS = TSx, emb = emb)
  q_y <- OPprob(TS = TSy, emb = emb)

  # Compute asymptotic means and variances

  switch(ent,

         "S" = {

           mu_x <- HShannon(q_x)
           sig2_x <- sigma2q(TS = TSx, emb = emb, ent = "S")

           mu_y <- HShannon(q_y)
           sig2_y <- sigma2q(TS = TSy, emb = emb, ent = "S")

         },


         "T" = {

           mu_x <- HTsallis(q_x, beta = beta)
           sig2_x <- sigma2q(TS = TSx, emb = emb, ent = "T", beta = beta)

           mu_y <- HTsallis(q_y, beta = beta)
           sig2_y <- sigma2q(TS = TSy, emb = emb, ent = "T", beta = beta)

         },

         "R" = {

           mu_x <- sum(q_x^beta)
           sig2_x <- sigma2q(TS = TSx, emb = emb, ent = "R", beta = beta)

           mu_y <- sum(q_y^beta)
           sig2_y <- sigma2q(TS = TSy, emb = emb, ent = "R", beta = beta)

         },

         "F" = {

           mu_x <- HFisher(q_x)
           sig2_x <- sigma2q(TS = TSx, emb = emb, ent = "F")

           mu_y <- HFisher(q_y)
           sig2_y <- sigma2q(TS = TSy, emb = emb, ent = "F")

         }

  )
  # Compute the statistic

  s <- (mu_x - mu_y) / sqrt(sig2_x / n_x + sig2_y / n_y)

  # Compute the p-value

  return(2 - 2 * pnorm(abs(s)))
}

