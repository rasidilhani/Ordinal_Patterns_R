% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma2q.R
\name{sigma2q}
\alias{sigma2q}
\title{Function that computes the asymptotic variance of ordinal patterns' entropies}
\usage{
sigma2q(TS, emb, ent = "S", beta = 1.5)
}
\arguments{
\item{TS}{time series}

\item{emb}{embedding dimension}

\item{ent}{type of entropy: S (Shannon, default), T (Tsallis, \eqn{\beta} is required), R (Rényi, \eqn{\beta} is required), or F (Fisher)}

\item{beta}{the parameter for the Tsallis and Rényi entropies (default 1.5)}
}
\value{
The asymptotic variance
}
\description{
Function that computes the asymptotic variance of ordinal patterns' entropies
}
\examples{
# Time series
set.seed(1234567890, kind="Mersenne-Twister")
x <- rnorm(1000) # white noise
y <- mov.av(x, order=51) # smoothed with moving averages 
# Asymptotic variances of their Shannon entropies
sigma2q(x, emb=4, ent="S")
sigma2q(y, emb=4, ent="S")

}
