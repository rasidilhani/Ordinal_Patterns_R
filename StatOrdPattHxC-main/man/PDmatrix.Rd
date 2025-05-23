% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDmatrix.R
\name{PDmatrix}
\alias{PDmatrix}
\title{Computes the matrix of partial derivatives
a given vector of probabilities without zero entries, and
a given entropy}
\usage{
PDmatrix(q, ent, beta = 0.9)
}
\arguments{
\item{q}{vector of probabilities}

\item{ent}{type of entropy: : S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)}

\item{beta}{(only required for Tsallis and Rényi; default value 0.9)}
}
\description{
Computes the matrix of partial derivatives
a given vector of probabilities without zero entries, and
a given entropy
}
\keyword{internal}
