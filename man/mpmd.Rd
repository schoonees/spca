% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpmd.R
\name{mpmd}
\alias{mpmd}
\title{Multifactor Penalized Matrix Decomposition}
\usage{
mpmd(
  Z,
  k = 2,
  c1 = 1,
  c2 = 1,
  maxit = 20,
  eps = sqrt(.Machine$double.eps),
  centre = FALSE
)
}
\arguments{
\item{Z}{Matrix to be decomposed}

\item{k}{Required rank of the result}

\item{c1}{L1-norm bound for U (greater than or equal to 1), either length-1, or 
with k entries (one for each component).}

\item{c2}{L1-norm bound for V (greater than or equal to 1), either length-1, or 
with k entries (one for each component)}

\item{maxit}{Maximum number of iterations}

\item{eps}{Stopping criterion, an absolute error tolerance on the mean squared reconstruction error}

\item{centre}{Logical indicating whether to centre the matrix Z using the overall mean 
before analysis}
}
\value{
A list with the penalized singular value decomposition (\code{d}, \code{U}, \code{V}).
}
\description{
Multifactor version of the rank-1 L1-penalized matrix decomposition in 
\code{\link{pmd}}. Components are found by applying \code{\link{pmd}} 
to the original matrix after deflation by deducting the components already 
found.
}
\examples{

## Simple random test matrix
set.seed(1)
Z <- matrix(rnorm(100), nrow = 20, ncol = 5)

## Ordinary SVD (equivalent up to changes in sign)
mpmd(Z, c1 = 5, c2 = 3, k = 5)
svd(Z)

## Test with constant c1 and c2
mpmd(Z, c1 = 2, c2 = 1.25, k = 5)

## Test with different c1 and c2 for different components
mpmd(Z, c1 = c(3, 1.25), c2 = c(2, 1.25), k = 2)

}
