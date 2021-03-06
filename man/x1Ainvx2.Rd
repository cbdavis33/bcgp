% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/miscFuncs.R
\name{x1Ainvx2}
\alias{x1Ainvx2}
\title{Multiply vector times matrix inverse times vector.}
\usage{
x1Ainvx2(x1, A, x2)
}
\arguments{
\item{x1}{A list of vectors of length \emph{n}. Generally, in the bcgp setting,
\emph{n} will be the number of training data locations.}

\item{A}{An \emph{n x n} matrix. In the \code{bcgp} setting, \emph{A} is a
covariance matrix.}
}
\value{
x2 A list of vectors of length \emph{n}. Generally, in the bcgp setting,
\emph{n} will be the number of training data locations.
}
\description{
\code{x1Ainvx2} returns a list where each element is the result of multiplying
a vector times a matrix inverse times another vector.
}
\details{
This returns a list where each element is the result of \eqn{x_1^{\top}A^{-1}x_2}
}
\examples{
n <- 10
x <- matrix(sort(runif(n)), ncol = 1)
A <- getCorMat(x, rho = 0.6)
x1 <- list(rep(1, n), rep(1,n))
x2 <- list(rep(1, n), MASS::mvrnorm(1, rep(0, n), A))
x1Ainvx2(x1, A, x2)
}
