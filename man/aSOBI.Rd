\name{aSOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{aSOBI}
\title{Alternative SOBI estimators}
\description{
The SOBI method solves the blind source separation problem in the case of second order stationary time series sources by jointly diagonalizing the covariance matrix and several autocovariance matrices. In the classical SOBI method, the sum of squares of the diagonal elements is used as diagonality measure. This function computes the SOBI estimate, when a choice from a family of alternative diagonality criteria is used.   
}
\usage{
aSOBI(X, k=12, a=4, eps=1e-06, maxiter=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix or a multivariate time series object of class \code{\link{ts}}.}
  \item{k}{if a single integer, then the lags 1:k are used, if an integer vector, then these are used as the lags.}
  \item{a}{numeric, determines the diagonality criterion, see details.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}
\details{ 
The classical SOBI estimator maximizes the sum of squares of the diagonal elements of the autocovariance matrices. A family of alternative SOBI estimators is obtained when the diagonality of \eqn{p\times p}{pxp} matrices \eqn{B_1,\dots,B_K}{B_1,...,B_K} is measured by
\deqn{\sum_{k=1}^K\sum_{i=1}^p |(B_k)_{ii}|^a} 
with \eqn{a>1}. The choice \eqn{a=2} gives the classical SOBI estimate and it can be computed more efficiently using function \code{\link[JADE]{SOBI}}.
}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{S}{estimated sources standardized to have mean 0 and unit variances.}
  \item{k}{lags used.}
  \item{a}{value of the diagonality criterion parameter used.}
}
\references{
\cite{Miettinen, J. (2015): Alternative diagonality criteria for SOBI.
In Nordhausen, K. and Taskinen, S. (editors), \emph{Modern Nonparametric, Robust and Multivariate methods, Festschrift in Honour of Hannu Oja}, 455--469, Springer.}
}
\author{Jari Miettinen}

\seealso{\link[JADE]{SOBI}}


\examples{
A <- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=0.6),1000)
s2 <- arima.sim(list(ma=c(0.2,0.3,-0.3)),1000)
s3 <- arima.sim(list(ar=-0.2,ma=c(0.5,-0.1,0.4)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

aSOBI(X, a=3)
}

\keyword{multivariate}
\keyword{ ts }
