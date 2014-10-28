\name{ASCOV_FastICAsym_est}
\alias{ASCOV_FastICAsym_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of symmetric FastICA estimate}
\description{
Symmetric FastICA solves the blind source separation problem in the case of independent components. This function computes an estimate of the covariance matrix of symmetric FastICA mixing and unmixing matrix estimates.
}
\usage{
ASCOV_FastICAsym_est(X, G, g, dg, mixed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix.}
  \item{G}{the integral function of the nonlinearity function.}
  \item{g}{the nonlinearity function.}
  \item{dg}{the first derivative function of the nonlinearity function.}
  \item{mixed}{logical, see details.}
}
\details{ If \code{mixed} is TRUE, then \code{X} will be transformed by the symmetric FastICA estimate. The option FALSE can be used, for example, to estimate the covariance when \code{X} are source estimates given by some other method than FastICA.}
\value{
 A list with the following components:
  \item{W}{estimated mean of the unmixing matrix estimate.}
  \item{COV_W}{estimated covariance matrix of the unmixing matrix estimate.}
  \item{A}{estimated mean of the mixing matrix estimate.}
  \item{COV_A}{estimated covariance matrix of the mixing matrix estimate.}
}
\references{
\cite{Hyv\"arinen, A. (1999), Fast and robust fixed-point algorithms for independent component analysis, \emph{{IEEE} Transactions on Neural Networks}, 
10, 626-634.}

\cite{Wei, T. (2014), The convergence and asymptotic analysis of the generalized symmetric FastICA algorithm, http://arxiv.org/abs/1408.0145.}

}
\author{Jari Miettinen}
\seealso{\link{ASCOV_FastICAsym}, \link[fICA]{fICA}}

\examples{
# source components have uniform- and exponential(1)- distribution  
n <- 1000
s1 <- runif(n,-sqrt(3),sqrt(3))
s2 <- rexp(n)
S <- cbind(s1,s2)

# couple of nonlinearities
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}
G_pow3 <- function(x){x^4/4}

g_gaus <- function(x){x*exp(-x^2/2)}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}
G_gaus <- function(x){-exp(-x^2/2)}

A<-matrix(rnorm(4),2,2)

X <- S \%*\% t(A)

round(n*ASCOV_FastICAsym_est(X,G=G_pow3,g=g_pow3,dg=dg_pow3)$COV_W,2)
round(n*ASCOV_FastICAsym_est(X,G=G_gaus,g=g_gaus,dg=dg_gaus)$COV_W,2)

}
\keyword{multivariate}
