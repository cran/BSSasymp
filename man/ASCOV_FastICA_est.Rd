\name{ASCOV_FastICAdefl_est}
\alias{ASCOV_FastICAdefl_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of adaptive FastICA estimate}
\description{
Adaptive FastICA solves the blind source separation problem in the case of independent components. This function computes an estimate of the covariance matrix of adaptive fastICA estimate for the mixing or the unmixing matrix.
}
\usage{
ASCOV_FastICAdefl_est(X, gs, dgs, mixed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix.}
  \item{gs}{the first derivative functions of the set of contrast functions.}
  \item{dgs}{the second derivative functions of the set of contrast functions.}
  \item{mixed}{logical, see details.}
}
\details{ If \code{mixed} is TRUE, then \code{X} will be transformed by the adaptive FastICA estimate. The option FALSE can be used, for example, to estimate the covariance when \code{X} are source estimates given by some other method than FastICA.}
\value{
 A list with the following components:
  \item{W}{estimated mean of the unmixing matrix estimate.}
  \item{COV_W}{estimated covariance matrix of the unmixing matrix estimate.}
  \item{A}{estimated mean of the mixing matrix estimate.}
  \item{COV_A}{estimated covariance matrix of the mixing matrix estimate.}
}
\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen S. (2013), Deflation-based FastICA with adaptive contrast functions, manuscript.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}

}
\author{Jari Miettinen}
\seealso{\link{ASCOV_FastICAdefl}, \link[fICA]{adapt_fICA}}

\examples{
# source components have uniform- and exponential(1)- distribution  
s1 <- runif(1000,-sqrt(3),sqrt(3))
s2 <- rexp(1000)
S <- cbind(s1,s2)

# couple of contrast functions
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}

g_gaus <- function(x){(x)*exp(-(x)^2/2)}
dg_gaus <- function(x){exp(-(x)^2/2)-(x)^2*exp(-(x)^2/2)}

gs <- c(g_pow3,g_gaus)
dgs <- c(dg_pow3,dg_gaus)

A<-matrix(rnorm(4),2,2)

X <- S \%*\% t(A)

round(1000*ASCOV_FastICAdefl_est(X, gs=gs, dgs=dgs)$COV_W,2)

}
\keyword{multivariate}
