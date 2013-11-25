\name{ASCOV_FastICAdefl}
\alias{ASCOV_FastICAdefl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of adaptive FastICA estimate}
\description{
Adaptive FastICA solves the blind source separation problem in the case of independent components. This function computes the asymptotic covariance matrix of adaptive fastICA estimate for the mixing or the unmixing matrix.
Note that reloaded FastICA is a special case of adaptive fastICA.
}
\usage{
ASCOV_FastICAdefl(sdf, gs, dgs, supp=NULL, A=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{a vector of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{gs}{the first derivative functions of the set of contrast functions.}
  \item{dgs}{the second derivative functions of the set of contrast functions.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}.}
  \item{A}{the mixing matrix, identity matrix as default.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{ The signs of the components are fixed so that the sum of the elements of each row of the unmixing matrix is positive.
}
\value{
 A list with the following components:
  \item{W}{mean of the unmixing matrix estimate.}
  \item{COV_W}{asymptotic covariance matrix of the unmixing matrix estimate.}
  \item{A}{mean of the mixing matrix estimate.}
  \item{COV_A}{asymptotic covariance matrix of the mixing matrix estimate.}
}
\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen S. (2013), Deflation-based FastICA with adaptive contrast functions, manuscript.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}


}
\author{Jari Miettinen}

\seealso{\link{ASCOV_FastICAdefl_est}, \link[fICA]{adapt_fICA}, \link[stats]{integrate}}

\examples{
# source components have uniform- and exponential(1)- distribution  
fu <- function(x){1/(sqrt(3)*2)}
fe <- function(x){exp(-x-1)}

supp <- matrix(c(-sqrt(3),sqrt(3),-1,Inf), nrow=2, ncol=2, byrow=TRUE)
# couple of contrast functions
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}

g_gaus <- function(x){(x)*exp(-(x)^2/2)}
dg_gaus <- function(x){exp(-(x)^2/2)-(x)^2*exp(-(x)^2/2)}

gs <- c(g_pow3,g_gaus)
dgs <- c(dg_pow3,dg_gaus)

A <- matrix(rnorm(4),2,2)

round(ASCOV_FastICAdefl(sdf=c(fu,fe), gs=gs, dgs=dgs, supp=supp, A=A)$COV_W, 2)

}

\keyword{multivariate}
