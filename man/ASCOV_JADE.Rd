\name{ASCOV_JADE}
\alias{ASCOV_JADE}
\alias{ASCOV_FOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of JADE and FOBI estimates}
\description{
JADE solves the blind source separation problem in the case of independent components with at most one component having kurtosis values zero, while FOBI requires distinct kurtosis values. The functions compute the asymptotic covariance matrices of JADE and FOBI estimates for the mixing or the unmixing matrices.
}
\usage{
ASCOV_JADE(sdf, supp=NULL, A=NULL, ...)

ASCOV_FOBI(sdf, supp=NULL, A=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{sdf}{a vector of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}.}
  \item{A}{the mixing matrix, identity matrix as default.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{ The order of the estimated components is fixed so that their kurtosis values are in a decreasing order. The signs of the components is fixed so that the sum of the elements of each row of the unmixing matrix is positive.

Since the unmixing matrix has asymptotic normal distribution, we have a connection between the asymptotic variances and the minimum distance index, which is defined as 
\deqn{MD(\hat{W},A)=\frac{1}{\sqrt{p-1}} \inf_{P D}{||PD \hat{W} A-I||,}}{MD(W.hat,A) = 1/sqrt(p-1) inf_(P D) ||P D W.hat A - I||,}
where \eqn{\hat{W}}{W.hat} is the unmixing matrix estimate, \eqn{A}{A} is the mixing matrix, \eqn{P}{P} is a permutation matrix and \eqn{D}{D} a diagonal matrix with nonzero diagonal entries. If \eqn{\hat{W}A}{W.hat A} converges to the identity matrix, the limiting expected value of \eqn{n(p-1)MD^2} is the sum of the asymptotic variances of the off-diagonal elements of \eqn{\hat{W}A}{W.hat A}. Here \eqn{n}{n} is the sample size and \eqn{p}{p} is the number of components.
}
\value{
 A list with the following components:
  \item{W}{mean of the unmixing matrix estimate.}
  \item{COV_W}{asymptotic covariance matrix of the unmixing matrix estimate.}
  \item{A}{mean of the mixing matrix estimate.}
  \item{COV_A}{asymptotic covariance matrix of the mixing matrix estimate.}
  \item{EMD}{The limiting expected value of \eqn{n(p-1)MD^2}, see details.}
}
\references{
\cite{Ilmonen, P., Nevalainen, J. and Oja, H. (2010), Characteristics of multivariate
distributions and the invariant coordinate system, \emph{Statistics and Probability Letters}, 80, 1844--1853.}

\cite{Ilmonen, P., Nordhausen, K., Oja, H. and Ollila, E. (2010): A New Performance Index for ICA: Properties, Computation and Asymptotic Analysis. 
In Vigneron, V., Zarzoso, V., Moreau, E., Gribonval, R. and Vincent, E. (editors) \emph{Latent Variable Analysis and Signal Separation}, 229--236, Springer.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen S. (2014), Fourth Moments and Independent Component Analysis, http://arxiv.org/abs/1406.4765.}
}
\author{Jari Miettinen}

\seealso{\link{ASCOV_JADE_est}, \link{ASCOV_FOBI_est}, \link[JADE]{JADE}, \link[JADE]{FOBI}}


\examples{
# source components have uniform- and exponential(1)- distribution  
fu <- function(x){1/(sqrt(3)*2)}
fe <- function(x){exp(-x-1)}

supp <- matrix(c(-sqrt(3),sqrt(3),-1,Inf), nrow=2, ncol=2, byrow=TRUE)
A<-matrix(rnorm(4),2,2)

jade <- ASCOV_JADE(sdf=c(fu,fe), supp=supp, A=A)
fobi <- ASCOV_FOBI(sdf=c(fu,fe), supp=supp, A=A)

round(jade$COV_W,2)
round(fobi$COV_W,2)

jade$EMD
fobi$EMD
}
\keyword{multivariate}
