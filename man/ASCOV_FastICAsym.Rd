\name{ASCOV_FastICAsym}
\alias{ASCOV_FastICAsym}
\alias{ASCOV_FastICAsym2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of symmetric FastICA estimates}
\description{
Symmetric FastICA estimators solves the blind source separation problem in the case of independent components. These functions compute the asymptotic covariance matrices of the regular and the squared symmetric FastICA mixing and unmixing matrix estimates.

}
\usage{
ASCOV_FastICAsym(sdf, G, g, dg, supp=NULL, A=NULL, ...)

ASCOV_FastICAsym2(sdf, G, g, dg, supp=NULL, A=NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{a list of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{G}{the integral function of the nonlinearity function, see details.}
  \item{g}{the nonlinearity function.}
  \item{dg}{the first derivative function of the nonlinearity function.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}. Can be \code{NULL} if each support is the whole real line.}
  \item{A}{the mixing matrix, identity matrix as default.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{ The signs of the components are fixed so that the sum of the elements of each row of the unmixing matrix is positive.

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
  \item{EMD}{the limiting expected value of \eqn{n(p-1)MD^2}, see details.}
}
\references{
\cite{Hyv\"arinen, A. (1999), Fast and robust fixed-point algorithms for independent component analysis, \emph{{IEEE} Transactions on Neural Networks}, 
10, 626-634.}

\cite{Wei, T. (2014), The convergence and asymptotic analysis of the generalized symmetric FastICA algorithm, http://arxiv.org/abs/1408.0145.}

\cite{Miettinen, J., Nordhausen, K., Oja, H., Taskinen, S. and Virta, J. (2015), The squared symmetric FastICA estimator, submitted. }

}
\author{Jari Miettinen}

\seealso{\link{ASCOV_FastICAsym_est}, \link[fICA]{fICA}, \link[stats]{integrate}}

\examples{
# source components have uniform- and exponential(1)- distribution  
fu <- function(x){1/(sqrt(3)*2)}
fe <- function(x){exp(-x-1)}

supp <- matrix(c(-sqrt(3),sqrt(3),-1,Inf), nrow=2, ncol=2, byrow=TRUE)

# couple of nonlinearities
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}
G_pow3 <- function(x){x^4/4}

g_gaus <- function(x){x*exp(-x^2/2)}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}
G_gaus <- function(x){-exp(-x^2/2)}

A <- matrix(rnorm(4),2,2)

res1 <- ASCOV_FastICAsym(sdf=c(fu,fe), G=G_pow3, g=g_pow3, dg=dg_pow3, 
supp=supp, A=A)
res2 <- ASCOV_FastICAsym(sdf=c(fu,fe), G=G_gaus, g=g_gaus, dg=dg_gaus, 
supp=supp, A=A)

round(res1$COV_W, 2)
res1$EMD

round(res2$COV_W, 2)
res2$EMD
}

\keyword{multivariate}
