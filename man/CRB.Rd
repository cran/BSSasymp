\name{CRB}
\alias{CRB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cramer-Rao bound for the unmixing matrix estimate in the independent component model.}
\description{
Cramer-Rao bound for the unmixing matrix estimate in the independent component model.}
\usage{
CRB(sdf,supp=NULL,A=NULL,eps=1e-04,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{a list of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}. Can be \code{NULL} if each support is the whole real line.}
  \item{A}{the mixing matrix, identity by default.}
  \item{eps}{a value which is used when the derivative functions of the density functions are approximated.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{
Let \eqn{\hat{W}} denote an unmixing matrix estimate. If the estimate is affine equivariant, then the matrix \eqn{\hat{G}=\hat{W}A} does not depend on the mixing matrix \eqn{A} and the estimated independent components are \eqn{\hat{S}=\hat{G}S}, where \eqn{S} is the matrix of the true independent components.
}
\value{
A list containing the following components:
\item{CRLB}{A matrix whose elements give the Cramer-Rao lower bounds for the asymptotic variances of the corresponding elements of \eqn{\hat{G}}.}
\item{FIM}{The Fisher information matrix.}
\item{EMD}{The sum of the Cramer-Rao lower bounds of the off-diagonal elements of \eqn{\hat{G}}.}
}
\references{
\cite{Ollila, E., Kim, H. J. and Koivunen, V. (2008), Compact Cramer-Rao bound expression for independent component analysis. \emph{IEEE Transactions on Signal Processing}, 56(4), 1421--1428.}
}
\author{Jari Miettinen}

\examples{
# source components have t(9)- and Gaussian distribution  

f1<-function(x)
{
 gamma(5)*(1+(x*sqrt(9/7))^2/9)^(-5)/
 (sqrt(9*pi/(9/7))*gamma(9/2))
}

f2<-function(x)
{
 exp(-(x)^2/2)/sqrt(2*pi)
}

CRB(sdf=c(f1,f2))

}
\keyword{multivariate}
