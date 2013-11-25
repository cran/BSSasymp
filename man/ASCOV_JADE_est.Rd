\name{ASCOV_JADE_est}
\alias{ASCOV_JADE_est}
\alias{ASCOV_FOBI_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of JADE and FOBI estimates}
\description{
JADE solves the blind source separation problem in the case of independent components with at most one component having kurtosis values zero, while FOBI requires distinct kurtosis values. The functions compute the asymptotic covariance matrices of JADE and FOBI estimates for the mixing or the unmixing matrix.
}
\usage{
ASCOV_JADE_est(X, mixed=TRUE)

ASCOV_FOBI_est(X, mixed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrixx.}
  \item{mixed}{logical, see details.}
}
\details{ If \code{mixed} is TRUE, then \code{X} will be transformed by the corresponding estimate. The option FALSE can be used, for example, to estimate the covariance when \code{X} are source estimates given by some other method than JADE or FOBI.
}

\value{
 A list with the following components:
  \item{W}{estimated mean of the unmixing matrix estimate.}
  \item{COV_W}{estimated covariance matrix of the unmixing matrix estimate.}
  \item{A}{estimated mean of the unmixing matrix estimate.}
  \item{COV_A}{estimated covariance matrix of the unmixing matrix estimate.}
}
\references{
\cite{Ilmonen, P., Nevalainen, J. and Oja, H. (2010), Characteristics of multivariate
distributions and the invariant coordinate system, \emph{Statistics and Probability Letters}, 80, 1844--1853.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen S. (2013), 4th Moments and Their Applications in ICA, manuscript.}
}
\author{Jari Miettinen}

\seealso{\link{ASCOV_JADE}, \link{ASCOV_FOBI}, \link[JADE]{JADE}, \link[JADE]{FOBI}}


\examples{
# source components have t-10-, uniform- and gaussian distribution 
s1 <- rt(1000,10)/sqrt(10/8)
s2 <- runif(1000,-sqrt(3),sqrt(3))
s3 <- rnorm(1000)
S <- cbind(s1,s2,s3)
A <- matrix(rnorm(9),3,3)

X <- S \%*\% t(A)
round(1000*ASCOV_JADE_est(X)$COV_W,2)
round(1000*ASCOV_FOBI_est(X)$COV_W,2)

}

\keyword{multivariate}
