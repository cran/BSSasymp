\name{eSOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{eSOBI}
\title{The efficient SOBI estimator}
\description{
The SOBI method solves the blind source separation problem in the case of second order stationary time series sources by jointly diagonalizing the covariance matrix and several autocovariance matrices. The separation performance depends on the lag set. The efficient SOBI estimator uses asymptotic results to estimate the variances of the elements of the SOBI unmixing matrices obtained by given lag sets. The unmixing matrix corresponding to the lag set which minimizes the sum of the variances is the efficient SOBI estimate.  
}
\usage{
eSOBI(X, taus=taus_def, M=200, fast=TRUE, eps=1e-06, maxiter=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix or a multivariate time series object of class \code{\link{ts}}.}
  \item{taus}{a list whose components are vectors of integers. The list gives the set of lag sets. The default set is \code{\link{taus_def}}.}
  \item{M}{the number of autocovariance matrices used for the estimation of the variance estimates, see \code{\link{ASCOV_SOBI_estN}}.}
  \item{fast}{logical, see details.}
  \item{eps}{convergence tolerance.}
  \item{maxiter}{maximum number of iterations.}
}
\details{ 
The function includes two versions of the efficient SOBI algorithm. The faster version uses only one SOBI estimate given by the first component in \code{taus} to estimate the sum of the limiting variances for all lag sets. In the other version, which is obtained by \code{fast=FALSE}, the sum of the limiting variances of each lag set is estimated using that particular lag set for the initial SOBI estimate. When the length of the time series is sufficient, say 5000 or more, the two versions yield equally good estimates and the use of the faster version is recommended. Otherwise we recommend the use of the slower version. The variance estimates are based on asymptotic results which assume that the time series are multivariate linear processes. Such processes include a wide class of stationary time series, for example all causal ARMA processes. It is also assumed that the innovations are Gaussian. This simplifies the computations and has practically no effect on which lag set is chosen. If the user does not want to make the Gaussianity assumption, the slower version of the function can be easily modified by replacing the function \code{\link{ASCOV_SOBI_estN}} in the code by \code{\link{ASCOV_SOBI_est}}. If the \code{\link[JADE]{SOBI}} algorithm fails to converge for some lag set, the corresponding value of the estimated sum of variances in \code{sum_var} is \code{Inf}.}
\value{
 A list with class 'bss' containing the following components:
  \item{W}{estimated unmixing matrix.}
  \item{S}{estimated sources standardized to have mean 0 and unit variances.}
  \item{taus_used}{the lag set which is considered best and used for the estimation of the unmixing matrix.}
  \item{sum_var}{estimated sum of variances of the unmixing matrix estimates for all lag sets given in \code{taus}.}
}
\references{
\cite{Miettinen, J., Illner, K., Nordhausen, K., Oja, H., Taskinen, S. and Theis, F. (2015),
Separation of uncorrelated stationary time series using autocovariance matrices, \emph{Journal of Time Series Analysis}, in print, http://arxiv.org/abs/1405.3388.}

\cite{Taskinen, S., Miettinen, J. and Nordhausen, K. (2016),
Efficient second order blind identification method for separation of uncorrelated stationary time series, submitted.}
}
\author{Jari Miettinen}

\seealso{\link{ASCOV_SOBI_est}, \link[JADE]{SOBI}}


\examples{
A <- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=0.6),1000)
s2 <- arima.sim(list(ma=c(0.2,0.3,-0.3)),1000)
s3 <- arima.sim(list(ar=-0.2,ma=c(0.5,-0.1,0.4)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

taus <- list(1,1:2,1:5,1:12)
eSOBI(X, taus=taus)
}

\keyword{multivariate}
\keyword{ ts }
