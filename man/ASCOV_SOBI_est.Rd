\name{ASCOV_SOBI_est}
\alias{ASCOV_SOBI_est}
\alias{ASCOV_SOBIdefl_est}
\alias{ASCOV_SOBI_estN}
\alias{ASCOV_SOBIdefl_estN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of symmetric and deflation-based SOBI estimates}
\description{
The symmetric and deflation-based SOBI methods solve the blind source separation problem in the case of second order stationary time series sources by jointly diagonalizing the covariance matrix and several autocovariance matrices at different lags. The functions compute an estimate of the covariance matrix of a SOBI estimate for the mixing or the unmixing matrix, under the assumption that the sources are \eqn{MA(\infty)} time series. Notice that, since AMUSE method is a special case of SOBI, also an estimate of the covariance matrix of an AMUSE estimate can be computed using these functions. 
}
\usage{
ASCOV_SOBI_estN(X, taus, mixed=TRUE, M=100)

ASCOV_SOBI_est(X, taus, arp=NULL, maq=NULL, mixed=TRUE, 
               M=100, ...)

ASCOV_SOBIdefl_estN(X, taus, mixed=TRUE,  M=100)

ASCOV_SOBIdefl_est(X, taus, arp=NULL, maq=NULL, mixed=TRUE, 
                   M=100, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix or a multivariate time series object of class \code{\link{ts}}.}
  \item{taus}{a vector of integers for the lags.}
  \item{arp}{a vector containing the AR orders used for the estimation of ARMA coefficients.}
  \item{maq}{a vector containing the MA orders used for the estimation of ARMA coefficients.}
  \item{mixed}{logical, see details.}
  \item{M}{the number of autocovariance matrices used for the estimation of the covariance matrices, see details.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{arima}}.}
}
\details{ Functions ASCOV_SOBI_estN and ASCOV_SOBIdefl_estN assume that the innovations of the components are gaussian. Therefore, they are faster than ASCOV_SOBI_est and  ASCOV_SOBI-defl_est, which estimate the fourth moments of the innovations by estimating the ARMA coefficients of the time series.
Fitting the univariate ARMA coefficients is done using the function \code{\link[stats]{arima}} based on the orders provided by \code{arp} and \code{maq}.

The estimation is mostly based on autocovariance matrices and all non-zero matrices should be included. On the other hand, too large value of \code{M} increases the computation time and it may even reduce the estimation accuracy.  

If \code{mixed} is TRUE, then \code{X} will be transformed by the corresponding SOBI estimate. The option FALSE can be used, for example, to estimate the covariance when \code{X} are source estimates given by some other method than SOBI.
}
\value{
 A list with the following components:
  \item{W}{estimated mean of the unmixing matrix estimate.}
  \item{COV_W}{estimated covariance matrix of the unmixing matrix estimate.}
  \item{A}{estimated mean of the mixing matrix estimate.}
  \item{COV_A}{estimated covariance matrix of the mixing matrix estimate.}
}
\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2012),
Statistical properties of a blind source separation estimator for stationary
time series, \emph{Statistics and Probability Letters}, 82, 1865--1873.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014),
Deflation-based separation of uncorrelated stationary time series, \emph{Journal of Multivariate Analysis}, 123, 214--227.}

\cite{Miettinen, J., Illner, K., Nordhausen, K., Oja, H., Taskinen, S. and Theis, F. (2014),
Separation of uncorrelated stationary time series using autocovariance matrices,
http://arxiv.org/abs/1405.3388.}

}
\author{Jari Miettinen}

\seealso{\link{ASCOV_SOBI}, \link[JADE]{SOBI}, \link[JADE]{AMUSE}, \link[stats]{arima}}


\examples{
A<- matrix(rnorm(9),3,3)
s1 <- arima.sim(list(ar=0.6),1000)
s2 <- arima.sim(list(ma=c(0.2,0.3,-0.3)),1000)
s3 <- arima.sim(list(ar=-0.2,ma=c(0.5,-0.1,0.4)),1000)

S <- cbind(s1,s2,s3)
X <- S \%*\% t(A)

round(1000*ASCOV_SOBI_estN(X, taus=1:10)$COV_W,2)
round(1000*ASCOV_SOBIdefl_estN(X, taus=1:10)$COV_W,2)

}

\keyword{multivariate}
\keyword{ ts }
