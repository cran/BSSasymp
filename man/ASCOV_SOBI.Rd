\name{ASCOV_SOBI}
\alias{ASCOV_SOBI}
\alias{ASCOV_SOBIdefl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrix of symmetric and deflation-based SOBI estimates}
\description{
The symmetric and deflation-based SOBI methods solve the blind source separation problem in the case of second order stationary time series sources by jointly diagonalizing the covariance matrix and several autocovariance matrices at different lags. The functions compute the asymptotic covariance matrices of a SOBI estimates for the mixing or the unmixing matrices, when the sources are \eqn{MA(\infty)} time series. Notice that, since AMUSE method is a special case of SOBI, also the asymptotic covariance matrix of an AMUSE estimate can be computed using these functions. 
}
\usage{
ASCOV_SOBI(psi, taus, Beta=NULL, A=NULL)

ASCOV_SOBIdefl(psi, taus, Beta=NULL, A=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{a numeric matrix containing the MA-coefficients of the time series, see details.}
  \item{taus}{a vector of integers for the lags.}
  \item{Beta}{a matrix of fourth moments of the innovations, see details.}
  \item{A}{the mixing matrix, identity matrix as default.}
}
\details{ 
Naturally, the function can deal with only finite number of coefficients. The larger is the number of the rows of psi, the more accurate is the result, but also the longer is the computation time. 

AR or ARMA coefficients can be transformed to MA coefficients by using \code{\link[stats]{ARMAtoMA}}.  

The \eqn{ij}th entry of Beta is \eqn{E[Z_{it}^2Z_{jt}^2]}, where \eqn{Z_{it}} are the innovations of \eqn{i}th source component such that \eqn{E[Z_{it}]=0} and \eqn{E[Z_{it}^2]=1}.      
 
The order of the estimated components is fixed so that the sums of their squared autocovariances over \code{taus} are in a decreasing order. The signs of the components are fixed so that the sum of the elements of each row of the unmixing matrix is positive.
}
\value{
 A list with the following components:
  \item{W}{mean of the unmixing matrix estimate.}
  \item{COV_W}{asymptotic covariance matrix of the unmixing matrix estimate.}
  \item{A}{mean of the mixing matrix estimate.}
  \item{COV_A}{asymptotic covariance matrix of the mixing matrix estimate.}
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

\seealso{\link{ASCOV_SOBI_est}, \link[JADE]{SOBI}, \link[JADE]{AMUSE}, \link[stats]{ARMAtoMA}}


\examples{

A<- matrix(rnorm(9),3,3)
psi1 <- ARMAtoMA(ar=0.6, lag.max=100)
psi2 <- ARMAtoMA(ar=c(0.2,0.3,-0.3),lag.max=100)
psi3 <- ARMAtoMA(ar=-0.2, ma=c(0.5,-0.1,0.4), lag.max=100)
psi <- cbind(psi1,psi2,psi3)

round(ASCOV_SOBI(psi=psi, taus=1:10, A=A)$COV_W,2)
round(ASCOV_SOBIdefl(psi=psi, taus=1:10, A=A)$COV_W,2)

}


\keyword{multivariate}
\keyword{ ts }
