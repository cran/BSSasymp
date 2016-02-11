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
ASCOV_SOBI(psi, taus, a=2, Beta=NULL, A=NULL)

ASCOV_SOBIdefl(psi, taus, Beta=NULL, A=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psi}{a numeric matrix containing the MA-coefficients of the time series, see details.}
  \item{taus}{a vector of integers for the lags.}
  \item{a}{numeric, see details.}
  \item{Beta}{a matrix of fourth moments of the innovations, see details.}
  \item{A}{the mixing matrix, identity matrix as default.}
}
\details{ 
Naturally, the function can deal with only finite number of coefficients. The larger is the number of the rows of psi, the more accurate is the result, but also the longer is the computation time. 

AR or ARMA coefficients can be transformed to MA coefficients by using \code{\link[stats]{ARMAtoMA}}.  

The \eqn{ij}th entry of Beta is \eqn{E[Z_{it}^2Z_{jt}^2]}, where \eqn{Z_{it}} are the innovations of \eqn{i}th source component such that \eqn{E[Z_{it}]=0} and \eqn{E[Z_{it}^2]=1}.      
 
The order of the estimated components is fixed so that the sums of their squared autocovariances over \code{taus} are in a decreasing order. The signs of the components are fixed so that the sum of the elements of each row of the unmixing matrix is positive.

Since the unmixing matrix has asymptotic normal distribution, we have a connection between the asymptotic variances and the minimum distance index, which is defined as 
\deqn{MD(\hat{W},A)=\frac{1}{\sqrt{p-1}} \inf_{P D}{||PD \hat{W} A-I||,}}{MD(W.hat,A) = 1/sqrt(p-1) inf_(P D) ||P D W.hat A - I||,}
where \eqn{\hat{W}}{W.hat} is the unmixing matrix estimate, \eqn{A}{A} is the mixing matrix, \eqn{P}{P} is a permutation matrix and \eqn{D}{D} a diagonal matrix with nonzero diagonal entries. If \eqn{\hat{W}A}{W.hat A} converges to the identity matrix, the limiting expected value of \eqn{n(p-1)MD^2} is the sum of the asymptotic variances of the off-diagonal elements of \eqn{\hat{W}A}{W.hat A}. Here \eqn{n}{n} is the sample size and \eqn{p}{p} is the number of components.

The symmetric SOBI estimator maximizes the sum of squares of the diagonal elements of the autocovariance matrices. Different SOBI estimators are obtained when the diagonality of \eqn{p\times p}{pxp} matrices \eqn{B_1,\dots,B_K}{B_1,...,B_K} is measured by
\deqn{\sum_{k=1}^K\sum_{i=1}^p |(B_k)_{ii}|^a} 
with \eqn{a>1}. The diagonality measure can be selected using the argument \code{a}. 
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
\cite{Ilmonen, P., Nordhausen, K., Oja, H. and Ollila, E. (2010): A New Performance Index for ICA: Properties, Computation and Asymptotic Analysis. 
In Vigneron, V., Zarzoso, V., Moreau, E., Gribonval, R. and Vincent, E. (editors) \emph{Latent Variable Analysis and Signal Separation}, 229--236, Springer.}

\cite{Miettinen, J. (2015): Alternative diagonality criteria for SOBI.
In Nordhausen, K. and Taskinen, S. (editors) \emph{Modern Nonparametric, Robust and Multivariate methods, Festschrift in Honour of Hannu Oja}, 455--469, Springer.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2012),
Statistical properties of a blind source separation estimator for stationary
time series, \emph{Statistics and Probability Letters}, 82, 1865--1873.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014),
Deflation-based separation of uncorrelated stationary time series, \emph{Journal of Multivariate Analysis}, 123, 214--227.}

\cite{Miettinen, J., Illner, K., Nordhausen, K., Oja, H., Taskinen, S. and Theis, F. (2015),
Separation of uncorrelated stationary time series using autocovariance matrices,
\emph{Journal of Time Series Analysis}, in print, DOI: 10.1111/jtsa.12159.}

}
\author{Jari Miettinen}

\seealso{\link{ASCOV_SOBI_est}, \link[JADE]{SOBI}, \link[JADE]{AMUSE}, \link[stats]{ARMAtoMA}}


\examples{

A<- matrix(rnorm(9),3,3)
psi1 <- ARMAtoMA(ar=0.6, lag.max=100)
psi2 <- ARMAtoMA(ar=c(0.2,0.3,-0.3),lag.max=100)
psi3 <- ARMAtoMA(ar=-0.2, ma=c(0.5,-0.1,0.4), lag.max=100)
psi <- cbind(psi1,psi2,psi3)

sym <- ASCOV_SOBI(psi=psi, taus=1:10, A=A)
defl <- ASCOV_SOBIdefl(psi=psi, taus=1:10, A=A)

round(sym$COV_W,2)
round(defl$COV_W,2)

sym$EMD
defl$EMD
}


\keyword{multivariate}
\keyword{ ts }
