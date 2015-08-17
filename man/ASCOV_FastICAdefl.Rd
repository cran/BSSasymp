\name{ASCOV_FastICAdefl}
\alias{ASCOV_FastICAdefl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrices of different deflation-based FastICA estimates}
\description{
The regular deflation-based FastICA finds the independent components one by one using a nonlinearity function. The adaptive deflation-based FastICA chooses, for each component separately, the best nonlinearity from a set of nonlinearities. This function computes asymptotic covariance matrices of the different deflation-based FastICA mixing and unmixing matrix estimates.
}
\usage{
ASCOV_FastICAdefl(sdf, gs, dgs, Gs=NULL, method="adapt", 
                  name=NULL, supp=NULL, A=NULL, ...)
}

\arguments{
  \item{sdf}{a list of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{gs}{a list of nonlinearity functions.}
  \item{dgs}{the first derivative functions of the nonlinearity functions.}
  \item{Gs}{the integral function of the nonlinearity function. Is needed only when \code{method="G"}.}
  \item{method}{"adapt", "G" or "regular", see details.}
  \item{name}{a list of strings, which labels the nonlinearities.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}. Can be \code{NULL} if each support is the whole real line.}
  \item{A}{the mixing matrix, identity matrix as default.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{
Depending on the argument \code{method}, the function computes the asymptotic covariance matrices for three different extraction orders of the independent components. The choice \code{method="adapt"} picks the adaptive deflation-based FastICA estimate, which extracts the components in asymptotically optimal order and uses the best nonlinearity from the set of nonlinearities \code{gs}. 
The other two methods use only one nonlinearity, and if \code{gs} and \code{dgs} contain more than one function, the first one is taken. When \code{method="G"}, the order is based on the deviance from normality measured by \code{Gs}. When \code{method="regular"}, the order is that of \code{sdf}.

The signs of the components are fixed so that the sum of the elements of each row of the unmixing matrix is positive.

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
  \item{used_gs}{indicates which nonlinearity is used in estimation of each rows of the unmixing matrix.}
}
\references{
\cite{Ilmonen, P., Nordhausen, K., Oja, H. and Ollila, E. (2010): A New Performance Index for ICA: Properties, Computation and Asymptotic Analysis. 
In Vigneron, V., Zarzoso, V., Moreau, E., Gribonval, R. and Vincent, E. (editors) \emph{Latent Variable Analysis and Signal Separation}, 229--236, Springer.}

\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014), Deflation-based FastICA with adaptive choices of nonlinearities, \emph{IEEE Transactions on Signal Processing}, 62(21), 5716--5724.}

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

# couple of nonlinearities
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}
G_pow3 <- function(x){x^4/4}

g_gaus <- function(x){x*exp(-x^2/2)}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}

gs <- c(g_pow3,g_gaus)
dgs <- c(dg_pow3,dg_gaus)

A <- matrix(rnorm(4),2,2)

res1 <- ASCOV_FastICAdefl(sdf=c(fu,fe), gs=gs, dgs=dgs, supp=supp, A=A)

round(res1$COV_W, 2)
res1$EMD
res1$used_gs

res2 <- ASCOV_FastICAdefl(sdf=c(fu,fe), gs=c(g_pow3), dgs=c(dg_pow3), 
Gs=c(G_pow3), method="G", supp=supp, A=A)

res2$EMD
}

\keyword{multivariate}
