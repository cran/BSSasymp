\name{ASCOV_FastICAdefl_est}
\alias{ASCOV_FastICAdefl_est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic covariance matrices of deflation-based FastICA estimates}
\description{
The regular deflation-based FastICA finds the independent components one by one using a nonlinearity function. The adaptive deflation-based FastICA chooses, for each component separately, the best nonlinearity from a set of nonlinearities. This function computes estimates of the covariance matrices of the different deflation-based FastICA mixing and unmixing matrix estimates.
}
\usage{
ASCOV_FastICAdefl_est(X, gs, dgs, Gs=NULL, method="adapt", 
                      name=NULL, mixed=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data matrix.}
  \item{gs}{a list of nonlinearity functions.}
  \item{dgs}{the first derivative functions of the nonlinearity functions.}
  \item{Gs}{the integral function of the nonlinearity function. Is needed only when \code{method="G"}.}
  \item{method}{"adapt" or "G", see details.}
\item{name}{a list of strings, which labels the nonlinearities.}
  \item{mixed}{logical, see details.}
}
\details{ 
Depending on the argument \code{method}, the function computes the asymptotic covariance matrices for two different extraction orders of the independent components. The choice \code{method="adapt"} picks the adaptive deflation-based FastICA estimate, which extracts the components in asymptotically optimal order and uses the best nonlinearity from the set of nonlinearities \code{gs}. 
When \code{method="G"}, the order is based on the deviance from normality measured by \code{Gs}. This method uses only one nonlinearity, and if \code{gs} and \code{dgs} contain more than one function, the first one is taken. 

If \code{mixed} is TRUE, then \code{X} will be transformed by the adaptive FastICA estimate. The option FALSE can be used, for example, to estimate the covariance when \code{X} are source estimates given by some other method than FastICA.}
\value{
 A list with the following components:
  \item{W}{estimated mean of the unmixing matrix estimate.}
  \item{COV_W}{estimated covariance matrix of the unmixing matrix estimate.}
  \item{A}{estimated mean of the mixing matrix estimate.}
  \item{COV_A}{estimated covariance matrix of the mixing matrix estimate.}
  \item{used_gs}{indicates which nonlinearity is used in estimation of each rows of the unmixing matrix.}
}
\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014), Deflation-based FastICA with adaptive choices of nonlinearities, \emph{IEEE Transactions on Signal Processing}, 62(21), 5716--5724.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}

}
\author{Jari Miettinen}
\seealso{\link{ASCOV_FastICAdefl}, \link[fICA]{adapt_fICA}}

\examples{
# source components have uniform- and exponential(1)- distribution  
s1 <- runif(1000,-sqrt(3),sqrt(3))
s2 <- rexp(1000)
S <- cbind(s1,s2)

# couple of nonlinearities
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}
G_pow3 <- function(x){x^4/4}

g_gaus <- function(x){x*exp(-x^2/2)}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}

gs <- c(g_pow3,g_gaus)
dgs <- c(dg_pow3,dg_gaus)

A<-matrix(rnorm(4),2,2)

X <- S \%*\% t(A)

round(1000*ASCOV_FastICAdefl_est(X, gs=gs, dgs=dgs)$COV_W,2)

round(1000*ASCOV_FastICAdefl_est(X, gs=c(g_pow3), dgs=c(dg_pow3),
Gs=c(G_pow3), method="G")$COV_W,2)

}
\keyword{multivariate}
