\name{alphas}
\alias{alphas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic variances of the deflation-based FastICA estimate}
\description{
Deflation-based FastICA solves the blind source separation problem in the case of independent components. These function computes some interesting theoretic quantities related to the deflation-based FastICA unmixing matrix estimate, see details.
}
\usage{
alphas(sdf, gs, dgs, name=NULL, supp=NULL,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf}{a list of density functions of the sources scaled so that the mean is 0 and variance is 1.}
  \item{gs}{a list of nonlinearity functions.}
  \item{dgs}{the first derivative functions of the nonlinearity functions.}
  \item{name}{a list of strings, which labels the nonlinearities.}
  \item{supp}{a two column matrix, where each row gives the lower and the upper limit used in numerical integration for the corresponding source component which is done using \code{\link[stats]{integrate}}. Can be \code{NULL} if each support is the whole real line.}
  \item{\dots}{arguments to be passed to \code{\link[stats]{integrate}}.}
}
\details{
When the mixing matrix is the identity matrix, the asymptotic variances of the first row elements of the deflation-based FastICA estimate depend only on the corresponding source component and the chosen nonlinearity function g. Furthermore, the asymptotic variances of the off-diagonal elements of the first row are equal, let us call this value alpha. Also the other asymptotic variances depend straightforwardly on alphas corresponding to different components and to (possibly) different nonlinearities. Alphas indicate which nonlinearities should be used and in which order the components should be separated. Reloaded (Nordhausen et al., 2011) and adaptive (Miettinen et al., 2014) deflation-based FastICA estimators are based on the estimation of alphas from the data. 
}
\value{
A matrix where the ith row gives the alphas for the ith nonlinearity and the jth column corresponds to the jth density in \code{sdf}.}
\references{
\cite{Miettinen, J., Nordhausen, K., Oja, H. and Taskinen, S. (2014), Deflation-based FastICA with adaptive choices of nonlinearities, \emph{IEEE Transactions on Signal Processing}, 62(21), 5716--5724.}

\cite{Nordhausen, K., Ilmonen, P., Mandal, A., Oja, H. and Ollila, E. (2011), Deflation-based FastICA reloaded,
in \emph{Proc. "19th European Signal Processing Conference 2011 (EUSIPCO 2011)",} Barcelona, 1854--1858.}

}
\author{Jari Miettinen}
\seealso{\link{ASCOV_FastICAdefl}, \link[fICA]{adapt_fICA}}

\examples{
# source components have uniform- and exponential(1)- distribution  
fu <- function(x){1/(sqrt(3)*2)}
fe <- function(x){exp(-x-1)}

supp <- matrix(c(-sqrt(3),sqrt(3),-1,Inf), nrow=2, ncol=2, byrow=TRUE)

# couple of nonlinearities
g_pow3 <- function(x){x^3}
dg_pow3 <- function(x){3*x^2}

g_gaus <- function(x){x*exp(-x^2/2)}
dg_gaus <- function(x){exp(-x^2/2)-x^2*exp(-x^2/2)}

gs <- c(g_pow3,g_gaus)
dgs <- c(dg_pow3,dg_gaus)

alphas(sdf=c(fu,fe), gs=c(g_pow3), dgs=c(dg_pow3), supp=supp)

alphas(sdf=c(fu,fe), gs=gs, dgs=dgs, supp=supp)

}
\keyword{multivariate}
