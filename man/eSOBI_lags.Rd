\name{eSOBI_lags}
%- Also NEED an '\alias' for EACH other topic documented here.
\alias{eSOBI_lags}
\alias{taus_def}
\title{The default set of lag sets for the efficient SOBI estimator}
\description{
The set of lag sets to be used as the default set in function \code{\link{eSOBI}.} 
}
\usage{
taus_def
}
\details{
The object \code{taus_def} contains the following lag sets:
\tabular{ll}{
taus_def[[1]  \tab 1,2,\dots,12 \cr
taus_def[[2]] \tab 1 \cr
taus_def[[3]] \tab 1,2 \cr
taus_def[[4]] \tab 1,2,3 \cr
taus_def[[5]] \tab 1,2,3,4,5 \cr
taus_def[[6]] \tab 1,2,\dots,8 \cr
taus_def[[7]] \tab 1,2,\dots,20 \cr
taus_def[[8]] \tab 1,2,\dots,50 \cr
taus_def[[9]] \tab 1,2,\dots,10,12,\dots,20 \cr
taus_def[[10]] \tab 5,6,\dots,10,12,\dots,20,25,\dots,50 \cr
taus_def[[11]] \tab 2,4,\dots,20 \cr
taus_def[[12]] \tab 1,2,3,5,7,11,13,17,19}
}
\references{
\cite{Taskinen, S., Miettinen, J. and Nordhausen, K. (2016),
Efficient second order blind identification method for separation of uncorrelated stationary time series, submitted.}
}
\author{Jari Miettinen}

\seealso{\link{eSOBI}}
