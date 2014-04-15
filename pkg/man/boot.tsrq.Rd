\name{boot.tsrq}
\alias{boot.tsrq2}
\alias{boot.nlrq2}
\alias{boot.tsrq}

\title{
Bootstrapping Tranformation Models
}

\description{
This function is used to obtain a bootstrap sample of a fitted tranformation quantile regression model.
}

\usage{
boot.tsrq(object, R = 50, seed = round(runif(1, 1, 10000)))
boot.tsrq2(object, R = 50, seed = round(runif(1, 1, 10000)))
boot.nlrq2(object, R = 50, seed = round(runif(1, 1, 10000)))
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{tsrq}, \code{tsrq2} or \code{nlrq2}.
}
  \item{R}{
number of bootstrap replications.
}
  \item{seed}{
optional random generator seed.
}
}


\value{
An array of dimension \code{c(R,p+q)} or \code{c(R,p+q,nt)}, where \code{R} is the number of bootstrap replications, \code{p} is the number of regression coefficients, \code{q} is the number of transformation parameters in the model and \code{nt = length(object$tau)}.

The following attributes are available:

\item{tau}{index of the quantile(s).}
\item{lambda}{the estimated transformation parameter(s) as given by \code{object}.}
\item{estimated}{the estimated regression coefficients and transformation parameter(s) as given by \code{object}.}
\item{R}{number of bootstrap replications.}
\item{seed}{the random number generator seed used to produce the bootstrap sample.}
\item{npars}{numer \code{p} of regression coefficients.}
\item{ntot}{total number \code{p+q} of parameters.}
\item{indices}{the bootstrap sample of data units.}
\item{rdf}{residual degrees of freedom.}
\item{term.labels}{names of coefficients.}
}

\author{
Marco Geraci
}

\keyword{ bootstrap }
