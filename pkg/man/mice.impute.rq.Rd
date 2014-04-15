\name{mice.impute.rq}
\alias{mice.impute.rrq}
\alias{mice.impute.rq}

\title{
QR-based Multiple Imputation
}

\description{
This function is used to multiply impute missing values using quantile regression imputation models.
}

\usage{
mice.impute.rq(y, ry, x, tsf = NULL, lambda = NULL, symm = TRUE,
	epsilon = 0.001, method.rq = "fn", ...)
mice.impute.rrq(y, ry, x, tsf = NULL, lambda = NULL, symm = TRUE,
	epsilon = 0.001, method.rq = "fn", ...)
}

\arguments{
  \item{y}{
numeric vector of length \code{n} with \code{nmis} missing values.
}
  \item{ry}{
missing data indicator. Logical vector of length \code{n}: \code{FALSE} if \code{y} is missing, \code{TRUE} if \code{y} is observed.
}
  \item{x}{
matrix \code{n x p} of completely observed covariates.
}
  \item{tsf}{
transformation to be used (see also \code{\link{ao}}). Possible options are \code{bc} for Box-Cox and \code{ao} for Aranda-Ordaz transformation models. For Proposal I transformation models (Geraci and Jones), \code{mcjI} (singly bounded outcome) or \code{mcjIb} (doubly bounded outcome).
}
  \item{lambda}{
transformation parameter.
}
  \item{symm}{
logical flag. If \code{TRUE} (default) a symmetric transformation is used.
}
  \item{epsilon}{
constant used to trim the values of the sample space.
}
  \item{method.rq}{
linear programming algorithm (see \code{\link{rq}}).
}
  \item{\dots}{
additional arguments.
}
}

\details{
This function implements the methods proposed by Geraci (2013) to impute missing values using quantile regression models. Uniform values are sampled from \emph{[epsilon, 1 - epsilon]}, therefore allowing the interval to be bounded away from 0 and 1 (default is 0.001). It is possible to specify a quantile regression transformation model with parameter \code{lambda} (Geraci and Jones). The function \code{mice.impute.rrq} performs imputation based on restricted
regression quantiles to avoid quantile crossing (see Geraci 2013 for details).
}

\value{

A vector of length \code{nmis} with imputations.

}

\references{
Geraci M. Estimation of regression quantiles in complex surveys with data missing at random: An application to birthweight determinants. Statistical Methods in Medical Research 2013.
doi:10.1177/0962280213484401

Geraci M and Jones MC. Improved transformation-based quantile regression. Unpublished Manuscript

}
\author{
Marco Geraci
}

\keyword{ multiple imputation}
\keyword{conditional quantiles}