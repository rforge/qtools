\name{predict.rq.counts}
\alias{predict.rq.counts}

\title{
Predictions from rq.counts Objects
}

\description{
This function computes predictions based on fitted linear quantile model.
}

\usage{
\method{predict}{rq.counts}(object, newdata,
	na.action = na.pass, ...)
}

\arguments{
  \item{object}{
an \code{rq.counts} object.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
 \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.}

  \item{\dots}{
further arguments passed to \code{boot.lqm}.
}

}

\value{
a vector or a matrix or an array of predictions.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{residuals.rq.counts}}, \code{\link{rq.counts}}, \code{\link{coef.rq.counts}}
}

\keyword{prediction}
