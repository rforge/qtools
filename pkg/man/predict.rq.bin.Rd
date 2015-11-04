\name{predict.rq.bin}
\alias{predict.rq.bin}

\title{
Predictions from rq.bin Objects
}

\description{
This function computes predictions based on fitted binary quantile models.
}

\usage{
\method{predict}{rq.bin}(object, newdata, na.action = na.pass,
	type = "latent", grid = TRUE, ...)
}

\arguments{
  \item{object}{
an \code{rq.counts} object.
}

  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}

 \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}

 \item{type}{
the type of prediction required. See details.
}

 \item{grid}{
either a vector of breakpoints or a logical flag. If \code{TRUE}, the breakpoints are set to \code{seq(0.05, 0.95, by = 0.05)}. If \code{FALSE}, the breakpoints are extracted from \code{object$tau}. At least 2 breakpoints must be provided.
}


  \item{\dots}{
not used.
}

}

\value{
a vector or a matrix of predictions.
}

\details{

If \code{type = "latent"} (default), the predictions are for the conditional quantiles of the latent response variable; if \code{type = "probability"}, the predictions are for the conditional probabilities of the observable response. Note that \code{type = "probability"} is possible only when \code{normalize = "last"} in \code{\link{rq.bin}} call. The probabilities are predicted over a \code{grid} of values.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{rq.bin}}, \code{\link{coef.rq.bin}}
}

\keyword{predict}
