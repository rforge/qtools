\name{predict.rqt}
\alias{predict.rqt}

\title{
Predictions from Quantile Regression Transformation Models
}

\description{
This function computes predictions based on fitted quantile regression transformation models.
}

\usage{
\method{predict}{rqt}(object, newdata, na.action = na.pass,
type = "response", ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{rqt}.
}

  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}

  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}

 \item{type}{
the type of prediction required. The default "response" is on the scale of the response variable, i.e. the values are back-transformed using the inverse of the transformation \emph{hinv(Xb)}; the alternative "link" is on the scale of the linear predictors \emph{h(y) = Xb}.
}

 \item{\dots}{
not used.
}
}


\value{
a vector or a matrix or an array of predictions.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{tsrq}}
}

\keyword{predict}

