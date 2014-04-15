\name{predict.tsrq}
\alias{predict.tsrq2}
\alias{predict.rcrq}
\alias{predict.nlrq2}
\alias{predict.tsrq}

\title{
Predictions from Quantile Regression Transformation Models
}

\description{
This function computes predictions based on fitted quantile regression transformation models.
}

\usage{
\method{predict}{tsrq}(object, newdata, na.action = na.pass, raw = TRUE, ...)
\method{predict}{tsrq2}(object, newdata, na.action = na.pass, raw = TRUE, ...)
\method{predict}{rcrq}(object, newdata, na.action = na.pass, raw = TRUE, ...)
\method{predict}{nlrq2}(object, newdata, na.action = na.pass, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{tsrq}, \code{tsrq2} or \code{nlrq2}.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
}
  \item{raw}{
logical flag. If \code{TRUE} (default) the predictions on the transformed scale \emph{Xb} are returned; otherwise, the values are back-transformed using the inverse of the transformation \emph{h} used to fit the model, i.e. \emph{hinv(Xb)}.
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


\keyword{ predict }

