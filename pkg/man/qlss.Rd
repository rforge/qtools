\name{qlss}
\alias{qlss.default}
\alias{qlss.formula}
\alias{qlss.numeric}
\alias{qlss}

\title{
Quantile-based Summary Statistics for Location, Scale and Shape
}

\description{
This function calculates quantile-based summary statistics for location, scale and shape of a distribution, unconditional or conditional.
}

\usage{
qlss(...)
\method{qlss}{default}(fun = "qnorm", probs = 0.1, ...)
\method{qlss}{formula}(formula, data, type = "rq", tsf = NULL, symm = TRUE, dbounded = FALSE,
lambda.p = NULL, delta.p = NULL, lambda.q = NULL, delta.q = NULL,
probs = 0.1, ci = FALSE, R = 500, predictLs = NULL, ...)
\method{qlss}{numeric}(x, probs = 0.1, ...)
}

\arguments{
  \item{fun}{
quantile function.
}
  \item{probs}{
a vector of probabilities.
}
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
}
  \item{data}{
a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
}
  \item{type}{
possible options are \code{rq} for quantile regression (default) or \code{rqt} for transformation-based quantile regression.
}
  \item{tsf}{
transformation to be used. Possible options are \code{mcjI} for Proposal I transformation models (default), \code{mcjI} for Proposal II transformation models, \code{bc} for Box-Cox and \code{ao} for Aranda-Ordaz transformation models.
}
  \item{symm}{
logical flag. If \code{TRUE} (default) a symmetric transformation is used.
}
  \item{dbounded}{
logical flag. If \code{TRUE} the response is assumed to be doubly bounded on [a,b]. If \code{FALSE} the response is assumed to be singly bounded (ie, strictly positive).
}
  \item{lambda.p, delta.p}{
vectors with values of transformation parameters to fit models for quantiles \code{probs}. Must be the same length as \code{probs}.
}
  \item{lambda.q, delta.q}{
vectors with values of transformation parameters to fit models for quantiles \code{1-probs}. Must be the same length as \code{probs}.
}
 \item{ci}{
logical flag. If \code{TRUE}, bootstrapped confidence intervals for the predictions are calculated.
}
 \item{R}{
number of bootstrap replications.
}
\item{x}{
a numeric vector.
}
  \item{predictLs}{
list of arguments for \code{predict.rq}.
}
  \item{\dots}{
other arguments for \code{fun}, \code{rq} or \code{tsrq}.
}
}

\details{
This function computes a number of quantile-based summary statistics for location (median), scale (inter-quartile range), and shape (Bowley skewness and shape index) of a distribution. These statistics can be computed for unconditional and conditional distributions. In the latter case, a \code{formula} specifies a linear quantile function, which is fitted with \code{\link{rq}}. The default \code{qlss} function computes the summary statistics of a standard normal distribution or of any other distribution via the argument \code{fun}. The latter must be a function with \code{p} as its probability argument (see for example \code{\link{qnorm}}, \code{\link{qt}}, \code{\link{qchisq}}, \code{\link{qgamma}}, etc.)

The argument \code{p} 
}

\value{
\code{qlss} returns an object of \code{\link{class}} \code{qlss}. This is a list that contains three elements:

\item{location}{summary statistic(s) for location.}
\item{scale}{summary statistic(s) for scale.}
\item{method}{summary statistic(s) for shape.}

}

\references{
Geraci M and Jones MC. Improved transformation-based quantile regression. Canadian Journal of Statistics 2015;43(1):118-132.

Gilchrist W. Statistical modelling with quantile functions. Chapman and Hall/CRC; 2000.
}

\author{
Marco Geraci
}

\examples{

# Compute summary statistics of a normal distribution
qlss()

# Compute summary statistics of a t distribution with 3 df
qlss(fun = "qt", df = 3, probs = 0.05)

# Compute summary statistics for a sample using a sequence of probabilities
x <- rnorm(1000)
qlss(x, probs = c(0.1, 0.2, 0.3, 0.4))

# Compute summary statistics for Volume conditional on Height
trees2 <- trees[order(trees$Height),]
fit <- qlss(Volume ~ Height, data = trees2)
plot(fit, z = trees2$Height, xlab = "height")

# Use a quadratic model for Height
fit <- qlss(Volume ~ poly(Height,2), data = trees2)
plot(fit, z = trees2$Height, xlab = "height")


}

\seealso{
\code{\link{rq}}, \code{\link{tsrq}}, \code{\link{tsrq2}}
}

\keyword{conditional quantiles}
\keyword{unconditional quantiles}
\keyword{location-scale-shape}
