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
\method{qlss}{default}(fun = "qnorm", p = 0.1, ...)
\method{qlss}{formula}(formula, p = 0.1, ..., predictLs = NULL)
\method{qlss}{numeric}(x, p = 0.1, ...)
}

\arguments{
  \item{p}{
a vector of probabilities.
}
  \item{fun}{
quantile function.
}
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted. The details of model specification are given under `Details'.
}
\item{x}{
a numeric vector.
}
  \item{\dots}{
other arguments for \code{fun} or \code{rq}.
}
  \item{predictLs}{
list of arguments for \code{predict.rq}.
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
Gilchrist W. Statistical modelling with quantile functions. Chapman and Hall/CRC; 2000.
}

\author{
Marco Geraci
}

\examples{

# Compute summary statistics of a normal distribution
qlss()

# Compute summary statistics of a t distribution with 3 df
qlss(fun = "qt", df = 3, p = 0.05)

# Compute summary statistics for a sample using a sequence of probabilities
x <- rnorm(1000)
qlss(x, p = c(0.1, 0.2, 0.3, 0.4))

# Compute summary statistics for Volume conditional on Height
trees2 <- trees[order(trees$Height),]
fit <- qlss(Volume ~ Height, data = trees2)
plot(fit, z = trees2$Height, xlab = "height")

# Use a quadratic model for Height
fit <- qlss(Volume ~ poly(Height,2), data = trees2, p = .1)
plot(fit, z = trees2$Height, xlab = "height")


}

\seealso{
\code{\link{rq}}
}

\keyword{ location-scale-shape }