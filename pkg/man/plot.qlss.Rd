\name{plot.qlss}
\alias{plot.qlss}

\title{
Quantile-based Summary Statistics for Location, Scale and Shape
}

\description{
This function plots location, scale and shape of a conditional distribution.
}

\usage{
\method{plot}{qlss}(x, z, which = 1, ci = FALSE, level = 0.95, type = "l", ...)
}

\arguments{
  \item{x}{
an object of class \code{qlss} as returned by \code{\link{qlss.formula}}.
}
  \item{z}{
numeric vector.
}
  \item{which}{
when \code{probs} in \code{qlss} is a vector, the argument \code{which} specifies which element of the vector should be used for plotting.
}
 \item{ci}{
logical flag. If \code{TRUE}, bootstrapped confidence intervals for the predictions are plotted.
}
 \item{level}{
if \code{ci = TRUE}, point-wise confidence intervals of a given \code{level} are plotted.
}
  \item{type}{
1-character string giving the type of plot desired. See \code{\link{plot.default}}.
}
\item{\dots}{
other arguments for \code{\link{plot.default}}.
}

}


\author{
Marco Geraci
}

\examples{

trees2 <- trees[order(trees$Height),]
fit <- qlss(Volume ~ Height, data = trees2, probs = c(.05, .1))
# Plot the results for probs = 0.1
plot(fit, z = trees2$Height, which = 2, xlab = "height")

}

\seealso{
\code{\link{qlss}}
}

\keyword{ plot }