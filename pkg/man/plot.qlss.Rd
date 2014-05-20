\name{plot.qlss}
\alias{plot.qlss}

\title{
Quantile-based Summary Statistics for Location, Scale and Shape
}

\description{
This function plots location, scale and shape of a conditional distribution.
}

\usage{
\method{plot}{qlss}(x, z, type = "l", ...)
}

\arguments{
  \item{x}{
an object of class \code{qlss} as returned by \code{\link{qlss.formula}}.
}
  \item{z}{
numeric vector.
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
fit <- qlss(Volume ~ Height, data = trees2)
plot(fit, z = trees2$Height, xlab = "height")

}

\seealso{
\code{\link{qlss}}
}

\keyword{ plot }