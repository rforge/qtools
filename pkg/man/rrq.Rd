\name{rrq}
\alias{rrq.fit}
\alias{rrq}

\title{
Restricted Regression Quantiles
}
\description{
This function fits a restricted quantile regression model to avoid crossing of quantile curves.
}

\usage{
rrq(formula, tau, data, method = "fn", model = TRUE, ...)
rrq.fit(x, y, tau, method = "fn")
}

\arguments{
  \item{formula}{
a formula object, with the response on the left of a \code{~} operator, and the terms, separated by \code{+} operators, on the right.
}
  \item{x}{
the design matrix.
}
  \item{y}{
the response variable.
}
  \item{tau}{
the quantile(s) to be estimated.
}
  \item{data}{
a data frame in which to interpret the variables named in the formula.
}
  \item{method}{
the algorithm used to compute the fit (see \code{\link{rq}}).
}
  \item{model}{
if \code{TRUE} then the model frame is returned. This is essential if one wants to call summary subsequently.
}
  \item{\dots}{
not used.
}

}

\references{
He X. Quantile curves without crossing. The American Statistician 1997;51(2):186-192.

Koenker R (2013). quantreg: Quantile Regression. R package version 5.05. URL \url{http://CRAN.R-project.org/package=quantreg}.
}

\author{
Marco Geraci
}

\examples{
data(esterase)

# Fit standard quantile regression
fit <- quantreg::rq(Count ~ Esterase, data = esterase, tau = c(.1,.25,.5,.75,.9))
yhat <- fit$fitted.values

# Fit restricted quantile regression
fitr <- rrq(Count ~ Esterase, data = esterase, tau = c(.1,.25,.5,.75,.9))
yhat2 <- predict(fitr)


# Plot results
par(mfrow = c(1, 2))

# Plot regression quantiles
with(esterase, plot(Count ~ Esterase, pch = 16, cex = .8))
apply(yhat, 2, function(y,x) lines(x,y,lwd = 1.5), x = esterase$Esterase)

# Plot restricted regression quantiles
with(esterase, plot(Count ~ Esterase, pch = 16, cex = .8))
apply(yhat2, 2, function(y,x) lines(x,y,lwd = 1.5), x = esterase$Esterase)

}

\keyword{ quantile crossing }
\keyword{ conditional quantiles }

