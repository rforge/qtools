\name{rq.counts}
\alias{rq.counts}

\title{
Quantile Regression for Counts}
\description{
This function is used to fit a quantile regression model when the response is a count variable.
}
\usage{
rq.counts(formula, tau = 0.5, data, tsf = "bc", symm = TRUE,
	lambda = 0,	weights = NULL, offset = NULL, contrasts = NULL,
	M = 50, zeta = 1e-05, B = 0.999, cn = NULL, alpha = 0.05,
	method = "fn")
}

\arguments{
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
}
  \item{tau}{
quantile to be estimated.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lqm is called.
}
  \item{tsf}{
transformation to be used. Possible options are \code{mcjI} for Proposal I and \code{bc} for Box-Cox. The latter is the default, with \code{lambda = 0}, i.e. the logarithmic transformation.
}
  \item{symm}{
logical flag. If \code{TRUE} and \code{tsf = 'mcj'}, the symmetric Proposal I transformation is used.
}
  \item{lambda}{
a numerical value for the transformation parameter. This is provided by the user or set to zero if not specified.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process.
}

  \item{offset}{
an optional offset to be included in the model frame.  
}

  \item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
}
  \item{M}{
number of dithered samples.
}
  \item{zeta}{
small constant (see References).
}

  \item{B}{
right boundary for uniform random noise U[0,B] to be added to the response variable (see References).
}

 \item{cn}{
small constant to be passed to \code{Fn} (see Theorem 3, Machado and Santos Silva).
}

 \item{alpha}{
significance level.
}
  \item{method}{
estimation method for the fitting process. See \code{\link{rq}}.
}

}

\details{A linear quantile regression model is fitted to the log--transformed response. Additional tranformation functions will be implemented. The notation used here follows closely that of Machado and Santos Silva (2005). This function is based on routines from package \code{quantreg} (Koenker, 2013). See also \code{lqm.counts} from package \code{lqmm} (Geraci, 2014) for Laplace gradient estimation.} 

\value{
a list of class \code{rq.counts} containing the following components

\item{call}{the matched call.}
\item{method}{the fitting algorithm for \code{rq}.}
\item{x}{the model matrix.}
\item{y}{the model response.}
\item{tau}{the order of the estimated quantile(s).}
\item{tsf}{tranformation used (see also \code{attributes(tsf)}).}
\item{coefficients}{regression quantile (on the log--scale).}
\item{fitted.values}{fitted values (on the response scale).}
\item{tTable}{coefficients, standard errors, etc.}
\item{offset}{offset.}
\item{M}{specified number of dithered samples for standard error estimation.}
\item{Mn}{actual number of dithered samples used for standard error estimation that gave an invertible D matrix (Machado and Santos Silva, 2005).}
\item{InitialPar}{starting values for coefficients.}
\item{terms}{the terms object used.}
\item{term.labels}{names of coefficients.}
\item{rdf}{the number of residual degrees of freedom.}
}


\references{

Geraci M and Jones MC. Improved transformation-based quantile regression. Canadian Journal of Statistics 2015;43(1):118-132.

Koenker R (2013). quantreg: Quantile Regression. R package version 5.05. URL \url{http://CRAN.R-project.org/package=quantreg}.

Machado JAF and Santos Silva JMC (2005). Quantiles for counts. Journal of the American Statistical Association, 100(472), 1226--1237.
}

\author{
Marco Geraci
}



\examples{

# Esterase data
data(esterase)

# Fit quantiles 0.25 and 0.75
fit1 <- rq.counts(Count ~ Esterase, tau = 0.25, data = esterase, M = 50)
coef(fit1)
fit2 <- rq.counts(Count ~ Esterase, tau = 0.75, data = esterase, M = 50)
coef(fit2)

# Plot
with(esterase, plot(Count ~ Esterase))
lines(esterase$Esterase, fit1$fitted.values, col = "blue")
lines(esterase$Esterase, fit2$fitted.values, col = "red")
legend(8, 1000, lty = c(1,1), col = c("blue", "red"), legend = c("tau = 0.25","tau = 0.75"))

}
\keyword{conditional quantiles}
\keyword{discrete}
