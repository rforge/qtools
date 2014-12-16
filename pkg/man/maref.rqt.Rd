\name{maref.rqt}
\alias{maref}
\alias{maref.rq.counts}
\alias{maref.rqt}

\title{
Marginal Effects
}

\description{
This function computes marginal effects for \code{rqt} and \code{rq.counts} objects.
}

\usage{
maref(object, newdata, index = 2, index.extra = NULL, ...)
\method{maref}{rqt}(object, newdata, index = 2, index.extra = NULL, ...)
\method{maref}{rq.counts}(object, newdata, index = 2, index.extra = NULL, ...)
}

\arguments{
  \item{object}{
an \code{rqt} object.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{index}{
a numeric value to specify the position of the effect among the model's coefficients.
}
  \item{index.extra}{
optional numeric value(s) to specify the position of additional effects to be included in the computation.
}
  \item{\dots}{
not used.
}
}

\details{
Given the general model \emph{Q(h(response)|X) = Xb}, where \emph{Q} is the conditional quantile function, \emph{X} a design matrix with \emph{p} columns, and \emph{h} is a one- or two-parameter transformation with inverse \emph{hinv}, \code{maref} allows computing the marginal effect:

d\emph{Q(response|X)}/d\emph{x[j]}

where \emph{j} specifies the covariate in the design matrix with respect to which the marginal effect is to be computed and is given in the argument \code{index}. Since the model may contain interactions with \emph{x[j]}, additional terms in the matrix \emph{X} to be included in the computation are given in the argument \code{index.extra}. See the examples below.
}


\value{
a vector for single quantiles or a matrix for multiple quantiles of marginal effects.
}


\author{
Marco Geraci
}

\examples{
# Box-Cox quantile regression model (dataset trees from package 'datasets')
fit <- tsrq(Volume ~ Height, data = trees, tsf = "bc", tau = 0.9)

# Coefficients (transformed scale)
coef(fit)

# Design matrix
head(fit$x)

# Marginal effect of 'Height'
maref(fit, index = 2)

# Plot marginal effect over grid of values (for fixed girth)
nd <- data.frame(Height = seq(min(trees$Height), max(trees$Height), length = 100),
	Girth = rep(mean(trees$Girth), 100))
x <- maref(fit, newdata = nd, index = 2)
plot(nd$Height, x, xlab = "height", ylab = "marginal effect on volume")

# Include interaction between 'Height' and 'Girth'
fit <- tsrq(Volume ~ Height * Girth, data = trees, tsf = "bc", tau = 0.5)
head(fit$x)

# Marginal effect of 'Height'
maref(fit, index = 2, index.extra = 4)

# Quantile regression for counts (log transformation)
data(esterase)
fit <- rq.counts(Count ~ Esterase, tau = 0.25, data = esterase, M = 50)
maref(fit, index = 2)

}


\seealso{
\code{\link{tsrq}}
}

\keyword{coefficients}
