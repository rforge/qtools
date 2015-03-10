\name{midquantile}
\alias{midecdf}
\alias{midquantile.ci}
\alias{midquantile}

\title{
Mid-distribution Functions
}

\description{
Compute mid-cumulative probabilities and mid-quantiles
}

\usage{
midecdf(x, na.rm = FALSE)
midquantile(x, probs = 1:3/4, na.rm = FALSE)
midquantile.ci(x, probs = 1:3/4, level = 0.95)
}

\arguments{
  \item{x}{
numeric vector of observations used to estimate the mid-cumulative distribution or the mid-quantiles.
}
  \item{probs}{
numeric vector of probabilities with values in [0,1].
}
  \item{level}{
level of the confidence intervals for mid-quantiles.
}

  \item{na.rm}{
logical value indicating whether NA values should be stripped before the computation proceeds.
}

}

\details{
Mid-distribution functions
}

\value{
Mid-cumulative probabilities and mid-quantiles
}

\references{
Ma Y., Genton M., and Parzen E. Asymptotic properties of sample quantiles of discrete distributions. Annals of the Institute of Statistical Mathematics 2011;63(2):227-243

Parzen E. Quantile probability and statistical data modeling. Statistical Science 2004;19(4):652-62.
}

\author{
Marco Geraci
}


\keyword{ unconditional quantiles }
\keyword{ discrete }
