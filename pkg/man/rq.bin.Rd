\name{rq.bin}
\alias{rqbin.fit}
\alias{rq.bin}

\title{
Binary Quantile Regression}

\description{
This function is used to fit a quantile regression model when the response is binary.
}

\usage{
rq.bin(formula, tau = 0.5, data, weights = NULL, contrasts = NULL,
	normalize = "last", control = NULL, fit = TRUE)
rqbin.fit(x, y, tau = 0.5, weights, control)
}

\arguments{
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
}
  \item{x}{
the design matrix.
}
  \item{y}{
the response variable.
}
  \item{tau}{
quantile to be estimated.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lqm is called.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector - not yet implemented.
}
  \item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
}

  \item{normalize}{
character specifying the type of normalization of the coefficients: if "last" (default), the last coefficient is set equal to 1; if "all", the vector of all coefficients except the intercept has norm equal to 1.
}

  \item{control}{
list of control parameters of the fitting process. See \code{\link{rqbinControl}}.
}

  \item{fit}{
logical flag. If \code{FALSE} the function returns a list of arguments to be passed to \code{\link{rqbin.fit}}.
}


}

\details{A binary quantile regression model is fitted as linear specification of the quantile function of a latent response variable (Manski 1975, 1985). The function \code{rqbin.fit} calls the Fortran routine \code{simann.f} implementing the simulated annealing algorithm of Goffe et al (1994) -- original code by William Goffe, modified by Gregory Kordas. Normalization is necessary for estimation to be possible. The normalization proposed by Horowitz (1992) assumes that there is a continuous regressor independent of the (latent) error and the corresponding regression coefficient is constrained to be equal to 1. Therefore, the user must ensure that the last term in \code{formula} or the last column in the matrix \code{x} corresponds to such regressor. If the argument \code{normalize = "all"}, then the normalization proposed by Manski (1975) is applied so that the norm of the vector with all the 'slopes' (i.e., excluding the intercept), is equal to 1.} 

\value{
a list of class \code{rq.bin} containing the following components

\item{coefficients}{a vector of coefficients. \code{coefficients} is a named matrix of coefficients when \code{tau} is a vector of values.}
\item{logLik}{the log--likelihood.}
\item{opt}{details on optimization.}
\item{call}{the matched call.}
\item{term.labels}{names for theta.}
\item{terms}{the terms object used.}
\item{nobs}{the number of observations.}
\item{edf}{the numer of parameters (minus 1 if normalize is \code{last}).}
\item{rdf}{the number of residual degrees of freedom.}
\item{tau}{the estimated quantile(s).}
\item{x}{the model matrix.}
\item{y}{the model response.}
\item{weights}{the weights used in the fitting process (a vector of 1's if \code{weights} = NULL).}
\item{levels}{factors levels.}
\item{control}{list of control parameters used for optimization (see \code{\link{rqbinControl}}).}
\item{normalize}{type of normalization.}
}


\references{
Goffe WL, Ferrier GD, Rogers J. Global optimization of statistical functions with simulated annealing. Journal of Econometrics 1994;60(1):65-99. Code retrieved from \url{http://EconPapers.repec.org/RePEc:wpa:wuwppr:9406001}.

Kordas G. Smoothed binary regression quantiles. Journal of Applied Econometrics 2006;21(3):387-407. Code retrieved from \url{http://qed.econ.queensu.ca/jae/2006-v21.3/kordas/}.

Horowitz JL. A Smoothed Maximum Score Estimator for the Binary Response Model. Econometrica 1992;60(3):505-531.

Manski CF. Maximum score estimation of the stochastic utility model of choice. Journal of Econometrics 1975;3(3):205-228.

Manski, CF. Semiparametric analysis of discrete response: Asymptotic properties of the maximum score estimator. Journal of Econometrics 1985;27(3):313-333. 
}

\author{
Marco Geraci
}



\keyword{conditional quantiles}
\keyword{discrete}
