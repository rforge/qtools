\name{rqbinControl}
\alias{rqbinControl}

\title{
Control parameters for rqbin estimation
}

\description{
A list of parameters for controlling the fitting process based on simulated annealing.
}

\usage{
rqbinControl(theta = NULL, lower = NULL, upper = NULL,
maximise = TRUE, rt = 0.15, tol = 1e-06, ns = 10, nt = 20, neps = 4,
maxiter = 1e+05, sl = NULL, vm = NULL, seed1 = 1, seed2 = 2, temp = 10,
sgn = 1)
}

\arguments{
  \item{theta}{
starting values for the regression coefficients.
}
  \item{lower}{
lower bound for the regression coefficients.
}
  \item{upper}{
upper bound for the regression coefficients.
}
  \item{maximise}{
logical flag. If \code{TRUE} the algorithm maximises the objective function.
}
  \item{rt}{
parameter.
}
  \item{tol}{
tolerance expressed as absolute change of the log-likelihood.
}
  \item{ns}{
parameter.
}
  \item{nt}{
parameter.
}
  \item{neps}{
parameter.
}
  \item{maxiter}{
maximum number of iterations.
}
  \item{sl}{
parameter.
}
  \item{vm}{
parameter.
}
  \item{seed1}{
random numbers seed.
}
  \item{seed2}{
random numbers seed.
}
  \item{temp}{
temperature of the algorithm.
}
  \item{sgn}{
sign of the objective function.
}
}

\details{
For more details on specific parameters, see Goffe et al (1994).
}

\references{
Goffe WL, Ferrier GD, Rogers J. Global optimization of statistical functions with simulated annealing. Journal of Econometrics 1994;60(1):65-99. Code retrieved from \url{http://EconPapers.repec.org/RePEc:wpa:wuwppr:9406001}.
}

\author{
Marco Geraci
}



\seealso{
\code{\link{rq.bin}}
}

\keyword{fitting}
\keyword{control}
