\name{boot.ci}
\alias{boot.ci}

\title{
Bootstrap Confidence Intervals
}

\description{
This function is used to calculate bootstrap percentile confidence intervals.
}

\usage{
boot.ci(object, conf = 0.95)
}

\arguments{
  \item{object}{
an object provided by \code{boot.tsrq}, \code{boot.tsrq2} or \code{boot.nlrq2}
}
  \item{conf}{
confidence level.
}
}

\details{
This is a modification of the function \code{boot.ci} from package \code{boot} (Canty and Ripley, 2013).
}

\author{
Marco Geraci
}

\references{

Angelo Canty and Brian Ripley (2013). boot: Bootstrap R (S-Plus) Functions. R package version 1.3-9.

}

\keyword{ bootstrap }
