\name{GOFTest}
\alias{GOFTest}

\title{
Goodness-of-Fit Test for Quantile Regression Models
}
\description{
This function calculate a goodness-of-fit test for quantile regression models based on a logistic likelihood ratio test.
}
\usage{
GOFTest(object, covariates = FALSE, nb = 10)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "rq","rqs", or "rqt".
}
  \item{covariates}{
logical flag. If \code{TRUE}, the covariates used in the fitted quantile regression model are included to assess the goodness of fit.
}
  \item{nb}{
number of intervals. This is an integer greater than or equal to 2.
}
}

\details{
This function computes a likelihood ratio test for the proportion of negative residuals \code{z < 0} from a quantile regression model (Ortaglia et al, 2015). A logistic regression model is fitted using the indicator function \code{I(z < 0)} as outcome and the quantile regression fitted values as predictor. A significant test provides evidence towards a poor fit.
}


\references{
Ortaglia, Geraci M, Bottai M. Simple and fast lack-of-fit tests for conditional quantile estimation. Unpublished manuscript 2015.

Redden DT, Fernandez JR and Allison DB. A simple significance test for quantile regression. Statistics in Medicine 2004;23: 2587-2597.
}

\author{
Marco Geraci
}

\examples{

data(barro, package = "quantreg")
fit <- quantreg::rq(y.net ~ lgdp2 + fse2 + gedy2 + Iy2 + gcony2, data = barro, tau = c(.1,.5,.9))
GOFTest(fit, covariates = FALSE, nb = 10)
# Some convergence problems
GOFTest(fit, covariates = TRUE, nb = 10)
# Reduce number of intervals
GOFTest(fit, covariates = TRUE, nb = 5)
}