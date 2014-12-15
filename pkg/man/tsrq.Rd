\name{tsrq}
\alias{tsrq2}
\alias{rcrq}
\alias{nlrq2}
\alias{tsrq}

\title{
Quantile Regression Tranformation Models
}
\description{
These functions are used to fit quantile regression transformation models.
}
\usage{
tsrq(formula, tsf = "mcjI", symm = TRUE, bounded = FALSE, lambda = NULL,
	conditional = FALSE, tau = 0.5, data, subset, weights, na.action,
	method = "fn", ...)
tsrq2(formula, bounded = FALSE, lambda = NULL, delta = NULL, tau = 0.5,
	conditional = FALSE, data, subset, weights, na.action, method = "fn", ...)
rcrq(formula, tsf = "mcjI", symm = TRUE, bounded = FALSE, lambda = NULL,
	tau = 0.5, data, subset, weights, na.action, method = "fn", ...)
nlrq2(formula, par = NULL, bounded = FALSE, tau = 0.5, data,
	subset, weights, na.action, ...)
}

\arguments{
  \item{formula}{
an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under `Details'.
}
  \item{tsf}{
transformation to be used. Possible options are \code{mcjI} for Proposal I transformation models (default) of Geraci and Jones, \code{bc} for Box-Cox and \code{ao} for Aranda-Ordaz transformation models (Dehbi, Cortina-Borja, Geraci).
}
  \item{symm}{
logical flag. If \code{TRUE} (default) a symmetric transformation is used.
}
  \item{bounded}{
logical flag. If \code{TRUE} the response is assumed to be doubly bounded on [a,b]. If \code{FALSE} the response is assumed to be singly bounded (ie, strictly positive).
}
  \item{lambda, delta}{
values of transformation parameters for grid search.
}
  \item{conditional}{
logical flag. If \code{TRUE}, the transformation parameter is assumed to be known and this must be provided via the arguments \code{lambda}, \code{delta} in vectors of the same length as \code{tau}.
}
  \item{par}{
vector of length \code{p + 2} of initial values for the parameters to be optimized over. The first \code{p} values are for the regression coefficients while the last 2 are for the transformation parameters \code{lambda} and \code{delta} in \code{\link{mcjII}}. These initial values are passed to \code{\link{optim}}.
}
  \item{tau}{
the quantile(s) to be estimated. See \code{\link{rq}}.
}
  \item{data}{
a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process. Should be NULL or a numeric vector. 
}
  \item{na.action}{
a function which indicates what should happen when the data contain \code{NA}s.
}
  \item{method}{
fitting algorithm for \code{\link{rq}} (default is Frisch-Newton interior point method "\code{fn}").
}
  \item{\dots}{
other arguments in \code{\link{summary.rq}}).
}
}

\details{
These functions implement quantile regression transformation models as discussed by Geraci and Jones (see references). The general model is assumed to be \emph{Q(h(response)) = xb}, where \emph{Q} is the conditional quantile function and \emph{h} is a one- or two-parameter transformation. A typical model specified in \code{formula} has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for the quantile of the transformed response. The \code{response}, which is singly or doubly bounded, i.e. \code{response > 0} or \code{0 <= response <= 1} respectively, undergoes the transformation specified in \code{tsf}. If the response is bounded in the generic [a,b] interval, the latter is automatically mapped to [0,1] and no further action is required. If, however, the response is singly bounded and contains negative values, it is left to the user to offset the response or the code will produce an error. The functions \code{tsrq} and \code{tsrq2} use a two-stage (TS) estimator (Fitzenberger et al, 2010) for, respectively, one- and two-parameter transformations. The function \code{rcrq} (one-parameter tranformations) is based on the residual cusum process estimator proposed by Mu and He (2007), while the function \code{nlrq2} (two-parameter tranformations) is based on Nelder-Mean optimization (Geraci and Jones).
}

\value{
\code{tsrq}, \code{tsrq2}, \code{rcrq}, \code{nlrq2} return an object of \code{\link{class}} \code{rqt}. This is a list that contains as typical components:

\item{}{the first \code{nt = length(tau)} elements of the list store the results from fitting linear quantile models on the tranformed scale of the response.}
\item{call}{the matched call.}
\item{method}{the fitting algorithm for \code{rq} or \code{optim}.}
\item{y}{the response -- untransformed scale.}
\item{theta}{if \code{bounded = TRUE}, the response mapped to the unit interval.}
\item{x}{the model matrix.}
\item{weights}{the weights used in the fitting process (a vector of 1's if \code{weights} is missing or \code{NULL}).}
\item{tau}{the order of the estimated quantile(s).}
\item{lambda}{the estimated parameter lambda.}
\item{eta}{the estimated parameters lambda and delta in the two-parameter Proposal II tranformation.}
\item{lambda.grid}{grid of lambda values used for estimation.}
\item{delta.grid}{grid of delta values used for estimation.}
\item{tsf}{tranformation used (see also \code{attributes(tsf)}).}
\item{objective}{values of the objective function minimised over the tranformation parameter(s). This is an array of dimension \code{c(nl,nt)} or \code{c(nl,nd,nt)}, where \code{nl = length(lambda.grid)}, \code{nd = length(delta.grid)} and \code{nt = length(tau)}.}
\item{optimum}{value of the objective function at solution.}
\item{coefficients}{quantile regression coefficients -- transformed scale.}
\item{fitted.values}{fitted values.}
\item{rejected}{proportion of inadmissible observations (Fitzenberger et al, 2010).}
\item{terms}{the \code{\link{terms}} used.}
\item{term.labels}{names of coefficients.}
\item{rdf}{residual degrees of freedom.}
}

\references{
Aranda-Ordaz FJ. On two families of transformations to additivity for binary response data. Biometrika 1981;68(2):357-363.

Box GEP, Cox DR. An analysis of transformations. Journal of the Royal Statistical Society Series B-Statistical Methodology 1964;26(2):211-252.

Dehbi H-M, Cortina-Borja M, Geraci M. Aranda-Ordaz quantile regression for student performance assessment. Unpublished Manuscript.

Fitzenberger B, Wilke R, Zhang X. Implementing Box-Cox quantile regression. Econometric Reviews 2010;29(2):158-181.

Geraci M and Jones MC. Improved transformation-based quantile regression. Canadian Journal of Statistics 2014 (in press).

Jones MC. Connecting distributions with power tails on the real line, the half line and the interval. International Statistical Review 2007;75(1):58-69.

Koenker R (2013). quantreg: Quantile Regression. R package version 5.05. URL \url{http://CRAN.R-project.org/package=quantreg}.

Mu YM, He XM. Power transformation toward a linear regression quantile. Journal of the American Statistical Association 2007;102(477):269-279.
}


\author{
Marco Geraci
}

\examples{
data(Chemistry)

Chemistry$gcse_gr <- cut(Chemistry$gcse, c(0,seq(4,8,by=0.5)))
with(Chemistry, plot(score ~ gcse_gr, xlab = "GCSE score", ylab = "A-level Chemistry score"))

\dontrun{

# The dataset has > 31000 observations and computation can be slow
set.seed(178)
chemsub <- Chemistry[sample(1:nrow(Chemistry), 2000), ]

# Fit symmetric Aranda-Ordaz quantile 0.9
tsrq(score ~ gcse, tsf = "ao", symm = TRUE, lambda = seq(0,2,by=0.01), tau = 0.9, data = chemsub)

# Fit symmetric Proposal I quantile 0.9
tsrq(score ~ gcse, tsf = "mcjI", symm = TRUE, bounded = TRUE, lambda = seq(0,2,by=0.01),
	tau = 0.9, data = chemsub)

# Fit Proposal II quantile 0.9 (Nelder-Mead)
nlrq2(score ~ gcse, bounded = TRUE, tau = 0.9, data = chemsub)

# Fit Proposal II quantile 0.9 (grid search)
# This is slower than nlrq2 but more stable numerically
tsrq2(score ~ gcse, bounded = TRUE, lambda = seq(0, 2, by = 0.1),
	delta = seq(0, 2, by = 0.1), tau = 0.9, data = chemsub)

}

}

\seealso{
\code{\link{predict.rqt}}, \code{\link{summary.rqt}}
}

\keyword{ transformations }
\keyword{conditional quantiles}