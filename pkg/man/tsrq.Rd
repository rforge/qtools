\name{tsrq}
\alias{tsrq2}
\alias{rcrq}
\alias{nlrq2}
\alias{tsrq}

\title{
Quantile Regression Transformation Models
}
\description{
These functions are used to fit quantile regression transformation models.
}
\usage{
tsrq(formula, tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = NULL,
	conditional = FALSE, tau = 0.5, data, subset, weights, na.action,
	method = "fn", ...)
tsrq2(formula, dbounded = FALSE, lambda = NULL, delta = NULL, conditional = FALSE,
	tau = 0.5, data, subset, weights, na.action, method = "fn", ...)
rcrq(formula, tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = NULL,
	tau = 0.5, data, subset, weights, na.action, method = "fn", ...)
nlrq2(formula, par = NULL, dbounded = FALSE, tau = 0.5, data,
	subset, weights, na.action, ...)
}

\arguments{
  \item{formula}{
an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under `Details'.
}
  \item{tsf}{
transformation to be used. Possible options are \code{mcjI} for Proposal I transformation models (default), \code{bc} for Box-Cox and \code{ao} for Aranda-Ordaz transformation models.
}
  \item{symm}{
logical flag. If \code{TRUE} (default) a symmetric transformation is used.
}
  \item{dbounded}{
logical flag. If \code{TRUE} the response is assumed to be doubly bounded on [a,b]. If \code{FALSE} (default) the response is assumed to be singly bounded (ie, strictly positive).
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
These functions implement quantile regression transformation models as discussed by Geraci and Jones (see references). The general model is assumed to be \emph{Q(h(response)) = xb}, where \emph{Q} is the conditional quantile function and \emph{h} is a one- or two-parameter transformation. A typical model specified in \code{formula} has the form \code{response ~ terms} where \code{response} is the (numeric) response vector and \code{terms} is a series of terms which specifies a linear predictor for the quantile of the transformed response. The \code{response}, which is singly or doubly bounded, i.e. \code{response > 0} or \code{0 <= response <= 1} respectively, undergoes the transformation specified in \code{tsf}. If the response is bounded in the generic [a,b] interval, the latter is automatically mapped to [0,1] and no further action is required. If, however, the response is singly bounded and contains negative values, it is left to the user to offset the response or the code will produce an error. The functions \code{tsrq} and \code{tsrq2} use a two-stage (TS) estimator (Fitzenberger et al, 2010) for, respectively, one- and two-parameter transformations. The function \code{rcrq} (one-parameter tranformations) is based on the residual cusum process estimator proposed by Mu and He (2007), while the function \code{nlrq2} (two-parameter tranformations) is based on Nelder-Mead optimization (Geraci and Jones).
}

\value{
\code{tsrq}, \code{tsrq2}, \code{rcrq}, \code{nlrq2} return an object of \code{\link{class}} \code{rqt}. This is a list that contains as typical components:

\item{}{the first \code{nt = length(tau)} elements of the list store the results from fitting linear quantile models on the tranformed scale of the response.}
\item{call}{the matched call.}
\item{method}{the fitting algorithm for \code{rq} or \code{optim}.}
\item{y}{the response -- untransformed scale.}
\item{theta}{if \code{dbounded = TRUE}, the response mapped to the unit interval.}
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

Dehbi H-M, Cortina-Borja M, Geraci M. Aranda-Ordaz quantile regression for student performance assessment. Journal of Applied Statistics 2015. \url{http://dx.doi.org/10.1080/02664763.2015.1025724}

Fitzenberger B, Wilke R, Zhang X. Implementing Box-Cox quantile regression. Econometric Reviews 2010;29(2):158-181.

Geraci M and Jones MC. Improved transformation-based quantile regression. Canadian Journal of Statistics 2015;43(1):118-132.

Jones MC. Connecting distributions with power tails on the real line, the half line and the interval. International Statistical Review 2007;75(1):58-69.

Koenker R (2013). quantreg: Quantile Regression. R package version 5.05. URL \url{http://CRAN.R-project.org/package=quantreg}.

Mu YM, He XM. Power transformation toward a linear regression quantile. Journal of the American Statistical Association 2007;102(477):269-279.
}


\author{
Marco Geraci
}

\examples{

###########################################################
## Example 1 - singly bounded (from Geraci and Jones, 2014)

data(trees)

\dontrun{
require(MASS)

dx <- 0.01

lambda0 <- boxcox(Volume ~ log(Height), data = trees,
	lambda = seq(-0.9, 0.5, by = dx))
lambda0 <- lambda0$x[which.max(lambda0$y)]
trees$z <- bc(trees$Volume,lambda0)
trees$y <- trees$Volume
trees$x <- log(trees$Height)
trees$x <- trees$x - mean(log(trees$Height))

fit.lm <- lm(z ~ x, data = trees)
newd <- data.frame(x = log(seq(min(trees$Height),
	max(trees$Height), by = 0.1)))
newd$x <- newd$x - mean(log(trees$Height))
ylm <- invbc(predict(fit.lm, newdata = newd), lambda0)

lambdas <- list(bc = seq(-10, 10, by=dx),
	mcjIs = seq(0,10,by = dx), mcjIa = seq(0,20,by = dx))

taus <- 1:3/4
fit0 <- tsrq(y ~ x, tsf = "bc", symm = FALSE, data = trees,
	lambda = lambdas$bc, tau = taus)
fit1 <- tsrq(y ~ x, tsf = "mcjI", symm = TRUE, dbounded = FALSE,
	data = trees, lambda = lambdas$mcjIs, tau = taus)
fit2 <- tsrq(y ~ x, tsf = "mcjI", symm = FALSE, dbounded = FALSE,
	data = trees, lambda = lambdas$mcjIa, tau = taus)


par(mfrow = c(1,3), mar = c(7.1, 7.1, 5.1, 2.1), mgp = c(5, 2, 0))

cx.lab <- 2.5
cx.ax <- 2
lw <- 2
cx <- 2
xb <- "log(Height)"
yb <- "Volume"
xl <- range(trees$x)
yl <- c(5,80)

yhat <- predict(fit0, newdata = newd)
plot(y ~ x,data = trees, xlim = xl, ylim = yl, main = "Box-Cox",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	cex = cx, xlab = xb, ylab = yb)
lines(newd$x, yhat[,1], lwd = lw)
lines(newd$x, yhat[,2], lwd = lw)
lines(newd$x, yhat[,3], lwd = lw)
lines(newd$x, ylm, lwd = lw, lty = 2)

yhat <- predict(fit1, newdata = newd)
plot(y ~ x,data = trees, xlim = xl, ylim = yl, main = "Proposal I (symmetric)",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	cex = cx, xlab = xb, ylab = yb)
lines(newd$x, yhat[,1], lwd = lw)
lines(newd$x, yhat[,2], lwd = lw)
lines(newd$x, yhat[,3], lwd = lw)
lines(newd$x, ylm, lwd = lw, lty = 2)

yhat <- predict(fit2, newdata = newd)
plot(y ~ x, data = trees, xlim = xl, ylim = yl, main = "Proposal I (asymmetric)",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	cex = cx, xlab = xb, ylab = yb)
lines(newd$x, yhat[,1], lwd = lw)
lines(newd$x, yhat[,2], lwd = lw)
lines(newd$x, yhat[,3], lwd = lw)
lines(newd$x, ylm, lwd = lw, lty = 2)
}

###########################################################
## Example 2 - doubly bounded

data(Chemistry)

Chemistry$gcse_gr <- cut(Chemistry$gcse, c(0,seq(4,8,by=0.5)))
with(Chemistry, plot(score ~ gcse_gr, xlab = "GCSE score",
	ylab = "A-level Chemistry score"))

\dontrun{

# The dataset has > 31000 observations and computation can be slow
set.seed(178)
chemsub <- Chemistry[sample(1:nrow(Chemistry), 2000), ]

# Fit symmetric Aranda-Ordaz quantile 0.9
tsrq(score ~ gcse, tsf = "ao", symm = TRUE, lambda = seq(0,2,by=0.01),
	tau = 0.9, data = chemsub)

# Fit symmetric Proposal I quantile 0.9
tsrq(score ~ gcse, tsf = "mcjI", symm = TRUE, dbounded = TRUE,
	lambda = seq(0,2,by=0.01), tau = 0.9, data = chemsub)

# Fit Proposal II quantile 0.9 (Nelder-Mead)
nlrq2(score ~ gcse, dbounded = TRUE, tau = 0.9, data = chemsub)

# Fit Proposal II quantile 0.9 (grid search)
# This is slower than nlrq2 but more stable numerically
tsrq2(score ~ gcse, dbounded = TRUE, lambda = seq(0, 2, by = 0.1),
	delta = seq(0, 2, by = 0.1), tau = 0.9, data = chemsub)

}

###########################################################
## Example 3 - doubly bounded

data(labor)

new <- labor
new$y <- new$pain
new$x <- (new$time-30)/30
new$x_gr <- as.factor(new$x)

par(mfrow = c(2,2), mai = c(1.3,1.3,0.5,0.5), mgp = c(5, 2, 0))

cx.lab <- 2.5
cx.ax <- 2.5
cx <- 2.5
yl <- c(0,0.06)

hist(new$y[new$treatment == 1], xlab = "Pain score", main = "Treated",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab, axes = TRUE,
	freq = FALSE, ylim = yl)

plot(y ~ x_gr, new, subset = new$treatment == 1, xlab = "Time (min)",
	ylab = "Pain score", cex.lab = cx.lab, cex.axis = cx.ax,
	cex.main = cx.lab, axes = FALSE, range = 0)
axis(1, at = 1:6, labels = c(0:5)*30 + 30, cex.axis = 2)
axis(2, cex.axis = cx.lab)
box()

hist(new$y[new$treatment == 0], xlab = "Pain score", main = "Placebo",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	axes = TRUE, freq = FALSE, ylim = yl)

plot(y ~ x_gr, new, subset = new$treatment == 0, xlab = "Time (min)",
	ylab = "Pain score", cex.lab = cx.lab, cex.axis = cx.ax,
	cex.main = cx.lab, axes = FALSE, range = 0)
axis(1, at = 1:6, labels = (0:5)*30 + 30, cex.axis = 2)
axis(2, cex.axis = cx.lab)
box()

#

\dontrun{

taus <- c(1:3/4)
ls <- seq(0,3.5,by=0.1)

fit.aos <- tsrq(y ~ x*treatment, tsf = "ao", symm = TRUE, dbounded = TRUE,
	data = new, tau = taus, lambda = ls)
fit.aoa <- tsrq(y ~ x*treatment, tsf = "ao", symm = FALSE, dbounded = TRUE,
	data = new, tau = taus, lambda = ls)
fit.mcjs <- tsrq(y ~ x*treatment, tsf = "mcjI", symm = TRUE, dbounded = TRUE,
	data = new, tau = taus, lambda = ls)
fit.mcja <- tsrq(y ~ x*treatment, tsf = "mcjI", symm = FALSE, dbounded = TRUE,
	data = new, tau = taus, lambda = ls)
fit.mcj2 <- tsrq2(y ~ x*treatment, dbounded = TRUE, data = new, tau = taus,
	lambda = seq(0,2,by=0.1), delta = seq(0,1.5,by=0.3))
fit.nlrq <- nlrq2(y ~ x*treatment, par = coef(fit.mcj2, all = TRUE)[,1],
	dbounded = TRUE, data = new, tau = taus)

sel <- 0
selx <- if(sel == 1) "a" else "b"
x <- new$x
nd <- data.frame(x = seq(min(x),max(x),length=200), treatment = sel)
xx <- nd$x+1

par(mfrow = c(2,2), mar = c(7.1, 7.1, 5.1, 2.1), mgp = c(5, 2, 0))
cx.lab <- 3
cx.ax <- 2
lw <- 3
cx <- 1

fit <- fit.aos
yhat <- predict(fit, newdata = nd)

plot(y ~ x_gr, new, subset = new$treatment == sel, xlab = "",
	ylab = "Pain score", cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	axes = FALSE, main = "Aranda-Ordaz (s)", range = 0, col = grey(4/5))
apply(yhat, 2, function(y,x) lines(x, y, lwd = lw), x = xx)
axis(1, at = 1:6, labels = (0:5)*30 + 30, cex.axis = cx.lab)
axis(2, at = c(0, 25, 50, 75, 100), cex.axis = cx.lab)
box()

fit <- fit.aoa
yhat <- predict(fit, newdata = nd)

plot(y ~ x_gr, new, subset = new$treatment == sel, xlab = "", ylab = "",
	cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab, axes = FALSE,
	main = "Aranda-Ordaz (a)", range = 0, col = grey(4/5))
apply(yhat, 2, function(y,x) lines(x, y, lwd = lw), x = xx)
axis(1, at = 1:6, labels = (0:5)*30 + 30, cex.axis = cx.lab)
axis(2, at = c(0, 25, 50, 75, 100), cex.axis = cx.lab)
box()

fit <- fit.mcjs
yhat <- predict(fit, newdata = nd)

plot(y ~ x_gr, new, subset = new$treatment == sel, xlab = "Time (min)",
	ylab = "Pain score", cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	axes = FALSE, main = "Proposal I (s)", range = 0, col = grey(4/5))
apply(yhat, 2, function(y,x) lines(x, y, lwd = lw), x = xx)
axis(1, at = 1:6, labels = (0:5)*30 + 30, cex.axis = cx.lab)
axis(2, at = c(0, 25, 50, 75, 100), cex.axis = cx.lab)
box()

fit <- fit.mcj2
yhat <- predict(fit, newdata = nd)

plot(y ~ x_gr, new, subset = new$treatment == sel, xlab = "Time (min)",
	ylab = "", cex.lab = cx.lab, cex.axis = cx.ax, cex.main = cx.lab,
	axes = FALSE, main = "Proposal II", range = 0, col = grey(4/5))
apply(yhat, 2, function(y,x) lines(x, y, lwd = lw), x = xx)
axis(1, at = 1:6, labels = (0:5)*30 + 30, cex.axis = cx.lab)
axis(2, at = c(0, 25, 50, 75, 100), cex.axis = cx.lab)
box()
}

}

\seealso{
\code{\link{predict.rqt}}, \code{\link{summary.rqt}}
}

\keyword{transformations}
\keyword{conditional quantiles}