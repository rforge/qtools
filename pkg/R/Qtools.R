###            Utilities for quantiles
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License or
#  any later version.
#
#  This program is distributed without any warranty,
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

".onAttach" <- function(lib, pkg) {
    if(interactive() || getOption("verbose"))
	packageStartupMessage(sprintf("Package %s (%s) loaded. Type citation(\"%s\") on how to cite this package\n", pkg,
		packageDescription(pkg)$Version, pkg))
}


######################################################################
### Sample quantiles
######################################################################

# mid-CDF

midecdf <- function(x, na.rm = FALSE){

if(!na.rm && any(is.na(x))) 
	return(NA)
if(na.rm && any(ii <- is.na(x))) 
	x <- x[!ii]
if(is.unsorted(x))
	x <- sort(x)

n <- length(x)
if (n < 1) 
	stop("'x' must have 1 or more non-missing values")
xo <- unique(x)

pmf <- as.numeric(table(x)/n)
val <- list()
val$x <- xo
val$y <- ecdf(x)(xo) - 0.5*pmf
return(val)

}

midquantile <- function(x, probs = NULL, na.rm = FALSE){

if(!na.rm && any(is.na(x))) 
	return(NA)
if(na.rm && any(ii <- is.na(x))) 
	x <- x[!ii]
if(!is.null(probs) && any(c(probs < 0, probs > 1)))
	stop("the probability must be between 0 and 1")
	
Gn <- midecdf(x)
qval <- approxfun(Gn$y, Gn$x, method = "linear", rule = 2)

if(is.null(probs))
	return(qval) else return(qval(probs))
}

######################################################################
### Q-based statistics
######################################################################

qlss <- function(...) UseMethod("qlss")

qlss.numeric <- function(x, p = 0.1, ...){

if(any(!p > 0) | any(!p < 0.5)) stop("Quantile index out of range: p must be > 0 and < 0.5")
if(any(duplicated(p))) p <- unique(p)
nq <- length(p)

vec3 <- as.numeric(do.call(what = quantile, args = list(x = x, probs = 1:3/4, ...)))
Me <- vec3[2]
IQR <- vec3[3] - vec3[1]

PD <- Ap <- Tp <- rep(NA, nq)
for(i in 1:nq){
	vecp <- as.numeric(do.call(what = quantile, args = list(x = x, probs = p[i], ...)))
	vecq <- as.numeric(do.call(what = quantile, args = list(x = x, probs = 1 - p[i], ...)))

	PD[i] <- vecq - vecp
	Ap[i] <- (vecq - 2*Me + vecp)/PD[i]
	Tp[i] <- PD[i]/IQR
}

val <- list(location = list(median = Me), scale = list(IQR = IQR), shape = list(Bowley = Ap, shape = Tp))
class(val) <- "qlss"
return(val)
}

qlss.default <- function(fun = "qnorm", p = 0.1, ...){

if(any(!p > 0) | any(!p < 0.5)) stop("Quantile index out of range: p must be > 0 and < 0.5")
if(any(duplicated(p))) p <- unique(p)
nq <- length(p)

vec3 <- as.numeric(do.call(what = match.fun(fun), args = list(p = 1:3/4, ...)))
Me <- vec3[2]
IQR <- vec3[3] - vec3[1]

PD <- Ap <- Tp <- rep(NA, nq)
for(i in 1:nq){
	vecp <- as.numeric(do.call(what = match.fun(fun), args = list(p = p[i], ...)))
	vecq <- as.numeric(do.call(what = match.fun(fun), args = list(p = 1 - p[i], ...)))

	PD[i] <- vecq - vecp
	Ap[i] <- (vecq - 2*Me + vecp)/PD[i]
	Tp[i] <- PD[i]/IQR
}

val <- list(location = list(median = Me), scale = list(IQR = IQR), shape = list(Bowley = Ap, shape = Tp))
class(val) <- "qlss"
return(val)
}

qlss.formula <- function(formula, p = 0.1, ..., predictLs = NULL){

if(any(!p > 0) | any(!p < 0.5)) stop("Quantile index out of range: p must be > 0 and < 0.5")
if(any(duplicated(p))) p <- unique(p)
nq <- length(p)

fit.list <- list(formula = formula, tau = p, method = "fn")
fit.list <- c(fit.list, list(...))

fit <- do.call(rq, args = fit.list)
predict.list <- c(list(object = fit), as.list(predictLs))
vecp <- as.matrix(do.call(predict, args = predict.list))

fit.list$tau <- 1-p
fit <- do.call(rq, args = fit.list)
predict.list$object <- fit
vecq <- as.matrix(do.call(predict, args = predict.list))

fit.list$tau <- 1:3/4
fit <- do.call(rq, args = fit.list)
predict.list$object <- fit
vec3 <- as.matrix(do.call(predict, args = predict.list))

Me <- as.matrix(vec3[,2])
IQR <- as.matrix(vec3[,3] - vec3[,1])
n <- length(Me)

PD <- Ap <- Tp <- matrix(NA, n, nq)
for(i in 1:nq){
	PD[,i] <- vecq[,i] - vecp[,i]
	Ap[,i] <- (vecq[,i] - 2*Me + vecp[,i])/PD[,i]
	Tp[,i] <- PD[,i]/IQR
}

val <- list(location = list(median = Me), scale = list(IQR = IQR), shape = list(Bowley = Ap, shape = Tp))
class(val) <- "qlss"
return(val)

}

plot.qlss <- function(x, z, type = "l", ...){

par(mfrow = c(2,2))
plot(z, x$location$median, ylab = "median", type = type, ...)
plot(z, x$scale$IQR, ylab = "IQR", type = type, ...)
plot(z, x$shape$Bowley, ylab = "Bowley", type = type, ...)
plot(z, x$shape$shape, ylab = "shape", type = type, ...)

}



######################################################################
### Transformation models
######################################################################

# Basis transformation functions

powerbasis <- function(x, lambda){
(x^(lambda) - 1)/lambda
}

invpowerbasis <- function(x, lambda){

val <- (lambda*x + 1)^(1/lambda)
return(val)
}

powrecbasis <- function(x, lambda){
1/(2*lambda) * (x^lambda - x^(-lambda))
}

# Logit transformation
logit <- function(theta, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1-omega
val <- log(theta/(1-theta))
return(val)
}

# Inverse logit transformation
invlogit <- function(x){
val <- exp(x)/(1 + exp(x))
val[val < 0] <- 0
val[val > 1] <- 1
return(val)
}

# c-log-log transformation
cloglog <- function(theta, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1-omega
val <- log(-log(1 - theta))
return(val)
}

# inverse c-log-log transformation
invcloglog <- function(x){
val <- 1 - exp(-exp(x))
return(val)
}

# Proposal I (one parameter)

mcjI <- function(x, lambda, symm = TRUE, bounded = FALSE, omega = 0.001){

if(bounded){
	if(any(x < 0) | any(x > 1)) stop("x outside interval")
	x[x == 0] <- omega
	x[x == 1] <- 1 - omega
} else {
	if(any(x <= 0)) stop("x must be strictly positive")
}

if(bounded){
	if(symm){
		x <- x/(1-x)
	} else {
		x <- -log(1-x)
	}
} else {
	if(!symm){
		x <- log(1+x)
	}
}

if(lambda != 0){
	val <-  powrecbasis(x, lambda)
	} else {val <- log(x)}

return(val)

}

# Inverse proposal I (one parameter)

invmcjI <- function(x, lambda, symm = TRUE, bounded = FALSE){

if(bounded){
	if(symm){
		if(lambda != 0){
				x <- lambda*x
				y <- (x + sqrt(1 + x^2))^(1/lambda)
				val <- y/(1+y)
			} else {
			val <- invlogit(x)
		}
	} else {
		if(lambda != 0){
				x <- lambda*x
				val <- (x + sqrt(1 + x^2))^(1/lambda)
				val <- 1 - exp(-val)
			} else {
			val <- invcloglog(x)
		}
	}
}
else {
	if(lambda != 0){
		x <- lambda*x
		val <- (x + sqrt(1 + x^2))^(1/lambda)
	} else {val <- exp(x)}

	if(!symm){
		val <- exp(val) - 1
	}
}

return(val)
}

# Proposal II (two parameters)

mcjII <- function(x, lambda, delta, bounded = FALSE, omega = 0.001){

if(bounded){
	if(any(x < 0) | any(x > 1)) stop("x outside interval")
	x[x == 0] <- omega
	x[x == 1] <- 1 - omega
} else {
	if(any(x <= 0)) stop("x must be strictly positive")
}

if(lambda == 0){
	lambda <- 1e-10
}
if(delta == 0){
	delta <- 1e-10
}

if(bounded){
	x <- ((1-x)^(-delta) - 1)/delta
	val <- powrecbasis(x, lambda)
}
else {
	x <- x/(1+x)
	x <- ((1-x)^(-delta) - 1)/delta
	val <- powrecbasis(x, lambda)
}


return(val)

}

# Inverse proposal II (two parameters)

invmcjII <- function(x, lambda, delta, bounded = FALSE){

if(lambda == 0){
	lambda <- 1e-10
}
if(delta == 0){
	delta <- 1e-10
}

if(bounded){
	x <- lambda*x
	val <- (x + sqrt(1 + x^2))^(1/lambda)
	val <- 1 - (val*delta + 1)^(-1/delta)
}
else {
	x <- lambda*x
	val <- (x + sqrt(1 + x^2))^(1/lambda)
	val <- 1 - (val*delta + 1)^(-1/delta)
	val <- val/(1-val)
}

return(val)
}

# Aranda-Ordaz transformation (symmetric and asymmetric)

ao <- function(theta, lambda, symm = TRUE, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1 - omega

if(symm){
	if(lambda != 0){
	val <- (2/lambda)* ((theta^lambda) - (1-theta)^lambda)/((theta^lambda) + (1-theta)^lambda)
	} else {
	val <- logit(theta)
	}
} else {
	if(lambda != 0){
	val <- log(((1-theta)^(-lambda) - 1)/lambda)
	} else {
	val <- cloglog(theta)
	}
}


return(val)

}

# Inverse Aranda-Ordaz transformation (symmetric and asymmetric)

invao <- function(x, lambda, symm  = TRUE){

if(symm){
	if(lambda != 0){
		y <- (lambda*x/2)
		a <- (1 + y)^(1/lambda)
		b <- (1 - y)^(1/lambda)
		val <- rep(1, length(x))
		val <- ifelse(abs(y) < 1, a/(a + b), val)
		val[y <= -1] <- 0
	} else {val <- invlogit(x)}
} else {
	if(lambda != 0){
		y <- lambda*exp(x)
		val <- ifelse(y > -1, 1 - (1 + y)^(-1/lambda), 1)
		
	} else {val <- invcloglog(x)}
}

return(as.numeric(val))
}

# Box-Cox transformation

bc <- function(x, lambda){
if(any(x <= 0)) stop("x must be strictly positive")
val <- if(lambda != 0) powerbasis(x, lambda) else log(x)
return(val)
}

# Inverse Box-Cox transformation

invbc <- function(x, lambda){

val <- if(lambda != 0) invpowerbasis(x, lambda) else exp(x)
return(val)
}

# Mapping from x.r[1],x.r[2] to 0,1

map <- function(x, x.r = NULL){
if(is.null(x.r)) x.r <- range(x, na.rm = TRUE)
theta <- (x - x.r[1])/(x.r[2] - x.r[1])
attr(theta, "range") <- x.r
return(theta)
}

# Mapping from 0,1 to x.r[1],x.r[2]

invmap <- function(x, x.r = NULL){

if(is.null(x.r)) x.r <- attr(x, "range")

(x.r[2] - x.r[1]) * x + x.r[1]

}

# L1-norm and residual cusum loss functions
l1Loss <- function(x, tau, weights){

ind <- ifelse(x < 0, 1, 0)
sum(weights * x * (tau - ind))/sum(weights)

}

rcLoss <- function(lambda, x, y, tsf, symm = TRUE, bounded = bounded, tau = 0.5, method.rq = "fn"){

if(length(tau) > 1) stop("One quantile at a time")
n <- length(y)
out <- rep(NA, n)

z <- switch(tsf,
	mcjI = mcjI(y, lambda, symm, bounded, omega = 0.001),
	bc = bc(y, lambda),
	ao = ao(y, lambda, symm, omega = 0.001)
	)

Rfun <- function(x, t, e) mean(apply(x, 1, function(xj,t) all(xj <= t), t = t) * e)

fit <- try(rq(z ~ x - 1, tau = tau, method = method.rq), silent = T)

if(class(fit)!="try-error"){
	e <- as.numeric(fit$residuals <= 0)
	#out <- apply(x, 1, function(t, z, e) Rfun(z, t, e), z = x, e = tau - e)
	for(i in 1:n){
	#out[i] <- mean(apply(x, 1, function(x,t) all(x <= t), t = x[i,]) * (tau - e))
	out[i] <- mean(apply(t(x) <= x[i,], 2, function(x) all(x)) * (tau - e))	
	}
}

return(mean(out^2))

}


######################################################################
# One-parameter transformations (MCJI, Box-Cox, Aranda-Ordaz)
######################################################################

# Two-stage estimator

tsrq <- function(formula, tsf = "mcjI", symm = TRUE, bounded = FALSE, lambda = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", ...){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(any(duplicated(tau))) tau <- unique(tau)

nq <- length(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
if(!is.data.frame(data)) stop("`data' must be a data frame")
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <-  y.old <- model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else weights

if(tsf == "mcjII") stop("For two-parameter transformations, see 'tsrq2' and 'nlrq2'")
isBounded <- (tsf == "mcjI" && bounded)
isBounded <- tsf == "ao" || isBounded
if(isBounded) y <- map(y)

if(is.null(lambda)){
	lambda <- if(tsf == "ao" & symm == FALSE) seq(-2, 2, by = 0.005)
		else seq(0, 2, by = 0.005)
}

n <- length(y)
p <- ncol(x)
nl <- length(lambda)

zhat <- res <- array(NA, dim = c(n, nq, nl))
matLoss <- rejected <- matrix(NA, nq, nl)
Ind <- array(NA, dim = c(n, nq, nl))

f.tmp <- update.formula(formula, newresponse ~ .)
data.tmp <- data
if(!missing(subset))
	data.tmp <- subset(data, subset)

for(i in 1:nl){

# transform response
data.tmp$newresponse <- switch(tsf,
	mcjI = mcjI(y, lambda[i], symm, bounded, omega = 0.001),
	bc = bc(y, lambda[i]),
	ao = ao(y, lambda[i], symm, omega = 0.001)
	)

# estimate linear QR for for sequence of lambdas
	for(j in 1:nq){
	fit <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)
		if(class(fit)!="try-error"){
		zhat[,j,i] <- predict(fit)
		
		Fitted <- switch(tsf,
			mcjI = invmcjI(zhat[,j,i], lambda[i], symm, bounded),
			bc = invbc(zhat[,j,i], lambda[i]),
			ao = invao(zhat[,j,i], lambda[i], symm),
		)
		
		res <- y - Fitted
		
		if(tsf == "bc"){
			FLAG <- lambda[i]*zhat[,j,i] + 1 > 0
			Ind[,j,i] <- FLAG
			rejected[j,i] <- mean(!FLAG)
			}

		if(tsf == "ao" & symm == TRUE){
			FLAG <- abs(lambda[i]*zhat[,j,i]/2) - 1 < 0
			Ind[,j,i] <- FLAG
			rejected[j,i] <- mean(!FLAG)
			}
		
		matLoss[j,i] <- l1Loss(res, tau = tau[j], weights = w)
		}
	}
}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)

betahat <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()

for(j in 1:nq){
# transform response with optimal lambda
data.tmp$newresponse <- switch(tsf,
	mcjI = mcjI(y, lambdahat[j], symm, bounded, omega = 0.001),
	bc = bc(y, lambdahat[j]),
	ao = ao(y, lambdahat[j], symm, omega = 0.001)
	)
fit[[j]] <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)

if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	tmp <- x%*%matrix(betahat[,j])
	Fitted[,j] <- switch(tsf,
		mcjI = invmcjI(tmp, lambdahat[j], symm, bounded),
		bc = invbc(tmp, lambdahat[j]),
		ao = invao(tmp, lambdahat[j], symm)
	)
	}

}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$mf <- mf
fit$y <- y.old
if(isBounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
attr(fit$tsf, "bounded") <- bounded
attr(fit$tsf, "isBounded") <- isBounded
attr(fit$tsf, "npar") <- 1
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
class(fit) <- "rqt"
return(fit)
}


# Cusum process estimator

rcrq <- function(formula, tsf = "mcjI", symm = TRUE, bounded = FALSE, lambda = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", ...){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(any(duplicated(tau))) tau <- unique(tau)

nq <- length(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
if(!is.data.frame(data)) stop("`data' must be a data frame")
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")
x <- model.matrix(mt, mf, contrasts)
y <-  y.old <- model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else weights

if(tsf == "mcjII") stop("For two-parameter transformations, see 'tsrq2' and 'nlrq2'")
isBounded <- (tsf == "mcjI" && bounded)
isBounded <- tsf == "ao" || isBounded
if(isBounded) y <- map(y)

if(is.null(lambda)){
	lambda <- if(tsf == "ao" & symm == FALSE) seq(-2, 2, by = 0.05)
		else seq(0, 2, by = 0.05)
}

n <- length(y)
p <- ncol(x)
nl <- length(lambda)

zhat <- res <- array(NA, dim = c(n, nq, nl))
matLoss <- rejected <- matrix(NA, nq, nl)
Ind <- array(NA, dim = c(n, nq, nl))

for(i in 1:nl){


# estimate linear QR for for sequence of lambdas

	for(j in 1:nq){
	matLoss[j,i] <- rcLoss(lambda[i], x, y, tsf, symm = symm, bounded = bounded, tau = tau[j], method.rq = method)
	}

}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)


betahat <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()

f.tmp <- update.formula(formula, newresponse ~ .)
data.tmp <- data
if(!missing(subset))
	data.tmp <- subset(data, subset)

for(j in 1:nq){
# transform response with optimal lambda
data.tmp$newresponse <- switch(tsf,
	mcjI = mcjI(y, lambdahat[j], symm, bounded, omega = 0.001),
	bc = bc(y, lambdahat[j]),
	ao = ao(y, lambdahat[j], symm, omega = 0.001)
	)

fit[[j]] <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)

	if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	tmp <- x%*%matrix(betahat[,j])
	Fitted[,j] <- switch(tsf,
		mcjI = invmcjI(tmp, lambdahat[j], symm, bounded),
		bc = invbc(tmp, lambdahat[j]),
		ao = invao(tmp, lambdahat[j], symm))
	}

}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$mf <- mf
fit$y <- y.old
if(isBounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
attr(fit$tsf, "bounded") <- bounded
attr(fit$tsf, "isBounded") <- isBounded
attr(fit$tsf, "npar") <- 1
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
class(fit) <- "rqt"
return(fit)

}


######################################################################
# Two-parameter transformations (MCJII)
######################################################################

# Two-stage estimator

tsrq2 <- function(formula, bounded = FALSE, lambda = NULL, delta = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", ...){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(any(duplicated(tau))) tau <- unique(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), 
	names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
if (method == "model.frame") 
	return(mf)
mt <- attr(mf, "terms")
y <- y.old <- model.response(mf)
w <- if (missing(weights)) rep(1, length(y)) else weights
x <- model.matrix(mt, mf, contrasts)
if(bounded) y <- map(y)

if(is.null(lambda)){
	lambda <- seq(0, 2, by = 0.005)
}
if(is.null(delta)){
	delta <- seq(0, 2, by = 0.005)
}
nl <- length(lambda)
nd <- length(delta)

n <- length(y)
p <- ncol(x)
nq <- length(tau)

matLoss <- array(NA, dim = c(nl, nd, nq), dimnames = list(lambda = 1:nl, delta = 1:nd, tau = tau))

f.tmp <- update.formula(formula, newresponse ~ .)
data.tmp <- data
if(!missing(subset))
	data.tmp <- subset(data, subset)

for(k in 1:nd){
	for(i in 1:nl){
	# transform response
	data.tmp$newresponse <- mcjII(y, lambda[i], delta[k], bounded, omega = 0.001)
		for(j in 1:nq){
		fit <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)
			if(class(fit)!="try-error"){
			Fitted <- invmcjII(predict(fit), lambda[i], delta[k], bounded)
			matLoss[i,k,j] <- l1Loss(y - Fitted, tau = tau[j], weights = w)
			}
		}

	}
}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
parhat <- apply(matLoss, 3, function(x, lambda, delta){
m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];
return(c(lambda[m[1]], delta[m[2]]))}, lambda = lambda, delta = delta)

betahat <- matrix(NA, p, nq)
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau

fit <- list()

for(j in 1:nq){
# transform response with optimal lambda
data.tmp$newresponse <- mcjII(y, parhat[1,j], parhat[2,j], bounded, omega = 0.001)
fit[[j]] <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)
	if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	Fitted[,j] <- invmcjII(x%*%matrix(betahat[,j]), parhat[1,j], parhat[2,j], bounded)
	}
}

if(bounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
dimnames(parhat) <- list(c("lambda","delta"), paste("tau =", format(round(tau, 3))))

fit$call <- call
fit$method <- method
fit$mf <- mf
fit$y <- y.old
if(bounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$lambda.grid <- lambda
fit$delta.grid <- delta
fit$tsf <- "mcjII"
attr(fit$tsf, "bounded") <- bounded
attr(fit$tsf, "isBounded") <- bounded
attr(fit$tsf, "npar") <- 2
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 3, function(x){m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];return(x[m[1],m[2]])})
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2
class(fit) <- "rqt"
return(fit)
}

# Nelder-Mead optimization (joint estimation)

nlrq2 <- function(formula, par = NULL, bounded = FALSE, tau = 0.5, data, subset, weights, na.action, ...){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
if(any(duplicated(tau))) tau <- unique(tau)

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), 
	names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
mt <- attr(mf, "terms")
y <- y.old <- model.response(mf)
x <- model.matrix(mt, mf, contrasts)
w <- if (missing(weights)) rep(1, length(y)) else weights
if(bounded) y <- map(y)

n <- length(y)
p <- ncol(x)
nq <- length(tau)

if(is.null(par)) par <- rep(0, p + 2)

f <- function(theta, data){
	beta <- matrix(theta[-c(1:2)])
	if(theta[2] < 0) return(Inf)
	Fitted <- invmcjII(data$x %*% beta, lambda = theta[1], delta = theta[2], bounded = data$bounded)
	return(l1Loss(data$y - Fitted, tau = data$tau, weights = data$weights))
}

fit <- list()
betahat <- matrix(NA, p, nq)
parhat <- matrix(NA, 2, nq)
Fitted <- matrix(NA, n, nq)

for(j in 1:nq){
	fit[[j]] <- try(optim(par = par, fn = f, method = "Nelder-Mead", data = list(x = x, y = y, bounded = bounded, tau = tau[j], weights = w)), silent = T)

	if(class(fit[[j]])!="try-error"){
		betahat[,j] <- fit[[j]]$par[-c(1:2)]
		parhat[,j] <- c(fit[[j]]$par[1], fit[[j]]$par[2])
		Fitted[,j] <- invmcjII(x %*% matrix(betahat[,j]), parhat[1,j], parhat[2,j], bounded)
	}
}

if(bounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y.old))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
dimnames(parhat) <- list(c("lambda","delta"), paste("tau =", format(round(tau, 3))))


fit$call <- call
fit$method <- "Nelder-Mead"
fit$mf <- mf
fit$y <- y.old
if(bounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$tsf <- "mcjII"
attr(fit$tsf, "bounded") <- bounded
attr(fit$tsf, "isBounded") <- bounded
attr(fit$tsf, "npar") <- 2
fit$coefficients <- betahat
fit$fitted.values <- drop(Fitted)
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2
class(fit) <- "rqt"
return(fit)
}

######################################################################
# Print, summary, bootstrap, predict, fitted for class rqt
######################################################################

print.rqt <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}
type <- if(attr(x$tsf, "isBounded")) "(bounded response)" else "(positive response)"
tsf <- switch(x$tsf,
	mcjI = "Proposal I",
	bc = "Box-Cox",
	ao = "Aranda-Ordaz",
	mcjII = "Proposal II")
if(x$tsf %in% c("mcjI", "ao")){
	tsf <- paste(tsf, if(attr(x$tsf, "symm"))
	"symmetric" else "asymmetric")
}
tsf <- paste(tsf, "transformation", type)

cat(tsf, "\n")
cat("\nOptimal transformation parameter:\n")
if(x$tsf == "mcjII") print(x$eta) else print(x$lambda)

coef <- x$coefficients
cat("\nCoefficients linear model (transformed scale):\n")
print(coef, ...)

nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
if (!is.null(attr(x, "na.message"))) 
	cat(attr(x, "na.message"), "\n")
invisible(x)
}

predict.rqt <- function(object, newdata, na.action = na.pass, raw = TRUE, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
bounded <- attributes(tsf)$bounded
isBounded <- attributes(tsf)$isBounded
if(tsf == "mcjII"){
	etahat <- object$eta
} else {
	lambdahat <- object$lambda
}
betahat <- object$coefficients

if(missing(newdata)) {x <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

linpred <- x %*% betahat

Fitted <- matrix(NA, nrow(linpred), ncol(linpred))
if(!raw){
	if(!missing(newdata)) attr(linpred, "x") <- x
	return(linpred)
}

if(tsf == "mcjII"){
	for(j in 1:nq){
		Fitted[,j] <- invmcjII(x = linpred[,j], lambda = etahat['lambda',j], delta = etahat['delta',j], bounded = bounded)
	}
} else {
	for(j in 1:nq){
		Fitted[,j] <- switch(tsf,
			mcjI = invmcjI(linpred[,j], lambdahat[j], symm, bounded),
			bc = invbc(linpred[,j], lambdahat[j]),
			ao = invao(linpred[,j], lambdahat[j], symm))
	}
}

if(isBounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
}

if(!missing(newdata)) attr(Fitted, "x") <- x

return(Fitted)

}

boot.fun <- function(data, inds, object){

tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
flag <- !object$tsf %in% c("mcjII")
nn <- if(flag) c(object$term.labels, "lambda") else c(object$term.labels, "lambda", "delta")

if(nq == 1){
	fit <- update(object, data = data[inds,])
	val <- fit$coefficients
	val <- if(flag) c(val, fit$lambda) else c(val, fit$eta)
} else {
	fit <- update(object, data = data[inds,])
	val <- fit$coefficients
	val <- if(flag) rbind(val, fit$lambda) else rbind(val, fit$eta)
}

val <- as.vector(val)
names(val) <- rep(nn, nq)
return(val)

}

summary.rqt <- function(object, alpha = 0.05, se = "boot", R = 50, sim = "ordinary", stype = "i", conditional = FALSE, ...){

call <- match.call(expand.dots = TRUE)

tau <- object$tau
nq <- length(tau)
mpar <- ncol(object$x)
ntot <- mpar + attr(object$tsf, "npar")
if(mpar == 1) object$mf$intercept <- 1

if(!conditional){
	Args <- list()
	Args$data <- object$mf
	Args$statistic <- boot.fun
	Args$object <- object
	Args$R <- R
	Args$sim <- sim
	Args$stype <- stype
	nn <- c("strata","L","m","weights","ran.gen","mle","simple","parallel","ncpus","cl")
	nn <- nn[pmatch(names(call), nn, duplicates.ok = FALSE)]
	nn <- nn[!is.na(nn)]
	if(length(nn) > 0) {tmp <- as.list(call[[nn]]); names(tmp) <- nn; Args <- c(Args, tmp)}
	B <- do.call(boot, args = Args)
	ci <- mapply(boot.ci, index = 1:(ntot*nq), MoreArgs = list(boot.out = B, conf = 1 - alpha, type = "perc"))[4,]
	ci <- t(sapply(ci, function(x) x[4:5]))

	S <- cov(B$t, use = "complete.obs")
	val <- cbind(B$t0, apply(B$t, 2L, mean, na.rm=TRUE) - B$t0, sqrt(diag(S)), ci)
	nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound")
	colnames(val) <- nn
	
	maxn <- seq(0, ntot*nq, by = ntot)[-1]
	minn <- seq(1, ntot*nq, by = ntot)
	ans <- list()
	for(j in 1:nq){
		ans[[j]] <- val[minn[j]:maxn[j], ]
	}
	names(ans) <- tau
	object$B <- B
} else {
	ans <- list()
	for(j in 1:nq){
		Args <- list()
		Args$object <- object[[j]]
		Args$se <- se
		if(se == "boot") Args$R <- R
		nn <- c("covariance","hs","bsmethod","mofn","iid")
		nn <- nn[pmatch(names(call), nn, duplicates.ok = FALSE)]
		nn <- nn[!is.na(nn)]
		if(length(nn) > 0) {tmp <- as.list(call[[nn]]); names(tmp) <- nn; Args <- c(Args, tmp)}
		ans[[j]] <- do.call(summary.rq, args = Args)$coefficients
		if(object$tsf == "mcjII") {
			tmp <- matrix(NA, nrow = 2, ncol = ncol(ans[[j]]))
			tmp[,1] <- object$eta[,j]
			rownames(tmp) <- c("lambda","delta")
		} else {
			tmp <- c(object$lambda[j], rep(NA, ncol(ans[[j]]) - 1))
		}
		ans[[j]] <- rbind(ans[[j]], tmp) 
	}
}

attr(ans, "conditional") <- conditional
object$coefficients <- ans
object$call <- call
class(object) <- c("summary.rqt", class(object))
return(object)

}

print.summary.rqt <- function(x, ...){

if (!is.null(cl <- x$call)) {
	cat("call:\n")
	dput(cl)
	cat("\n")
}

tau <- x$tau
nq <- length(tau)
mpar <- ncol(x$x)

type <- if(attr(x$tsf, "isBounded")) "(bounded response)" else "(positive response)"
tsf <- switch(x$tsf,
	mcjI = "Proposal I",
	bc = "Box-Cox",
	ao = "Aranda-Ordaz",
	mcjII = "Proposal II")
if (x$tsf %in% c("mcjI", "ao")){ 
	tsf <- paste(tsf, if (attr(x$tsf, "symm")) 
		"symmetric"
	else "asymmetric")
}
tsf <- paste(tsf, "transformation", type)

cat(tsf, "\n")

conditional <- if(attr(x$coefficients, "conditional")) "conditional" else "unconditional"
cat("\nSummary for", conditional, "inference\n")

for(i in 1:nq){
cat("\ntau = ", tau[i], "\n")

cat("\nOptimal transformation parameter:\n")
print(x$coefficients[[i]][-c(1:mpar),], ...)

cat("\nCoefficients linear model (transformed scale):\n")
print(x$coefficients[[i]][1:mpar,], ...)
}

nobs <- length(x$y)
p <- ncol(x$x)
rdf <- nobs - p
cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
if (!is.null(attr(x, "na.message"))) 
	cat(attr(x, "na.message"), "\n")
invisible(x)


}

fitted.rqt <- fitted.values.rqt <- function(object, ...){

return(object$fitted.values)

}


coef.rqt <- coefficients.rqt <- function(object, all = FALSE, ...){

if(!class(object) %in% c("rqt")) stop("Class 'rqt' only")

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
nn <- object$term.labels

if(all){
	nn <- if(tsf %in% "mcjII") c(nn, "lambda", "delta") else c(nn, "lambda")
}

if(!all){
	ans <- object$coefficients
} else {
	ans <- if(nq == 1) c(object$coefficients, object$lambda)
		else rbind(object$coefficients, object$lambda)
}

if(nq == 1){
	names(ans) <- nn
} else {
	rownames(ans) <- nn
	colnames(ans) <- paste("tau =", tau)
}

return(ans)

}

##################################################
### Restricted quantiles
##################################################

rrq <- function(formula, tau, data, method = "fn", model = TRUE, ...){

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0)
mf <- mf[c(1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name("model.frame")
mf <- eval.parent(mf)
if (method == "model.frame")
	return(mf)
mt <- attr(mf, "terms")
weights <- as.vector(model.weights(mf))
y <- model.response(mf)
x <- model.matrix(mt, mf, contrasts)

eps <- .Machine$double.eps^(2/3)

if (length(tau) > 1) {
        if (any(tau < 0) || any(tau > 1)) 
            stop("invalid tau:  taus should be >= 0 and <= 1")
        if (any(tau == 0)) 
            tau[tau == 0] <- eps
        if (any(tau == 1)) 
            tau[tau == 1] <- 1 - eps
}

nq <- length(tau)

fit.lad <- rq(formula, tau = 0.5, data = data, method = method)
data$r.lad <- fit.lad$residuals
data$r.abs <- abs(fit.lad$residuals)
beta <- fit.lad$coefficients

fit.lad <- rq(update.formula(formula, r.abs ~ .), tau = 0.5, data = data, method = method)
data$s.lad <- fit.lad$fitted
gamma <- fit.lad$coefficients

zeta <- rq(r.lad ~ s.lad - 1, tau = tau, data = data, method = method)$coefficients

val <- if (length(tau) > 1)
apply(outer(gamma, zeta, "*"), 3, function(x, b) x + b, b = beta)
	else beta + zeta * gamma

fit <- list(coefficients = val, c = zeta, beta = beta, gamma = gamma, tau = tau)
fit$na.action <- attr(mf, "na.action")
fit$formula <- formula
fit$terms <- mt
fit$xlevels <- .getXlevels(mt, mf)
fit$call <- call
fit$weights <- weights
fit$residuals <- drop(fit$residuals)
fit$method <- method
fit$x <- x
fit$y <- y
fit$fitted.values <- x %*% val
attr(fit, "na.message") <- attr(m, "na.message")
if (model)
	fit$model <- mf
class(fit) <- "rrq"
return(fit)
}

rrq.fit <- function(x, y, tau, method = "fn"){

if(length(tau) > 1) stop("only one quantile")

fit.lad <- rq.fit(x, y, tau = 0.5, method = method)
r.lad <- fit.lad$residuals
r.abs <- abs(fit.lad$residuals)
beta <- fit.lad$coefficients

fit.lad <- rq.fit(x, r.abs, tau = 0.5, method = method)
s.lad <- fit.lad$fitted.values
gamma <- fit.lad$coefficients

zeta <- rq.fit(s.lad, r.lad, tau = tau, method = method)$coefficients

val <- beta + zeta * gamma

return(list(coef = val, c = zeta, beta = beta, gamma = gamma, tau = tau))
}

predict.rrq <- function(object, newdata, na.action = na.pass, ...){


tau <- object$tau
nq <- length(tau)
betahat <- object$coefficients

if(missing(newdata)) {x <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

return(x %*% betahat)

}


##################################################
### Marginal effects
##################################################

# Generic

maref <- function(...) UseMethod("maref")

# rqt object

maref.rqt <- function(object, newdata, index = 2, index.extra = NULL, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
bounded <- attributes(tsf)$bounded
isBounded <- attributes(tsf)$isBounded
if (tsf == "mcjII") {
	etahat <- object$eta
}
else {
	lambdahat <- object$lambda
}

betahat <- as.matrix(object$coefficients)
linpred <- predict(object, newdata, raw = FALSE)

if(missing(newdata)) {x <- object$x}
	else {x <- attr(linpred, "x")}
n <- nrow(x)

if(identical(x[,index], rep(1,n))) return(kronecker(matrix(betahat[index,],nrow = 1), rep(1,n)))

if(!is.null(index.extra)){
	if(index %in% index.extra) warning(paste("Index", index, "appears twice as 'index' and in 'index.extra'."))
	param <- matrix(NA, n, nq)
	for(i in 1:nq){
		tmp <- sweep(as.matrix(x[,index.extra]), 2, betahat[index.extra,i], "*")
		param[,i] <- betahat[index,i] + rowSums(tmp)
	}
} else {param <- matrix(betahat[index,], nrow = 1)}

val <- matrix(NA, n, nq)

for(j in 1:nq)(
val[,j] <- switch(tsf,
	mcjI = marefmcjI(linpred[,j], param[,j], lambdahat[j], symm, bounded),
	bc = marefbc(linpred[,j], param[,j], lambdahat[j]),
	ao = marefao(linpred[,j], param[,j], lambdahat[j], symm),
	)
)

return(val)

}

marefbc <- function(x, param, lambda){

if(lambda != 0){
	val <- param * (lambda * x + 1)^(1/lambda - 1)
} else {
	val <- param * exp(x)
}

return(val)
}

marefao <- function(x, param, lambda, symm){

if(symm){
	if(lambda != 0){
		y <- (lambda * x/2)
		a <- (1 + y)^(1/lambda)
		b <- (1 - y)^(1/lambda)
		aprime <- 1/2*param * (1 + y)^(1/lambda - 1)
		bprime <- -1/2*param * (1 - y)^(1/lambda - 1)
		val <- (b/aprime^2 - bprime/a)/(1 + b/a)^2
	} else {
		val <- (param * exp(-x))/(1 + exp(-x))^2
	}
} else {
	if(lambda != 0){
		y <- lambda * exp(x)
		val <- (param * exp(x)) * (1 + y)^(-1/lambda - 1)
	} else {
		val <- 1 + param * exp(x - exp(x))
	}
}

return(val)
}

marefmcjI <- function(x, param, lambda, symm, bounded){

if(bounded){
	if(symm){
		if(lambda != 0){
			y <- lambda*x + sqrt(1 + (lambda*x)^2)
			a <- 1/lambda * (y^(-1/lambda - 1)) * (lambda * param + 1/2*(1 + (lambda*x)^2)^(-1/2) * (2*x*param*lambda^2))
			b <- (1 + y^(-1/lambda))^2
			val <- a/b
		} else {
			val <- (param * exp(-x))/(1 + exp(-x))^2
		}
	} else {
		if(lambda != 0){
			y <- lambda*x + sqrt(1 + (lambda*x)^2)
			val <- 1/lambda * (y^(1/lambda - 1)) * (lambda * param + 1/2*(1 + (lambda*x)^2)^(-1/2) * (2*x*param*lambda^2)) * (exp(-y^(1/lambda)))
		} else {
			val <- 1 + param * exp(x - exp(x))	
		}
	}
} else {
	if (lambda != 0) {
		y <- lambda*x + sqrt(1 + (lambda*x)^2)
		val <- 1/lambda * (y^(1/lambda - 1)) * (lambda*param + 1/2*(1 + (lambda*x)^2)^(-1/2) * (2*x*param*lambda^2))
	} else {
		val <- param*exp(x)
	}
	if (!symm) {
		val <- val * exp(val)
	}
}

return(val)

}


##################################################
### Multiple imputation
##################################################

mice.impute.rq <- function (y, ry, x, tsf = "none", symm = TRUE, bounded = FALSE, lambda = NULL, epsilon = 0.001, method.rq = "fn", ...) 
{
    x <- cbind(1, as.matrix(x))
	y.old <- y
	isBounded <- (tsf == "mcjI" && bounded)
	isBounded <- tsf == "ao" || isBounded

	if(is.null(lambda))
		lambda <- 0

	if(tsf %in% c("mcjI","bc","ao")){
		if(isBounded) y <- map(y)
		z <- switch(tsf,
			mcjI = mcjI(y, lambda, symm, bounded, omega = 0.001),
			bc = bc(y, lambda),
			ao = ao(y, lambda, symm, omega = 0.001)
			)
	} else {z <- y}

	n <- sum(!ry)
	p <- ncol(x)
	u <- round(runif(n, epsilon, 1 - epsilon)*1e3)
	u <- ifelse(u %in% c(1:4,996:999), u/1e3, (u - u %% 5)/1e3)
	taus <- unique(u)
	nt <- length(taus)

	xobs <- x[ry, ]
	yobs <- z[ry]
	xmis <- x[!ry,]
	fit <- matrix(NA, p, nt)
	for(j in 1:nt){
		fit[,j] <- as.numeric(rq.fit(xobs, yobs, tau = taus[j], method = method.rq)$coefficients)
	}
	# n times nt matrix
	ypred <- xmis%*%fit
	# diagonal of n times n matrix
	ypred <- diag(ypred[,match(u, taus)])

	if(tsf %in% c("mcjI","bc","ao")){
		val <- switch(tsf,
			mcjI = invmcjI(ypred, lambda, symm, bounded),
			bc = invbc(ypred, lambda),
			ao = invao(ypred, lambda, symm));
		if(isBounded) val <- invmap(val, range(y.old))
	} else {val <- ypred}
	
    return(val)
}

# Impute using restricted quantiles

mice.impute.rrq <- function (y, ry, x, tsf = "none", symm = TRUE, bounded = FALSE, lambda = NULL, epsilon = 0.001, method.rq = "fn", ...) 
{
    x <- cbind(1, as.matrix(x))
	y.old <- y
	isBounded <- (tsf == "mcjI" && bounded)
	isBounded <- tsf == "ao" || isBounded

	if(is.null(lambda))
		lambda <- 0

	if(tsf %in% c("mcjI","bc","ao")){
		if(isBounded) y <- map(y)
		z <- switch(tsf,
			mcjI = mcjI(y, lambda, symm, bounded, omega = 0.001),
			bc = bc(y, lambda),
			ao = ao(y, lambda, symm, omega = 0.001)
			)
	} else {z <- y}

	n <- sum(!ry)
	p <- ncol(x)
	u <- round(runif(n, epsilon, 1 - epsilon)*1e3)
	u <- ifelse(u %in% c(1:4,996:999), u/1e3, (u - u %% 5)/1e3)
	taus <- unique(u)
	nt <- length(taus)

	xobs <- x[ry, ]
	yobs <- z[ry]
	xmis <- x[!ry,]
	fit <- matrix(NA, p, nt)
	for(j in 1:nt){
		fit[,j] <- as.numeric(rrq.fit(xobs, yobs, tau = taus[j], method = method.rq)$coef)
	}
	# n times nt matrix
	ypred <- xmis%*%fit
	# diagonal of n times n matrix
	ypred <- diag(ypred[,match(u, taus)])
	
	if(tsf %in% c("mcjI","bc","ao")){
		val <- switch(tsf,
			mcjI = invmcjI(ypred, lambda, symm, bounded),
			bc = invbc(ypred, lambda),
			ao = invao(ypred, lambda, symm));
		if(isBounded) val <- invmap(val, range(y.old))
	} else {val <- ypred}

    return(val)
}


##################################################
### QR for counts
##################################################

rq.counts <- function (formula, data, tsf = "bc", symm = TRUE, lambda = 0, weights = NULL, offset = NULL, contrasts = NULL, tau = 0.5, M = 50, zeta = 1e-5, B = 0.999, cn = NULL, alpha = 0.05, method = "fn") 
{
nq <- length(tau)
if (nq > 1) 
	stop("One quantile at a time")

call <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "weights"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame())
mt <- attr(mf, "terms")

y <- model.response(mf, "numeric")
w <- as.vector(model.weights(mf))
if (!is.null(w) && !is.numeric(w)) 
  stop("'weights' must be a numeric vector")
if(is.null(w))
  w <- rep(1, length(y))
x <- model.matrix(mt, mf, contrasts)
p <- ncol(x)
n <- nrow(x)
term.labels <- colnames(x)

if (is.null(offset)) 
	offset <- rep(0, n)


Fn <- function(x, cn){
	xf <- floor(x)
	df <- x - xf
	if(df < cn & x >= 1){
		val <- xf - 0.5 + df/(2*cn)
	}
	if(any(cn <= df & df < (1 - cn), x < 1)){
		val <- xf
	}

	if(df >= (1 - cn)){
		val <- xf + 0.5 + (df - 1)/(2*cn)
	}

	return(val)
}

Fvec <- Vectorize(Fn)

# Add noise
Z <- replicate(M, addnoise(y, centered = FALSE, B = B))

# Transform Z
TZ <- apply(Z, 2, function(x, off, tsf, symm, lambda, tau, zeta){
	z <- ifelse((x - tau) > zeta, x - tau, zeta);
	switch(tsf,
		mcjI = mcjI(z, lambda, symm, bounded = FALSE, omega = 0.001),
		bc = bc(z, lambda)) - off
	}, off = offset, tsf = tsf, symm = symm, lambda = lambda, tau = tau, zeta = zeta)

# Fit linear QR on TZ
fit <- apply(TZ, 2, function(y, x, weights, tau, method) 
	rq.wfit(x = x, y = y, tau = tau, weights = weights, method = method), x = x, tau = tau, weights = w, method = method)
	
# predicted values
yhat <- sapply(fit, function(obj, x) x %*% obj$coefficients, x = x)
yhat <- as.matrix(yhat)

# sweep offset back in
linpred <- sweep(yhat, 1, offset, "+")

# back-transform + offset tau
zhat <- matrix(NA, n, M)
for(i in 1:M){
zhat[,i] <- tau + switch(tsf,
	mcjI = invmcjI(linpred[,i], lambda, symm, bounded = FALSE),
	bc = invbc(linpred[,i], lambda))
}
	
# covariance matrix
if(is.null(cn)) cn <- 0.5 * log(log(n))/sqrt(n)
F <- apply(zhat, 2, Fvec, cn = cn)
Fp <- apply(zhat + 1, 2, Fvec, cn = cn)

multiplier <- (tau - (TZ <= yhat))^2
a <- array(NA, dim = c(p, p, M))
for (i in 1:M) a[, , i] <- t(x * multiplier[, i]) %*% x/n

multiplier <- tau^2 + (1 - 2 * tau) * (y <= (zhat - 1)) + 
	((zhat - y) * (zhat - 1 < y & y <= zhat)) * (zhat - y - 
		2 * tau)
b <- array(NA, dim = c(p, p, M))
for (i in 1:M) b[, , i] <- t(x * multiplier[, i]) %*% x/n

multiplier <- (zhat - tau) * (F <= Z & Z < Fp)
d <- array(NA, dim = c(p, p, M))
sel <- rep(TRUE, M)
for (i in 1:M) {
	tmpInv <- try(solve(t(x * multiplier[, i]) %*% x/n), 
		silent = TRUE)
	if (class(tmpInv) != "try-error") 
		{d[, , i] <- tmpInv}
	else {sel[i] <- FALSE}
}
    
dad <- 0
dbd <- 0
for (i in (1:M)[sel]) {
	dad <- dad + d[, , i] %*% a[, , i] %*% d[, , i]
	dbd <- dbd + d[, , i] %*% b[, , i] %*% d[, , i]
}
    
m.n <- sum(sel)
if (m.n != 0) {
	V <- dad/(m.n^2) + (1 - 1/m.n) * dbd * 1/m.n
	V <- V/n ## CHECK V AND WEIGHTS
	stds <- sqrt(diag(V))
	} else {
	stds <- NA
	warning("Standard error not available")
	}

betahat <- sapply(fit, function(x) x$coefficients)
betahat <- if (p == 1) mean(betahat) else rowMeans(betahat)

linpred <- if (p == 1) {
	mean(linpred[1, ])
} else {
	rowMeans(linpred)
}

Fitted <- tau + switch(tsf,
	mcjI = invmcjI(linpred, lambda, symm, bounded = FALSE),
	bc = invbc(linpred, lambda))

lower <- betahat + qt(alpha/2, n - p) * stds
upper <- betahat + qt(1 - alpha/2, n - p) * stds
tP <- 2 * pt(-abs(betahat/stds), n - p)

ans <- cbind(betahat, stds, lower, upper, tP)
colnames(ans) <- c("Value", "Std. Error", "lower bound", "upper bound", "Pr(>|t|)")
rownames(ans) <- names(betahat) <- term.labels

fit <- list()
fit$call <- call
fit$method <- method
fit$mf <- mf
fit$x <- x
fit$weights <- w
fit$offset <- offset
fit$tau <- tau
fit$lambda <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
fit$coefficients <- betahat
fit$M <- M
fit$Mn <- m.n
fit$fitted.values <- Fitted
fit$tTable <- ans
fit$Cov <- V
fit$levels <- .getXlevels(mt, mf)
fit$terms <- mt
fit$term.labels <- term.labels
fit$rdf <- n - p

class(fit) <- "rq.counts"
	
return(fit)
}

coef.rq.counts <- coefficients.rq.counts <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$coefficients

if(nq == 1){
  names(ans) <- object$term.labels
}

return(ans)

}

fitted.rq.counts <- fitted.values.rq.counts <- function(object, ...){

return(object$fitted.values)

}

predict.rq.counts <- function(object, newdata, offset, na.action = na.pass, raw = TRUE, ...) 
{

tsf <- object$tsf
symm <- attributes(tsf)$symm
lambda <- object$lambda

if(missing(newdata)){
	linpred <- drop(object$x %*% object$coefficients) + object$offset
} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
	x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
	if(missing(offset)) offset <- rep(0, nrow(x))
	linpred <- drop(x %*% object$coefficients) + offset
}

if (!raw) {
	if(!missing(newdata)) attr(linpred, "x") <- x
	return(linpred)
}

Fitted <- object$tau + switch(tsf,
	mcjI = invmcjI(linpred, lambda, symm, bounded = FALSE),
	bc = invbc(linpred, lambda))
    
return(Fitted)
}

residuals.rq.counts <- function(object, ...){

ans <- drop(object$y) - predict(object, raw = TRUE, ...)
return(ans)

}

print.rq.counts <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    tau <- x$tau
    nq <- length(tau)
    cat("Call: ")
    dput(x$call)
    cat("\n")
    if (nq == 1) {
        cat(paste("Quantile", tau, "\n"))
        cat("\n")
        cat("Fixed effects:\n")
        printCoefmat(x$tTable, signif.stars = TRUE, P.values = TRUE)
    }
    else {
    NULL
	}
}

maref.rq.counts <- function(object, newdata, index = 2, index.extra = NULL, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
lambda <- object$lambda

betahat <- as.matrix(object$coefficients)
linpred <- as.matrix(predict(object, newdata, raw = FALSE))

if(missing(newdata)) {x <- object$x}
	else {x <- attr(linpred, "x")}
n <- nrow(x)

if(identical(x[,index], rep(1,n))) return(kronecker(matrix(betahat[index,],nrow = 1), rep(1,n)))

if(!is.null(index.extra)){
	if(index %in% index.extra) warning(paste("Index", index, "appears twice as 'index' and in 'index.extra'."))
	param <- matrix(NA, n, nq)
	for(i in 1:nq){
		tmp <- sweep(as.matrix(x[,index.extra]), 2, betahat[index.extra,i], "*")
		param[,i] <- betahat[index,i] + rowSums(tmp)
	}
} else {param <- matrix(betahat[index,], nrow = 1)}

val <- matrix(NA, n, nq)

for(j in 1:nq)(
val[,j] <- switch(tsf,
	mcjI = marefmcjI(linpred[,j], param[,j], lambda, symm, bounded = FALSE),
	bc = marefbc(linpred[,j], param[,j], lambda))
)

return(val)

}

addnoise <- function(x, centered = TRUE, B = 0.999) 
{

	n <- length(x)
    if (centered) 
        z <- x + runif(n, -B/2, B/2)
    else z <- x + runif(n, 0, B)
	
    return(z)
}


##################################################
### Khmaladze test
##################################################


pvalue.KT <- function(object, epsilon){

if(class(object) != "KhmaladzeTest") stop("class(object) must be 'KhmaladzeTest'")
tt <- get("KhmaladzeTable")
if(!(epsilon %in% unique(tt$epsilon))) stop("'epsilon' must be in c(0.05,0.10,0.15,0.20,0.25,0.30)")

p <- length(object$THn)
ans <- matrix(NA, p + 1, 2)
colnames(ans) <- c("Value", "Pr")

sel <- tt[tt$p == p & tt$epsilon == epsilon,3:5]
alpha <- c(0.01,0.05,0.1)

ans[1,1] <- object$Tn
ans[1,2] <- min(c(1,alpha[object$Tn > sel]))

ans[2:(p+1),1] <- object$THn

if(p==1){val <- min(c(1,alpha[object$THn > sel]))} else
{val <- rep(0,p); for(i in 1:p) val[i] <- min(c(1,alpha[object$THn[i] > sel]))}
ans[2:(p+1),2] <- val

return(ans)

}


KhmaladzeTable <- structure(list(p = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 
11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 
4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 
18L, 19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 
12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 4L, 
5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 
19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 
13L, 14L, 15L, 16L, 17L, 18L, 19L, 20L, 1L, 2L, 3L, 4L, 5L, 6L, 
7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 19L, 
20L), epsilon = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 
0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 
0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3), alpha01 = c(2.721, 4.119, 
5.35, 6.548, 7.644, 8.736, 9.876, 10.79, 11.81, 12.91, 14.03, 
15, 15.93, 16.92, 17.93, 18.85, 19.68, 20.63, 21.59, 22.54, 2.64, 
4.034, 5.267, 6.34, 7.421, 8.559, 9.573, 10.53, 11.55, 12.54, 
13.58, 14.65, 15.59, 16.52, 17.53, 18.46, 19.24, 20.21, 21.06, 
22.02, 2.573, 3.908, 5.074, 6.148, 7.247, 8.355, 9.335, 10.35, 
11.22, 12.19, 13.27, 14.26, 15.22, 16.12, 17.01, 17.88, 18.78, 
19.7, 20.53, 21.42, 2.483, 3.742, 4.893, 6.023, 6.985, 8.147, 
9.094, 10.03, 10.9, 11.89, 12.85, 13.95, 14.86, 15.69, 16.55, 
17.41, 18.19, 19.05, 19.96, 20.81, 2.42, 3.633, 4.737, 5.818, 
6.791, 7.922, 8.856, 9.685, 10.61, 11.48, 12.48, 13.54, 14.34, 
15.26, 16, 16.81, 17.59, 18.49, 19.4, 20.14, 2.32, 3.529, 4.599, 
5.599, 6.577, 7.579, 8.542, 9.413, 10.27, 11.15, 12.06, 12.96, 
13.82, 14.64, 15.46, 16.25, 17.04, 17.85, 18.78, 19.48), alpha05 = c(2.14, 
3.393, 4.523, 5.56, 6.642, 7.624, 8.578, 9.552, 10.53, 11.46, 
12.41, 13.34, 14.32, 15.14, 16.11, 16.98, 17.9, 18.83, 19.72, 
20.58, 2.102, 3.287, 4.384, 5.43, 6.465, 7.412, 8.368, 9.287, 
10.26, 11.17, 12.1, 13, 13.9, 14.73, 15.67, 16.56, 17.44, 18.32, 
19.24, 20.11, 2.048, 3.199, 4.269, 5.284, 6.264, 7.197, 8.125, 
9.044, 9.963, 10.85, 11.77, 12.61, 13.48, 14.34, 15.24, 16.06, 
16.93, 17.8, 18.68, 19.52, 1.986, 3.1, 4.133, 5.091, 6.07, 6.985, 
7.887, 8.775, 9.672, 10.52, 11.35, 12.22, 13.09, 13.92, 14.77, 
15.58, 16.43, 17.3, 18.09, 18.95, 1.923, 3, 4.018, 4.948, 5.853, 
6.76, 7.611, 8.51, 9.346, 10.17, 10.99, 11.82, 12.66, 13.46, 
14.33, 15.09, 15.95, 16.78, 17.5, 18.3, 1.849, 2.904, 3.883, 
4.807, 5.654, 6.539, 7.357, 8.211, 9.007, 9.832, 10.62, 11.43, 
12.24, 13.03, 13.85, 14.61, 15.39, 16.14, 16.94, 17.74), alpha1 = c(1.872, 
3.011, 4.091, 5.104, 6.089, 7.047, 7.95, 8.89, 9.82, 10.72, 11.59, 
12.52, 13.37, 14.28, 15.19, 16.06, 16.97, 17.84, 18.73, 19.62, 
1.833, 2.946, 3.984, 4.971, 5.931, 6.852, 7.77, 8.662, 9.571, 
10.43, 11.29, 12.2, 13.03, 13.89, 14.76, 15.65, 16.53, 17.38, 
18.24, 19.11, 1.772, 2.866, 3.871, 4.838, 5.758, 6.673, 7.536, 
8.412, 9.303, 10.14, 10.98, 11.86, 12.69, 13.48, 14.36, 15.22, 
16.02, 16.86, 17.7, 18.52, 1.73, 2.781, 3.749, 4.684, 5.594, 
6.464, 7.299, 8.169, 9.018, 9.843, 10.66, 11.48, 12.31, 13.11, 
13.91, 14.74, 15.58, 16.37, 17.17, 17.97, 1.664, 2.693, 3.632, 
4.525, 5.406, 6.241, 7.064, 7.894, 8.737, 9.517, 10.28, 11.11, 
11.93, 12.67, 13.47, 14.26, 15.06, 15.83, 16.64, 17.38, 1.602, 
2.602, 3.529, 4.365, 5.217, 6.024, 6.832, 7.633, 8.4, 9.192, 
9.929, 10.74, 11.51, 12.28, 13.05, 13.78, 14.54, 15.3, 16.05, 
16.79)), .Names = c("p", "epsilon", "alpha01", "alpha05", "alpha1"
), class = "data.frame", row.names = c(NA, -120L))


