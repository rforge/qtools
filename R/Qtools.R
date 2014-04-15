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


##################################################
### Sample quantiles
##################################################

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

##################################################
### Transformation models (Geraci and Jones, 2014)
##################################################


# Modified update function

my_update <- function(mod, formula = NULL, weights, data = NULL) {
  call <- getCall(mod)
  if (is.null(call)) {
    stop("Model object does not support updating (no call)", call. = FALSE)
  }
  term <- terms(mod)
  if (is.null(term)) {
    stop("Model object does not support updating (no terms)", call. = FALSE)
  }

  if (!is.null(data)) call$data <- data
  if (!is.null(weights)) call$weights <- weights
  if (!is.null(formula)) call$formula <- update.formula(call$formula, formula)
  env <- attr(term, ".Environment")

  eval(call, env, parent.frame())
}

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

# Proposal I (half line - symmetric and asymmetric)

mcjI <- function(x, lambda, symm = TRUE){
if(any(x <= 0)) stop("x must be positive")

if(!symm){
	x <- log(1+x)
}

if(lambda != 0){
	val <-  powrecbasis(x, lambda)}
	else {val <- log(x)}

return(val)
}

# Inverse proposal I (half line - symmetric and asymmetric)

invmcjI <- function(x, lambda, symm = TRUE){

if(lambda != 0){
	x <- lambda*x
	val <- (x + sqrt(1 + x^2))^(1/lambda)
} else {val <- exp(x)}

if(!symm){
	val <- exp(val) - 1
}

return(val)
}


# Proposal I (bounded - symmetric and asymmetric)
mcjIb <- function(theta, lambda, symm = TRUE, omega = 0.001){
if(any(theta < 0) | any(theta > 1)) stop("theta outside interval")
theta[theta == 0] <- omega
theta[theta == 1] <- 1 - omega

if(symm){
	x <- theta/(1-theta)
} else {x <- -log(1-theta)}

if(lambda != 0){
	val <- powrecbasis(x, lambda)} else {val <- log(x)}

return(val)
}

# Inverse proposal I (bounded - symmetric and asymmetric)
invmcjIb <- function(x, lambda, symm = TRUE){

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

return(val)
}

# Proposal II (two parameters)

mcjII <- function(x, lambda, delta, bounded = FALSE, omega = 0.001){

if(bounded){
	if(any(x < 0) | any(x > 1)) stop("x outside interval")
	x[x == 0] <- omega
	x[x == 1] <- 1 - omega
} else {
	if(any(x <= 0)) stop("x must be positive")
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
if(any(x <= 0)) stop("x must be positive")
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


# Two-stage estimator (Chamberlain 1994, Buchinsky 1995). The proportion of rejected observations for which the loss function is not defined (and hence to be discarded according to Fiztenberger et al's (2010) method) is also reported.

# L1-norm loss function
L1loss <- function(x, tau, weights){

ind <- ifelse(x < 0, 1, 0)
sum(weights * x * (tau - ind))/sum(weights)

}

# 1-parameter transformations (all transformations except mcjII)

tsrq <- function(formula, tsf, symm = TRUE, lambda = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", se = "nid", ...){

if(any(tau <= 0) | any(tau >= 1)) stop("Quantile index out of range")
tau <- round(sort(tau), digits = 4)
if(any(duplicated(tau))) stop("Quantile indices duplicated or too close")

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
y <-  model.response(mf, "numeric")
w <- if (missing(weights)) rep(1, length(y)) else weights

theta <- map(y)

if(se == "rank") stop("Not implemented")
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
	bc = bc(y, lambda[i]),
	mcjI = mcjI(y, lambda[i], symm),
	ao = ao(theta, lambda[i], symm, omega = 0.001),
	mcjIb = mcjIb(theta, lambda[i], symm, omega = 0.001)
	)

# estimate linear QR for for sequence of lambdas
	for(j in 1:nq){
	fit <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)
		if(class(fit)!="try-error"){
		zhat[,j,i] <- predict(fit)
		
		Fitted <- switch(tsf,
			bc = invbc(zhat[,j,i], lambda[i]),
			mcjI = invmcjI(zhat[,j,i], lambda[i], symm),
			ao = invao(zhat[,j,i], lambda[i], symm),
			mcjIb = invmcjIb(zhat[,j,i], lambda[i], symm)
		)
		
		res <- switch(tsf,
			bc = y - Fitted,
			mcjI = y - Fitted,
			ao = theta - Fitted,
			mcjIb = theta - Fitted
		)
		
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
		
		matLoss[j,i] <- L1loss(res, tau = tau[j], weights = w)
		}
	}
}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)


betahat <- ses <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()

for(j in 1:nq){
# transform response with optimal lambda
data.tmp$newresponse <- switch(tsf,
	bc = bc(y, lambdahat[j]),
	mcjI = mcjI(y, lambdahat[j], symm),
	ao = ao(theta, lambdahat[j], symm, omega = 0.001),
	mcjIb = mcjIb(theta, lambdahat[j], symm, omega = 0.001)
	)
fit[[j]] <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)

if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	ses[,j] <- summary(fit[[j]], se = se, ...)$coefficients[,2]
	tmp <- x%*%matrix(betahat[,j])
	Fitted[,j] <- switch(tsf,
		bc = invbc(tmp, lambdahat[j]),
		mcjI = invmcjI(tmp, lambdahat[j], symm),
		ao = invao(tmp, lambdahat[j], symm),
		mcjIb = invmcjIb(tmp, lambdahat[j], symm)
	)
	}

}

if(tsf %in% c("ao","mcjIb")){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$y <- y
if(tsf %in% c("ao","mcjIb")) fit$theta <- theta
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$std.error <- ses
fit$Fitted <- Fitted
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
class(fit) <- "tsrq"
return(fit)
}

predict.tsrq <- function(object, newdata, na.action = na.pass, raw = TRUE, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
lambdahat <- object$lambda
betahat <- object$coefficients

if(missing(newdata)) {X <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

tmp <- X %*% betahat
Fitted <- matrix(NA, nrow(tmp), ncol(tmp))

if(!raw) return(tmp)

for(j in 1:nq){
	Fitted[,j] <- switch(tsf,
		bc = invbc(tmp[,j], lambdahat[j]),
		mcjI = invmcjI(tmp[,j], lambdahat[j], symm),
		ao = invao(tmp[,j], lambdahat[j], symm),
		mcjIb = invmcjIb(tmp[,j], lambdahat[j], symm))
}

if(tsf %in% c("ao","mcjIb")){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
}

return(Fitted)

}


boot.tsrq <- function(object, R = 50, seed = round(runif(1, 1, 10000))){

set.seed(seed)
tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
obsS <- replicate(R, sample(all.obs, size = n, replace = TRUE))
npars <- ncol(object$x)
ntot <- npars + if(class(object) == "tsrq") 1 else if(class(object) == "tsrq2") 2
nn <- if(class(object) == "tsrq") c(object$term.labels, "lambda") else if(class(object) == "tsrq2") c(object$term.labels, "lambda", "delta")

if(nq == 1){
  bootmat <- matrix(NA, R, ntot);
  colnames(bootmat) <- nn
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i])
    w[match(names(a), all.obs)] <- as.numeric(a)
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars] <- fit$coefficients
		if(class(object) == "tsrq") bootmat[i,(npars+1):ntot] <- fit$lambda
		if(class(object) == "tsrq2") bootmat[i,(npars+1):ntot] <- fit$eta
	}
  }
} else {
  bootmat <- array(NA, dim = c(R, ntot, nq), dimnames = list(NULL, nn, paste("tau = ", format(tau, digits = 4), sep ="")));
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i]);
    w[match(names(a), all.obs)] <- as.numeric(a);
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars,] <- fit$coefficients
		if(class(object) == "tsrq") bootmat[i,(npars+1):ntot,] <- fit$lambda
		if(class(object) == "tsrq2") bootmat[i,(npars+1):ntot,] <- fit$eta
	}
  }
}

class(bootmat) <- "boot.tsrq"
attr(bootmat, "tau") <- tau
if(class(object) == "tsrq") lambda <- object$lambda
if(class(object) == "tsrq2") lambda <- object$eta
attr(bootmat, "estimated") <- rbind(object$coefficients, lambda)
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "ntot") <- ntot
attr(bootmat, "indices") <- obsS
attr(bootmat, "rdf") <- object$rdf
attr(bootmat, "term.labels") <- nn

return(bootmat)

}


summary.boot.tsrq <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), which = NULL, ...){

tau <- attr(object, "tau")
nq <- length(tau)

est <- attr(object, "estimated")
ntot <- attr(object, "ntot")
rdf <- attr(object, "rdf")
R <- attr(object, "R")


nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound", "Pr(>|t|)")
		
if(nq == 1){
  sel <- complete.cases(object)
  R <- sum(sel)
  if(R < 2) stop("Insufficient boostrap sample")
  object <- as.matrix(object[sel,])
  bias <- est - apply(object, 2, mean)
  Cov <- cov(as.matrix(object))
  stds <- sqrt(diag(Cov))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2 * pt(-abs(est/stds), R - 1)
  ans <- cbind(est, bias, stds, lower, upper, tP)
  colnames(ans) <- nn
  rownames(ans) <- attr(object, "term.labels")
  printCoefmat(ans, digits = digits, signif.stars = TRUE, P.values = TRUE)
}
else {
  ans <- list()
  bias <- est - apply(object, 3, colMeans, na.rm = TRUE)
  Cov <- apply(object, 3, function(x) cov(as.matrix(x), use = "complete.obs"))
  if(ntot == 1) Cov <- matrix(Cov, nrow = 1)
  stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = ntot))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2*pt(-abs(est/stds), R - 1)
  sel <- (1:nq)
  if(!is.null(which)) sel <- sel[which]

  for(i in 1:nq){
    if(ntot == 1){
    ans[[i]] <- c(est[i], bias[i], stds[i], lower[i], upper[i], tP[i]);
    ans[[i]] <- matrix(ans[[i]], nrow = 1)
    } else {ans[[i]] <- cbind(est[,i], bias[,i], stds[,i], lower[,i], upper[,i], tP[,i])}
    colnames(ans[[i]]) <- nn
    rownames(ans[[i]]) <- attr(object, "term.labels")
    if(i %in% sel) cat(paste("tau = ", tau[i], "\n", sep =""))
    if(i %in% sel) printCoefmat(ans[[i]], digits = digits, signif.stars = TRUE, P.values = TRUE)
    if(i %in% sel) cat("\n")
  }

}
invisible(ans)
}


# 2-parameter transformations (mcjII)

tsrq2 <- function(formula, bounded = TRUE, lambda = NULL, delta = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", se = "nid", ...){

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

if(se == "rank") stop("Not implemented")
if(is.null(lambda)){
	lambda <- seq(0, 2, by = 0.005)
}
if(is.null(delta)){
	delta <- seq(0, 2, by = 0.005)
}

n <- length(y)
p <- ncol(x)
nq <- length(tau)
nl <- length(lambda)
nd <- length(delta)

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
			matLoss[i,k,j] <- L1loss(y - Fitted, tau = tau[j], weights = w)
			}
		}

	}
}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
parhat <- apply(matLoss, 3, function(x, lambda, delta){
m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];
return(c(lambda[m[1]], delta[m[2]]))}, lambda = lambda, delta = delta)

betahat <- ses <- matrix(NA, p, nq)
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau

fit <- list()

for(j in 1:nq){
# transform response with optimal lambda
data.tmp$newresponse <- mcjII(y, parhat[1,j], parhat[2,j], bounded, omega = 0.001)
fit[[j]] <- try(do.call(rq, args = list(formula = f.tmp, data = data.tmp, tau = tau[j], method = method, weights = w)), silent = TRUE)
	if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	ses[,j] <- summary(fit[[j]], se = se, ...)$coefficients[,2]
	Fitted[,j] <- invmcjII(x%*%matrix(betahat[,j]), parhat[1,j], parhat[2,j], bounded)
	}
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
dimnames(parhat) <- list(c("lambda","delta"), paste("tau =", format(round(tau, 3))))

fit$call <- call
fit$y <- y.old
if(bounded) fit$theta <- y
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$lambda.grid <- lambda
fit$delta.grid <- delta
fit$tsf <- "mcjII"
fit$bounded <- bounded
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 3, function(x){m <- which(x == min(x, na.rm = TRUE), arr.ind = TRUE)[1,];return(x[m[1],m[2]])})
fit$coefficients <- betahat
fit$std.error <- ses
fit$Fitted <- Fitted
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2
class(fit) <- "tsrq2"
return(fit)
}

predict.tsrq2 <- function(object, newdata, na.action = na.pass, raw = TRUE, ...){


tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
bounded <- object$bounded
etahat <- object$eta
betahat <- object$coefficients

if(missing(newdata)) {X <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

tmp <- X %*% betahat
Fitted <- matrix(NA, nrow(tmp), ncol(tmp))

if(!raw) return(tmp)

for(j in 1:nq){
	Fitted[,j] <- invmcjII(x = tmp[,j], lambda = etahat[1,j], delta = etahat[2,j], bounded = bounded)
}

if(bounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
}

return(Fitted)

}

boot.tsrq2 <- function(object, R = 50, seed = round(runif(1, 1, 10000))){

set.seed(seed)
tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
obsS <- replicate(R, sample(all.obs, size = n, replace = TRUE))
npars <- ncol(object$x)
ntot <- npars + 2
nn <- c(object$term.labels, "lambda", "delta")

if(nq == 1){
  bootmat <- matrix(NA, R, ntot);
  colnames(bootmat) <- nn
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i])
    w[match(names(a), all.obs)] <- as.numeric(a)
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars] <- fit$coefficients
		bootmat[i,(npars+1):ntot] <- fit$eta
	}
  }
} else {
  bootmat <- array(NA, dim = c(R, ntot, nq), dimnames = list(NULL, nn, paste("tau = ", format(tau, digits = 4), sep ="")));
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i]);
    w[match(names(a), all.obs)] <- as.numeric(a);
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars,] <- fit$coefficients
		bootmat[i,(npars+1):ntot,] <- fit$eta
	}
  }
}

class(bootmat) <- "boot.tsrq2"
attr(bootmat, "tau") <- tau
eta <- object$eta
attr(bootmat, "estimated") <- rbind(object$coefficients, object$eta)
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "ntot") <- ntot
attr(bootmat, "indices") <- obsS
attr(bootmat, "rdf") <- object$rdf
attr(bootmat, "term.labels") <- nn

return(bootmat)

}

summary.boot.tsrq2 <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), which = NULL, ...){

tau <- attr(object, "tau")
nq <- length(tau)

est <- attr(object, "estimated")
ntot <- attr(object, "ntot")
rdf <- attr(object, "rdf")
R <- attr(object, "R")


nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound", "Pr(>|t|)")
		
if(nq == 1){
  sel <- complete.cases(object)
  R <- sum(sel)
  if(R < 2) stop("Insufficient boostrap sample")
  object <- as.matrix(object[sel,])
  bias <- est - apply(object, 2, mean)
  Cov <- cov(as.matrix(object))
  stds <- sqrt(diag(Cov))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2 * pt(-abs(est/stds), R - 1)
  ans <- cbind(est, bias, stds, lower, upper, tP)
  colnames(ans) <- nn
  rownames(ans) <- attr(object, "term.labels")
  printCoefmat(ans, digits = digits, signif.stars = TRUE, P.values = TRUE)
}
else {
  ans <- list()
  bias <- est - apply(object, 3, colMeans, na.rm = TRUE)
  Cov <- apply(object, 3, function(x) cov(as.matrix(x), use = "complete.obs"))
  if(ntot == 1) Cov <- matrix(Cov, nrow = 1)
  stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = ntot))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2*pt(-abs(est/stds), R - 1)
  sel <- (1:nq)
  if(!is.null(which)) sel <- sel[which]

  for(i in 1:nq){
    if(ntot == 1){
    ans[[i]] <- c(est[i], bias[i], stds[i], lower[i], upper[i], tP[i]);
    ans[[i]] <- matrix(ans[[i]], nrow = 1)
    } else {ans[[i]] <- cbind(est[,i], bias[,i], stds[,i], lower[,i], upper[,i], tP[,i])}
    colnames(ans[[i]]) <- nn
    rownames(ans[[i]]) <- attr(object, "term.labels")
    if(i %in% sel) cat(paste("tau = ", tau[i], "\n", sep =""))
    if(i %in% sel) printCoefmat(ans[[i]], digits = digits, signif.stars = TRUE, P.values = TRUE)
    if(i %in% sel) cat("\n")
  }

}
invisible(ans)
}


## Nelder-Mead optimization for mcjII (joint estimation)

nlrq2 <- function(formula, par = NULL, bounded = TRUE, tau = 0.5, data, subset, weights, na.action, method = "fn", se = "nid", ...){

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
x <- model.matrix(mt, mf, contrasts)
w <- if (missing(weights)) rep(1, length(y)) else weights

if(bounded) y <- map(y)
if(se == "rank") stop("Not implemented")

n <- length(y)
p <- ncol(x)
nq <- length(tau)

if(is.null(par)) par <- rep(0, p + 2)

f <- function(theta, data){
	beta <- matrix(theta[-c(1:2)])
	if(theta[2] < 0) return(Inf)
	Fitted <- invmcjII(data$x %*% beta, lambda = theta[1], delta = theta[2], bounded = data$bounded)
	return(L1loss(data$y - Fitted, tau = data$tau, weights = data$weights))
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
fit$y <- y.old
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$eta <- parhat
fit$tsf <- "mcjII"
fit$bounded <- bounded
fit$coefficients <- betahat
fit$Fitted <- Fitted
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 2

class(fit) <- "nlrq2"
return(fit)
}


boot.nlrq2 <- function(object, R = 50, seed = round(runif(1, 1, 10000))){

set.seed(seed)
tau <- object$tau
nq <- length(tau)
all.obs <- rownames(object$x)
n <- length(all.obs)
obsS <- replicate(R, sample(all.obs, size = n, replace = TRUE))
npars <- ncol(object$x)
ntot <- npars + 2
nn <- c(object$term.labels, "lambda", "delta")

if(nq == 1){
  bootmat <- matrix(NA, R, ntot);
  colnames(bootmat) <- nn
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i])
    w[match(names(a), all.obs)] <- as.numeric(a)
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars] <- fit$coefficients
		bootmat[i,(npars+1):ntot] <- fit$eta
	}
  }
} else {
  bootmat <- array(NA, dim = c(R, ntot, nq), dimnames = list(NULL, nn, paste("tau = ", format(tau, digits = 4), sep ="")));
  for(i in 1:R){
	w <- rep(0, n)
    a <- table(obsS[,i]);
    w[match(names(a), all.obs)] <- as.numeric(a);
	fit <- try(my_update(object, weights = w), silent = TRUE)
    if(class(fit)!="try-error"){
		bootmat[i,1:npars,] <- fit$coefficients
		bootmat[i,(npars+1):ntot,] <- fit$eta
	}
  }
}

class(bootmat) <- "boot.nlrq2"
attr(bootmat, "tau") <- tau
eta <- object$eta
attr(bootmat, "estimated") <- rbind(object$coefficients, object$eta)
attr(bootmat, "R") <- R
attr(bootmat, "seed") <- seed
attr(bootmat, "npars") <- npars
attr(bootmat, "ntot") <- ntot
attr(bootmat, "indices") <- obsS
attr(bootmat, "rdf") <- object$rdf
attr(bootmat, "term.labels") <- nn

return(bootmat)

}


summary.boot.nlrq2 <- function(object, alpha = 0.05, digits = max(3, getOption("digits") - 3), which = NULL, ...){

tau <- attr(object, "tau")
nq <- length(tau)

est <- attr(object, "estimated")
ntot <- attr(object, "ntot")
rdf <- attr(object, "rdf")
R <- attr(object, "R")


nn <- c("Value", "Bias", "Std. Error", "Lower bound", "Upper bound", "Pr(>|t|)")
		
if(nq == 1){
  sel <- complete.cases(object)
  R <- sum(sel)
  if(R < 2) stop("Insufficient boostrap sample")
  object <- as.matrix(object[sel,])
  bias <- est - apply(object, 2, mean)
  Cov <- cov(as.matrix(object))
  stds <- sqrt(diag(Cov))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2 * pt(-abs(est/stds), R - 1)
  ans <- cbind(est, bias, stds, lower, upper, tP)
  colnames(ans) <- nn
  rownames(ans) <- attr(object, "term.labels")
  printCoefmat(ans, digits = digits, signif.stars = TRUE, P.values = TRUE)
}
else {
  ans <- list()
  bias <- est - apply(object, 3, colMeans, na.rm = TRUE)
  Cov <- apply(object, 3, function(x) cov(as.matrix(x), use = "complete.obs"))
  if(ntot == 1) Cov <- matrix(Cov, nrow = 1)
  stds <- sqrt(apply(Cov, 2, function(x, n) diag(matrix(x, n, n, byrow = TRUE)), n = ntot))
  lower <- est + qt(alpha/2, R - 1)*stds
  upper <- est + qt(1 - alpha/2, R - 1)*stds
  tP <- 2*pt(-abs(est/stds), R - 1)
  sel <- (1:nq)
  if(!is.null(which)) sel <- sel[which]

  for(i in 1:nq){
    if(ntot == 1){
    ans[[i]] <- c(est[i], bias[i], stds[i], lower[i], upper[i], tP[i]);
    ans[[i]] <- matrix(ans[[i]], nrow = 1)
    } else {ans[[i]] <- cbind(est[,i], bias[,i], stds[,i], lower[,i], upper[,i], tP[,i])}
    colnames(ans[[i]]) <- nn
    rownames(ans[[i]]) <- attr(object, "term.labels")
    if(i %in% sel) cat(paste("tau = ", tau[i], "\n", sep =""))
    if(i %in% sel) printCoefmat(ans[[i]], digits = digits, signif.stars = TRUE, P.values = TRUE)
    if(i %in% sel) cat("\n")
  }

}
invisible(ans)
}


predict.nlrq2 <- function(object, newdata, na.action = na.pass, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
bounded <- object$bounded
etahat <- object$eta
betahat <- object$coefficients

if(missing(newdata)) {X <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

tmp <- X %*% betahat
Fitted <- matrix(NA, nrow(tmp), ncol(tmp))

for(j in 1:nq){
	Fitted[,j] <- invmcjII(tmp[,j], lambda = etahat['lambda',j], delta = etahat['delta',j], bounded)
}

if(bounded){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
}

return(Fitted)

}


# Cusum process estimator (Mu and He, 2007) 1-parameter.

rcLoss <- function(lambda, x, y, tsf, symm = TRUE, tau = 0.5, method.rq = "fn"){

if(length(tau) > 1) stop("One quantile at a time")
n <- length(y)
theta <- map(y)
out <- rep(NA, n)

z <- switch(tsf,
	bc = bc(y, lambda),
	mcjI = mcjI(y, lambda, symm),
	ao = ao(theta, lambda, symm, omega = 0.001),
	mcjIb = mcjIb(theta, lambda, symm, omega = 0.001)
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

rcrq <- function(formula, tsf, symm = TRUE, lambda = NULL, tau = 0.5, data, subset, weights, na.action, method = "fn", se = "nid", ...){

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
w <- as.vector(model.weights(mf))
y <- model.response(mf)
x <- model.matrix(mt, mf, contrasts)
theta <- map(y)

if(se == "rank") stop("Not implemented")
if(is.null(lambda)){
	lambda <- if(tsf == "ao" & symm == FALSE) seq(-2, 2, by = 0.05)
		else seq(0, 2, by = 0.05)
}

n <- length(y)
p <- ncol(x)
nq <- length(tau)
nl <- length(lambda)

zhat <- res <- array(NA, dim = c(n, nq, nl))

matLoss <- rejected <- matrix(NA, nq, nl)
Ind <- array(NA, dim = c(n, nq, nl))

for(i in 1:nl){


# estimate linear QR for for sequence of lambdas

	for(j in 1:nq){
	matLoss[j,i] <- rcLoss(lambda[i], x, y, tsf, symm = symm, tau = tau[j], method.rq = method)
	}

}

if(all(is.na(matLoss))) return(list(call = call, y = y, x = x))

# minimise for lambda
lambdahat <- apply(matLoss, 1, function(x, lambda) lambda[which.min(x)], lambda = lambda)


betahat <- ses <- matrix(NA, p, nq)
colnames(betahat) <- tau
Fitted <- matrix(NA, n, nq)
colnames(Fitted) <- tau
fit <- list()

for(j in 1:nq){
# transform response with optimal lambda
newresponse <- switch(tsf,
	bc = bc(y, lambdahat[j]),
	mcjI = mcjI(y, lambdahat[j], symm),
	ao = ao(theta, lambdahat[j], symm, omega = 0.001),
	mcjIb = mcjIb(theta, lambdahat[j], symm, omega = 0.001)
	)

fit[[j]] <- try(rq(newresponse ~ x - 1, tau = tau[j], method = method, ...), silent = T)
	if(class(fit[[j]])!="try-error"){
	betahat[,j] <- coefficients(fit[[j]])
	ses[,j] <- summary(fit[[j]], se = se, ...)$coefficients[,2]
	tmp <- x%*%matrix(betahat[,j])
	Fitted[,j] <- switch(tsf,
		bc = invbc(tmp, lambdahat[j]),
		mcjI = invmcjI(tmp, lambdahat[j], symm),
		ao = invao(tmp, lambdahat[j], symm),
		mcjIb = invmcjIb(tmp, lambdahat[j], symm)
	)
	}

}

if(tsf %in% c("ao","mcjIb")){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(y))
}

dimnames(betahat) <- list(dimnames(x)[[2]], paste("tau =", format(round(tau, 3))))
names(lambdahat) <- paste("tau =", format(round(tau, 3)))

fit$call <- call
fit$method <- method
fit$y <- y
if(tsf %in% c("ao","mcjIb")) fit$theta <- theta
fit$x <- x
fit$weights <- w
fit$tau <- tau
fit$lambda <- lambdahat
fit$lambda.grid <- lambda
fit$tsf <- tsf
attr(fit$tsf, "symm") <- symm
fit$objective <- matLoss
fit$optimum <- apply(matLoss, 1, function(x) x[which.min(x)])
fit$coefficients <- betahat
fit$std.error <- ses
fit$Fitted <- Fitted
fit$rejected <- rejected
fit$terms <- mt
fit$term.labels <- colnames(x)
fit$rdf <- n - p - 1
class(fit) <- "rcrq"
return(fit)

}

predict.rcrq <- function(object, newdata, na.action = na.pass, raw = TRUE, ...){

tau <- object$tau
nq <- length(tau)
tsf <- object$tsf
symm <- attributes(tsf)$symm
lambdahat <- object$lambda
betahat <- object$coefficients

if(missing(newdata)) {X <- object$x} else {
	objt <- terms(object)
	Terms <- delete.response(objt)
	m <- model.frame(Terms, newdata, na.action = na.action, 
		xlev = object$levels)
	if (!is.null(cl <- attr(Terms, "dataClasses"))) 
		.checkMFClasses(cl, m)
	X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
}

tmp <- X %*% betahat
Fitted <- matrix(NA, nrow(tmp), ncol(tmp))

if(!raw) return(tmp)

for(j in 1:nq){
	Fitted[,j] <- switch(tsf,
		bc = invbc(tmp[,j], lambdahat[j]),
		mcjI = invmcjI(tmp[,j], lambdahat[j], symm),
		ao = invao(tmp[,j], lambdahat[j], symm),
		mcjIb = invmcjIb(tmp[,j], lambdahat[j], symm))
}

if(tsf %in% c("ao","mcjIb")){
	Fitted <- apply(Fitted, 2, function(x,x.r) invmap(x, x.r), x.r = range(object$y))
}

return(Fitted)

}


# Print functions for class tsrq, tsrq2, rcrq

print.tsrq <- print.rcrq <- function(x, ...){

    if (!is.null(cl <- x$call)) {
        cat("call:\n")
        dput(cl)
		cat("\n")
    }
	tsf <- switch(x$tsf,
		bc = "Box-Cox",
		mcjI = "Proposal I (positive)",
		ao = "Aranda-Ordaz",
		mcjIb = "Proposal I (bounded)")
	if(x$tsf %in% c("mcjI", "ao","mcjIb")) tsf <- paste(tsf, if(attr(x$tsf, "symm")) "symmetric" else "asymmetric")
	cat(tsf, "transformation\n")
    cat("\nOptimal transformation parameter:\n")
	print(x$lambda)
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

print.tsrq2 <- function(x, ...){

    if (!is.null(cl <- x$call)) {
        cat("call:\n")
        dput(cl)
		cat("\n")
    }
	tsf <- "Jones proposal II"
	tsf <- paste(tsf, if(x$bounded) "(bounded)" else "(positive)")
	cat(tsf, "transformation\n")
    cat("\nOptimal transformation parameter:\n")
	print(x$eta)
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

print.nlrq2 <- function(x, ...){

    if (!is.null(cl <- x$call)) {
        cat("call:\n")
        dput(cl)
		cat("\n")
    }
	tsf <- "Jones proposal II"
	tsf <- paste(tsf, if(x$bounded) "(bounded)" else "(positive)")
	cat(tsf, "transformation\n")
    cat("\nOptimal transformation parameter:\n")
	print(x$eta)
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

# From package boot


perc.ci <- function(x, conf = 0.95, hinv = function(x) x)
#
#  Bootstrap Percentile Confidence Interval Method
#
{
    alpha <- (1+c(-conf,conf))/2
    qq <- norm.inter(x,alpha)
    matrix(hinv(qq[,2]),ncol=2L)
}


norm.inter <- function(x,alpha)
#
#  Interpolation on the normal quantile scale.  For a non-integer
#  order statistic this function interpolates between the surrounding
#  order statistics using the normal quantile scale.  See equation
#  5.8 of Davison and Hinkley (1997)
#
{
    x <- x[is.finite(x)]
    R <- length(x)
    rk <- (R+1)*alpha
    if (!all(rk>1 & rk<R))
        warning("extreme order statistics used as endpoints")
    k <- trunc(rk)
    inds <- seq_along(k)
    out <- inds
    kvs <- k[k>0 & k<R]
    tstar <- sort(x, partial = sort(union(c(1, R), c(kvs, kvs+1))))
    ints <- (k == rk)
    if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
    out[k == 0] <- tstar[1L]
    out[k == R] <- tstar[R]
    not <- function(v) xor(rep(TRUE,length(v)),v)
    temp <- inds[not(ints) & k != 0 & k != R]
    temp1 <- qnorm(alpha[temp])
    temp2 <- qnorm(k[temp]/(R+1))
    temp3 <- qnorm((k[temp]+1)/(R+1))
    tk <- tstar[k[temp]]
    tk1 <- tstar[k[temp]+1L]
    out[temp] <- tk + (temp1-temp2)/(temp3-temp2)*(tk1 - tk)
    cbind(round(rk, 2), out)
}


boot.ci <- function(object, conf = 0.95){

n <- dim(object)
x <- attributes(object)$estimated

	if(length(n) == 2){
		val <- apply(object, 2, function(x, conf) perc.ci(x, conf = conf), conf = conf)
		val <- cbind(x, t(val))
	} else {
		val <- list()
		for(j in 1:n[3]){
			val[[j]] <- apply(object[,,j], 2, function(x, conf) perc.ci(x, conf = conf), conf = conf)
			val[[j]] <- cbind(x[,j], t(val[[j]]))
		}
			
	}

	return(val)

}

##################################################
### Restricted quantiles (He, 1997, AmStat; Zhao, 2000, JMA)
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
### Multiple imputation (Geraci, 2013, SMMR)
##################################################

mice.impute.rq <- function (y, ry, x, tsf = NULL, lambda = NULL, symm = TRUE, epsilon = 0.001, method.rq = "fn", ...) 
{
    x <- cbind(1, as.matrix(x))
	theta <- map(y)
	if(!is.null(tsf) && is.null(lambda))
		lambda <- 0

	z <- if(!is.null(tsf)){
		switch(tsf,
		bc = bc(y, lambda),
		mcjI = mcjI(y, lambda, symm),
		ao = ao(theta, lambda, symm, omega = 0.001),
		mcjIb = mcjIb(theta, lambda, symm, omega = 0.001)
		)
	} else y

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

	val <- if(!is.null(tsf)){
		switch(tsf,
		bc = invbc(ypred, lambda),
		mcjI = invmcjI(ypred, lambda, symm),
		ao = invao(ypred, lambda, symm),
		mcjIb = invmcjIb(ypred, lambda, symm)
		)
	} else ypred
	
	if(tsf %in% c("ao","mcjIb")){
		val <- invmap(val, range(y))
	}

    return(val)
}

# Impute using restricted quantiles

mice.impute.rrq <- function (y, ry, x, tsf = NULL, lambda = NULL, symm = TRUE, epsilon = 0.001, method.rq = "fn", ...) 
{
    x <- cbind(1, as.matrix(x))
	theta <- map(y)
	if(!is.null(tsf) & is.null(lambda))	lambda <- 0

	z <- if(!is.null(tsf)){
		switch(tsf,
		bc = bc(y, lambda),
		mcjI = mcjI(y, lambda, symm),
		ao = ao(theta, lambda, symm, omega = 0.001),
		mcjIb = mcjIb(theta, lambda, symm, omega = 0.001)
		)
	} else y

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
	
	val <- if(!is.null(tsf)){
		switch(tsf,
		bc = invbc(ypred, lambda),
		mcjI = invmcjI(ypred, lambda, symm),
		ao = invao(ypred, lambda, symm),
		mcjIb = invmcjIb(ypred, lambda, symm)
		)
	} else ypred
	
	if(tsf %in% c("ao","mcjIb")){
		val <- invmap(val, range(y))
	}

    return(val)
}


##################################################
### QR for counts (Machado and Santos Silva, 2005, JASA)
##################################################

addnoise <- function(x, centered = TRUE, B = 0.999) 
{

	n <- length(x)
    if (centered) 
        z <- x + runif(n, -B/2, B/2)
    else z <- x + runif(n, 0, B)
	
    return(z)
}

F.rq <- function(x, cn){

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

rq.counts <- function (formula, data, weights = NULL, offset = NULL, contrasts = NULL, tau = 0.5, M = 50, zeta = 1e-5, B = 0.999, cn = NULL, alpha = 0.05, method = "fn") 
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

    y <-  model.response(mf, "numeric")
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

	# Add noise
	Z <- replicate(M, addnoise(y, centered = FALSE, B = B))
	# Transform Z
    TZ <- apply(Z, 2, function(x, off, tau, zeta) log(ifelse((x - 
        tau) > zeta, x - tau, zeta)) - off, off = offset, tau = tau, zeta = zeta)
	# Fit linear QR on TZ
    fit <- apply(TZ, 2, function(y, x, weights, tau, method) 
        rq.wfit(x = x, y = y, tau = tau, weights = weights, method = method), x = x, tau = tau, weights = w, method = method)
	# Trasform back
    yhat <- sapply(fit, function(obj, x) x %*% obj$coefficients, x = x)
	yhat <- as.matrix(yhat)
    eta <- sweep(yhat, 1, offset, "+")
    zhat <- tau + exp(eta)
	#
    Fvec <- Vectorize(F.rq)
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
    
    multiplier <- exp(eta) * (F <= Z & Z < Fp)
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
		V <- V/n
		stds <- sqrt(diag(V))
		} else {
		stds <- NA
        warning("Standard error not available")
		}

    est <- sapply(fit, function(x) x$coefficients)
    est <- if (p == 1) mean(est) else rowMeans(est)

    qfit <- if (p == 1) {
        tau + exp(mean(eta[1, ]))
    } else {
        tau + exp(rowMeans(eta))
    }

	lower <- est + qt(alpha/2, n - p) * stds
	upper <- est + qt(1 - alpha/2, n - p) * stds
	tP <- 2 * pt(-abs(est/stds), n - p)

	ans <- cbind(est, stds, lower, upper, tP)
	colnames(ans) <- c("Value", "Std. Error", "lower bound", "upper bound", 
        "Pr(>|t|)")
	rownames(ans) <- names(est) <- term.labels
	
	fit <- list()
	fit$call <- call
	fit$na.action <- attr(mf, "na.action")
	fit$contrasts <- attr(x, "contrasts")
	fit$term.labels <- term.labels
	fit$terms <- mt

	fit$coefficients <- est
	fit$tau <- tau
	fit$nobs <- n
	fit$M <- M
	fit$Mn <- m.n
	fit$rdf <- n - p
	fit$x <- x
	fit$y <- y
	fit$fitted <- qfit
	fit$offset <- offset
	fit$Cov <- V
	fit$tTable <- ans
	fit$levels <- .getXlevels(mt, mf)
	fit$method <- method

	class(fit) <- "rq.counts"
	
    return(fit)
}

coef.rq.counts <- function(object, ...){

tau <- object$tau
nq <- length(tau)
ans <- object$coefficients

if(nq == 1){
  names(ans) <- object$term.labels
}

return(ans)

}

predict.rq.counts <- function(object, newdata, na.action = na.pass, ...) 
{

tau <- object$tau

	if(missing(newdata)){
		yhat <- drop(object$x %*% object$coefficients)
	}
	else {
		objt <- terms(object)
		Terms <- delete.response(objt)
		m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$levels)
		if (!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
		x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
		yhat <- drop(x %*% object$coefficients)
		
	}

return(yhat)
}

residuals.rq.counts <- function(object, ...){

ans <- as.numeric(object$y) - predict(object)
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


##################################################
### Khmaladze test (Koenker, 2005, Appendix)
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


