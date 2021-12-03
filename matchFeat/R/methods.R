print.matchFeat <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	cat("Dimensions: units = ", ncol(x$sigma),
		", classes = ", ncol(x$mu), 
		", variables = ", nrow(x$mu), "\n", sep = "")
	cat("Objective value =", x$objective, "\n")
	invisible(NULL)	
}



summary.matchFeat <- function(object, ...)
{
	cat("Call:\n")
	print(object$call)
	n <- ncol(object$sigma)
	cat("Dimensions: units = ", ncol(object$sigma),
		", classes = ", ncol(object$mu), 
		", variables = ", nrow(object$mu), "\n", sep = "")
	cat("Objective value =", object$objective, "\n")
	mu <- object$mu
	ssb <- n * sum((mu-rowMeans(mu))^2)
	ssw <- n * sum(apply(object$V,3,diag))
	cat("Sum of squares between clusters =",ssb,"\n")
	cat("Sum of squares within clusters =",ssw,"\n")
	invisible(NULL)	
}

predict.matchFeat <- function(object, newdata, unit = NULL, ...)
{
	stopifnot(is(object,"matchFeat"))

	## Prediction
	if (is.element("match.gaussmix", as.character(object$call))) {
		out <- match.gaussmix(x = newdata, unit = unit,
			mu = object$mu, V = object$V, fixed = TRUE)
	} else {
		out <- match.template(x = newdata, template = object$mu, 
			unit = unit)
		out$mu <- object$mu
		out$V <- object$V
	}
	out$call <- sys.call()
	return(out)
	
}


