print.matchFeats <- function(obj)
{
	cat("Call:\n")
	print(obj$call)
	cat("Dimensions: units =",ncol(obj$sigma),
		"classes =",ncol(obj$mu),"variables =",nrow(obj$mu),"\n")
	cat("Objective value =",obj$cost,"\n")
	invisible(NULL)	
}



summary.matchFeats <- function(obj)
{
	cat("Call:\n")
	print(obj$call)
	n <- ncol(obj$sigma)
	cat("Dimensions: units =",n,
		"classes =",ncol(obj$mu),"variables =",nrow(obj$mu),"\n")
	cat("Objective value =",obj$cost,"\n")
	mu <- obj$mu
	ssb.matched <- n * sum((mu-rowMeans(mu))^2)
	ssw.matched <- n * sum(apply(obj$V,3,diag))
	tab <- matrix(c(obj$ss.between.unmatched, obj$ss.within.unmatched, 
		ssb.matched, ssw.matched),2,2)
	rownames(tab) <- c("Sum of Squares Between", "Sum of Squares Within")
	colnames(tab) <- c("Unmatched","Matched")
	print(tab)
	invisible(NULL)	
}


