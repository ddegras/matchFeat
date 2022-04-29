#########################
# Preprocessing function 
# for input arguments
#########################


preprocess <- function(x, unit, mu = NULL, V = NULL, w = NULL)
{

	out <- list(x = NULL, m = NULL, n = NULL, p = NULL, mu = mu, V = V,
		R = NULL)	
	
	if (any(is.infinite(x) | is.na(x)))
		stop(paste("Please make sure that 'x' does not contain missing",
		"or infinite values"))

	dimx <- dim(x)
	
	if (length(dimx) < 3) {
	## Case: x specified as matrix
	if (is.vector(x)) dim(x) <- c(length(x),1)
	# if (!is.matrix(x)) x <- as.matrix(x)
		# Check that 'unit' is not NULL 
		if (is.null(unit)) 
			stop(paste("Please provide 'x' as an array of dimensions",
				"(p,m,n) OR provide 'x' as a matrix of dimensions (m*n,p)\n",
				"together with the vector 'unit' (length m*n) where m,n,p",
				"are the respective numbers of feature vectors, units,",
				"and variables"))
		unit <- as.integer(unit)
		if (length(unit) == 1) {
			n <- unit
			m <- nrow(x) / n
			if (m != floor(m)) 
				stop(paste("Please make sure that the number of rows in 'x'",
					"is a multiple of 'n'"))			
		} else {
			# Check that dimensions of x and 'unit' are compatible  
			if (length(unit) != nrow(x))
				stop(paste("Please make sure that the number of rows in 'x'",
					"and the length  of 'unit' are equal"))
			# Check that units are balanced 
			# (i.e. have same number of feature vectors)
			if (!all(diff(unit) >= 0))
				stop(paste("Please specify 'unit' as a vector of nondecreasing",
					"integers and sort the rows of 'x' accordingly"))
			freq <- table(unit)
			if (!all(freq == freq[1]))
				stop(paste("Unbalanced units not supported yet.",
					"Make sure that all unique values in 'unit'",
					"have the same frequency"))
			# Check dimensions of mu and V if specified
			m <- as.integer(freq[1])
			n <- length(freq)
		}
		p <- ncol(x)
		out$x <- t(x)
	} else {		
	## Case: x specified as array	
		if (length(dimx) > 3)
			stop(paste("Please provide 'x' as an array of dimensions (p,m,n)",
				"OR provide 'x' as a matrix of dimensions (n*p,m)\n",
				"together with the vector 'unit' (length n*p) where m,n,p",
				"are the respective numbers of feature vectors, units,",
				"and variables"))
		p <- dimx[1]
		m <- dimx[2]
		n <- dimx[3]
		# dim(x) <- c(p,m*n)		
		# if (any(is.infinite(x) | is.na(x)))
			# stop(paste("Please make sure that 'x' does not contain",
			 	# "missing or infinite values"))
	}
	out$m <- m; out$n <- n; out$p <- p

	## Check mu and V if provided
	if (!is.null(mu) || !is.null(V)) {
		if (is.null(mu) || is.null(V))	
			stop("Please provide both arguments 'mu' and 'V' or neither")
		if (!((is.vector(mu) && length(mu) == p) || 
			(is.matrix(mu) && dim(mu)[1] == p && dim(mu)[2] %in% c(1,m))))
			stop(paste("Please specify 'mu' as a vector of length p",
				"or as a matrix of dimensions (p,m)\n where p is the number", 
				"of variables and m the number of feature vectors"))
		# Reshape and/or recycle mu if needed
		if (!identical(dim(mu),c(p,m))) mu <- matrix(mu,p,m)
		if (!((length(V) == 1 && p == 1) || 
			(is.matrix(V) && all(dim(V) == c(p,p))) || 
			(length(dim(V)) == 3 && all(dim(V) == c(p,p,m)))))
			stop(paste("Please specify 'V' as a matrix of dimensions (p,p)",
				"or as an array of dimensions (p,p,m)\n where p is the number", 
				"of variables and m the number of feature vectors"))
		# Reshape and/or recycle V if needed 
		if (!(length(dim(V) == 3 && all(dim(V) == c(p,p,m))))) 
			V <- array(V,c(p,p,m))
	}	
	
	## Check w if provided
	if (!is.null(w)) {
		err <- paste("Please make sure that 'w' is a vector of", 
			"positive numbers or a positive definite matrix")
		if (is.vector(w)) {
			if (any(w <= 0) || !is.element(length(w),c(1,p))) stop(err)
		} else if (is.matrix(w)) {
			out$R <- tryCatch(chol(w), error=function(e) stop(err))
		} else stop(err)
	}
		
	return(out)
}






################################
# Function for feature matching 
# in trivial cases
################################
 

trivial <- function(x, m, n, p, w, R, equal.variance, syscall)
{

	stopifnot(m == 1 || n == 1 || p == 1)
	
	## Case: only one unit
	if (n == 1) {
		dim(x) <- c(p,m)
		xbar <- rowMeans(x)
		if (is.null(w)) {
			ssb <- sum((x-xbar)^2)
		} else {
			ssb <- if (is.vector(w)) {
				rowSums((x-xbar)^2) * w
			} else { sum((R %*% (x-xbar))^2) }
		} 		
		out <- list(sigma=matrix(1:m,m,1), objective=0,
			mu=x, V=array(0,c(p,p,m)), call=syscall)
			# ss.between.unmatched=ssb, 
			# ss.within.unmatched=0,

	## Case: only one class/feature 
	} else if (m == 1) {
		dim(x) <- c(p,n)
		xbar <- rowMeans(x)
		if (is.null(w)) {
			ssw <- sum((x-xbar)^2)
		} else {
			ssw <- if (is.vector(w)) {
				rowSums((x-xbar)^2) * w
			} else { sum((R %*% (x-xbar))^2) }
		} 
		out <- list(sigma = matrix(1,1,n), cost = ssw/(n-1),
			mu = xbar, V = array(tcrossprod(x-xbar)/n, c(p,p,1)), call=syscall)
						
			
	## Case: only one variable 
	} else {				
		dim(x) <- c(m,n)
		xbar <- as.vector(rowMeans(x))
		if (is.null(w)) w <- 1
		sigma <- apply(x,2,order)
		x <- x[cbind(as.vector(sigma),rep(1:n,each=m))]
		dim(x) <- c(m,n)
		mu <- rowMeans(x)
		V <- rowMeans(x^2) - mu^2
		cost <- n/(n-1) * sum(V) 
		if (equal.variance) V <- rep(mean(V),m)
		out <- list(sigma=sigma, objective=cost, mu=mu,
			V=V, call=syscall)
	}	

	class(out) <- "matchFeat"	
	return(out)		
}




#####################
# Objective function 
# (balanced data)
#####################


objective.fun <- function(x, sigma = NULL, unit = NULL, w = NULL)
{
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x=x, unit=unit, w=w)
	m <- pre$m; n <- pre$n; p <- pre$p
	if (!is.null(w))
		{ w <- pre$w; R <- pre$R }
	if (!is.null(pre$x)) 
		{ x <- pre$x } else { dim(x) <- c(p,m*n) }
	rm(pre)	
	if (is.null(sigma)) 
		sigma <- matrix(1:m,m,n)
 
	## Rescale data if required
	if (!is.null(w)) {
		if (is.vector(w)) {
			x <- sqrt(w) * x
		} else {
			x <- R %*% x
		}
	}

	## Calculate cost of assignment	
	mu <- matrix(0,p,m)
	shift <- seq.int(0, by=m, length.out=n)
	for (l in 1:m)
		mu[,l] <- rowMeans(x[,sigma[l,]+shift,drop=FALSE])
	objective <- (sum(x^2) - n * sum(mu^2)) / (n-1)
	return(objective)
	
}





#######################################
# Objective function for general case:
# balanced or unbalanced data
#######################################


objective.gen.fun <- function(x, unit, cluster)
{
	n <- length(unique(unit))
	u <- unique(cluster)
	nclass <- length(u)
	obj <- numeric(nclass)
	nrmx2 <- rowSums(x^2)
	for (k in 1:nclass) {
		idx <- which(cluster == u[k])
		nk <- length(idx)
		if (nk > 0) 
			obj[k] <- nk * sum(nrmx2[idx]) - 
				sum(colSums(x[idx,,drop=FALSE])^2)
	}
	obj <- sum(obj) / (2 * n * (n-1))
	return(obj)		
}





##############################
# Rand index between two 
# partitions of a set
##############################


Rand.index <- function(x,y)
{
	stopifnot(length(x) == length(y))
	n <- length(x)
	freq <- table(x,y)
	nx <- rowSums(freq)
	ny <- colSums(freq)
	1 - sum(nx^2, ny^2,-2*freq^2) / (n * (n-1))
	# 1 - sum(nx*(nx-1), ny*(ny-1),-2*freq * (freq-1)) / (n * (n-1))
}

