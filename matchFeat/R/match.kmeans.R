match.kmeans <- function(x, unit = NULL, w = NULL, 
	method=c("hungarian", "bruteforce"), control = list())
{
	
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x, unit)
	m <- pre$m; n <- pre$n; p <- pre$p
	R <- pre$R # Cholesky decomposition of w or null
	if (!is.null(pre$x)) 
		{ x <- pre$x; dim(x) <- c(p,m,n) }
	rm(pre)
	syscall <- sys.call()
			
	## Tuning parameters
	sigma <- matrix(1:m,m,n)
	maxit <- 1000L
	equal.variance <- FALSE
	if (is.list(control)) {
		if (!is.null(control$sigma)) sigma <- control$sigma
		if (!is.null(control$maxit)) maxit <- control$maxit
		if (!is.null(control$equal.variance)) 
			equal.variance <- control$equal.variance
	}
	method <- match.arg(method)	
		
	## Trivial cases 
	if (m == 1 || n == 1 || p == 1)
		return(trivial(x,m,n,p,w,R,equal.variance,syscall))

	## Rescale data if required
	if (!is.null(w)) {
		if (is.vector(w)) {
			x <- sqrt(w) * x
		} else {
			dim(x) <- c(p,m*n)
			x <- R %*% x
			dim(x) <- c(p,m,n)
		}
	}
	
	## Ensure that data values are non-negative
	xmin <- min(x)
	if (xmin < 0) x <- x - xmin

	# Assignment method
	method <- match.arg(method)
	if (method == "bruteforce") {
		# Generate permutations
		perms <- permn(m)
		# Convert to list
		perms <- lapply(1:NROW(perms), 
			function(i) perms[i,])
	}
	
	## Initial centers 
	xbar <- matrix(0,p,m)
	for (i in 1:n) 
		xbar <- xbar + x[,sigma[,i],i]
	xbar <- xbar / n

	## Initial cost
	sumx2 <- sum(x^2)
	cost <- (sumx2 - n * sum(xbar^2)) / (n-1)
	
	for (count in 1:maxit)
	{
		# Store previous cost
		cost.old <- cost
		
		# New assignment
		for (i in 1:n)
			sigma[,i] <- if (method == "brute") {
				brute(x[,,i],xbar,perms)
		} else { solve_LSAP(crossprod(xbar,x[,,i]), 
			maximum = TRUE)}
			
		# New centers
		xbar <- matrix(0,p,m)
		for (i in 1:n) 
			xbar <- xbar + x[,sigma[,i],i]
		xbar <- xbar / 	n
		
		# New cost 
		cost <- (sumx2 - n * sum(xbar^2)) / (n-1)
		
		# Terminate if no reduction in cost
		if (cost == cost.old) break
	}	

	## Sample means and covariances of matched feature vectors
	mu <- xbar
	V <- array(,c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0,by=m,length.out=n) + sigma[l,]
		V[,,l] <- tcrossprod(x[,idx,drop=F])/n - 
			tcrossprod(mu[,l])	
	}	
	if (xmin < 0) mu <- mu + xmin
	if (equal.variance) V <- rowMeans(V,dims=2L) # apply(V,1:2,mean)
 	if (!is.null(w)) {
 		if (is.vector(w)) {
 			mu <- mu / sqrt(w)
 			V <- V / w
 		} else {
 			mu <- backsolve(R,mu)
 			dim(V) <- c(p,m*p)
 			V <- t(backsolve(R,V))
			dim(V) <- c(p,p,m)
		}
 	}

 	## Cluster assignment
 	cluster <- matrix(,m,n)
 	for (i in 1:n)
 		cluster[sigma[,i],i] <- 1:m

	out <- list(sigma = sigma, cluster = cluster, 
	objective = cost, mu = mu, V = V, call = syscall)
	class(out) <- "matchFeat"
	return(out)
	
}
