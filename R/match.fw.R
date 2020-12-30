match.fw <- function(x, unit = NULL, w = NULL, control = list())
{

	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- matchFeat:::preprocess(x=x, unit=unit, w=w)
	m <- pre$m; n <- pre$n; p <- pre$p
	if (!is.null(w)) w <- pre$w
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
 
	## Trivial cases 
	if (m == 1 || n == 1 || p == 1)
		return(matchFeat:::trivial(x,m,n,p,w,R,equal.variance,syscall))
	
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

	## Ensure that data are non-negative
	xmin <- min(x)
	if (xmin < 0) x <- x - xmin

	## Change notations	
	P <- sigma

	## Initialize objective	
	sumxP <- matrix(0,p,m)
	for (i in 1:n)
		sumxP <- sumxP + x[,P[,i],i]
	objective <- sum(sumxP^2)
			
	# Search direction
	Q <- matrix(,m,n)
	
	for (it in 1:maxit)
	{	
		# Store previous objective
		objective.old <- objective
					
		# Find next candidate solution 
		sumxQ <- matrix(0,p,m)
		for (i in 1:n) {
			grad <- crossprod(sumxP,x[,,i])
			Q[,i] <- solve_LSAP(grad, maximum=TRUE)
			sumxQ <- sumxQ + x[,Q[,i],i] 
		}

		# Compare to previous solution
		objective <- sum(sumxQ^2)
		if (objective > objective.old) {
			P <- Q
			sumxP <- sumxQ
		} else break
				
	}

	cost <- (sum(x^2) - (objective/n)) / (n-1)

	## Sample means and covariances of matched vectors
	mu <- sumxP/n
	V <- array(,c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0,by=m,len=n) + sigma[l,]
		V[,,l] <- tcrossprod(x[,idx])/n - tcrossprod(mu[,l])	
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
 		cluster[P[,i],i] <- 1:m
		
	out <- list(sigma=P, cluster=cluster, 
		objective=cost, mu=mu, V=V, call=syscall)
	class(out) <- "matchFeat"
	return(out)
	
}
