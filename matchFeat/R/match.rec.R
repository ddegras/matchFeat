match.rec <- function(x, unit = NULL, w = NULL, control = list())
{

	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x, unit, w=w)
	m <- pre$m; n <- pre$n; p <- pre$p
	R <- pre$R # Cholesky decomposition of w or null
	if (!is.null(pre$x))
	{ x <- pre$x; dim(x) <- c(p,m,n) }
	rm(pre)
	syscall <- sys.call()
	
	## Tuning parameters
	equal.variance <- FALSE
	if (is.list(control) && !is.null(control$equal.variance)) 
		equal.variance <- control$equal.variance		

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

	## Ensure that data are non-negative
	xmin <- min(x)
	if (xmin < 0) x <- x - xmin
	
	## Match feature vectors	
	sigma <- matrix(,m,n)
	sigma[,1] <- 1:m
	sumxP <- x[,sigma[,1],1]	
	for (i in 2:n) {
		sigma[,i] <- solve_LSAP(crossprod(sumxP,x[,,i]), 
			maximum = TRUE)
		sumxP <- sumxP + x[,sigma[,i],i]
	}			
	
	## Objective value
	objective <- (sum(x^2) - (sum(sumxP^2)/n)) / (n-1)

	## Means and covariances of matched vectors
	# mu <- matrix(,p,m)
	mu <- sumxP / n
	V <- array(dim=c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0, by = m, length.out = n) + sigma[l,]
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
 			V <- t(backsolve(R,V))
			dim(V) <- c(p,p,m)	
		}		
 	}
 	
 	## Cluster assignment
 	cluster <- matrix(,m,n)
 	for (i in 1:n)
 		cluster[sigma[,i],i] <- 1:m
 		
	out <- list(sigma = sigma, cluster = cluster, 
		objective = objective, 
		mu = mu, V = V, call = syscall)
	class(out) <- "matchFeat"
	return(out)
	
}
