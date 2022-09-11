match.bca <- function(x, unit = NULL, w = NULL, 
	method = c("cyclical","random"), control = list())
{

	## Check inputs
	if (is.matrix(x)) {
		if (is.null(unit))
			stop("If 'x' is a matrix, 'unit' should also be provided.")
		freq <- table(unit)
		if (any(freq != freq[1])) 
			stop(paste("The data are unbalanced.",
				"Please use function 'match.bca.gen' instead"))
	}
	
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x = x, unit = unit, w = w)
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

	## Ensure that data are non-negative
	xmin <- min(x)
	if (xmin < 0) x <- x - xmin

	## Sweeping method
	method <- match.arg(method)			
	
	## Initialize objective	
	sumxP <- matrix(0,p,m)
	for (i in 1:n)
		sumxP <- sumxP + x[,sigma[,i],i]
	objective <- sum(sumxP^2)
		
	## Define sweeping order if method = cyclical
	if (method == "cyclical") 
		sweep <- 1:n 
	
	for (count in 1:maxit)
	{
		## Store previous objective
		objective.old <- objective
		
		## Randomize sweeping order if method = random
		if (method == "random")
			sweep <- sample(1:n)
					
		for (i in sweep)
		{
			sumxPi <- sumxP - x[,sigma[,i],i]  
			sigma[,i] <- solve_LSAP(crossprod(sumxPi,x[,,i]), 
					maximum = TRUE)
			sumxP <- sumxPi + x[,sigma[,i],i]
		}			
		
		## Periodically recalculate X1 P1 + ... + Xn Pn
		## to contain roundoff errors
		if (count %% 10 == 0) {
			sumxP <- matrix(0,p,m)
			for (i in 1:n)
				sumxP <- sumxP + x[,sigma[,i],i]
		}
	
		## Update objective 
		objective <- sum(sumxP^2)
		
		## Terminate if no improvement in objective
		if (objective <= objective.old) break
		
	}	
	
	cost <- (sum(x^2) - (objective/n)) / (n-1)

	## Sample means and covariances of matched vectors
	# mu <- matrix(,p,m)
	mu <- sumxP/n
	V <- array(,c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0,by=m,length.out=n) + sigma[l,]
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
 		cluster[sigma[,i],i] <- 1:m
		
	out <- list(sigma = sigma, cluster = cluster, 
		objective = cost, mu = mu, V = V, call = syscall)
	class(out) <- "matchFeat"
	return(out)
	
}
