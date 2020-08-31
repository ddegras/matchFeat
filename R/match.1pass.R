match.1pass <- function(x, unit = NULL, w = NULL, control = list())
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

	# ## Trivial case p == 1 
	# if (p == 1) {
		# dim(x) <- c(m,n)
		# sigma <- apply(x,2,order)
		# x <- apply(x,2,sort)
		# mu <- rowMeans(x)
		# V <- rowMeans(x^2) - mu^2
		# cost <- if (n>1) { n/(n-1) * sum(V) } else 0
		# if (equal.variance) V <- rep(mean(V),m)
		# out <- list(sigma=sigma, cost=cost, mu=mu,
			# V=V, ss.between.unmatched=ssb, 
			# ss.within.unmatched=ssw, call=syscall)
		# class(out) <- "matchFeats"
		# return(out)
	# }	
	
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

 	## Sums of squares for unmatched data
 	mu <- rowMeans(x, dims=2L)
 	ssw <- sum((x-as.vector(mu))^2)
 	ssb <- n * sum((mu-rowMeans(mu))^2)

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
	cost <- (sum(x^2) - (sum(sumxP^2)/n)) / (n-1)

	## Sample means and covariances of matched vectors
	# mu <- matrix(,p,m)
	mu <- sumxP / n
	V <- array(dim=c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0,by=m,len=n) + sigma[l,]
		# mu[,l] <- rowMeans(x[,idx,drop=F])
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
 		
	out <- list(sigma=sigma, cost=cost, mu=mu, V=V, 
		ss.between.unmatched=ssb, ss.within.unmatched=ssw,
		call=syscall)
	class(out) <- "matchFeats"
	return(out)
	
}
