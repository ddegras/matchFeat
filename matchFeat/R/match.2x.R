match.2x <- function(x, sigma=NULL, unit=NULL, w=NULL, control=list())
{
	test <- requireNamespace("gurobi", quietly = TRUE)
	if (!test) stop(paste(
		"The software 'gurobi' and its R interface package",
		"are required to run the function 'match.2x'"))
	pre <- preprocess(x,unit)
	m <- pre$m; n <- pre$n; p <- pre$p
	if (!is.null(w)) w <- pre$w
	R <- pre$R
	if (!is.null(pre$x)) 
		{ x <- pre$x } else { dim(x) <- c(p,m*n) }
	rm(pre)
	syscall <- sys.call()
	
	## Tuning parameters
	if (is.null(sigma)) sigma <- matrix(1:m,m,n)
	maxit <- 1000L
	equal.variance <- FALSE
	timeout <- NULL
	if (is.list(control)) {
		if (!is.null(control$sigma)) sigma <- control$sigma
		if (!is.null(control$maxit)) maxit <- control$maxit
		if (!is.null(control$equal.variance)) 
			equal.variance <- control$equal.variance	
		if (!is.null(control$timeout))
			timeout <- control$timeout
	}
 
	## Trivial cases 
	if (m == 1 || n == 1 || p == 1)
		return(trivial(x,m,n,p,w,R,equal.variance,syscall))
	
	## Rescale data if required
	if (!is.null(w)) {
		if (is.vector(w)) {
			x <- sqrt(w) * x
		} else {
			x <- R %*% x
		}
	}


	## Model parameters	
	model <- list()
	model$modelsense <- "max"
	model$vtype <- "B"
	params <- list(OutputFlag=0)
	params$timeLimit <- timeout
	model$A <- matrix(0,1,n)
	active <- 1:m
	a1 <- NULL
	shift <- seq(0,by=m,length.out=n)
	while(length(active) > 1) {
		a1 <- sample(setdiff(active,a1),1) 
		# prevent choosing same candidate twice in a row
		cand <- setdiff(active,a1) # candidate set
		ncand <- length(cand)
		bestgain <- 0
		for (k in 1:ncand) {
			a2 <- cand[k]
			dx <- x[,sigma[a1,]+shift] - x[,sigma[a2,]+shift]
			dbar <- rowMeans(dx)
			model$Q <- crossprod(dx)
			model$obj <- -n*crossprod(dx,dbar)
			result <- gurobi::gurobi(model,params)
			nnz <- sum(result$x) # number of 1's in (binary) solution
			if (result$objval > bestgain && nnz > 0 && nnz < n) {
				cc <- result$x
				bestgain <- result$objval
				a2s <- a2 # best candidate
			}	
		}
		if (bestgain > 0) {
			active <- 1:m
			cc <- as.logical(cc)
			tmp <- sigma[a1,cc]
			sigma[a1,cc] <- sigma[a2s,cc]
			sigma[a2s,cc] <- tmp
		} else { 
			active <- setdiff(active,a1)
		}
	}
	
	## Calculate objective
	mu <- matrix(,p,m)
	for (l in 1:m) 
		mu[,l] <- rowMeans(x[,sigma[l,]+shift])
	cost <- (sum(x^2) - n * sum(mu^2)) / (n-1)

	## Sample means and covariances of matched vectors
	V <- array(,c(p,p,m))
	for (l in 1:m) {
		idx <- seq.int(0,by=m,length.out=n) + sigma[l,]
		V[,,l] <- tcrossprod(x[,idx])/n - tcrossprod(mu[,l])	
	}	
	if (equal.variance) 
		V <- rowMeans(V,dims=2L) # apply(V,1:2,mean)
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

