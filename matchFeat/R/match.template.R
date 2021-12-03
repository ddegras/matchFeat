match.template <- function(x, template = 1L, unit = NULL, w = NULL, 
	method = c("hungarian","bruteforce"), equal.variance = FALSE)
{
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x, unit, w=w)
	m <- pre$m; n <- pre$n; p <- pre$p
	R <- pre$R # Cholesky decomposition of w or null
	if (!is.null(pre$x)) 
		{ x <- pre$x; dim(x) <- c(p,m,n) }
	tmpl <- template
	if (length(tmpl) == 1) { 
		tmpl <- x[,,tmpl] 
	} else {
		stopifnot(NROW(tmpl) == p && NCOL(tmpl) == m)	
	}
	rm(pre)
	syscall <- sys.call()
	
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

	## Ensure that data and template are non-negative
	xtmin <- min(x,tmpl)
	if (xtmin < 0) {
		x <- x - xtmin
		tmpl <- tmpl - xtmin }
		
	## Assignment method
	method <- match.arg(method)
	if (method == "bruteforce") {
		perms <- permn(m)
		perms <- lapply(1:NROW(perms), 
			function(i) perms[i,]) }
		
	## Match feature vectors 
	sigma <- matrix(,m,n)
	if (method == "bruteforce") {
		for (i in 1:n)
			sigma[,i] <- brute(x[,,i],tmpl,perms)
	} else { 
		for (i in 1:n)
			sigma[,i] <- solve_LSAP(crossprod(tmpl,x[,,i]), 
				maximum = TRUE) 
	}
				
	## Sample means and covariances of matched vectors
	mu <- matrix(,p,m)
	V <- array(,c(p,p,m))
	dim(x) <- c(p,m*n)
	cost <- numeric(m)
	for (l in 1:m) {
		idx <- seq.int(0 ,by = m, length.out = n) + sigma[l,]
		mu[,l] <- rowMeans(x[,idx,drop=F])
		V[,,l] <- tcrossprod(x[,idx,drop=F])/n - 
			tcrossprod(mu[,l])	
		cost[l] <- sum(diag(V[,,l]))
	}	
	objective <- sum(cost) / (n-1) # sum of within-cluster variances
	if (xtmin < 0) mu <- mu + xtmin
	if (equal.variance) V <- apply(V,1:2,mean)
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
		objective = objective, mu = mu, V = V, call = syscall)
	class(out) <- "matchFeat"
	return(out)	
}
