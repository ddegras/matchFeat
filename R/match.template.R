match.template <- function(x, y = 1, unit = NULL, w = NULL, 
	equal.variance = FALSE, method = c("hungarian","bruteforce"))
{
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x,unit,w=w)
	m <- pre$m; n <- pre$n; p <- pre$p
	R <- pre$R # Cholesky decomposition of w or null
	x <- pre$x; dim(x) <- c(p,m,n)
	if (length(y) == 1) { 
		y <- x[,,y] 
	} else {
		stopifnot(NROW(y) == p && NCOL(y) == m)	
	}
	rm(pre)
	syscall <- sys.call()
	
	## Trivial case p == 1 
	if (p == 1) {
		dim(x) <- c(m,n)
		sigma <- apply(x,2,order)
		x <- apply(x,2,sort)
		mu <- rowMeans(x)
		V <- rowMeans(x^2) - mu^2
		cost <- n * sum(V)
		if (equal.variance) V <- rep(mean(V),m)
		return(list(sigma=sigma, cost=cost, mu=mu, V=V))
	}	
	
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
	xymin <- min(x,y)
	if (xymin < 0) {
		x <- x - xymin
		y <- y - xymin }
		
	## Assignment method
	method <- match.arg(method)
	if (method == "bruteforce") {
		perms <- permn(k)
		perms <- lapply(1:NROW(perms), 
			function(i) perms[i,]) }
		
	## Match feature vectors 
	sigma <- matrix(,m,n)
	for (i in 1:n)
		sigma[,i] <- if (method == "bruteforce") {
			brute(x[,,i],y,perms)
	} else { solve_LSAP(crossprod(y,x[,,i]), 
		maximum = TRUE) }
			
	## Calculate cost 
	cost <- 0
	for (i in 1:n)
		cost <- cost + sum((x[,sigma[,i],i] - y)^2)	
	
	## Sample means and covariances of matched vectors
	mu <- matrix(,p,m)
	V <- array(,c(p,p,m))
	dim(x) <- c(p,m*n)
	for (l in 1:m) {
		idx <- seq.int(0,by=m,len=n) + sigma[l,]
		mu[,l] <- rowMeans(x[,idx,drop=F])
		V[,,l] <- tcrossprod(x[,idx,drop=F])/n - 
			tcrossprod(mu[,l])	
	}	
	if (xymin < 0) mu <- mu + xymin
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
 	
 	## Sums of squares for unmatched data
 	xbarl <- matrix(,p,m)
 	for (l in 1:m) 
 		xbarl[,l] <- rowMeans(x[,seq.int(l,by=m,len=n),drop=FALSE])
 	xbar <- rowMeans(xbarl)
 	ssb <- n * sum((xbarl-xbar)^2)
 	sst <- sum((x-xbar)^2)

	out <- list(sigma=sigma, cost=cost, mu=mu, V=V, 
		ss.between.unmatched=ssb, ss.within.unmatched=sst-ssb,
		call=syscall)
	class(out) <- "matchFeats"
	return(out)	
}
