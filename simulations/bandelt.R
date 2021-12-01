##########################################################
# ADAPTATION OF LAGRANGIAN RELAXATION METHOD OF BANDELT 
# ET AL (2004) TO OBTAIN LOWER BOUND IN MULTIDIMENSIONAL 
# ASSIGNMENT PROBLEM WITH DECOMPOSABLE COSTS (MDADC) 
##########################################################


# References: 
# Bandelt,  Maas, and Spieksma (2004) "Local search heuristics for 
# multi-index assignment problems with decomposable costs"
# Degras (2021) "Scalable Feature Matching Across Large Data Collections"





bandelt.lb <- function(x, unit, hub = 1L, rho = 0.1, maxit = 1000L, 
	verbose = FALSE)
{
	n <- length(unique(unit))
	m <- nrow(x) / n
	p <- ncol(x)
	x <- t(x)
	dim(x) <- c(p,m,n)
	npairs <- choose(n,2)	
	
	z <- array(dim=c(m,m,npairs))
	
	## Calculate distances
	d <- array(dim=c(m,m,npairs))
	nrmx2 <- colSums(x^2)
	pairs <- combn(n,2)
	for (k in 1:npairs) 
		d[,,k] <- matrix(nrmx2[,pairs[1,k]],m,m) + 
			matrix(nrmx2[,pairs[2,k]],m,m, byrow=TRUE) + 
			- 2 * crossprod(x[,,pairs[1,k]],x[,,pairs[2,k]])	
	
	## Logical: does given pair (i,j) contain reference dimension?
	h <- hub
	has.h <- (pairs[1,] == h | pairs[2,] == h)
	pairs.h <- which(has.h)
	pairs.noh <- (1:npairs)[-pairs.h]

	## Index of {h,r} in list of pairs
	hr.idx <- c(head(pairs.h,h-1),NA,tail(pairs.h,n-h)) 
	
	## Lagrange multipliers
	lambda <- array(0, dim=c(m,m,m,choose(n,2)))
	
	## Indices of pairs containing dimension r but not h
	## (r=1,...,n, r!=h) with other index smaller than r 
	## (left) or larger than r (right)
	locs <- vector("list",n)
	for (r in (1:n)[-h]) {
		locs[[r]]$left <- setdiff(which(pairs[2,] == r), pairs.h)
		locs[[r]]$right <- setdiff(which(pairs[1,] == r), pairs.h) 
	}
	
	## Lower bound 
	lb <- numeric(maxit)

	## Loop
	for (it in 1:maxit)
	{
		## Update z for pairs that do not contain the hub  
		z[,,pairs.noh] <- 
			(d[,, pairs.noh] < colSums(lambda[,,, pairs.noh]))	

		## Update z for (n-1) pairs that contain the hub
		for (k in pairs.h) {
			# Is h the smaller index in the pair?
			h1 <- pairs[1,k] == h
			# Other member of the pair 
			r <- if (h1) pairs[2,k] else pairs[1,k]
			# Contribution of cost
			# Make sure that dimension 1 = h, dimension 2 = r 
			A1 <- if (h1) d[,,k] else t(d[,,k])
			# Contribution of Lagrange multipliers associated with
			# triplet (h,r,s) such that s < r 
			# lambda(u,v,w,(r,s)): first permute dimensions 2 & 3
			# so that dimension 2 corresponds to r, then sum over 
			# dimensions 3 & 4
			idx <- locs[[r]]$left
			A2 <- rowSums(aperm(lambda[,,,idx,drop=FALSE],
				c(1,3,2,4)), dims=2L)
			idx <- locs[[r]]$right
			A3 <- rowSums(lambda[,,,idx,drop=FALSE], dims=2L)
			A <- A1 + A2 + A3	
			sigma <- clue::solve_LSAP(A)
			ztmp <- matrix(0,m,m)	
			if (h1) { ztmp[cbind(1:m,sigma)] <- 1	
				} else { ztmp[cbind(sigma,1:m)] <- 1 }
			z[,,k] <- ztmp
		}

		## Gradient mu
		# mu[u,v,w,rs] = z[u,v,hr] + z[u,w,hs] - z[v,w,rs] - 1
		# zexp <- array(z,c(dim(z),m)) # (u,v,rs,C)
		# This array is constant along dimension C
		mu <- array(0, dim=c(m,m,m,npairs))
		for (k in pairs.noh) {
			r <- pairs[1,k]
			s <- pairs[2,k]
			mu[,,,k] <- if (h < r && h < s) {
				rep(z[,,hr.idx[r]],m) + rep(z[,,hr.idx[s]],m) - 
				rep(z[,,k], each=m) - 1
			} else if (h > r && h < s) {
				rep(t(z[,,hr.idx[r]]),m) + rep(z[,,hr.idx[s]],m) - 
				rep(z[,,k], each=m) - 1
			} else if (h < r && h > s) {
				rep(z[,,hr.idx[r]],m) + rep(t(z[,,hr.idx[s]]),m) - 
				rep(z[,,k], each=m) - 1
			} else {
				rep(t(z[,,hr.idx[r]]),m) + rep(t(z[,,hr.idx[s]]),m) - 
				rep(z[,,k], each=m) - 1	
			}
		}
				
		## Update lambda
		lambda <- pmax(lambda + rho * mu, 0)
		
		## Calculate associated lower bound
		idx <- which(z == 1)
		lb[it] <- (sum(z[idx] * d[idx]) + sum(lambda * mu)) / n / (n-1)
		if (verbose && it %% 10 == 0) 
			cat("Iteration",it,"Lower bound",lb[it],"\n")	
	}
	return(max(lb))
	
}


