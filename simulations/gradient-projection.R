########################################################
# CONVEX RELAXATION OF THE FEATURE MATCHING PROBLEM (7) 
# OF DEGRAS (2021) "SCALABLE FEATURE MATCHING ACROSS 
# LARGE DATA COLLECTIONS"
########################################################


# (7) min F(P1,...,Pn) = sum(i,j) || Xi Pi - Xj Pj ||_F^2 / n / (n-1)
# where i,j vary in 1,...,n, X1,...,Xn are p x m matrices, 
# and P1,...,Pn are m x m permutation matrices
# Convex relaxation: P1,...,Pn as doubly stochastic matrices
# The function grad.proj implements an accelerated gradient
# projection method for the convex relaxation of (7)


grad.proj <- function(x, D = NULL, maxit = 1000, tol = 1e-8)
	# stepsize = c("fixed","backtrack"), # eta = 1.25, 
{
	p <- dim(x)[1]
	m <- dim(x)[2]
	n <- dim(x)[3]
	
	## Initialize stepsize
	# stepsize <- match.arg(stepsize)
	nrmx <- apply(x,3,norm,type="2")
	imax <- which.max(nrmx)
	vmax <- svd(x[,,imax],nu=0,nv=1)$v
	xv <- matrix(0,p,n)
	for (i in 1:n)
		xv[,i] <- x[,,i] %*% vmax
	L <- nrmx[imax]^2 - ifelse(n>=100,0.9,0.5)/n * sum(xv[,imax]*xv)
	
	## Initialize solution D if needed 
	if (is.null(D))
		{ D <- rep(diag(m),n); dim(D) <- c(m,m,n) } 
		
	## Initial objective
	xD <- array(dim=c(p,m,n))
	for (i in 1:n)
		xD[,,i] <- x[,,i] %*% D[,,i]
	xDbar <- as.vector(rowMeans(xD, dims=2L))
	objective <- (sum(xD^2) - n*sum(xDbar^2))/(n-1)
	objective.best <- objective
	D.best <- D 
	
	## Misc
	grad <- array(0,dim=c(m,m,n))
	D.adj <- D
	tt <- 1
	counter <- 0

	for (count in 1:maxit)
	{
		objective.old <- objective
		D.old <- D
		tt.old <- tt
		
		## Gradient step
		for (i in 1:n)
			xD[,,i] <- x[,,i] %*% D.adj[,,i]
		xDbar <- as.vector(rowMeans(xD, dims=2L))
		xD <- xD - xDbar
		for (i in 2:n)
			grad[,,i] <- crossprod(x[,,i], xD[,,i])
		D <- D.adj - grad / L
		
		## Projection step
		D <- apply(D, 3, matchFeat:::proj.ds)
		dim(D) <- c(m,m,n)
		
		## Acceleration step
		tt <- (1+sqrt(4*tt.old^2))/2
		D.adj <- D + ((tt.old-1)/tt) * (D - D.old)
		
		## Objective 
		for (i in 1:n)
			xD[,,i] <- x[,,i] %*% D[,,i]
		xDbar <- as.vector(rowMeans(xD, dims=2L))
		objective <- (sum(xD^2) - n*sum(xDbar^2))/(n-1)
		if (objective < objective.best) {
			objective.best <- objective
			D.best <- D
		}
		
		## Monitor progress
		progress <- (objective < (1-tol) * objective.best)
		counter <- ifelse(progress, 0, counter + 1)	
		if (counter == 5) break
	}
	
	return(list(D=D.best, objective=objective.best))
	
}



############################################
# RELAXATION METHOD OF BANDELT ET AL (2004)
# TO OBTAIN LOWER BOUND ON MDADC PROBLEM
############################################

match.lb <- function(x, unit, hub = 1L, rho = 0.1, maxit = 1000L, 
	verbose = FALSE)
{
	n <- length(unique(unit))
	m <- nrow(x) / n
	p <- ncol(x)
	x <- t(x)
	dim(x) <- c(p,m,n)
	npairs <- choose(n,2)	
	
	# if (is.null(z)) {
		# z <- rep(diag(m),choose(n,2))
		# dim(z) <- c(m,m,choose(n,2))
	# }
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
			sigma <- solve_LSAP(A)
			ztmp <- matrix(0,m,m)	
			if (h1) { ztmp[cbind(1:m,sigma)] <- 1	
				} else { ztmp[cbind(sigma,1:m)] <- 1 }
			z[,,k] <- ztmp
		}

		## Update lambda

		# Gradient mu
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
				
		lambda <- pmax(lambda + rho * mu, 0)
		
		## Calculate associated lower bound
		idx <- which(z == 1)
		lb[it] <- sum(z[idx] * d[idx]) + sum(lambda * mu)
		if (verbose && it %% 10 == 0) 
			cat("Iteration",it,"Lower bound",lb[it],"\n")	
	}
	return(lb)
	
}


