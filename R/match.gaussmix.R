

#################################
# Regularize covariance matrices
#################################


regularize.cov <- function(V, cond = 1e8)
{
	p <- nrow(V)
	eigV <- 	eigen(V,TRUE)
	lambda <- eigV$values		
	test <- (lambda < lambda[1] / cond)
	if (any(test)) {
		lambda[test] <- lambda[1] / cond
		P <- eigV$vectors
		V <- tcrossprod(P * rep(sqrt(lambda), each=p))
	}	
	return(V)
}



####################################
# Calculate multivariate normal 
# density values phi(x_ik,mu_l,V_l) 
####################################


norm.dens <- function(x,mu,V,log=FALSE)
{
	m <- ncol(mu)
	n <- ncol(x)/m
	phi <- array(dim=c(m,m,n))
	for (l in 1:m) {
		R <- tryCatch(chol(V[,,l]), error=function(e) NULL)
		if (!is.null(R)) {
			# Case: non-degenerate covariance
			e <- backsolve(R, x-mu[,l], transpose=TRUE)
			log.det <- 2 * sum(log(diag(R))) + p * log(2*pi)
		} else {
			# Case: degenerate covariance
			eigV <- eigen(V[,,l], TRUE)
			lambda <- eigV$values
			nlam <- sum(lambda > 1e-8 * lambda[1])
			if (nlam == 0) { 
				# Case: zero matrix
				e <- matrix(Inf,1,1)
				log.det <- Inf
			} else {
				lambda <- lambda[1:nlam]
				P <- t(eigV$vectors[,1:nlam]) / sqrt(lambda)
				e <- P %*% (x-mu[,l]) 
				log.det <- sum(log(lambda)) + p * log(2*pi)	
			}
		}
		phi[,l,] <- -0.5 * (colSums(e^2) + log.det)
	}
	if (!log) phi <- exp(phi)

	return(phi)
}



####################################################################


#################################
# Iteratively normalize row sums 
# and column sums of a matrix 
#################################


scale.rc <- function(x, niter = NULL)
{
	m <- nrow(x)
	n <- ncol(x)
	if (is.null(niter)) niter <- max(m,n)
	for (i in 1:niter) {
		x <- x / rowSums(x)
		x <- x / rep(colSums(x), each = m) 
	}
	return(x)
}




######################################
# Function to calculate the permanent 
# of a matrix using Ryser formula
######################################


ryser.perm <- function(A)
{
	# if (length(A) == 1) return(A)
	n <- ncol(A)
	rsum <- rowSums(A)
	csum <- colSums(A)
	if(any(rsum == 0 | csum == 0)) return(0)
	per <- numeric(n)
	per[1] <- prod(rsum)
	for (k in 1:(n-1)) {
		Ak <- array(A[,combn(n,n-k)], c(n,n-k,choose(n,k)))
		Ak <- apply(Ak,c(1,3),sum)
		Ak <- apply(Ak,2,prod)
		per[k+1] <- (-1)^k * sum(Ak)
	}
	# Prevent negative values caused by underflow
	ub1 <- per[1]
	per <- sum(per)
	if (per > 0) return(per) else return(min(ub1,prod(csum)))
}



####################################################################




###################################################
# Function to project a square matrix orthogonally 
# on the set of doubly stochastic matrices
###################################################



proj.ds <- function(x, maxit=1000, eps=1e-8)
{
	n <- ncol(x)
	for (i in 1:maxit)
	{
		rsum <- rowSums(x)
		if (all(abs(rsum-1) <= eps) && 
			all(abs(colSums(x)-1) <= eps) && 
			all(x >= 0 & x <= 1)) break		
		x <- x - (rsum - 1)/n
		x <- x - matrix((colSums(x)-1)/n,n,n,byrow=TRUE)	
		x[x > 1] <- 1
		x[x < 0] <- 0
	}
	return(x)
}




####################################################################



###################################
# Calculate the class 
# probabilities P(I_ikl = 1 | X_i)
###################################

class.probs <- function(phi, method = c("exact","approx"), maxit, eps)
{

	m <- dim(phi)[1]
	n <- dim(phi)[3]
	eps <- 1e-10
	method <- match.arg(method)
	out <- array(0,dim=c(m,m,n))
	
	if (method == "exact") {
		for (i in 1:n) {
			A <- phi[,,i]
			rmax <- apply(A,1,max)
			cmax <- apply(A,2,max)
			nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
			if (all(rowSums(nz) == 1 & colSums(nz) == 1))
				{ A[nz] <- 1; A[!nz] <- 0; out[,,i] <- A; next }
			nzr <- (rmax > 0)
			nzc <- (cmax > 0)
			if (all(!nzr)) { out[,,i] <- 1/m; next }
			A[nzr,nzc] <- scale.rc(A[nzr,nzc,drop=F])
			for (k in 1:m) {
				for (l in 1:m) {
					if (A[k,l] > 0) 
						out[k,l,i] <- A[k,l] * ryser.perm(A[-k,-l,drop=F])
				}
			}
		}
	} else if (method == "approx") {
		for (i in 1:n) {
			A <- phi[,,i]
			rmax <- apply(A,1,max)
			cmax <- apply(A,2,max)
			nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
			if (all(rowSums(nz) == 1 & colSums(nz) == 1))
				{ A[nz] <- 1; A[!nz] <- 0; out[,,i] <- A; next }
			if (all(!nzr)) { out[,,i] <- 1/m; next }
			A[nzr,nzc] <- scale.rc(A[nzr,nzc,drop=F])
			for (k in 1:m) {
				for (l in 1:m) {
					if (A[k,l] == 0) next
					P <- A[-k,-l,drop=F]
					logP <- log(P)
					logP[is.infinite(logP)] <- NA
					logP <- logP + m * max(logP,na.rm=TRUE) - 
						(m+1) * min(logP,na.rm=TRUE)
					logP[is.na(logP)] <- 0					
					sigma <- clue::solve_LSAP(logP,TRUE)			
					out[k,l,i] <- A[k,l] * prod(P[cbind(1:(m-1),sigma)])
				}
			}
		}
	} 

	## Project on set of doubly stochastic matrices
	for (i in 1:n) {
		sumout <- sum(out[,,i])/m
		out[,,i] <- if (sumout > 0) 
			{ proj.ds(out[,,i]/sumout,maxit,eps) } else { 1/m }
	}
	out[out < .Machine$double.eps] <- 0
	return(out)	
}




####################################################################




############################
# Calculate log-likelihood 
############################


logLik <- function(phi)
{
	m <- dim(phi)[1]
	n <- dim(phi)[3]
	eps <- 1e-10
	
	## Mixture density values and log-likelihood 
	logL <- numeric(n)
	for (i in 1:n) {
		# sigma <- clue::solve_LSAP(phi[,,i]-min(phi[,,i]),TRUE)
		# maxval <- sum(phi[cbind(1:m,sigma,rep(i,m))]) 
		# z <- exp(phi[,,i] - maxval/m)
		# logL[i] <- log(ryser.perm(z)) + maxval  
		A <- phi[,,i]
		rmax <- apply(A,1,max)
		cmax <- apply(A,2,max)
		if (any(rmax == 0 | cmax == 0))
			{logL[i] <- -Inf; next}
		nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
		if (all(rowSums(nz) == 1 & colSums(nz) == 1))
			{logL[i] <- sum(log(A[nz])); next}
		rsum <- rowSums(A)
		logL[i] <- log(ryser.perm(A/rsum)) + sum(log(rsum))
	}

	test <- is.infinite(logL)
	if (any(test)) logL[test] <- min(logL[!test]) 
	logL <- sum(logL) - n*sum(log(1:m))
	return(logL)
}




#####################################################################





##########################################
# EM algorithm for Gaussian mixture model
##########################################


match.gaussmix <- function(x, unit=NULL, mu=NULL, V=NULL, 
	equal.variance=FALSE, method=c("exact","approx"), 
	control=list())
{
	
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	pre <- preprocess(x,unit,mu,V)
	m <- pre$m; n <- pre$n; p <- pre$p
	x <- pre$x; mu <- pre$mu; V <- pre$V
	rm(pre)
			
	## Tuning parameters
	con <- list(maxit=1e4, eps=1e-8, verbose=FALSE)
	if (length(control) > 0) {
		name <- intersect(names(control), names(con))
		con[name] <- control[name]	
	}
	for (name in names(con))
		assign(name, con[[name]])
	method <- match.arg(method)
	maxit.proj <- 1000 
	eps.proj <- 1e-6
		
	## Initialize parameter estimates
	if (is.null(mu) || is.null(V)) {
		dim(x) <- c(p,m,n)
		init <- match.bca(x = x, unit = NULL, method = "c", 
			w = NULL, equal.variance = equal.variance)
		mu <- init$mu
		V <- init$V
		dim(x) <- c(p,m*n)
	} 
	if (equal.variance) 
		V <- array(V[1:(p^2)],c(p,p,m))
		
	## Regularize covariance matrices if needed
	for (l in 1:m) 
		V[,,l] <- regularize.cov(V[,,l])
		
	logL <- NA
	
	for (count in 1:maxit) {
				
		## Log-likelihood	
		logL.old <- logL
		logphi <- norm.dens(x, mu, V, log=TRUE)
		phi <- exp(logphi)
		logL <- logLik(phi)
		if (verbose) 
			cat("Iteration:",count,"Log-likelihood:",logL,"\n")
			
		## Check progress
		if (	(count > 1 && (logL-logL.old) <= eps * abs(logL.old)) || 
			count == maxit) break		

		## E step
		P <- class.probs(phi,method,maxit.proj,eps.proj)

		## M step
		V.tmp <- array(dim=c(p,p,m))
		for (l in 1:m) {
			Pl <- P[,l,] / sum(P[,l,])
			dim(Pl) <- NULL
			mu[,l] <- x %*% Pl 	
			V.tmp[,,l] <- tcrossprod(x * matrix(sqrt(Pl),p,m*n,byrow=T)) - 
				tcrossprod(mu[,l])
			V.tmp[,,l] <- regularize.cov(V.tmp[,,l])
		}
		if (equal.variance) 
			V.tmp <- array(apply(V.tmp,1:2,mean),dim(V))
			
		# Check that covariance update increases Q function
		logphi.tmp <- norm.dens(x, mu, V.tmp, log=TRUE)
		if (equal.variance) {
			if (sum(logphi * P) < sum(logphi.tmp * P))
				V <- V.tmp	
		} else {
			for (l in 1:m) {
				if (sum(logphi[,,l] * P[,,l]) < 
					sum(logphi.tmp[,,l] * P[,,l]))
					V[,,l] <- V.tmp[,,l]		
			}
		}
		
	}
	
	## Determine the most likely class allocation for each unit
	sigma <- matrix(,m,n)
	for (i in 1:n) {
		nz <- (P[,,i] > 0)
		if (all(rowSums(nz) == 1) && all(colSums(nz) == 1)) { 
			P[,1,i] <- 1
			sigma[,i] <- which(nz, arr.ind=TRUE)[,1]
		} else {
			logP <- log(P[,,i])
			logP[is.infinite(logP)] <- NA
			logP <- logP + m * max(logP,na.rm=TRUE) - 
				(m+1) * min(logP,na.rm=TRUE)
			logP[is.na(logP)] <- 0
			sigma[,i] <- clue::solve_LSAP(logP,maximum=TRUE)  
			P[,1,i] <- P[cbind(1:m,sigma[,i],rep(i,m))]
		}
	}
	if (equal.variance) V <- V[,,1]
	
	return(list(sigma=sigma, P=P[,1,], mu=mu, V=V, loglik=logL))
	
}






