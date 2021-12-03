

#################################
# Regularize covariance matrices
#################################


regularize.cov <- function(V,cond)
{
	p <- nrow(V)
	eigV <- 	eigen(V,TRUE)
	lambda <- eigV$values		
	a <- lambda[1] / cond
	test <- (lambda < a)
	if (any(test)) {
		lambda[test] <- a
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
	p <- nrow(x)
	logphi <- array(dim=c(m,m,n))
	for (l in 1:m) {
		R <- tryCatch(chol(V[,,l]), error=function(e) NULL)
		if (!is.null(R)) {
			# Case: non-degenerate covariance
			e <- backsolve(R, x-mu[,l], transpose=TRUE)
			log.det <- 2 * sum(log(diag(R)))
		} else {
			# Case: degenerate covariance
			eigV <- eigen(V[,,l], TRUE)
			lambda <- eigV$values
			nlam <- sum(lambda > 0)
			if (nlam == 0) { 
				# Case: zero matrix
				test <- colSums(x == mu[,l])
				test[test == p] <- Inf
				test[test < p] <- -Inf
				logphi[,l,] <- test
				next
			} else {
				lambda <- lambda[1:nlam]
				P <- t(eigV$vectors[,1:nlam]) / sqrt(lambda)
				e <- P %*% (x-mu[,l]) 
				log.det <- sum(log(lambda)) 
			}
		}
		logphi[,l,] <- -0.5 * (colSums(e^2) + log.det + p * log(2*pi))
	}
	if (log) return(logphi) else return(exp(logphi))

}



####################################################################


#################################
# Iteratively normalize row sums 
# and column sums of a matrix 
#################################


scale.rc <- function(x, niter = NULL)
{
	if (length(x) == 1) return(1)
	m <- nrow(x)
	n <- ncol(x)
	if (is.null(niter)) niter <- max(m,n)
	for (i in 1:niter) {
		x <- x / rowSums(x)
		x <- x / rep(colSums(x), each = m) 
	}
	return(x)
}



scale.rc2 <- function(x, niter = NULL)
{
	if (length(x) == 1) return(list(x=1, c=1))
	m <- nrow(x)
	n <- ncol(x)
	if (is.null(niter)) niter <- max(m,n)
	trace.rsum <- matrix(,m,niter)
	trace.csum <- matrix(,n,niter)
	for (i in 1:niter) {
		trace.rsum[,i] <- rowSums(x)
		x <- x / trace.rsum[,i]
		trace.csum[,i] <- colSums(x)
		x <- x / rep(trace.csum[,i], each = m) 
	}
	return(list(x=x, c=c(trace.rsum,trace.csum)))
}




######################################
# Function to calculate the permanent 
# of a matrix using Ryser formula
######################################


ryser.perm <- function(A)
{
	if (length(A) == 1) return(as.numeric(A))
	n <- ncol(A)
	hn <- floor(n/2)
	even <- (n == 2 * hn)
	rsum <- rowSums(A)
	csum <- colSums(A)
	if(any(rsum == 0 | csum == 0)) return(0)
	per <- numeric(n)
	per[1] <- prod(rsum)
	f <- function(idx) rowSums(A[,idx])
	kmax <- if(even) (hn-1) else hn
	if (n > 2) {
		for (k in 1:kmax) {
			Ak <- apply(combn(n,n-k), 2, f)
			Ak2 <- rsum - Ak
			Ak <- apply(Ak,2,prod)
			per[k+1] <- (-1)^k * sum(Ak)
			Ak2 <- apply(Ak2,2,prod)
			per[n-k+1] <- (-1)^(n-k) * sum(Ak2)
		}
	}
	if (even) {
		g <- function(idx) prod(rowSums(A[,idx]))
		Ak <- apply(combn(n,hn),2,g)
		per[hn+1] <- (-1)^hn * sum(Ak)	
	}
	
	per <- sum(per)	
	if (per >= 0) return(per) else return(NA)

}


all.perm <- function(A, eps = 1e-10) 
{
	m <- NCOL(A)
	if (m == 1) return(A)
	if (all(A == 0)) return(matrix(1/m,m,m))
	nzr <- (rowSums(A) > 0)
	nzc <- (colSums(A) > 0)
	A[nzr,nzc] <- scale.rc(A[nzr,nzc])
	nz <- (A > eps)
	# rmax <- apply(A,1,max)
	# cmax <- apply(A,2,max)
	# nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
	if (all(rowSums(nz) == 1 & colSums(nz) == 1))
		{ A[nz] <- 1; A[!nz] <- 0; return(A) }
	out <- matrix(0,m,m)
	for (k in 1:m) {
		for (l in 1:m) {
			if (A[k,l] > 0) 
				out[k,l] <- A[k,l] * ryser.perm(A[-k,-l])
		}
	}
	out[is.na(out)] <- 0
	return(out)
}


####################################################################




#############################################
# Orthogonal projection of a square matrix 
# onto the set of doubly stochastic matrices
#############################################



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

class.probs <- function(phi, method = c("exact","approx"), 
	maxit, eps, beta, parallel)
{

	m <- dim(phi)[1]
	n <- dim(phi)[3]
	eps0 <- 0 # 1e-10
	method <- match.arg(method)
	out <- array(0,dim=c(m,m,n))
	
	if (method == "exact") {
		if (parallel) {
			out <- foreach(i = 1:n, 
				.packages = "matchFeats") %dopar% all.perm(phi[,,i], eps0)
			out <- unlist(out)
		} else {
			out <- apply(phi, 3, all.perm, eps=eps0)
		}
		dim(out) <- c(m,m,n)

	} else if (method == "approx") {
		for (i in 1:n) {
			A <- phi[,,i]
			rmax <- apply(A,1,max)
			cmax <- apply(A,2,max)
			nz <- (A >= (eps0*rmax)) | (A >= matrix(eps0*cmax,m,m,TRUE))
			if (all(rowSums(nz) == 1 & colSums(nz) == 1))
				{ A[nz] <- 1; A[!nz] <- 0; out[,,i] <- A; next }
			nzr <- (rowSums(A) > 0)
			nzc <- (colSums(A) > 0)
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
					sigma <- solve_LSAP(logP,TRUE)			
					out[k,l,i] <- A[k,l] * prod(P[cbind(1:(m-1),sigma)])
				}
			}
		}
	} 

	## Project on set of doubly stochastic matrices
	if (beta < 1) 
		out <- out^beta
	for (i in 1:n) {
		sumout <- sum(out[,,i])/m
		out[,,i] <- if (sumout > 0) 
			{ proj.ds(out[,,i]/sumout,maxit,eps) } else { 1/m }
	}
	# out[out < .Machine$double.eps] <- 0
	return(out)	
}




####################################################################




############################
# Calculate log-likelihood 
############################


logLik <- function(phi, parallel)
{
	m <- dim(phi)[1]
	n <- dim(phi)[3]
	eps <- 0 # 1e-10
	
	## Mixture density values and log-likelihood 
	if (parallel) {
		logL <- 	foreach(i = 1:n, .packages="matchFeats", 
			.combine=c) %dopar% {
		A <- phi[,,i]
		rmax <- apply(A,1,max)
		cmax <- apply(A,2,max)
		if (any(rmax == 0 | cmax == 0))
			return(-Inf)
		tmp <- scale.rc2(A)
		A <- tmp$x
		c <- tmp$c
		nz <- (A >= eps)
		# nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
		if (all(rowSums(nz) == 1 & colSums(nz) == 1))
			return(sum(log(A[nz]),log(c)))
		return(log(ryser.perm(A)) + sum(log(c))) }
		# rsum <- rowSums(A)
			# return(log(ryser.perm(A/rsum)) + sum(log(rsum))) }
	} else {
		logL <- numeric(n)
		for (i in 1:n) {
			A <- phi[,,i]
			rmax <- apply(A,1,max)
			cmax <- apply(A,2,max)
			if (any(rmax == 0 | cmax == 0))
				{logL[i] <- -Inf; next}
			tmp <- scale.rc2(A)
			A <- tmp$x
			c <- tmp$c
			nz <- (A >= eps)
			# nz <- (A >= (eps*rmax)) | (A >= matrix(eps*cmax,m,m,TRUE))
			logL[i] <- if (all(rowSums(nz) == 1 & colSums(nz) == 1))
				{ sum(log(A[nz]),log(c)) 
				} else { sum(log(ryser.perm(A)),log(c)) }
			# rsum <- rowSums(A)
			# logL[i] <- log(ryser.perm(A/rsum)) + sum(log(rsum))
		}
	}
	
	test <- is.infinite(logL) | is.na(logL)
	if (any(test)) logL[test] <- min(logL[!test]) 
	logL <- sum(logL) - n*sum(log(1:m))
	return(logL)
}




#####################################################################




###################################
# Check for trivial solutions to 
# the one-to-one matching problem
###################################



# This function checks whether all sets of vectors are identical
# The reasons for this verification are: 1) the verification 
# is very quick (negligible time in comparison to running the 
# Gaussian mixture EM) and 2) it may avoid issues of degeneracy
# (covariance) and/or lack of identifiability in the EM



check.set.equality <- function(x,n,cond)
{	
	p <- nrow(x)
	mn <- ncol(x)
	m <- mn/n
	out <- list(sigma=NULL, P=matrix(1,m,n), mu=NULL, 
		V=array(diag(1/cond,p,p),c(p,p,m)), loglik=Inf)
	
	## Test if all vectors are identical
	if (all(x[,1:(mn-1)] == x[,2:mn])) {	
		out$mu <- x[,1:m]
		out$sigma <- matrix(1:m,m,n)
		return(out)		
	}
	
	## Test if all vectors are setwise identical (set = unit) 
	test <- matrix(,p,n)
	idx <- matrix(1:mn,m,n)
	flag <- FALSE
	for (r in 1:p) {
		for (i in 1:n) {
			test[r,i] <- length(unique(x[r,idx[,i]])) == m
		}
		ntrue <- sum(test[r,])
		if (ntrue > 0 && ntrue < n) return(NULL)
		if (ntrue == n) {flag <- TRUE; break}
	}	
	## Case: at least one row has only unique values
	if (flag) {
		sigma <- matrix(,m,n)
		o1 <- order(x[r,1:m])
		sigma[,1] <- o1
		for (i in 2:n) {
			o2 <- order(x[r,idx[,i]])
			sigma[,i] <- o2 
			if (any(x[,o1] != x[,idx[o2,i]])) 
				return(NULL)
		}
		out$mu <- x[,o1]
		out$sigma <- sigma
		return(out)		
	}
	## Other cases
	sigma <- matrix(,m,n)
	sigma[,1] <- 1:m
	for (i in 2:n) {
		active <- idx[,i]
		for (k in 1:m) {
			matched <- FALSE
			for (l in active) {
				if (all(x[,k] == x[,l])) {
					active <- setdiff(active,l)
					sigma[l,i] <- k
					matched <- TRUE 
					break }
			if (!matched) return(NULL)
			}			
		}
	}
	out$mu <- x[,1:m]
	out$sigma <- sigma
	return(out)		
	
}



#####################################################################





##########################################
# EM algorithm for Gaussian mixture model
##########################################



match.gaussmix <- function(x, unit = NULL, mu = NULL, V = NULL, 
	equal.variance = FALSE, method = c("exact","approx"), 
	fixed = FALSE, control = list())
{
	
	## Preprocess input arguments: check dimensions, 
	## reshape and recycle as needed	
	if (fixed && (is.null(mu) || is.null(V))) 
		stop("If 'fixed' is TRUE, please specify 'mu' and 'V'")
	pre <- preprocess(x,unit,mu,V)
	m <- pre$m; n <- pre$n; p <- pre$p
	mu <- pre$mu; V <- pre$V
	if (!is.null(pre$x)) x <- pre$x

	rm(pre)
	syscall <- sys.call()
				
	## Tuning parameters
	con <- list(maxit = 1e4, eps = 1e-8, cond = 1e8, beta = 1,
		betarate = 1, parallel = FALSE, verbose = FALSE)
	if (length(control) > 0) {
		name <- intersect(names(control), names(con))
		con[name] <- control[name]	
	}
	maxit <- con$maxit
	eps <- con$eps
	cond <- con$cond
	beta <- con$beta
	betarate <- con$betarate
	parallel <- con$parallel
	verbose <- con$verbose
	method <- match.arg(method)
	maxit.proj <- 1000 
	eps.proj <- 1e-6
		
 	## Sums of squares for unmatched data
 	dim(x) <- c(p,m,n)
 	mu <- rowMeans(x,dims=2)
 	ssw <- sum((x-as.vector(mu))^2)
 	ssb <- n * sum((mu-rowMeans(mu))^2)
	dim(x) <- c(p,m*n)
 	
 	## Check set-wise equality of feature vectors between units  	
	test <- check.set.equality(x,n,cond)
	if (!is.null(test)) {
		test$call <- syscall
		class(test) <- "matchFeat"
		return(test)	
	}
	
	## Initialize parameter estimates
	if (is.null(mu) || is.null(V)) {
		dim(x) <- c(p,m,n)
		init <- match.bca(x, 
			control=list(equal.variance=equal.variance))
		mu <- init$mu
		V <- init$V
		dim(x) <- c(p,m*n)
	}
	if (equal.variance) 
		V <- array(V[1:(p^2)],c(p,p,m))
		
	## Regularize covariance matrices if needed
	for (l in 1:m) 
		V[,,l] <- regularize.cov(V[,,l], cond)
		
	logL.best <- -Inf
	counter <- 0
	max.counter <- 2 # for simulations, set. back to 3 later
	
	for (it in 1:maxit) {
				
		## Log-likelihood	
		logphi <- norm.dens(x, mu, V, log=TRUE)
		phi <- exp(logphi)
		logL <- logLik(phi, parallel)
		if (verbose) 
			cat("Iteration:",it,"Log-likelihood:",logL,"\n")
			
		## Monitor progress
		progress <- (it == 1 || (logL-logL.best) > eps * abs(logL.best))
		counter <- if (progress) 0 else counter + 1 

		## Best estimate
		if (logL > logL.best) {
			logL.best <- logL; mu.best <- mu; V.best <- V; 
			P.best <- if (it == 1) NULL else P }
		if (counter >= max.counter) break

		## E step
		P <- class.probs(phi,method,maxit.proj,eps.proj,beta,parallel)
		beta <- min(beta * betarate, 1)

		## M step
		if (!fixed) {
			Q.old <- apply(logphi*P,2,sum)
			V.tmp <- array(,c(p,p,m))
			for (l in 1:m) {
				w <- P[,l,] / sum(P[,l,])
				dim(w) <- NULL
				mu[,l] <- x %*% w
				V.tmp[,,l] <- tcrossprod(x * matrix(sqrt(w),p,m*n,byrow=T)) - 
					tcrossprod(mu[,l])
				V.tmp[,,l] <- regularize.cov(V.tmp[,,l], cond)
			}
			if (equal.variance) {
				V.tmp <- apply(V.tmp,1:2,mean)
				V.tmp <- array(V.tmp,c(p,p,m))
			}			
			logphi <- norm.dens(x,mu,V.tmp,log=TRUE)
			Q <- apply(logphi*P,2,sum)
			idx <- (Q > Q.old)
			V[,,idx] <- V.tmp[,,idx]
		}		
	}
	
	if (is.null(P.best)) {
		phi <- norm.dens(x,mu.best,V.best)
		P.best <- class.probs(phi,method,maxit.proj,eps.proj,beta,parallel)
	}
	
	## Determine the most likely class allocation for each unit
	sigma <- matrix(,m,n)
	for (i in 1:n) {
		nz <- (P.best[,,i] > 0)
		if (all(rowSums(nz) == 1) && all(colSums(nz) == 1)) { 
			P.best[,1,i] <- 1
			sigma[,i] <- which(nz, arr.ind=TRUE)[,1]
		} else {
			logP <- log(P.best[,,i])
			logP[is.infinite(logP)] <- NA
			logP <- logP + m * max(logP,na.rm=TRUE) - 
				(m+1) * min(logP,na.rm=TRUE)
			logP[is.na(logP)] <- 0
			sigma[,i] <- solve_LSAP(logP,maximum=TRUE)  
			P.best[,1,i] <- P.best[cbind(1:m,sigma[,i],rep(i,m))]
		}
	}
	if (equal.variance) V.best <- V.best[,,1]
	
 	## Cluster assignment
 	cluster <- matrix(,m,n)
 	for (i in 1:n)
 		cluster[sigma[,i],i] <- 1:m
 		
	out <- list(sigma = sigma, cluster = cluster, P = P.best[,1,], 
		mu = mu.best, V = V.best, loglik = logL.best, call = syscall) 
	class(out) <- "matchFeat"
	return(out)
	
}






