######################################
# METHODS OF KUROKI AND MATSUI (2009)
# FOR THE MULTIDIMENSIONAL ASSIGNMENT 
# PROBLEM WITH DECOMPOSABLE COSTS 
######################################





###################################
# Method 1: Integer Linear Program
###################################



kuroki.matsui.ilp <- function(x, parallel = FALSE)
{
	require(clue)
	p <- dim(x)[1]
	m <- dim(x)[2]
	n <- dim(x)[3]
	
	xmin <- min(x)
	if (xmin < 0) x <- x - xmin
		
	# Solve LAP for each pair (i,j) (1 <= i < j <= n)
	combs <- combn(n,2)
	ncombs <- ncol(combs)
	if (parallel) { 
		nwork <- getDoParWorkers()
		start <- round(seq.int(0,ncombs,len=nwork+1))
		match2 <- foreach(i = 1:nwork, .combine = cbind, 
			.packages = "clue") %dopar% {
			idx <- (start[i]+1):start[i+1]
			out <- matrix(,m,length(idx))
			for (k in 1:length(idx)) {
				i <- combs[1,idx[k]]
				j <- combs[2,idx[k]]
				cp <- crossprod(x[,,i],x[,,j])
				out[,k] <- clue::solve_LSAP(cp, maximum=TRUE)	
			}
			return(out)
		}
	} else {
		match2 <- matrix(,m,ncombs)
		for (k in 1:ncombs) {
			i <- combs[1,k]
			j <- combs[2,k]
			cp <- crossprod(x[,,i],x[,,j])
			match2[,k] <- clue::solve_LSAP(cp, maximum=TRUE)	
		}
	}

	## Recover feasible solutions w.r.t. each vertex subset (1:n)
	dim(x) <- c(p,m*n)
	if (parallel) {
		start <- round(seq.int(0,n,len=nwork+1))
		result <- foreach(i = 1:nwork) %dopar% {
			hash <- (start[i]+1):start[i+1]
			sigma <- matrix(,m,n)
			ncostbest <- -Inf
			shift <- seq(0,by=m,len=n)
			for (i in hash) {
				sigma[,i] <- 1:m
				idx <- which(combs[1,] == i)
				sigma[,combs[2,idx]] <- match2[,idx]
				idx <- which(combs[2,] == i)
				sigma[cbind(as.vector(match2[,idx]),
					rep(combs[1,idx],each=m))] <- 1:m
				## Calculate associated negative cost (for maximization problem)
				ncost <- 0
				for (k in 1:m)
					ncost <- ncost + sum(rowSums(x[,shift+sigma[k,]])^2)
				if (ncost > ncostbest) {
					ibest <- i
					sigmabest <- sigma
					ncostbest <- ncost
				}
			}
			return(list(sigma=sigmabest, ncost=ncostbest, hub=ibest))
		}
		ibest <- which.max(sapply(result,"[[","ncost"))
		sigmabest <- result[[ibest]]$sigma
		ncostbest <- result[[ibest]]$ncost
		ibest <- result[[ibest]]$hub				
	} else {
		sigma <- matrix(0L,m,n)
		ncostbest <- -Inf
		shift <- seq(0,by=m,len=n)
		dim(x) <- c(p,m*n)
		for (i in 1:n) {
			sigma[,i] <- 1:m
			idx <- which(combs[1,] == i)
			sigma[,combs[2,idx]] <- match2[,idx]
			idx <- which(combs[2,] == i)
			sigma[cbind(as.vector(match2[,idx]),
				rep(combs[1,idx],each=m))] <- 1:m
			## Calculate associated negative cost (for maximization problem)
			ncost <- 0
			for (k in 1:m)
				ncost <- ncost + sum(rowSums(x[,shift+sigma[k,]])^2)
			if (ncost > ncostbest) {
				ibest <- i
				sigmabest <- sigma
				ncostbest <- ncost
			}
		}	
	}
		
	## Cost of best solution (for minimization problem)
	objective <- (sum(x^2) - ncostbest / n) / (n-1)
	
	return(list(sigma = sigmabest, objective = objective, hub = ibest))
		
}



############################################################################





########################################
# Method 2: Integer Quadratic Program
########################################



kuroki.matsui.iqp <- function(x, hub = 1L, parallel = FALSE, 
	timeLimit = NULL)
{
	require(gurobi)
	require(slam)
	p <- dim(x)[1]
	m <- dim(x)[2]
	n <- dim(x)[3]
	model <- list()
	
	## Construct components of objective function
	nvar <- (n-1) * m^2
	## Edge weights
	o <- matrix(,n*(n-1)/2,2)
	# o[,1] <- rep(1:(n-1),(n-1):1)
	# o[,2] <- unlist(lapply(2:n, seq.int, to=n))
	w <- array(dim=c(m,m,n*(n-1)/2))
	count <- 0
	nrmx <- colSums(x^2)
	for (i in 1:(n-1))
	for (j in (i+1):n)
	{
		count <- count + 1
		o[count,] <- c(i,j)
		w[,,count] <- matrix(nrmx[,i],m,m) + 
			matrix(nrmx[,j],m,m,byrow=TRUE) - 
			2 * crossprod(x[,,i],x[,,j])
	}
	## Quadratic objective matrix
	idx <- which(o[,1] != hub & o[,2] != hub)
	Qmat <- matrix(0,nvar,nvar)
	for (count in idx) {
		i <- o[count,1] 
		if (i > hub) i <- i - 1
		j <- o[count,2]
		if (j > hub) j <- j - 1
		rowidx <- seq.int((i-1)*m^2+1,i*m^2)
		colidx <- seq.int((j-1)*m^2+1,j*m^2)
		val <- kronecker(w[,,count], diag(0.5,m,m))	
		Qmat[rowidx,colidx] <- val
		Qmat[colidx,rowidx] <- t(val)
	}
	model$Q <- as.simple_triplet_matrix(Qmat)
	rm(Qmat)
	
	## Linear objective vector
	cvec <- numeric(nvar)
	idx <- 	which(o[,1] == hub | o[,2] == hub)
	for (count in idx) {
		i <- setdiff(o[count,],hub)
		if (i < hub) {
			cvec[seq.int((i-1)*m^2+1,i*m^2)] <- as.vector(t(w[,,count]))
		} else {
			i <- i - 1
			cvec[seq.int((i-1)*m^2+1,i*m^2)] <- as.vector(w[,,count])
		} 
	}
	model$obj <- cvec
	rm(w,cvec)

	## Linear constraints
	triplet <- matrix(,2*nvar,3)
	triplet[,1] <- rep(1:(2*m*(n-1)),each=m)
	triplet[,3] <- 1
	triplet[1:nvar,2] <- 1:nvar # sum(u) w(u,v) = 1
	start <- rep(seq.int(0,by=m^2,len=(n-1)), each=m) + (1:m)
	triplet[(nvar+1):(2*nvar),2] <- 
		sapply(start, seq.int, by=m, len=m) # sum(v) w(u,v) = 1
	model$A <- simple_triplet_matrix(triplet[,1],triplet[,2],triplet[,3])
	rm(triplet)
	
	## Set up and solve quadratic problem
	model$rhs <- 1
	model$sense <- "="
	model$vtype <- "B"
	params <- list(NonConvex=2, OutputFlag=0)
	if (!is.null(timeLimit)) 
		params$timeLimit <- timeLimit
	result <- gurobi(model,params)
	sigma <- array(dim=c(m,m,n))
	sigma[,,hub] <- diag(m)
	sigma[,,-hub] <- result$x
	sigma[sigma < 1e-10] <- 0
	sigma[sigma > 0] <- 1
	sigma <- which(aperm(sigma,c(2,1,3)) == 1)
	sigma <- sigma %% m
	sigma[sigma == 0] <- m
	dim(sigma) <- c(m,n)
	out <- list(sigma=sigma, objective=result$objval/n/(n-1))
	return(out)
}




############################################################################



########
# Tests 
########


# library(matchFeats)
# library(slam)
# library(gurobi)

# data(optdigits)
# x <- optdigits$x
# x <- t(x)
# dim(x) <- c(64,10,100)

# m <- 5
# n <- 4
# p <- 10
# x <- array(rnorm(m*n*p),c(p,m,n))


# test1 <- match.bca(x)
# test2 <- kuroki.matsui.ilp(x) 
# test3 <- kuroki.matsui.nqp(x) 



# library(microbenchmark)
# microbenchmark(kuroki.matsui.ilp(x), times=10)

# costfun <- function(x,sigma)
# {
	# n <- ncol(x)
	# p <- dim(x)
	# sumxsigma <- matrix(0,p,m)
	# for (i in 1:n) sumxsigma <- sumxsigma + x[,sigma[,i],i]
	# (sum(x^2) - sum(sumxsigma^2) / n) / (n-1)
# }

# result <- numeric(100)
# for (i in 1:100) 
	# result[i] <- match.bca(x, control=list(sigma=replicate(n,sample(m))))$cost
	



