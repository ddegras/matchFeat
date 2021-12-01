####################################
# Utility functions for simulations 
####################################





########################
# Random starting point
# (permutation)
########################


rand.start <- function(m,n) replicate(n,sample(m))


################## 
# Multi-hub start
################## 

multihub.start <- function(x, unit=NULL)
{
	dimx <- dim(x)
	n <- if (is.null(unit)) dimx[3] else length(unique(unit))
	m <- if (length(dimx) == 3) dimx[2] else dimx[1] / n
	objective.best <- Inf
	for (i in 1:n) {
		out <- matchFeat::match.template(x,i,unit)
		if (out$objective < objective.best) {
			objective.best <- out$objective
			sigma <- out$sigma
		}
	}
	return(sigma)
}


#######################################
# Conversion from cluster indicator
# vector to permutation and vice-versa
#######################################

cluster2perm <- function(x,n)
{
	m <- length(x)/n
	dim(x) <- c(m,n)
	apply(x,2,order)	
}


perm2cluster <- function(sigma)
{
	m <- nrow(sigma)
	n <- ncol(sigma)
	cluster <- matrix(,m,n)
	for (i in 1:n) cluster[sigma[,i],i] <- 1:m
	as.vector(cluster)
}


#######################################
# Convert matrix of permutations 
# to 3D array of permutation matrices
#######################################

# Input: 
# sigma - matrix of permutation functions (each column = 1 permutation)

# Output: 
# 3D array of permutation matrices (clases by classes by units)

perm2arr <- function(sigma)
{
	m <- nrow(sigma)
	n <- ncol(sigma)
	out <- array(0, dim = c(m,m,n))
	idx <- cbind(as.vector(sigma), rep(1:m,n), rep(1:n, each = m))
	out[idx] <- 1
	out
}


##############################################
# Calculate objective function with 3D arrays 
# of data and permutation matrices
##############################################

objective.arr <- function(x,P)
{
	n <- dim(x)[3]
	p <- dim(x)[1] 
	m <- dim(x)[2]
	obj <- 0
	sumxP <- matrix(0,p,m)
	for (i in 1:n)
	{
		xP <- x[,,i] %*% P[,,i]
		obj <- obj + sum(xP^2)
		sumxP <- sumxP + xP
	}
	(obj - sum(sumxP^2)/n) / (n-1)
}





