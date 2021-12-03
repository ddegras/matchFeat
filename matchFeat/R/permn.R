## Recursive function to create all permutations of {1,...,k}
permn <- function(k)
{
	k <- as.integer(k)
	if (k <= 0) return(NULL)
	factkm1 <- factorial(k-1)
	factk <- k*factkm1
	out <- matrix(k,factk,k)
	tmp <- permn(k-1)
	for (i in 1:k) 
		out[seq.int((i-1)*factkm1+1,length.out=factkm1),-i] <- tmp	
	return(out)	
}


## Brute force determination of the best permutation 
## for matching two sets of feature vectors
brute <- function(x, y, perms)
{
	k <- NCOL(x)
	cost <- matrix(,k,k)
	for (j in 1:k)
		cost[,j] <- colSums((x-y[,j])^2)
	sumcost <- lapply(perms, 
		function(idx) sum(cost[cbind(idx,1:k)]))
	return(perms[[which.min(sumcost)]])	
}
