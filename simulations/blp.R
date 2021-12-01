

##################################################
# BINARY LINEAR PROGRAMMING IMPLEMENTATION
# OF MULTIDIMENSIONAL ASSIGNMENT PROBLEM WITH 
# DECOMPOSABLE COSTS AKA MINIMUM WEIGHT CLIQUE 
# COVER IN COMPLETE N-PARTITE GRAPH
##################################################



 


############
# Version 1
############


# n x (n-1) / 2 x m^2 variables 
# edges e(i,j,k,l) with 1 <= i < j <= m, k=1:m, l=1:m

# 2 x (n-1) x m equalities (one-to-one map between clusters 1 and j for all j>1)
# sum(k) e(1,j,k,l) = 1 for all j,l (j=2:n)
# sum(l) e(1,j,k,l) = 1 for all j,k (j=2:n)

# n x (n-1) / 2 x m^3 inequalities (clique property)
# If node k from cluster i and node l from cluster j 
# are both connected to node r from cluster 1 (with 1,i,j all distinct)
# then these two nodes are also connected to each other


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
###### THERE SEEM TO BE NEW BUGS ########
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



blp.solve1 <- function(x, timeout = NULL, solver=c("gurobi","glpk"))
{
	
	## Model dimensions
	m <- dim(x)[2]
	n <- dim(x)[3]
	p <- dim(x)[1]
	
	## Solver
	solver <- match.arg(solver)
	
	
	#########################
	# Objective coefficients 
	#########################
	
	# Squared Euclidean distances between vectors
	obj <- array(dim=c(m,m,n,n))
	for (i in 1:(n-1))
	for (j in (i+1):n)
	for (l in 1:m)
		obj[,l,j,i] <- colSums((x[,,i]-x[,l,j])^2)
	dim(obj) <- NULL
	obj <- na.omit(obj)
	nvar <- length(obj)
	# Indexing: (k,l,j,i), k moves fastest, then l, then j, then i	
	
	
	#######################
	# Equality constraints 
	#######################

	# For every node in every cluster, 
	# exactly one edge with every other cluster
	
	neq <- n*(n-1)*m # number of equality constraints
	a <- matrix(,neq*m,3) # triplet matrix
	# Each equality constraint involves m coefficients --> m*neq rows in 'a' 
	a[,1] <- rep(1:neq, each=m) # row indices (1 row = 1 constraint)
	# 1:nvar --> constraint sum(k) e(k,l,j,i) = 1 for all (i,j,l)
	a22 <- array(1:nvar,dim=c(m,m,n*(n-1)/2)) 
	# constraint sum(l) e(k,l,j,i) = 1 for all (i,j,k)
	a22 <- aperm(a22,c(2,1,3))
	a[,2] <- c(1:nvar,a22) # column indices 
	# Note: each edge coefficient appears exactly once in all equality constraints
	a[,3] <- 1 # constraint coefficients
	a <- slam::simple_triplet_matrix(a[,1],a[,2],a[,3],neq,nvar)

		
		
	#########################	
	# Inequality constraints
	#########################	

	## Given three nodes in three distinct clusters, 
	## if there are (at least) two edges between these 
	## nodes, then there are three (transitivity/clique property)
	nineq <- 3 * choose(n,3) * m^3
	b <- matrix(,3*nineq,3)
	b[,1] <- rep(1:nineq,each=3)
	b[,3] <- rep(c(1,1,-1),nineq)
	unitslow <- rep.int(1:(n-1),(n-1):1)
	unitfast <- unlist(lapply(2:n, function(i) i:n))
	blocksize <- 9 * m^3
	
	# Node triplets
	nodes <- matrix(,m^3,3)
	nodes[,1] <- rep(1:m,m^2) # k
	nodes[,2] <- rep(1:m,each=m) # l 
	nodes[,3] <- rep(1:m,each=m^2) # r
	
	# Corresponding edge triplets (indices in {1,...,m^2})
	edges <- matrix(,m^3,3)
	edges[,1] <- m * (nodes[,2]-1) + nodes[,1] # (k,l)
	edges[,2] <- m * (nodes[,3]-1) + nodes[,1] # (k,r)
	edges[,3] <- m * (nodes[,3]-1) + nodes[,2] # (l,r)
	colnames(edges) <- c("kl","kr","lr")	
	rm(nodes)
		
	combs <- combn(n,3) # combinations of three units among n
	for (comb in 1:ncol(combs))
	{
		# Unit indices
		i <- combs[1,comb]
		j <- combs[2,comb]
		q <- combs[3,comb]
		
		# Column indices for pairs of units/clusters
		pair <- c(ij = which(unitslow == i & unitfast == j),
			iq = which(unitslow == i & unitfast == q),
			jq = which(unitslow == j & unitfast == q))
			
		# If edges between (i,k)-(j,l) and between (j,l)-(q,r),
		# then edge between (i,k)-(q,r)	
		mat1 <- rbind(
			(pair["ij"]-1) * m^2 + edges[,"kl"], # coef 1
			(pair["jq"]-1) * m^2 + edges[,"lr"], # coef 1 
			(pair["iq"]-1) * m^2 + edges[,"kr"]) # coef -1
	
		# If edges between (i,k)-(j,l) and between (i,k)-(q,r),
		# then edge between (i,k)-(q,r)	
		mat2 <- rbind(
			(pair["ij"]-1) * m^2 + edges[,"kl"], # coef 1
			(pair["iq"]-1) * m^2 + edges[,"kr"], # coef 1 
			(pair["jq"]-1) * m^2 + edges[,"lr"]) # coef -1
	
		# If edges between (i,k)-(q,r) and between (j,l)-(q,r),
		# then edge between (i,k)-(j,l)	
		mat3 <- rbind(
			(pair["iq"]-1) * m^2 + edges[,"kr"], # coef 1
			(pair["jq"]-1) * m^2 + edges[,"lr"], # coef 1 
			(pair["ij"]-1) * m^2 + edges[,"kl"]) # coef -1
	
		idx <- seq.int((comb-1)*blocksize+1,comb*blocksize)
		b[idx,2] <- c(mat1,mat2,mat3)
	}
	rm(mat1,mat2,mat3)
	b <- slam::simple_triplet_matrix(b[,1],b[,2],b[,3],nineq,nvar)
	
	
	########################
	# Binary linear program
	########################

	constraints <- rbind(a,b)
	rhs <- rep(1,neq+nineq)
	
	if (solver == "glpk") {
		result <- Rglpk_solve_LP(obj=obj, mat=constraints, 	
			dir=rep(c("==","<="), c(neq,nineq)), rhs=rhs, 
			types="B", control=list(tm_limit=1000*timeout))
	} else {	
		model <- list(obj=obj, A=constraints, rhs=rhs, 
			sense=rep(c("=","<="), c(neq,nineq)), 
			modelsense="min", vtype="B")
		params <- list(OutputFlag=0)
		if (!is.null(timeout))
			params$timeLimit <- timeout
		result <- gurobi(model,params)
		idx <- match(c("x","objval"), names(result))
		names(result)[idx] <- c("solution","optimum")
		
	}
	
	
	#######################################
	# Recover solution as clique partition 
	# and permutations
	#######################################
	
	# The set of edges starting from unit 1 (or any unit) 
	# is sufficient to reconstruct the m cliques or equivalently, 
	# the m permutations of feature vectors

	cluster <- integer(m*n)	
	## All edges startinig from unit/subgraph 1
	idx <- which(head(result$solution,(n-1)*m^2) == 1)
	for (k in 1:m) {
		# Indices (in edge list) of all nodes 
		# linked with node k 
		idx2 <- idx[edges[idx,1] == k]
		# Cluster assignment
		cluster[c(k,edges[idx2,2])] <- k
	}	
	sigma <- apply(matrix(cluster,m,n),2,order)

	return(list(sigma=sigma, cluster=cluster, 
		objective=result$optimum/n/(n-1)))
		
		
}	




##################################################################################




############
# VERSION 2
############


# n x (n-1) / 2 x m^2 variables 
# edges e(i,j,k,l) with 1 <= i < j <= n, k=1:m, l=1:m

# 2 x (n-1) x m equalities (one-to-one map between clusters 1 and j for all j>1)
# sum(k) e(1,j,k,l) = 1 for all j,l (j=2:n)
# sum(l) e(1,j,k,l) = 1 for all j,k (j=2:n)

# n x (n-1) / 2 x m^3 inequalities (clique property)
# If node k from cluster i and node l from cluster j 
# are both connected to node r from cluster 1 (with 1,i,j all distinct)
# then these two nodes are also connected to each other



blp.solve2 <- function(x)
{
	
	## Model dimensions
	m <- dim(x)[2]
	n <- dim(x)[3]
	p <- dim(x)[1]
	
	
	#########################
	# Objective coefficients 
	#########################
	
	# Squared Euclidean distances between vectors
	obj <- array(dim=c(m,m,n,n))
	for (i in 1:(n-1))
	for (j in (i+1):n)
	for (l in 1:m)
		obj[,l,j,i] <- colSums((x[,,i]-x[,l,j])^2)
	dim(obj) <- NULL
	obj <- na.omit(obj)
	nvar <- length(obj)
	# Indexing: (k,l,j,i), k moves fastest, then l, then j, then i	
	
	
	#######################
	# Equality constraints 
	#######################

	# For every node in cluster 1, 
	# exactly one edge with every other cluster
	neq <- 2*(n-1)*m # number of equality constraints
	a <- matrix(,neq*m,3) # triplet matrix
	# Each equality constraint involves m coefficients --> m*neq rows in 'a' 
	a[,1] <- rep(1:neq, each=m) # row indices (1 row = 1 constraint)
	idx <- seq.int(1,(n-1)*m^2)
	a[1:(m*neq/2),2] <- idx
	dim(idx) <- c(m,m,n-1) # index: (k,l,j)
	# constraint sum(l) e(k,l,j,1) = 1 for all (j,k)
	a[(m*neq/2+1):(m*neq),2] <- aperm(idx,c(2,1,3)) # index: (l,k,j)
	a[,3] <- 1 # constraint coefficients
	a <- slam::simple_triplet_matrix(a[,1],a[,2],a[,3],neq,nvar)
		
		
	#########################	
	# Inequality constraints
	#########################	

	## If two nodes from clusters i and j (1,i,j are distinct) are connected 
	## to the same node in cluster 1, then these two nodes are connected to 
	## each other (transitivity/clique property)
	nineq <- choose(n-1,2) * m^3
	b <- matrix(,3*nineq,3)
	b[,1] <- rep(1:nineq,each=3)
	b[,3] <- rep(c(1,1,-1),nineq)
	unitslow <- rep.int(1:(n-1),(n-1):1)
	unitfast <- unlist(lapply(2:n, function(i) i:n))
	blocksize <- 3 * m^3
	
	# Node triplets
	nodes <- matrix(,m^3,3)
	colnames(nodes) <- c("k","l","r")
	nodes[,"k"] <- rep(1:m,m^2) # k
	nodes[,"l"] <- rep(1:m,each=m) # l 
	nodes[,"r"] <- rep(1:m,each=m^2) # r
	
	# Corresponding edge triplets (indices in {1,...,m^2})
	edges <- matrix(,m^3,3)
	colnames(edges) <- c("kl","kr","lr")
	tmp <- matrix(1:(m^2),m,m)
	edges[,"kl"] <- tmp[nodes[,c("k","l")]]
	edges[,"kr"] <- tmp[nodes[,c("k","r")]]
	edges[,"lr"] <- tmp[nodes[,c("l","r")]]
	rm(nodes)
	count <- 0 # block counter in b 		
	for (i in 2:(n-1))
	for (j in (i+1):n)
	{
		# Column indices for pairs of units/clusters
		colidx <- which(unitslow == i & unitfast == j)
			
		# If edges between (1,k)-(i,l) and between (1,k)-(j,r),
		# then edge between (i,k)-(j,r)	
		mat <- rbind((i-2) * m^2 + edges[,"kl"], # coef 1
			(j-2) * m^2 + edges[,"kr"], # coef 1 
			(colidx-1) * m^2 + edges[,"lr"]) # coef -1
		
		idx <- seq.int(count*blocksize+1,(count+1)*blocksize)
		b[idx,2] <- mat
		count <- count + 1
	}
	rm(mat)
	b <- slam::simple_triplet_matrix(b[,1],b[,2],b[,3],nineq,nvar)

	
	
	########################
	# Binary linear program
	########################

	constraints <- rbind(a,b)
	rhs <- rep(1,neq+nineq)
	direction <- rep(c("==","<="),c(neq,nineq))
	result <- Rglpk::Rglpk_solve_LP(obj=obj, mat=constraints, 
		dir=direction, rhs=rhs, types="B")
	
	# The set of edges starting from unit 1 (or any unit) 
	# is sufficient to reconstruct the m cliques or equivalently, 
	# the m permutations of feature vectors
	edge1 <- which(result$solution[1:((n-1)*m^2)] == 1)
	idx1 <- arrayInd(edge1,c(m,m,n-1))
	o <- order(idx1[,3],idx1[,1])
	sigma <- matrix(c(1:m,idx1[o,2]),m,n)
	objective <- result$optimum / (2*n)
	return(list(sigma=sigma, objective = objective))
}	





############
# VERSION 3
############


# n * (n-1) / 2 * m^2 variables (edge indicator variables)
# edges e(i,j,k,l) with 1 <= i < j <= n, k=1:m, l=1:m

# n * (n-1) * m linear equality constraints
# Each node of each cluster/subgraph is connected to 
# exactly one node in each other subgraph/cluster

# n * (n-1) * (n-2) / 6 * m^3 linear inequality constraints 
# (clique property)
# If node q from cluster i and node r from cluster j 
# are both connected to node s from cluster k (i,j,k distinct)
# then (q,i) and (r,j) are connected to each other

# This is the "full" version of the BLP formulation of the MDADC
# No cluster is used as a hub, therefore the number of clique 
# constraints is much higher (O(n^3 m^3) vs O(n^2 m^3))

# The code below is much more transparent than blp.solve1 and blp.solve2
# It uses (the R interface to) Gurobi but can easily be ported to 
# RGplk if needed


blp.solve3 <- function(x, timeout=NULL, solver=c("gurobi","glpk"))
{
	
	## Problem dimensions
	m <- dim(x)[2]
	n <- dim(x)[3]
	p <- dim(x)[1]
	
	## Solver
	solver <- match.arg(solver)
	
	## Nodes and edges
	nc2 <- choose(n,2)
	mc2 <- choose(m,2)
	nedges <- nc2 * m^2
	nodes <- matrix(1:(m*n),m,n)
	edges <- matrix(,nedges,2)
	count <- 0
	for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		idx <- (count+1):(count+m^2)
		edges[idx,1] <- rep(nodes[,i],each=m)
		edges[idx,2] <- rep(nodes[,j],m)
		count <- count + m^2
	}}
	
	## Distances aka costs
	cost <- numeric(nedges)
	sqnrmx <- colSums(x^2)
	dim(x) <- c(p,m*n)
	for (e in 1:nedges)		
		cost[e] <- sqnrmx[edges[e,1]] + sqnrmx[edges[e,2]] - 
			2 * sum(x[,edges[e,1]] * x[,edges[e,2]])
	
	
	#######################
	# Equality constraints 
	#######################

	# Each node of each subgraph connected to exactly one node 
	# in every other subgraph --> C(n,2) * 2m = n * (n-1) * m
	neq <- nc2 * 2 * m 
	
	## Lists of edges between pairs of subgraphs 	
	edges2 <- vector("list",nc2)
	
	A <- matrix(,m,neq) # locations of 1's for each constraint
	count.eq <- 0
	count.edge <- 1
	for (i in 1:(n-1)) {
	for (j in (i+1):n) {
		## Corresponding edges 
		idx <- ((count.edge-1)*m^2+1):(count.edge*m^2)
		## Each node in s(i) connected to exactly one node in s(j) 
		A[,(count.eq+1):(count.eq+m)] <- idx[order(edges[idx,1])]
		## Each node in s(j) connected to exactly one node in s(i) 
		A[,(count.eq+m+1):(count.eq+2*m)] <- idx[order(edges[idx,2])]

		edges2[[count.edge]] <- idx
		count.eq <- count.eq + 2*m
		count.edge <- count.edge + 1 
	}}	

	
	
	#########################
	# Inequality constraints
	#########################
	
	# Clique constraints on triangles
	# If node i in subgraph q (i,q) is connected to node j in subgraph r
	# (j,r) and to node k in subgraph s (k,s), then (j,r) is connected to 
	# to (k,s) 
	
	## Matrix to map pairs of indices (i,j) to linear index in 'edges2'
	subgpairs <- matrix(,n,n)
	subgpairs[lower.tri(subgpairs)] <- 1:nc2
	subgpairs <- t(subgpairs)
	
	nineq <- choose(n,3) * m^3
	A2 <- matrix(,3,nineq)
	count <- 0
	for (i in 1:(n-2)) {
	for (j in (i+1):(n-1)) {
	for (k in (j+1):n) {
		
		idx <- (count+1):(count+m^3)	
		idx.ij <- edges2[[subgpairs[i,j]]]
		idx.ik <- edges2[[subgpairs[i,k]]]
		idx.jk <- edges2[[subgpairs[j,k]]]
		idx.qr <- rep(idx.ij,each=m)
		tmp <- matrix(idx.ik,m,m)
		idx.qs <- as.vector(matrix(t(tmp),m^2,m,byrow=TRUE))
		idx.rs <- rep(idx.jk,m)
		A2[1,idx] <- idx.qr
		A2[2,idx] <- idx.qs
		A2[3,idx] <- idx.rs		
		count <- count + m^3		
	}}}
		
	# Acopy <- A
	## Combine equality and inequality constraints
	A <- slam::simple_triplet_matrix(
		i=c(rep(1:neq, each=m),rep((neq+1):(neq+nineq), each=3)),
		j=c(as.vector(A),as.vector(A2)), 
		v=c(rep(1,neq*m),rep(c(1,1,-1),nineq)),
		nrow=neq+nineq, ncol=nedges)
	

	
	########################
	# Binary linear program 
	########################

	if (solver == "glpk") {
		result <- Rglpk_solve_LP(obj=cost, mat=A, 	
			dir=rep(c("==","<="), c(neq,nineq)), 
			rhs=rep(1,neq+nineq), types="B",
			control=list(tm_limit=1000*timeout))
	} else {	
		model <- list(obj=cost, A=A, rhs=1, 
			sense=rep(c("=","<"), c(neq,nineq)),
			modelsense="min", vtype="B")
		params <- list(OutputFlag=0)
		if (!is.null(timeout))
			params$timeLimit <- timeout
		result <- gurobi(model,params)
		idx <- match(c("x","objval"), names(result))
		names(result)[idx] <- c("solution","optimum")
		
	}
	
	#######################################
	# Recover solution as clique partition 
	# and permutations
	#######################################
	
	# The set of edges starting from unit 1 (or any unit) 
	# is sufficient to reconstruct the m cliques or equivalently, 
	# the m permutations of feature vectors
	
	cluster <- integer(m*n)	
	## All edges startinig from unit/subgraph 1
	idx <- which(head(result$solution,(n-1)*m^2) == 1)
	for (k in 1:m) {
		# Indices (in edge list) of all nodes 
		# linked with node k 
		idx2 <- idx[edges[idx,1] == k]
		# Cluster assignment
		cluster[c(k,edges[idx2,2])] <- k
	}	
	sigma <- apply(matrix(cluster,m,n),2,order)

	return(list(sigma=sigma, cluster=cluster, 
		objective=result$optimum/n/(n-1)))

}


##################################################################################




########
# Tests 
########


# ## Required packages
# library(slam)
# library(Rglpk)
# library(gurobi)
# library(matchFeat)

# Note: Gurobi is much faster than GLPK

## Generate random data
# m <- 5
# n <- 10
# p <- 5
# x <- runif(m*n*p)
# dim(x) <- c(p,m,n)

# ## Timing
# system.time(res1 <- match.bca(x))  # fastest 
# system.time(res2 <- blp.solve1(x, timeout=10)) # too slow as soon as m and/or n > 10 or so
# system.time(res3 <- blp.solve2(x)) # even slower than blp.solve1
# system.time(res4 <- blp.solve3(x, timeout=30, solver="gurobi")) 



## Precise timing
# library(microbenchmark)
# microbenchmark(match.bca(x), blp.solve1(x), blp.solve2(x), times=10)

# ## Objective
# res1$objective
# res2$objective 
# res3$objective 
# res4$objective 

# ## Permutation
# res1$sigma
# res2$sigma
# res3$sigma

# ## Use multiple random starts for BCA method
# system.time(res1 <- replicate(100, match.bca(x, method="random"), simplify=FALSE))
# best <- which.min(sapply(res1,"[[","cost"))
# res1 <- res1[[best]]

# ## Check objective
# # sigma <- res2$sigma
# # xmatched <- array(dim=c(p,m,n))
# # for (i in 1:n) xmatched[,,i] <- x[,sigma[,i],i]
# # xbar <- rowMeans(xmatched,dims=2)
# # sum((xmatched-as.vector(xbar))^2)



