########################################################
# CONVEX RELAXATION OF THE FEATURE MATCHING PROBLEM (7) 
# OF DEGRAS (2021) "SCALABLE FEATURE MATCHING ACROSS 
# LARGE DATA COLLECTIONS"
########################################################


# (7) min F(P1,...,Pn) = sum(i,j) || Xi Pi - Xj Pj ||_F^2 / n / (n-1)
# where i,j vary in 1,...,n, X1,...,Xn are p x m matrices, 
# and P1,...,Pn are m x m permutation matrices
# Convex relaxation: P1,...,Pn as doubly stochastic matrices
# The function lb.fw implements the Frank-Wolfe algorithm 
# for the convex relaxation of (7)


# Inputs: 
# x - 3D array (variables by classes by units)
# P - optional 3D array (of permutation matrices classes by classes by units) 
# maxit - maximum number of iterations 
# tol - tolerance for relative error in convergence
# verbose - if TRUE, display algorithm progress

lb.fw <- function(x, P = NULL, maxit = 1000, tol = 1e-6, verbose = FALSE)
{

	## Data dimensions
	m <- dim(x)[2]
	n <- dim(x)[3]
	p <- dim(x)[1]
	
	## Initialize solution	
	if (is.null(P)) 
		P <- array(diag(m), dim = c(m,m,n))

	## Initialize objective	
	sumxP <- matrix(0,p,m)
	objective <- 0
	for (i in 1:n) {
		xPi <- x[,,i] %*% P[,,i]
		sumxP <- sumxP + xPi
		objective <- objective + sum(xPi^2)
	}
	
	objective <- (objective - sum(sumxP^2)/n) / (n-1)
	objective.best <- objective
			
	## Search direction
	Q <- matrix(,m,n)
	
	for (it in 1:maxit)
	{	
		## Store previous objective
		objective.old <- objective
		
		## Find search direction
		sumxP <- matrix(0,p,m)
		for (i in 1:n) sumxP <- sumxP + x[,,i] %*% P[,,i]
		sumxQ <- matrix(0,p,m)
		sum.nrm.xQP <- 0
		sum.cp.xP.xQP <- 0
		for (i in 1:n) {
			xPi <- x[,,i] %*% P[,,i]
			tgrad <- crossprod(xPi - sumxP/n, x[,,i])
			Q[,i] <- clue::solve_LSAP(tgrad - min(tgrad))
			sumxQ <- sumxQ + x[,Q[,i],i]
			xQPi <- x[,Q[,i],i] - xPi 
			sum.cp.xP.xQP <- sum.cp.xP.xQP + sum(xPi*xQPi)
			sum.nrm.xQP <- sum.nrm.xQP + sum(xQPi^2)	 
		}

		## Calculate step size by line search
		a <- (n * sum.nrm.xQP - sum((sumxQ-sumxP)^2)) / n / (n-1)
		b <- (n * sum.cp.xP.xQP - sum(sumxP * (sumxQ-sumxP))) / n / (n-1)
		alpha <- min(max(-b/a,0),1)
				
		## Update solution
		P <- (1-alpha) * P
		idx <- cbind(as.vector(Q),rep(1:m,n),rep(1:n,each=m))
		P[idx] <- P[idx] + alpha

		## Calculate objective 
		objective <- objective - b^2/a
		if (verbose) 
			cat("Iteration",it,"Objective",objective,"\n")
		
		## Monitor progress
		if (objective >= (1-tol) * objective.old) {
			count <- count + 1
		} else { count <- 0 }
		if (count == 10) break 
				
	}

 	## Round the doubly stochastic solutions to recover permutations 
	sigma <- apply(P, 3, clue::solve_LSAP, maximum = TRUE)
		
	list(sigma = sigma, objective = objective.fun(x, sigma), 
		P = P, lb = objective)

}


