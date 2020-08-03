preprocess <- function(x, unit, mu = NULL, V = NULL, w = NULL)
{

	if (any)
	if (length(dim(x)) < 3) {
	## Case: x specified as matrix
	if (is.vector(x)) dim(x) <- c(length(x),1)
	if (!is.matrix(x)) x <- as.matrix(x)
	if (any(is.infinite(x) | is.na(x)))
		stop(paste("Please make sure that 'x' does not contain missing",
		"or infinite values"))
		# Check that 'unit' is not NULL 
		if (is.null(unit)) 
			stop(paste("Please provide 'x' as an array of dimensions",
				"(p,m,n) OR provide 'x' as a matrix of dimensions (n*p,m)\n",
				"together with the vector 'unit' (length n*p) where m,n,p",
				"are the respective numbers of feature vectors, units,",
				"and variables"))
		unit <- as.integer(unit)
		if (length(unit) == 1) {
			n <- unit
			m <- nrow(x) / n
			if (m != floor(m)) 
				stop(paste("Please make sure that the number of rows in 'x'",
					"is a multiple of 'n'"))			
		} else {
			# Check that dimensions of x and 'unit' are compatible  
			if (length(unit) != nrow(x))
				stop(paste("Please make sure that the number of rows in 'x'",
					"and the length  of 'unit' are equal"))
			# Check that units are balanced 
			# (i.e. have same number of feature vectors)
			if (!all(diff(unit) >= 0))
				stop(paste("Please specify 'unit' as a vector of nondecreasing",
					"integers and sort the rows of 'x' accordingly"))
			freq <- table(unit)
			if (!all(diff(freq) == 0))
				stop(paste("Unbalanced units not supported yet.",
					"Make sure that all unique values in 'unit'",
					"have the same frequency"))
			# Check dimensions of mu and V if specified
			m <- as.integer(freq[1])
			n <- length(freq)
		}
		p <- ncol(x)
		x <- t(x)
	} else {		
	## Case: x specified as array	
		dim.x <- dim(x)
		if (length(dim.x) > 3)
			stop(paste("Please provide 'x' as an array of dimensions (p,m,n)",
				"OR provide 'x' as a matrix of dimensions (n*p,m)\n",
				"together with the vector 'unit' (length n*p) where m,n,p",
				"are the respective numbers of feature vectors, units,",
				"and variables"))
		p <- dim.x[1]
		m <- dim.x[2]
		n <- dim.x[3]
		dim(x) <- c(p,m*n)		
		if (any(is.infinite(x) | is.na(x)))
			stop(paste("Please make sure that 'x' does not contain",
			 	"missing or infinite values"))
	}

	## Check mu and V if provided
	if (!is.null(mu) || !is.null(V)) {
		if (is.null(mu) || is.null(V))	
			stop("Please provide both arguments 'mu' and 'V' or neither")
		if (!((is.vector(mu) && length(mu) == p) || 
			(is.matrix(mu) && dim(mu)[1] == p && dim(mu)[2] %in% c(1,m))))
			stop(paste("Please specify 'mu' as a vector of length p",
				"or as a matrix of dimensions (p,m)\n where p is the number", 
				"of variables and m the number of feature vectors"))
		# Reshape and/or recycle mu if needed
		if (!identical(dim(mu),c(p,m))) mu <- matrix(mu,p,m)
		if (!((length(V) == 1 && p == 1) || 
			(is.matrix(V) && all(dim(V) == c(p,p))) || 
			(length(dim(V)) == 3 && all(dim(V) == c(p,p,m)))))
			stop(paste("Please specify 'V' as a matrix of dimensions (p,p)",
				"or as an array of dimensions (p,p,m)\n where p is the number", 
				"of variables and m the number of feature vectors"))
		# Reshape and/or recycle V if needed 
		if (!(length(dim(V) == 3 && all(dim(V) == c(p,p,m))))) 
			V <- array(V,c(p,p,m))
	}	
	
	## Check w if provided
	R <- NULL
	if (!is.null(w)) {
		err <- paste("Please make sure that 'w' is a vector of", 
			"positive numbers or a positive definite matrix")
		if (is.vector(w) && any(w <= 0)) {
			stop(err)
		} else if (is.matrix(w)) {
			R <- tryCatch(chol(w), error=function(e) stop(err))
		} else stop(err)
	}
		
	return(list(x=x, m=m, n=n, p=p, mu=mu, V=V, R=R))
}