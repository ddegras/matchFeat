match.bca.gen <- function (x, unit = NULL, cluster = NULL, w = NULL, 
	method = c("cyclical", "random"), control = list()) 
{
	## Check and reformat input arguments
	stopifnot(length(unit) == nrow(x))
	stopifnot(is.null(w) || length(w) == 1 || length(w) == ncol(x))
	if (!is.null(w)) 
		stopifnot(all(w >= 0))
	if (!is.integer(unit)) 
		unit <- unclass(factor(unit))
	m <- table(unit)
	n <- length(m)
	unit <- split(1:nrow(x), unit)
	names(unit) <- NULL	
	err <- paste("'cluster' must be a single integer or an",
		"integer vector of length equal to the number of rows of 'x'")
	if (!is.numeric(cluster)) stop(err)
	cluster <- as.integer(cluster)
	if (length(cluster) == 1) {
		ncluster <- cluster
		cluster[unlist(unit)] <- rep_len(1:ncluster, nrow(x))
	} else if (length(cluster) == nrow(x)) {
		ncluster <- length(unique(cluster))
		if (!all(cluster %in% 1:ncluster)) 
		stop(paste("The values in 'cluster' must be integers between 1",
			"and the number of clusters"))
	} else stop(err)
	if (ncluster < max(m))
		stop(paste0("The number of clusters (", ncluster, ") is too small ",
			"compared to the largest number of observations per subject (",
			max(m), "). Please set the input 'cluster' to an integer ",
			"at least ", max(m), " or to an integer vector with at least ", 
			max(m), " different values."))
	p <- ncol(x)
	x <- t(x)
	syscall <- sys.call()
	maxit <- 1000L
	if (is.list(control)) {
		if (!is.null(control$maxit)) 
			maxit <- control$maxit
	}
	method <- match.arg(method)
	if (!is.null(w)) {
		if (is.vector(w)) {
			x <- sqrt(w) * x
		}
		else {
			R <- chol(w)
			x <- R %*% x
		}
	}
	nrmx2 <- colSums(x^2)
	sumxP <- matrix(0, p, ncluster)
	sumnrmx2 <- numeric(ncluster)
	size <- integer(ncluster)
	for (k in 1:ncluster) {
		idx <- which(cluster == k)
		size[k] <- length(idx)
		if (size[k] == 0) 
			next
		sumxP[, k] <- rowSums(x[, idx, drop = FALSE])
		sumnrmx2[k] <- sum(nrmx2[idx])
	}
	objective <- 2/(n * (n - 1)) * (sum(size * sumnrmx2) - sum(sumxP^2))
	if (method == "cyclical") 
		sweep <- 1:n
	for (count in 1:maxit) {
		objective.old <- objective
		objective.after <- NULL
		if (method == "random") 
			sweep <- sample(1:n)
		for (i in sweep) {
			si <- cluster[unit[[i]]]
			assigned <- which(si > 0)
			si <- si[assigned]
			sumnrmx2i <- sumnrmx2
			sumnrmx2i[si] <- sumnrmx2i[si] - nrmx2[unit[[i]][assigned]]
			sizei <- size
			sizei[-si] <- sizei[-si] + 1L
			sumxPi <- sumxP
			sumxPi[, si] <- sumxPi[, si] - x[, unit[[i]][assigned]]			
 			ai <- if (m[i] <= ncluster) { sumnrmx2i } else {
 				sizei * sumnrmx2i - colSums(sumxPi^2) }
			A <- matrix(ai, m[i], ncluster, byrow=TRUE) -
				2 * crossprod(x[,unit[[i]],drop=FALSE], sumxPi) + 
				tcrossprod(nrmx2[unit[[i]]], sizei-1) 
 			if (min(A) < 0) 
				A <- A - min(A)
			assigned.old <- assigned
			si.old <- si
			if (m[i] <= ncluster) {
				si <- solve_LSAP(A)
			} else {
				si <- numeric(m[i])
				map <- solve_LSAP(t(A))
				si[map] <- 1:ncluster
			}
			cluster[unit[[i]]] <- si
			assigned <- which(si > 0)
			si <- si[assigned]
			if (any(assigned.old != assigned) || any(si != si.old)) {
				if (m[i] <= ncluster) {
				  sumxP[,-si] <- sumxPi[,-si]
				  sumxP[, si] <- sumxPi[, si] + x[, unit[[i]]]
				}
				else {
				  sumxP <- sumxPi + x[, unit[[i]][map]]
				}
				sumnrmx2 <- sumnrmx2i
				sumnrmx2[si] <- sumnrmx2[si] + nrmx2[unit[[i]][assigned]]
				size[si.old] <- size[si.old] - 1
				size[si] <- size[si] + 1
			}
		}
		if (count%%10 == 0) {
			for (k in 1:ncluster) 
				sumxP[,k] <- rowSums(x[, cluster == k, drop=FALSE])
		}
		objective <- 2/(n*(n-1)) * (sum(size*sumnrmx2) - sum(sumxP^2))
		if (objective >= objective.old) 
			break
	}
	if (!is.null(w)) {
		if (is.vector(w)) {
			invsqrtw <- 1/sqrt(w)
			invsqrtw[is.infinite(invsqrtw)] <- 0
			x <- x * invsqrtw
		}
		else {
			x <- backsolve(R, x)
		}
	}
	mu <- matrix(0, p, ncluster)
	V <- array(0, c(p, p, ncluster))
	# withinss <- numeric(ncluster)
	size <- integer(ncluster)
	for (k in 1:ncluster) {
		idx <- which(cluster == k)
		size[k] <- length(idx)
		if (size[k] > 0) {
			mu[, k] <- rowMeans(x[, idx, drop = FALSE])
			if (size[k] > 1)
			V[,,k] <- tcrossprod(x[, idx] - mu[,k]) / size[k]
			# withinss[k] <- sum(x[, idx]^2) - size[k] * sum(centers[, 
				# k]^2)
		}
	}
	# xbar <- as.vector(centers %*% (size/sum(size)))
	# betweenss <- sum(size * colSums((centers - xbar)^2))
	# tot.withinss <- sum(withinss)
	# totss <- betweenss + tot.withinss
	out <- list(cluster = cluster, objective = objective, 
		mu = mu, V = V, size = size, call = syscall)
		# totss = totss, withinss = withinss, 
		# tot.withinss = sum(withinss), betweenss = betweenss, 
	class(out) <- "matchFeat"
	return(out)
}
