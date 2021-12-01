match.bca.gen <- function (x, unit = NULL, nclass = NULL, w = NULL, 
	method = c("cyclical", "random"), control = list()) 
{
	stopifnot(!is.null(unit))
	stopifnot(length(unit) == nrow(x))
	stopifnot(is.null(w) || length(w) == 1 || length(w) == ncol(x))
	if (!is.null(w)) 
		stopifnot(all(w >= 0))
	if (!is.integer(unit)) 
		unit <- as.integer(factor(unit))
	if (any(diff(unit) < 0)) 
		stop("'unit' should be a nondecreasing vector.")
	m <- table(unit)
	n <- length(m)
	nfeat <- nrow(x)
	if (is.null(nclass)) 
		nclass <- max(m)
	p <- ncol(x)
	x <- t(x)
	syscall <- sys.call()
	start <- cumsum(c(1, m[-n]))
	len <- pmin(m, nclass)
	pos <- unlist(mapply(seq.int, from = start, length.out = len))
	val <- rep_len(1:nclass, length(pos))
	cluster <- numeric(nfeat)
	cluster[pos] <- val
	rm(start, len, pos, val)
	unit <- split(1:nfeat, unit)
	names(unit) <- NULL
	maxit <- 1000L
	if (is.list(control)) {
		if (!is.null(control$sigma)) 
			cluster <- control$sigma
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
	method <- match.arg(method)
	nrmx2 <- colSums(x^2)
	sumxP <- matrix(0, p, nclass)
	sumnrmx2 <- numeric(nclass)
	size <- integer(nclass)
	for (k in 1:nclass) {
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
 			ai <- if (m[i] <= nclass) { sumnrmx2i } else {
 				sizei * sumnrmx2i - colSums(sumxPi^2) }
			A <- matrix(ai, m[i], nclass, byrow=TRUE) -
				2 * crossprod(x[,unit[[i]],drop=FALSE], sumxPi) + 
				tcrossprod(nrmx2[unit[[i]]], sizei-1) 
 			if (min(A) < 0) 
				A <- A - min(A)
			assigned.old <- assigned
			si.old <- si
			if (m[i] <= nclass) {
				si <- solve_LSAP(A)
			} else {
				si <- numeric(m[i])
				map <- solve_LSAP(t(A))
				si[map] <- 1:nclass
			}
			cluster[unit[[i]]] <- si
			assigned <- which(si > 0)
			si <- si[assigned]
			if (any(assigned.old != assigned) || any(si != si.old)) {
				if (m[i] <= nclass) {
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
			for (k in 1:nclass) 
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
	centers <- matrix(0, p, nclass)
	withinss <- numeric(nclass)
	size <- integer(nclass)
	for (k in 1:nclass) {
		idx <- which(cluster == k)
		size[k] <- length(idx)
		if (size[k] > 0) {
			centers[, k] <- rowMeans(x[, idx, drop = FALSE])
			withinss[k] <- sum(x[, idx]^2) - size[k] * sum(centers[, 
				k]^2)
		}
	}
	xbar <- as.vector(centers %*% (size/sum(size)))
	betweenss <- sum(size * colSums((centers - xbar)^2))
	tot.withinss <- sum(withinss)
	totss <- betweenss + tot.withinss
	out <- list(cluster = cluster, objective = objective, centers = centers, 
		totss = totss, withinss = withinss, tot.withinss = sum(withinss), 
		betweenss = betweenss, size = size, iter = count, call = syscall)
	return(out)
}
