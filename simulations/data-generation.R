

##############
# Import data
##############

train <- read.csv("optdigits.tes", header = FALSE)
x <- as.matrix(train[,1:64])
label <- as.vector(train[,65])
rm(train)


#############################
# Perform PCA for each class
#############################

ulabel <- sort(unique(label))
nlabel <- length(ulabel)
pca <- vector("list",nlabel)
for (k in 1:nlabel) 
	pca[[k]] <- prcomp(x[label == ulabel[k],])
rm(ulabel,nlabel)


########################
## Simulation function
########################

# pca: 	list of PCA objects (one per class), 
# typically obtained by a call to 'prcomp'
# n: number of units to simulate
# k: number of PCs to use for simulation
# sigma: standard deviation of additional white noise
# shuffle: should instances be randomly shuffled? 

simulate.digits <- function(pca, n, k = 25, sigma = 2.5, 
	shuffle = TRUE)
{
	m <- length(pca)
	p <- length(pca[[1]]$sdev)
	x <- matrix(,p,m*n)
	unit <- rep(1:n, each=m)
	label <- rep(0:(m-1),n)
	for (r in 1:m) {
		idx <- seq(r, by=m, len=n)
		rand <- matrix(rnorm(k*n) * pca[[r]]$sdev[1:k],
			k, n)
		x[,idx] <- pca[[r]]$center + 
			pca[[r]]$rotation[,1:k] %*% rand
	}
	if (sigma > 0) 
		x <- x + rnorm(m*n*p, sd=sigma)		
	if (shuffle) {
		dim(x) <- c(p,m,n)
		dim(label) <- c(m,n)
		for (i in 1:n) {
			perm <- sample(m)
			x[,,i] <- x[,perm,i]
			label[,i] <- label[perm,i]
		}
		dim(x) <- c(p,m*n)
		dim(label) <- NULL
	}
	return(list(x=t(x), unit=unit, label=label))	
}


