############################################################
# CODE FOR THE SIMULATIONS OF "SCALABLE FEATURE MATCHING 
# ACROSS LARGE DATA COLLECTIONS" (DEGRAS, 2021) 
############################################################


## Clear workspace
rm(list=ls())


#########################
# Load required packages
#########################

## Install CRAN packages as needed
required.pkgs <- c("conclust", "devtools", "doParallel",
	"foreach", "R.utils", "Rglpk", "slam")
pkgs.to.install <- setdiff(required.pkgs, installed.packages()[,1])
if (length(pkgs.to.install) > 0) install.packages(pkgs.to.install)

## Load CRAN packages
library(conclust)
library(devtools)
library(foreach)
library(R.utils, quietly = TRUE, warn.conflicts = FALSE, verbose = FALSE)
library(doParallel)
registerDoParallel() # change the number of cores as needed	

## Install package 'matchFeat' if needed and load it  
test <- require(matchFeat)
if (!test) devtools::install_github("ddegras/matchFeat/matchFeat")
library(matchFeat)

## Try loading 'Rglpk' and 'gurobi'
gurobi.flag <- require(gurobi)
if (!gurobi.flag) warning(
	paste("It is strongly recommended to install Gurobi and its R interface", 
	"before running the simulations. Without it, several methods (pairwise",
	"interchange, integer quadratic programming) will not run. Gurobi can be", 
	"downloaded for free after obtaining an academic license at", 
	"https://www.gurobi.com"))
glpk.flag <- require(Rglpk)
if (!glpk.flag) warning(
	paste("Package 'Rglpk' failed to load. Please make sure that GLPK is installed.", 
	"GLPK can be downloaded at https://www.gnu.org/software/glpk"))
lp.pkg <- if (gurobi.flag) "gurobi" else if (glpk.flag) "glpk" else NULL


##########################
# Source simulation files
##########################

source("data-generation.R")
source("kuroki-matsui.R")
source("blp.R")
source("utils.R")


########################
# Simulation parameters
########################

## Number of units
nvals <- c(5,10,20,30,40,50,75,100,200,500,1000)
## Number of classes/clusters
m <- 10
## Additive noise level 
stdev <- 2.5
## Number of PCs to use in simulation
k <- 25
## Number of replications
nrep <- 100


##############
# Simulations 
##############

## Method names 
start <- c("id","r100","rec","hub","label")
method <- c(start, "2x", paste0(start,".kmm"), 
	paste0(start,".bca"),  paste0(start,".kmm.2x"), 
	paste0(start,".bca.2x"), "gaussmix", "ilp", "kur.ilp",  
	"kur.iqp", "ccls", "copkm", "lcvqe", "mpckm")
if (!gurobi.flag) 
	method <- grep("2x|iqp", method, value = TRUE, invert = TRUE)
nmethod <- length(method)

for (n in nvals)
{
	cat("\nn =",n,"\n")
	## "Must" and "cannot" links for constrained clustering methods  
	mustlink <- t(c(1,1))
	mC2 <- choose(m,2)
	cantlink <- matrix(,n*mC2,2)
	for (i in 1:n) {
		idx1 <- seq.int((i-1)*mC2+1,i*mC2)
		idx2 <- seq.int((i-1)*m+1,i*m)
		cantlink[idx1,] <- t(combn(idx2,2))
	}
	
	#### MAIN COMPUTATIONS 	
	result <- foreach(i = 1:nrep, .errorhandling = "pass") %dopar% 
	# For serial computations, comment out the line above  
	# and uncomment the 2 lines below (and at the end of the loop)
	# result <- vector("list",nrep)  
	# for (i in 1:nrep)
	{
	
		out <- matrix(,nmethod,3)
		rownames(out) <- method
		colnames(out) <- c("objective","Rand","time")
		timeout <- 300
		
		# Simulate data 
		set.seed(i)
		sim <- simulate.digits(pca[1:m], n, k, stdev)
		x <- sim$x
		label <-sim$label
		unit <- sim$unit
		p <- ncol(x)
		xx <- t(x)
		dim(xx) <- c(p, m, n)
		
		# Starting points
		sigma <- list(
			id = matrix(1:m, m, n), 
			r100 = replicate(100, rand.start(m, n)), 
			rec = match.rec(xx)$sigma, 
			hub = multihub.start(xx),
			label = cluster2perm(label, n)) 
		for (name in names(sigma)) {
			if (name == "r100") {
				objective.best <- Inf
				for (j in 1:100)	 {
					objective <- objective.fun(xx, sigma$r100[,,j])
					if (objective < objective.best) {
						objective.best <- objective
						idx.best <- j
					}
				}
				out[name,] <- c(objective.best, Rand.index(label,
							perm2cluster(sigma$r100[,,idx.best])), 0)
			} else {
			out[name,] <- c(objective.fun(xx, sigma[[name]]), 
				Rand.index(label, perm2cluster(sigma[[name]])), 0)
			}
		}
		 	
		## Pairwise interchange (too slow when starting from scratch)
		if (n <= 30 && gurobi.flag) {
			dt <- system.time(res <- match.2x(xx, 
				control = list(timeout = timeout)))
			out["2x",] <- if (!is.null(dt)) 
				c(res$objective, Rand.index(label, res$cluster), dt[3])	
		}
		 	
		# K-means matching + BCA
		for (l in 1:2) { 
			name0 <- switch(l, "kmm", "bca")
			match.fun <- switch(l, match.kmeans, match.bca)	
			for (name in names(sigma)) {
				if (name == "r100") {
					objective.best <- Inf
					dt <- system.time({
						for (j in 1:100) {
							res <- match.fun(xx, 
								control=list(sigma=sigma$r100[,,j]))
							if (res$objective < objective.best) {
								objective.best <- res$objective
								sigma.best <- res$sigma
							}
						}
					})
					out[paste0("r100.",name0),] <- c(objective.best, 
						Rand.index(label, perm2cluster(sigma.best)), dt[3])
					if (n <= 100 && gurobi.flag) {
						dt <- system.time(res <- match.2x(xx, sigma.best,
							control = list(timeout = timeout))) 
					if (!is.null(dt))
						out[paste0("r100.",name0,".2x"),] <- c(res$objective, 
						Rand.index(label, perm2cluster(res$sigma)), dt[3]) }	
				} else {
					dt <- system.time(res <- match.fun(xx, 
							control=list(sigma=sigma[[name]])))
					out[paste0(name,".",name0),] <- c(res$objective, 
						Rand.index(label, res$cluster), dt[3])
					if (n < 100 && gurobi.flag) {
						dt <- system.time(res <- match.2x(xx, res$sigma,
							control = list(timeout = timeout)))
						if (!is.null(dt))
							out[paste0(name,".",name0,".2x"),] <- c(res$objective, 
								Rand.index(label,res$cluster), dt[3]) }
				}
			}
		}
	
		## Gaussian mixture with constraints (EM)
		dt <- withTimeout(system.time(res <- match.gaussmix(xx, 
			control=list(maxit=2))),
			timeout = timeout, onTimeout = "silent") 
		if (!is.null(dt)) 
			out["gaussmix",] <- c(objective.fun(xx,res$sigma), 
					Rand.index(label,res$cluster), dt[3])
			
		## Standard ILP
		if (n <= 10) {
		dt <- system.time(res <- blp.solve3(xx, timeout, solver = lp.pkg))
		out["ilp",] <- c(res$objective, 
			Rand.index(label, perm2cluster(res$sigma)), dt[3]) }
		
		## ILP relaxation and IQP of Kuroki and Matsui (2009)
		dt <- withTimeout(system.time(res <- kuroki.matsui.ilp(xx)),
			timeout = timeout, onTimeout = "silent")
		out["kur.ilp",] <- c(res$objective, 
			Rand.index(label, perm2cluster(res$sigma)), dt[3])		
		if (n <= 7 && gurobi.flag) {
			dt <- system.time(res <- kuroki.matsui.iqp(xx, timeLimit=timeout))
			out["kur.iqp",] <- c(res$objective,
				Rand.index(label, perm2cluster(res$sigma)), dt[3]) }
		
		
		# Constrained cluster methods
		if (n <= 100) {
			dt <- withTimeout(
			system.time(res <- conclust::ckmeans(x,m,mustlink,cantlink)),
				timeout = timeout, onTimeout = "silent")
		if (!is.null(dt))
			out["copkm",] <- c(objective.fun(xx, cluster2perm(res,n), unit), 
			Rand.index(label,res), dt[3]) }
		if (n <= 200) {
		dt <- withTimeout(
			system.time(res <- conclust::ccls(x,m,mustlink,cantlink)),
				timeout = timeout, onTimeout = "silent")
		if (!is.null(dt))
			out["ccls",] <- c(objective.fun(xx, cluster2perm(res,n), unit), 
			Rand.index(label,res), dt[3]) }
		dt <- withTimeout(
			system.time(res <- conclust::lcvqe(x,m,mustlink,cantlink)),
			timeout = timeout, onTimeout = "silent")
		if (!is.null(dt))
			out["lcvqe",] <- c(objective.fun(xx, cluster2perm(res,n), unit), 
			Rand.index(label,res), dt[3])
		if (n <= 100) {
		dt <- withTimeout(
			system.time(res <- conclust::mpckm(x,m,mustlink,cantlink)),
			timeout = timeout, onTimeout = "silent")
		if (!is.null(dt))
			out["mpckm",] <- c(objective.fun(xx, cluster2perm(res,n), unit), 
			Rand.index(label,res), dt[3]) }

		# For serial computations, uncomment the line below 
		# and comment out the 2 following lines	
		# result[[i]] <- out
		out[is.na(out[,3]),3] <- timeout
		return(out)
	} # END i loop 

	## Remove outputs from loop iterations that produced errors - if any	
	test <- sapply(result, function(x) all(is.numeric(x)))
	if (any(!test)) {
		pb <- result[!test] 
		result <- result[test] }
	
	## Reshape the results as 3D array and save them  
	result <- do.call(cbind, result)
	dim(result) <- c(nmethod,3,sum(test))
	dimnames(result) <- list(method=method, 
		stat=c("objective","Rand","time"), NULL)
	save(n,result,stdev, file=paste0("sim_results_n",n,".rda"))
	
} # END n loop 


