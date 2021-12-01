###############################################
# ADDITIONAL EXPERIMENTS FOR "SCALABLE FEATURE 
# MATCHING ACROSS LARGE DATA COLLECTIONS" 
# (DEGRAS, 2021)
###############################################


# This file runs additional experiments on relaxation methods 
# to obtain lower bounds in problem (1) of Degras (2021). 
# Two methods are considered:  
# (i) Convex relaxation of problem (7) of Degras (2021)
# (ii) Lagrangian relaxation of Bandelt et al (2004)



########
# Setup
########

## Clear workspace
rm(list=ls())

## Load packages 
# install.packages("devtools")
# devtools::install_github("ddegras/matchFeat/matchFeat")
library(matchFeat)
library(clue)

## Source files
source("bandelt.R")
source("lbfw.R")
source("utils.R")



##############################
# Load handwritten digit data 
##############################

data(optdigits)
m <- 10
n <- 100
p <- 64
x <- optdigits$x
unit <- optdigits$unit



#############################
# Fit BCA matching algorithm
#############################

# On this small example, it has been checked that the solution found 
# is a global minimizer 
bca <- match.bca(x, unit)
cat("BCA algorithm\nAttained (=optimum) objective value:",bca$objective,"\n")



####################
# Convex relaxation 
####################

## Reshape data matrix to 3D array
xx <- t(x)
dim(xx) <- c(p,m,n)

## Run Frank-Wolfe algorithm 
fw <- lb.fw(xx, tol = 1e-7)
cat("\nConvex relaxation (Frank-Wolfe algorithm)\nAttained objective value:",
	fw$lb,"\n")

P <- fw$P # solution in 3D array form 
P0 <- array(1/m, dim=c(m,m,n)) # theoretical global optimizer
cat("Difference between empirical and theoretical solution:\n",
all.equal(P,P0),"\n\n") # the two should be nearly identical 
cat("Global minimum for convex relaxation:",objective.arr(xx,P0),"\n\n") 



########################
# Lagrangian relaxation
########################

# Calculate lower bounds by Lagrangian relaxation. This procedure is SLOW. 
# Decrease 'maxit' for faster (or use a subset of hubs) but less accurate results

lb <- numeric(n)
for (i in 1:n) {
	lb[i] <- bandelt.lb(x, unit, hub = i, maxit = 25L) 
	if (i %% 10 == 0) cat("Current best lower bound:", max(lb), "\n")
}
cat("Best lower bound:", max(lb), "\n")

# Compare the best lower bound for problem (1) to its global minimum
# The lower bound should be far from the global minimum
# According to Bandelt et al (2004), this situation is caused 
# by the fact that the dissimilarity used (*squared* Euclidean distance)
# in problem (1) is not a Euclidean distance 

