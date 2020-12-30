\name{match.2x}
\alias{match.2x}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pairwise Interchange Heuristic (2-Assignment-Exchange)
}

\description{
This function implements the Pairwise Interchange Heuristic for the multidimensional assignment problem with decomposable costs (MDADC).
}

\usage{
match.2x(x, sigma = NULL, unit = NULL, w = NULL, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{sigma}{
permutations: matrix of dimensions \eqn{(m,n)}
}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for loss function: single positive number, 
	\eqn{p}-vector of length, or \eqn{(p,p)} positive definite matrix}

\item{method}{sweeping method for block coordinate ascent: \code{cyclical} or \code{random} (simple random sampling without replacement)}
  
\item{control}{tuning parameters}
}

\details{
Use of this function requires to have the GUROBI software and its R package interface installed. Both can be downloaded from \url{https://www.gurobi.com} after obtaining a free academic license.  
}

\value{
A list of class \code{matchFeats} with components
\item{sigma}{best assignment as set of permutations (\eqn{(m,n)} matrix)}
\item{cluster}{best assignment as a cluster membership vector}
\item{objective}{minimum objective value}
\item{mu}{mean vector for each class/label (\eqn{(p,m)} matrix)}
\item{V}{covariance matrix for each class/label (\eqn{(p,p,m)} array)}
\item{call}{function call}
}

\references{
Degras (2020). One-to-one feature matching with applications to multilevel analysis.}

\seealso{
\code{\link{match.1pass}}, \code{\link{match.bca}}, 
\code{\link{match.gaussmix}}, \code{\link{match.kmeans}}, 
\code{\link{match.template}}
}

\examples{
## Generate small example
  m <- 3  # number of classes
  n <- 10 # number of statistical units 
  p <- 5  # number of variables
  mu <- matrix(rnorm(p*m),p,m) # mean vectors
  sigma <- 0.1
  x <- array(as.vector(mu) + rnorm(p*m*n,sigma), c(p,m,n))

## Match all feature vectors
  result <- match.2x(x)

## Display results 
  result$cost  # objective value = assignment cost
  result$sigma # solution permutations
  xmatched <- array(dim=dim(x)) 
  
## Matched feature vectors
  for (i in 1:n)
	  xmatched[,,i] <- x[,result$sigma[,i],i]
	
}
