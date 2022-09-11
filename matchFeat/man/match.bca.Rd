\name{match.bca}
\alias{match.bca}

\title{
Block Coordinate Ascent Method
}

\description{
This function solves the multidimensional assignment problem with decomposable costs (MDADC) by block coordinate ascent. The dissimilarity function is the squared Euclidean distance. 
}

\usage{
match.bca(x, unit = NULL, w = NULL, 
	method = c("cyclical", "random"), control = list())
}


\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}
\item{w}{weights for loss function: single positive number, 
	\eqn{p}-vector of length, or \eqn{(p,p)} positive definite matrix}

\item{method}{sweeping method for block coordinate ascent: \code{cyclical} or \code{random} (simple random sampling without replacement)}
  
\item{control}{tuning parameters}
}

\details{
Given a set of \eqn{n} statistical units, each having \eqn{m} possibly mislabeled feature vectors, the one-to-one matching problem is to find a set of \eqn{n} label permutations that produce the best match of feature vectors across units. The objective function to minimize is the sum of (weighted) squared Euclidean distances between all pairs of feature vectors having the same (new) label. This amounts to minimizing the sum of the within-label variances.  
The sample means and sample covariances of the matched feature vectors are calculated as a post-processing step.  
	
The block-coordinate ascent (BCA) algorithm successively sweeps through the statistical units (=blocks), each time relabeling the \eqn{m} feature vectors of a unit to best match those of the other \eqn{n-1} units. 

If \code{x} is a matrix, the rows should be sorted by increasing unit label and  \code{unit} should be a nondecreasing sequence of integers, for example \eqn{(1,...,1,2,...,2,...,n,...,n)} with each integer \eqn{1,...,n} replicated \eqn{m} times. 

The argument \code{w} can be specified as a vector of positive numbers (will be recycled to length \eqn{p} if needed) or as a positive definite matrix of size \eqn{(p,p)}.

The optional argument \code{control} is a list with three fields: \code{sigma}, starting point for the optimization (\eqn{(m,n)} matrix of permutations; \code{maxit}, maximum number of iterations; and \code{equal.variance}, logical value that specifies whether the returned sample covariance matrices \code{V} for matched features should be equal between labels/classes (TRUE) or label-specific (FALSE, default).  
}

\value{
A list of class \code{matchFeat} with components
\describe{
\item{\code{sigma}}{best set of permutations for feature vectors (\eqn{(m,n)} matrix)}
\item{\code{cluster}}{associated clusters (=inverse permutations)}
\item{\code{objective}}{minimum objective value}
\item{\code{mu}}{sample mean for each class/label (\eqn{(p,m)} matrix)}
\item{\code{V}}{sample covariance for each class/label (\eqn{(p,m)} matrix}
\item{\code{call}}{function call}
}
}

\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429} \cr
Wright (2015). Coordinate descent algorithms. \url{https://arxiv.org/abs/1502.04759}
}

\seealso{
\code{\link{match.2x}},
\code{\link{match.bca.gen}}, \code{\link{match.gaussmix}}, 
\code{\link{match.kmeans}}, \code{\link{match.rec}}, \code{\link{match.template}}
}

\examples{
data(optdigits)
m <- length(unique(optdigits$label)) # number of classes
n <- nrow(optdigits$x) / m # number of units

## Use function with data in matrix form
fit1 <- match.bca(optdigits$x, unit=n)

## Use function with data in array form
p <- ncol(optdigits$x)
x <- t(optdigits$x)
dim(x) <- c(p,m,n)
fit2 <- match.bca(x)

}

