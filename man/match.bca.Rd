\name{match.bca}
\alias{match.bca}

\title{
Block Coordinate Ascent Method
}

\description{
This function solves the one-to-one feature matching problem by block coordinate ascent.
}

\usage{
function(x, unit = NULL, w = NULL, equal.variance = FALSE,
	method = c("cyclical","random"), control = list())
}


\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for loss function: single positive number, 
	\eqn{p}-vector of length, or \eqn{(p,p)} positive definite matrix}

\item{method}{sweeping method for block coordinate ascent: \code{cyclical} or \code{random} (simple random sampling without replacement)}

\item{equal.variance}{logical; if TRUE, resp. FALSE, calculate common, resp. label-specific, covariance of matched features}
  
\item{control}{optional tuning parameters}
}

\details{
Given a set of \eqn{n} statistical units, each having \eqn{m} possibly mislabeled feature vectors, the one-to-one matching problem is to find a set of \eqn{n} label permutations that produce the best match of feature vectors across units. The objective function to minimize is the sum of (weighted) squared Euclidean distances between all pairs of feature vectors having the same (new) label. This amounts to minimizing the sum of the within-label variances.  
The sample means and sample covariances of the matched feature vectors are calculated as a post-processing step.  
	
The block-coordinate ascent (BCA) algorithm successively sweeps through the statistical units (=blocks), each time relabeling the \eqn{m} feature vectors of a unit to best match those of the other \eqn{n-1} units. 

If \code{x} is a matrix, the rows should be sorted by increasing unit label and  \code{unit} should be a nondecreasing sequence of integers, for example \eqn{(1,...,1,2,...,2,...,n,...,n)} with each integer \eqn{1,...,n} replicated \eqn{m} times. 

The argument \code{w} can be specified as a vector of positive numbers (will be recycled to length \eqn{p} if needed) or as a positive definite matrix of size \eqn{(p,p)}.

To calculate a common sample covariance for all matched feature vectors, set \code{equal.variance} to TRUE; otherwise (FALSE) one sample covariance matrix is calculated for each label/class. 

The optional argument \code{control} is a list with two fields: \code{sigma}, starting point for the optimization (\eqn{(m,n)} matrix of permutations; and \code{maxit}, maximum number of iterations/sweeps. 
}

\value{
A  list with components
\item{sigma}{best set of permutations for feature vectors (\eqn{(m,n)} matrix)}
\item{cost}{minimum objective value}
\item{mu}{sample mean for each class/label (\eqn{(p,m)} matrix)}
\item{V}{sample covariance for each class/label (\eqn{(p,m)} matrix}

}
\references{
Degras and Chen (2020). One-to-one feature matching under association ambiguity. \cr
Wright (2015). Coordinate descent algorithms. \url{https://arxiv.org/abs/1502.04759}
}




\seealso{
\code{\link{match.gaussmix}}, \code{\link{match.kmeans}}, \code{\link{match.template}}
}

\examples{

}

