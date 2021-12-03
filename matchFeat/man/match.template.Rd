\name{match.template}
\alias{match.template}

\title{Template Matching}

\description{
This function solves the multidimensional assignment problem with decomposable costs (MDADC) by matching the data to a pre-specified set of vectors (the template). The dissimilarity function is the squared Euclidean distance. 
}

\usage{
match.template(x, template = 1L, unit = NULL, w = NULL, 
	method = c("hungarian","bruteforce"), equal.variance = FALSE)
}

\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{template}{integer (= which sample unit to take as template) or \eqn{(p,m)} matrix}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for the loss function. Can be specified as a \eqn{p}-vector of diagonal weights or as a full \eqn{p \times p} (positive definite) matrix}

\item{method}{method for the linear assignment problem: \code{hungarian} algorithm or \code{bruteforce}}

\item{equal.variance}{logical; if TRUE, resp. FALSE, return common, resp. label-specific, covariance of matched features}
}

\details{
Given a sample of \eqn{n} statistical units, each having \eqn{m} possibly mislabeled feature vectors, the one-to-one matching problem is to find a set of \eqn{n} label permutations that produce the best match of feature vectors across units. The objective function to minimize is the sum of squared (Euclidean) distances between all feature vectors having the same (new) label. This amounts to minimizing the sum of the within-label variances.  
	
The template-based method consists in relabeling successively each sample unit to best match a template matrix of feature vectors. This method is very fast but its optimization performance is only as good as the template. A template that is typical or representative of the collected data. 

If \code{x} is a matrix, the rows should be sorted by increasing unit label and  \code{unit} should be a nondecreasing sequence of integers, for example \eqn{(1,...,1,2,...,2,...,n,...,n)} with each integer \eqn{1,...,n} replicated \eqn{m} times. 

The argument \code{w} can be specified as a vector of positive numbers (will be recycled to length \eqn{p} if needed) or as a positive definite matrix of size \eqn{(p,p)}.
}

\value{
A list of class \code{matchFeat} with fields 
\describe{
\item{\code{sigma}}{best assignement as set of permutations (\eqn{(m,n)} matrix)}
\item{\code{cluster}}{best assignement as cluster indicators (\eqn{m*n}-vector)}
\item{\code{objective}}{minimum objective value}
\item{\code{mu}}{mean vector for each class/label (\eqn{(p,m)} matrix)}
\item{\code{V}}{covariance matrix for each class/label (\eqn{(p,p,m)} array if \code{equal.variance} is FALSE, \eqn{(p,p)} matrix otherwise}
\item{\code{call}}{function call}
}
}

\references{
Degras (2021) "Scalable feature matching across large data collections". 
\url{https://arxiv.org/abs/2101.02035}\cr
\url{https://en.wikipedia.org/wiki/Assignment_problem}\cr
\url{https://en.wikipedia.org/wiki/Hungarian_algorithm}
}


\seealso{
 \code{\link{match.2x}}, \code{\link{match.bca}}, 
\code{\link{match.bca.gen}}, \code{\link{match.gaussmix}}, \code{\link{match.kmeans}}, \code{\link{match.rec}}
}

\examples{
## Generate data
n <- 10
k <- 3
d <- 5
mu <- matrix(1:k, nrow=d, ncol=k, byrow=TRUE)
sigma <- 0.3
x <- array(mu, c(d,k,n)) + rnorm(d*k*n,sigma)
## Match all feature vectors with first case as template
result <- match.template(x,1)
## Display results 
result$cost # cost function
xmatched <- array(dim=dim(x)) 
# re-arranged (matched) feature vectors
for (i in 1:n)
	xmatched[,,i] <- x[,result$sigma[,i],i]
}
