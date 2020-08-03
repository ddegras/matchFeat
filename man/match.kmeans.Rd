\name{match.kmeans}
\alias{match.kmeans}
\title{
K-Means-Like Method
}
\description{
This function matches collections of feature vectors in a one-to-one fashion 
using a K-means-like method.
}
\usage{
match.kmeans(x, unit = NULL, w = NULL, equal.variance = FALSE, 
	method = c("hungarian", "bruteforce"), control = list())
}

\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{unit}{integer (= number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for the (squared Euclidean) loss function. Can be specified as single positive number, \eqn{p}-vector, or \eqn{p \times p} positive definite matrix}

\item{equal.variance}{logical; if TRUE, resp. FALSE, calculate common, resp. label-specific, covariance after matching features}

\item{method}{method for linear assignment problem: \code{hungarian} algorithm or \code{bruteforce}}

\item{control}{optional tuning parameters}

}

\details{
Given a set of \eqn{n} statistical units, each having \eqn{m} unlabeled feature vectors, the one-to-one matching problem is to find a set of \eqn{n} labels that produce the best match of feature vectors across units. The objective function to minimize is the sum of (weighted) squared Euclidean distances between all pairs of feature vectors having the same label. This amounts to minimizing the sum of the within-label variances.  
The sample means and sample covariances of the matched feature vectors are calculated as a post-processing step.  

If \code{x} is a matrix, the rows should be sorted by increasing unit label and  \code{unit} should be a nondecreasing sequence of integers, for example \eqn{(1,...,1,2,...,2,...,n,...,n)} with each integer \eqn{1,...,n} replicated \eqn{m} times. 

The argument \code{w} can be specified as a vector of positive numbers (will be recycled to length \eqn{p} if needed) or as a positive definite matrix of size \eqn{(p,p)}.

To calculate a common sample covariance for all matched feature vectors, set \code{equal.variance} to TRUE; otherwise (FALSE) one sample covariance matrix is calculated for each label/class. 

The optional argument \code{control} is a list with two fields: \code{sigma}, starting point for the optimization (\eqn{(m,n)} matrix of permutations; and \code{maxit}, maximum number of iterations. 

}


\value{
A  list with components
\item{sigma}{best set of permutations for feature vectors (\eqn{(m,n)} matrix)}
\item{cost}{minimum objective value}
\item{mu}{sample mean for each class/label (\eqn{(p,m)} matrix)}
\item{V}{sample covariance for each class/label (\eqn{(p,m)} matrix}
}

\references{
Degras and Chen (2020). One-to-one feature matching under association ambiguity.\cr
\url{https://en.wikipedia.org/wiki/K-means_clustering}\cr
\url{https://en.wikipedia.org/wiki/Assignment_problem}\cr
\url{https://en.wikipedia.org/wiki/Hungarian_algorithm}
%% Hartigan and Wong (1979). Algorithm AS 136: A K-means clustering algorithm.
}

\seealso{
\code{\link{match.bca}}, \code{\link{match.gaussmix}}, \code{\link{match.template}}
}

\examples{
## Generate data
m <- 3
n <- 10
p <- 5
mu <- matrix(rnorm(p*m),p,m)
sigma <- 0.1
x <- array(mu + rnorm(p*m*n,sigma), c(p,m,n))

## Match all feature vectors
result <- match.kmeans(x)

## Display results 
result$cost # cost function
xmatched <- array(dim=dim(x)) 
  
## re-arranged (matched) feature vectors
for (i in 1:n){
  xmatched[,,i] <- x[,result$sigma[,i],i]}
}
