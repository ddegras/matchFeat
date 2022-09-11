\name{match.bca.gen}
\alias{match.bca.gen}

\title{
Block Coordinate Ascent Method for General (Balanced or Unbalanced) Data
}

\description{
Solve a feature matching problem by block coordinate ascent
}

\usage{
match.bca.gen(x, unit = NULL, cluster = NULL, w = NULL, 
	method = c("cyclical", "random"), control = list())
}

\arguments{
\item{x}{data matrix (rows=instances, columns=features)}

\item{unit}{vector of unit labels (length = number of rows of \code{x})}

\item{cluster}{integer specifying the number of classes/clusters to assign the feature vectors to OR integer vector specifiying the initial cluster assignment. }

\item{w}{feature weights in loss function. Can be specified as single positive number, vector, or positive definite matrix}

\item{method}{sweeping method for block coordinate ascent: \code{cyclical} or \code{random} (simple random sampling without replacement)}
  
\item{control}{optional list of tuning parameters}

}

\details{
If \code{cluster} is an integer vector, it must have the same length as \code{unit} 
and its values must range between 1 and the number of clusters. 

The list \code{control} can contain a field \code{maxit}, 
an integer that fixes the maximum number of algorithm iterations.
}

\value{
A list of class \code{matchFeat} with components
\describe{
\item{\code{cluster}}{integer vector of cluster assignments (length = \code{now(x)})}
\item{\code{objective}}{minimum objective value}
\item{\code{mu}}{sample mean for each cluster/class (feature-by-cluster matrix)}
\item{\code{V}}{sample covariance for each cluster/class (feature-by-feature-by-cluster 3D array)}
\item{\code{size}}{integer vector of cluster sizes}	
\item{\code{call}}{function call}
}
}

\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429} \cr
Wright (2015). Coordinate descent algorithms. 
\url{https://arxiv.org/abs/1502.04759}
}


\seealso{
\code{\link{match.2x}}, \code{\link{match.bca}}, 
\code{\link{match.bca.gen}}, \code{\link{match.gaussmix}}, 
\code{\link{match.kmeans}}, \code{\link{match.rec}}, \code{\link{match.template}}
}


\examples{
data(optdigits)
nobs <- nrow(optdigits$x) # total number of observations
n <- length(unique(optdigits$unit)) # number of statistical units
rmv <- sample.int(nobs, n-1) # remove (n-1) observations to make data unbalanced
min.m <- max(table(optdigits$unit[-rmv])) # smallest possible number of clusters
# lower values will result in an error message 
m <- min.m
result <- match.bca.gen(optdigits$x[-rmv,], optdigits$unit[-rmv], m)



}

