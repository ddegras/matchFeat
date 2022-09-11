\name{predict.matchFeat}
\alias{predict.matchFeat}

\title{
Match New Feature Vectors To Existing Clusters}

\description{
\code{predict} method for class \code{"matchFeat"}
}

\usage{\method{predict}{matchFeat}(object, newdata, unit = NULL, ...)}

\arguments{
  \item{object}{an object of class \code{"matchFeat"}.}
  \item{newdata}{new dataset of feature vectors}
  \item{unit}{unit labels for new data. Only necessary if \code{newdata} is a matrix}
  \item{...}{for compatibility with the generic \code{predict} method; argument not currently used.}
}

\details{
The function \code{predict.matchFeat} finds the best matching for new feature vectors relative to an existing set of cluster/class centers. If code{object} results from a call to \code{\link{match.gaussmix}}, the same function is used for prediction (with fixed mean vectors and covariance matrices). In other cases, the function \code{\link{match.template}} is used for prediction. 
}

\value{A list of class \code{matchFeat} with fields
\describe{
\item{\code{sigma}}{best matching as set of permutations (\eqn{(m,n)} matrix)}
\item{\code{cluster}}{best matching as cluster indicators (\eqn{(m,n)}-matrix)}
\item{\code{objective}}{minimum objective value}
\item{\code{mu}}{mean vector for each class/label (\eqn{(p,m)} matrix)}
\item{\code{V}}{covariance matrix for each class/label (\eqn{(p,p,m)} array if \code{equal.variance} is FALSE, \eqn{(p,p)} matrix otherwise}
\item{\code{call}}{function call}
}
}

\seealso{
\code{\link{print.matchFeat}}, \code{\link{summary.matchFeat}}
}

\examples{
data(optdigits)
train.result <- match.bca(optdigits$x[1:900,], optdigits$unit[1:900])  
test.result <- predict(train.result, optdigits$x[901:1000,], optdigits$unit[901:1000])
test.result
}

\keyword{predict}
\keyword{methods}
