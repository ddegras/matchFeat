\name{objective.gen.fun}
\alias{objective.gen.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Objective Value in One-To-One Feature Matching with Balanced or Unbalanced Data
}

\description{
Calculates the objective value in the multidimensional assignment problem with decomposable costs (MDADC). The dissimilarity function used in this problem is the squared Euclidean distance. The data can be balanced OR unbalanced.}

\usage{
objective.gen.fun(x, unit, cluster)
}

\arguments{
  \item{x}{
data matrix with feature vectors in rows
}
  \item{unit}{
vector of unit labels (length should equal number of rows in \code{x})
}
  \item{cluster}{
vector of cluster labels (length should equal number of rows in \code{x})}
}

\details{
See equation (2) in Degras (2022). This function gives the same value as \code{\link{objective.fun}} when the data are balanced.   
}

\value{Objective value}

\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{objective.fun}}
}

\examples{
data(optdigits)
m <- 10
n <- 100

## Balanced example: both 'objective.fun' and 'objective.gen.fun' work
sigma <- matrix(1:m,m,n)
cluster <- rep(1:m,n)
objective.fun(optdigits$x, sigma, optdigits$unit)
objective.gen.fun(optdigits$x, optdigits$unit, cluster)

## Unbalanced example
idx <- 1:999
objective.gen.fun(optdigits$x[idx,], optdigits$unit[idx], cluster[idx])


}
