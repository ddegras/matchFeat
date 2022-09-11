\name{match.rec}
\alias{match.rec}

\title{
Recursive Initialization Method
}

\description{
RECUR1 algorithm of Bandelt et al (2004) to find starting point in the multidimensional assignment problem with decomposable costs (MDADC)  
}

\usage{
match.rec(x, unit = NULL, w = NULL, control = list())
}

\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for loss function: single positive number, 
	\eqn{p}-vector of length, or \eqn{(p,p)} positive definite matrix}
  
\item{control}{tuning parameters}
}



\value{
A list of class \code{matchFeat} with components
\describe{
\item{\code{sigma}}{best set of permutations for feature vectors (\eqn{(m,n)} matrix)}
\item{\code{cluster}}{associated clusters (= inverse permutations)}
\item{\code{cost}}{minimum objective value}
\item{\code{mu}}{sample mean for each class/label (\eqn{(p,m)} matrix)}
\item{\code{V}}{sample covariance for each class/label (\eqn{(p,m)} matrix}
\item{\code{call}}{function call}
}
}
\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429} \cr
Bandelt, Maas, and Spieksma (2004), "Local search heuristics for multi-index assignment problems with decomposable costs." \doi{10.1057/palgrave.jors.2601723}
}



\seealso{
\code{\link{match.2x}}, \code{\link{match.bca}},
\code{\link{match.gaussmix}}, \code{\link{match.template}}, 
\code{\link{match.kmeans}}
}

\examples{
data(optdigits)
m <- length(unique(optdigits$label)) # number of classes
n <- nrow(optdigits$x) / m # number of units

## Use function with data in matrix form
fit1 <- match.rec(optdigits$x, unit=n)

## Use function with data in array form
p <- ncol(optdigits$x)
x <- t(optdigits$x)
dim(x) <- c(p,m,n)
fit2 <- match.rec(x)
}

