\name{objective.fun}
\alias{objective.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Cost of Multidimensional Assignment
}
\description{
Calculates the objective value in the multidimensional assignment problem with decomposable costs (MDADC). The dissimilarity function used in this problem is the squared Euclidean distance.}
\usage{
objective.fun(x, sigma = NULL, unit = NULL, w = NULL)
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

}

\details{
	Given \eqn{n} datasets having each \eqn{m} vectors of same size, 
	say \eqn{{x_{11},...,x_{1m}},...,x_{n1},...,x_{nm}}, and permutations 
	\eqn{\sigma_1,...,\sigma_n} of \eqn{{1,...,m}}, the function calculates 
 \eqn{1/(n(n-1)) sum_{i,j} sum_{k} || x_{i,sigma_i(k)- x_{j,\sigma_j(k) \|^2}}} where \eqn{i} and \eqn{n} run from 1 to \eqn{n} and \eqn{k} runs from 1 to \eqn{m}. This is the objective value (1) of Degras (2021), up to the factor \eqn{1/(n(n-1))}.
}

\value{Objective value}

\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429}
 }

\seealso{\code{\link{objective.gen.fun}}}


\examples{
data(optdigits)
m <- 10
n <- 100
sigma <- matrix(1:m,m,n) # identity permutations
objective.fun(optdigits$x, sigma, optdigits$unit)
}

