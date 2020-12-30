\name{match.1pass}
\alias{match.1pass}

\title{
One-Pass Method
}

\description{
Bla
}

\usage{
match.1pass(x, unit = NULL, w = NULL, control = list())
}

\arguments{
\item{x}{data: matrix of dimensions \eqn{(mn,p)} or 3D array of dimensions \eqn{(p,m,n)} with \eqn{m} = number of labels/classes, \eqn{n} = number of sample units, and \eqn{p} = number of variables)}

\item{unit}{integer (=number of units) or vector mapping rows of \code{x} to sample units (length \eqn{mn}). Must be specified only if \code{x} is a matrix.}

\item{w}{weights for loss function: single positive number, 
	\eqn{p}-vector of length, or \eqn{(p,p)} positive definite matrix}
  
\item{control}{tuning parameters}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list of class \code{matchFeats} with components
\item{sigma}{best set of permutations for feature vectors (\eqn{(m,n)} matrix)}
\item{cost}{minimum objective value}
\item{mu}{sample mean for each class/label (\eqn{(p,m)} matrix)}
\item{V}{sample covariance for each class/label (\eqn{(p,m)} matrix}
\item{ss.between.unmatched}{sum of squares between labels for unmatched data}
\item{ss.within.unmatched}{sum of squares within labels for unmatched data}
\item{call}{function call}
}

\references{
%% ~put references to the literature/web site here ~
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
fit1 <- match.1pass(optdigits$x, unit=n)

## Use function with data in array form
p <- ncol(optdigits$x)
x <- t(optdigits$x)
dim(x) <- c(p,m,n)
fit2 <- match.1pass(x)
stopifnot(identical(fit1,fit2))
}

