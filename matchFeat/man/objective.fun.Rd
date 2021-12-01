\name{objective.fun}
\alias{objective.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Cost of Multidimensional Assignment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
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
%%  ~~ If necessary, more details than the description above ~~
}

\value{
The assignment cost
}

\references{
Degras (2020). One-to-one feature matching.
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
}

