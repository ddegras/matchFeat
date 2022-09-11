\name{summary.matchFeat}
\alias{summary.matchFeat}

\title{Summarize a matchFeat Object}

\description{\code{summary} method for class \code{"matchFeat"} }

\usage{\method{summary}{matchFeat}(object, ...)}

\arguments{\item{object}{an object of class \code{"matchFeat"}}
\item{...}{additional arguments; not currently used.}
}

\details{
The function \code{summary.matchFeat} displays all fields of a \code{matchFeat} object at the exception of \code{x}, \code{y}, \code{par}, \code{nonpar}, \code{normscb}, and \code{bootscb} which are potentially big. It provides information on the function call, data, local polynomial fit, SCB, and statistical tests.
}

\value{No return value, called for side effects}

\seealso{
\code{\link{predict.matchFeat}},
\code{\link{print.matchFeat}} 
}

\examples{
data(optdigits)
result <- match.bca(optdigits$x, optdigits$unit)
summary(result)
}


\keyword{methods}
