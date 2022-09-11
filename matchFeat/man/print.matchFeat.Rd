\name{print.matchFeat}
\alias{print.matchFeat}

\title{
Print a matchFeat Object}
\description{
\code{print} method for class \code{"matchFeat"}.
}
\usage{\method{print}{matchFeat}(x,...)}

\arguments{
  \item{x}{an object of class \code{"matchFeat"}.}
  \item{...}{for compatibility with the generic \code{print} method; argument not currently used.}
}

\details{
The function \code{print.matchFeat} concisely displays the information of an object of class \code{"matchFeat"}. More precisely it shows the  
data range, bandwidth used in local polynomial estimation, and key information on SCB and statistical tests. 
}

\value{No return value, called for side effects}

\seealso{
\code{\link{predict.matchFeat}}, \code{\link{summary.matchFeat}}
}

\examples{
data(optdigits)
result <- match.bca(optdigits$x, optdigits$unit)
print(result)
}

\keyword{print}
\keyword{methods}
