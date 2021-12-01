\name{match.bca.gen}
\alias{match.bca.gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Block Coordinate Ascent Method for General (Balanced or Unbalanced) Data
}
\description{
Solve a feature matching problem by block coordinate ascent.
}
\usage{
match.bca.gen(x, unit = NULL, nclass = NULL, w = NULL, method = c("cyclical", "random"), control = list())
}

\arguments{
\item{x}{data matrix (rows=instances, columns=features)}

\item{unit}{vector mapping rows of \code{x} to observation units}

\item{w}{feature weights in loss function. Can be specified as single positive number, vector, or positive definite matrix}

\item{method}{sweeping method for block coordinate ascent: \code{cyclical} or \code{random} (simple random sampling without replacement)}
  
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
Degras (2020). One-to-one feature matching with applications to multilevel analysis. \cr
Wright (2015). Coordinate descent algorithms. \url{https://arxiv.org/abs/1502.04759}
}




\seealso{
\code{\link{match.1pass}}, \code{\link{match.2x}},
\code{\link{match.bca}}, \code{\link{match.gaussmix}}, 
\code{\link{match.kmeans}}, \code{\link{match.template}}
}

\examples{
}

