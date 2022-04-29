\name{Rand.index}
\alias{Rand.index}

\title{
Rand Index of Agreement Between Two Partitions}

\description{
Calculates the Rand Index between two partitions of a set}

\usage{
Rand.index(x, y)
}

\arguments{
  \item{x}{
first partition vector}
  \item{y}{
second partition vector}
}

\details{
The two vectors \code{x} and \code{y} must have equal length. Given a set \eqn{S} and two partitions \eqn{X} and \eqn{Y} of \eqn{S}, the Rand index is the proportion of pairs of elements in \eqn{S} (out of all pairs) that are either concordant in both \eqn{X} and \eqn{Y} (i.e., they belong to the same member of \eqn{X} and to the same member of \eqn{Y}) or discordant (i.e., not concordant) in both \eqn{X} and {Y}.  
}

\value{The Rand index (not adjusted for chance)}

\references{
W. M. Rand (1971). "Objective criteria for the evaluation of clustering methods"\cr
\url{https://en.wikipedia.org/wiki/Rand_index}
}



\examples{
## Example 1
x <- sample.int(3, 20, replace = TRUE)
y <- sample.int(3, 20, replace = TRUE)
table(x,y)
Rand.index(x,y)

## Example 2
data(optdigits)
label <- optdigits$label 
m <- length(unique(label)) # 10 
n <- length(unique(optdigits$unit)) # 100
dim(label) <- c(m,n)
p <- ncol(optdigits$x) # 64
x <- array(t(optdigits$x),c(p,m,n))
## Permute data and labels to make problem harder
for (i in 1:n) {
	sigma <- sample.int(m)
	x[,,i] <- x[,sigma,i]
	label[,i] <- label[sigma,i]
}
## Compare Rand indices of matching methods
Rand.index(match.bca(x)$cluster, label)
Rand.index(match.rec(x)$cluster, label)
Rand.index(match.template(x)$cluster, label)
Rand.index(match.kmeans(x)$cluster, label)


}
