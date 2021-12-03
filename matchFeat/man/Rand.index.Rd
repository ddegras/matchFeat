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
This is the Rand index not adjusted for chance}

\references{
W. M. Rand (1971). "Objective criteria for the evaluation of clustering methods"\cr
\url{https://en.wikipedia.org/wiki/Rand_index}
}



\examples{
x <- sample.int(3, 20, replace = TRUE)
y <- sample.int(3, 20, replace = TRUE)
table(x,y)
Rand.index(x,y)
}
