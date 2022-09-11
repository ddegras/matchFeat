\name{matchFeat-package}
\alias{matchFeat-package}

\docType{package}
\title{
\packageTitle{matchFeat}
}

\description{
\packageDescription{matchFeat}
}

\details{
This package serves to match feature vectors across a collection of datasets 
in a one-to-one fashion. This task is formulated as a multidimensional assignment problem with decomposable costs (MDADC). 
We propose fast algorithms with time complexity roughly linear in the number \eqn{n} of datasets and space complexity a small fraction of the data size. 

\itemize{
\item Initialization methods: \code{\link{match.rec}} (recursive) and \code{\link{match.template}} (template-based). \cr
\item Main matching algorithms: \code{\link{match.bca}}, \code{\link{match.bca.gen}} (for unbalanced data), and \code{\link{match.kmeans}} (\eqn{k}-means matching). \cr
\item Refinement methods (post-processing): \code{\link{match.2x}} (pairwise interchange) and \code{\link{match.gaussmix}} (Gaussian mixture model with permutation constraints).
}
}

\author{
Author: \packageAuthor{matchFeat}\cr
Maintainer: \packageMaintainer{matchFeat}
}

\references{
Degras (2022) "Scalable feature matching across large data collections."  
\doi{10.1080/10618600.2022.2074429}
 \cr
Wright (2015). Coordinate descent algorithms.
\url{https://arxiv.org/abs/1502.04759} \cr
McLachlan and Krishnan (2008). \emph{The EM Algorithm and Extensions}
}
