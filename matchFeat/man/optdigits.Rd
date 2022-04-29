\name{optdigits}
\alias{optdigits}
\docType{data}
\title{
Handwritten Digits Data
}
\description{
Digitized images of handwritten digits used in optical recognition tasks}

\usage{data("optdigits")}

\format{
  The format is:\cr
List of 2\cr
 $ x    : int [1:1000, 1:64] 0 0 0 0 0 0 0 0 0 0 ...\cr
 $ label: int [1:1000] 0 1 2 3 4 5 6 7 8 9 ...
}

\details{
This is a subset of a larger dataset containing handwritten digits contributed by 30 people  on a preprinted form. The forms were converted to normalized bitmaps of size 32x32 which were divided into nonoverlapping blocks of size 4x4. The number of pixels was counted in each block, producing a matrix of size 8x8 with integer coefficients ranging in 0..16.  
These matrix are vectorized in the rows of \code{optdigits$x}. The corresponding digits are in \code{optdigits$label}. 100 examples are available for each digit 0..9.   
}

\source{
	UCI Machine Learning Repository. 
\url{https://archive.ics.uci.edu/ml/datasets/Optical+Recognition+of+Handwritten+Digits
}}
\references{
Alpaydin and Kaynak (1998). Cascading Classifiers. 
\url{ftp://ftp.icsi.berkeley.edu/pub/ai/ethem/kyb.ps.Z}
}

\examples{
data(optdigits)
## Quick visualization
oldpar <- par(no.readonly = TRUE)
par(mfrow=c(2,5))
for (i in 1:10) {
	mat <- matrix(optdigits$x[i,],8,8)
	image(mat[,8:1], xaxt="n", yaxt="n")
	title(optdigits$label[i])
}
par(oldpar)
}
\keyword{datasets}
