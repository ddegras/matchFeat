\name{match.fw}
\alias{match.fw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Franck-Wolfe Algorithm for Matching Problem
}
\description{
}

\usage{
match.fw(x, unit = NULL, w = NULL, control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
data matrix or array. Dimensions \eqn{(mn,p)} for matrix, \eqn{(p,m,n)} for array with m: number of classes/clusters, n: number of observation units, p: number of variables (size of observation vector)
}
  \item{unit}{
integer vector of unit labels (length \eqn{n}) or integer: number of units. Not needed id \code{x} is an array 
}
  \item{w}{
vector or matrix weights for distance between observations. Default: \code{NULL}
}
  \item{control}{
optional list of tuning parameters for algorithm
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{matchFeat}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Degras (2020). Fast methods for the feature matching problem
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, unit = NULL, w = NULL, control = list()) 
{
    pre <- matchFeat:::preprocess(x = x, unit = unit, w = w)
    m <- pre$m
    n <- pre$n
    p <- pre$p
    if (!is.null(w)) 
        w <- pre$w
    R <- pre$R
    if (!is.null(pre$x)) {
        x <- pre$x
        dim(x) <- c(p, m, n)
    }
    rm(pre)
    syscall <- sys.call()
    sigma <- matrix(1:m, m, n)
    maxit <- 1000L
    equal.variance <- FALSE
    if (is.list(control)) {
        if (!is.null(control$sigma)) 
            sigma <- control$sigma
        if (!is.null(control$maxit)) 
            maxit <- control$maxit
        if (!is.null(control$equal.variance)) 
            equal.variance <- control$equal.variance
    }
    if (m == 1 || n == 1 || p == 1) 
        return(matchFeat:::trivial(x, m, n, p, w, R, equal.variance, 
            syscall))
    if (!is.null(w)) {
        if (is.vector(w)) {
            x <- sqrt(w) * x
        }
        else {
            dim(x) <- c(p, m * n)
            x <- R \%*\% x
            dim(x) <- c(p, m, n)
        }
    }
    xmin <- min(x)
    if (xmin < 0) 
        x <- x - xmin
    P <- sigma
    sumxP <- matrix(0, p, m)
    for (i in 1:n) sumxP <- sumxP + x[, P[, i], i]
    objective <- sum(sumxP^2)
    Q <- matrix(, m, n)
    for (it in 1:maxit) {
        objective.old <- objective
        sumxQ <- matrix(0, p, m)
        for (i in 1:n) {
            grad <- crossprod(sumxP, x[, , i])
            Q[, i] <- solve_LSAP(grad, maximum = TRUE)
            sumxQ <- sumxQ + x[, Q[, i], i]
        }
        objective <- sum(sumxQ^2)
        if (objective > objective.old) {
            P <- Q
            sumxP <- sumxQ
        }
        else break
    }
    cost <- (sum(x^2) - (objective/n))/(n - 1)
    mu <- sumxP/n
    V <- array(, c(p, p, m))
    dim(x) <- c(p, m * n)
    for (l in 1:m) {
        idx <- seq.int(0, by = m, len = n) + sigma[l, ]
        V[, , l] <- tcrossprod(x[, idx])/n - tcrossprod(mu[, 
            l])
    }
    if (xmin < 0) 
        mu <- mu + xmin
    if (equal.variance) 
        V <- rowMeans(V, dims = 2L)
    if (!is.null(w)) {
        if (is.vector(w)) {
            mu <- mu/sqrt(w)
            V <- V/w
        }
        else {
            mu <- backsolve(R, mu)
            dim(V) <- c(p, m * p)
            V <- t(backsolve(R, V))
            dim(V) <- c(p, p, m)
        }
    }
    cluster <- matrix(, m, n)
    for (i in 1:n) cluster[P[, i], i] <- 1:m
    out <- list(sigma = P, cluster = cluster, objective = cost, 
        mu = mu, V = V, call = syscall)
    class(out) <- "matchFeat"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
