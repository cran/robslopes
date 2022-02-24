\name{PassingBablok}
\alias{PassingBablok}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Passing-Bablok slope and intercept estimator.
}
\description{
Computes the equivariant Passing-Bablok regression. The implemented algorithm was proposed by Raymaekers and Dufey (2022) and runs in an expected \eqn{O(n log n)} time while requiring \eqn{O(n)} storage.
}
\usage{
PassingBablok(x, y, alpha = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of predictor values.
}
\item{y}{
A vector of response values.
}
\item{alpha}{
Determines the order statistic of the target slope, which is equal to \eqn{[alpha*n*(n-1)]}, where \eqn{n} denotes the sample size. Defaults to \code{NULL}, which corresponds
with the (upper) median. 
}
\item{verbose}{
Whether or not to print out the progress of the algorithm. Defaults to \code{TRUE}.
}
}
\details{
Given two input vectors \code{x} and \code{y} of length \eqn{n}, the equivariant Passing-Bablok estimator is computed as \eqn{med_{ij} |(y_i - y_j)/(x_i-x_j)|}.  By default, the median in this experssion is the upper median, defined as \eqn{\lfloor (n +2) / 2 \rfloor}. 
By changing \code{alpha}, other order statistics of the slopes can be computed.
}
\value{
A list with elements:
\item{intecept }{The estimate of the intercept.}
\item{slope }{The Theil-Sen estimate of the slope.}
}
\references{

Passing, H., Bablok, W. (1983). A new biometrical procedure for testing the equality of measurements from two different analytical methods. Application of linear regression procedures for method comparison studies in clinical chemistry, Part I,  \emph{Journal of clinical chemistry and clinical biochemistry}, \bold{21},709-720.

Bablok, W., Passing, H., Bender, R., Schneider, B. (1988). A general regression procedure for method transformation. Application of linear
regression procedures for method comparison studies in clinical chemistry, Part III.  \emph{Journal of clinical chemistry and clinical biochemistry}, \bold{26},783-790.

Raymaekers J., Dufey F. (2022). Equivariant Passing-Bablok regression in quasilinear time. \href{https://arxiv.org/abs/2202.08060}{(link to open access pdf)}

}
\author{
Jakob Raymaekers
}

\examples{
# We compare the implemented algorithm against a naive brute-force approach.

bruteForcePB <- function(x, y) {
  
  n <- length(x)
  medind1 <- floor(((n * (n - 1)) / 2 + 2) / 2) # upper median
  medind2 <- floor((n + 2) / 2)
  temp <-  t(sapply(1:n, function(z)  apply(cbind(x, y), 1 ,
                                            function(k) (k[2] - y[z]) /
                                              (k[1] - x[z]))))
  PBslope <- sort(abs(as.vector(temp[lower.tri(temp)])))[medind1]
  PBintercept <- sort(y - x * PBslope)[medind2]
  return(list(intercept = PBintercept, slope = PBslope))
}


n = 100
set.seed(2)
x = rnorm(n)
y = x + rnorm(n)

t0 <- proc.time()
PB.fast <- PassingBablok(x, y, NULL, FALSE)
t1 <- proc.time()
t1 - t0

t0 <- proc.time()
PB.naive <- bruteForcePB(x, y)
t1 <- proc.time()
t1 - t0

PB.fast$slope - PB.naive$slope

}
