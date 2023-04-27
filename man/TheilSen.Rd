\name{TheilSen}
\alias{TheilSen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theil-Sen slope and intercept estimator.
}
\description{
Computes the Theil-Sen median slope estimator by Theil (1950) and Sen (1968).
The implemented algorithm was proposed by Dillencourt et. al (1992) and runs in an expected \eqn{O(n log n)} time while requiring \eqn{O(n)} storage.
}
\usage{
TheilSen(x, y, alpha = NULL, verbose = TRUE)
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
Given two input vectors \code{x} and \code{y} of length \eqn{n}, the Theil-Sen estimator is computed as \eqn{med_{ij} (y_i - y_j)/(x_i-x_j)}.  By default, the median in this experssion is the upper median, defined as \eqn{\lfloor (n +2) / 2 \rfloor}. 
By changing \code{alpha}, other order statistics of the slopes can be computed.
}
\value{
A list with elements:
\item{intecept }{The estimate of the intercept.}
\item{slope }{The Theil-Sen estimate of the slope.}
}
\references{
Theil, H. (1950), A rank-invariant method of linear and polynomial regression
analysis (Parts 1-3),  \emph{Ned. Akad. Wetensch. Proc. Ser. A}, \bold{53}, 386-392,
521-525, 1397-1412.

Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau.  \emph{Journal of the American statistical association}, \bold{63}(324), 1379-1389.

Dillencourt, M. B., Mount, D. M., & Netanyahu, N. S. (1992). A randomized algorithm for slope selection. \emph{International Journal of Computational Geometry & Applications}, \bold{2}(01), 1-27.

Raymaekers (2023). "The R Journal: robslopes: Efficient Computation of the (Repeated) Median Slope", The R Journal. \href{https://journal.r-project.org/articles/RJ-2023-012/}{(link to open access pdf)}
}
\author{
Jakob Raymaekers
}

\examples{
# We compare the implemented algorithm against a naive brute-force approach.

bruteForceTS <- function(x, y) {
  
  n <- length(x)
  medind1 <- floor(((n * (n - 1)) / 2 + 2) / 2)
  medind2 <- floor((n + 2) / 2)
  temp <-  t(sapply(1:n, function(z)  apply(cbind(x, y), 1 ,
                                                  function(k) (k[2] - y[z]) /
                                                    (k[1] - x[z]))))
  TSslope <- sort(as.vector(temp[lower.tri(temp)]))[medind1]
  TSintercept <- sort(y - x * TSslope)[medind2]
  return(list(intercept = TSintercept, slope = TSslope))
}


n = 100
set.seed(2)
x = rnorm(n)
y = x + rnorm(n)

t0 <- proc.time()
TS.fast <- TheilSen(x, y, NULL, FALSE)
t1 <- proc.time()
t1 - t0

t0 <- proc.time()
TS.naive <- bruteForceTS(x, y)
t1 <- proc.time()
t1 - t0

TS.fast$slope - TS.naive$slope

}
