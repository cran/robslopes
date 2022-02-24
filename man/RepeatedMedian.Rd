\name{RepeatedMedian}
\alias{RepeatedMedian}
\title{
Siegel's repeated median slope and intercept estimator.
}
\description{
Computes the repeated median slope proposed by Siegel (1982) using the algorithm by 
Matousek et. al (1998). The algorithm runs in an expected \eqn{O(n (log n)^2)} time,
which is typically significantly faster than the \eqn{O(n^2)} computational cost of 
the naive algorithm, and requires \eqn{O(n)} storage.
}
\usage{
RepeatedMedian(x, y, alpha = NULL, beta = NULL, verbose = TRUE)
}

\arguments{
  \item{x}{
A vector of predictor values.
}
\item{y}{
A vector of response values.
}
\item{alpha}{
Determines the outer order statistic, which is equal to \eqn{[alpha*n]}, where \eqn{n} denotes the sample size. Defaults to \code{NULL}, which corresponds
with the (upper) median. 
}
\item{beta}{
Determines the inner order statistic, which is equal to \eqn{[beta*(n-1)]}, where \eqn{n} denotes the sample size. Defaults to \code{NULL}, which corresponds
with the (upper) median. 
}
\item{verbose}{
Whether or not to print out the progress of the algorithm. Defaults to \code{TRUE}.
}
}
\details{
Given two input vectors \code{x} and \code{y} of length \eqn{n}, the repeated median is computed as \eqn{med_i med_j (y_i - y_j)/(x_i-x_j)}.  The default "outer'' median is the \eqn{\lfloor (n + 2) / 2 \rfloor} largest element in the ordered median slopes. The inner median, which for each observation is calculated as the median of the slopes connected to this observation, is the \eqn{\lfloor (n +1) / 2 \rfloor} largest element in the ordered slopes. By changing \code{alpha} and  \code{beta}, other repeated order statistics of the slopes can be computed.
}
\value{
A list with elements:
\item{intecept}{The estimate of the intercept.}
\item{slope}{The Theil-Sen estimate of the slope.}
}

\references{
Siegel, A. F. (1982). Robust regression using repeated medians. \emph{Biometrika}, \bold{69}(1), 242-244.  

Matousek, J., Mount, D. M., & Netanyahu, N. S. (1998). Efficient randomized algorithms for the repeated median line estimator. \emph{Algorithmica}, \bold{20}(2), 136-150.
}
\author{
Jakob Raymaekers
}

\seealso{
\code{\link{TheilSen}}
}
\examples{

# We compare the implemented algorithm against a naive brute-force approach.

bruteForceRM <- function(x, y) {
  
  n <- length(x)
  medind1 <- floor((n+2) / 2)
  medind2 <- floor((n+1) / 2)
  temp <-  t(sapply(1:n, function(z)  sort(apply(cbind(x, y), 1 ,
                                                  function(k) (k[2] - y[z]) /
                                                    (k[1] - x[z])))))
  RMslope <- sort(temp[, medind2])[medind1]
  RMintercept <- sort(y - x * RMslope)[medind1]
  return(list(intercept = RMintercept, slope = RMslope))
}

n = 100
set.seed(2)
x = rnorm(n)
y = x + rnorm(n)

t0 <- proc.time()
RM.fast <- RepeatedMedian(x, y, NULL, NULL, FALSE)
t1 <- proc.time()
t1 - t0

t0 <- proc.time()
RM.naive <- bruteForceRM(x, y)
t1 <- proc.time()
t1 - t0

RM.fast$slope - RM.naive$slope
}
