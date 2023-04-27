\name{robslope.fit}
\alias{robslope.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust slope estimator
}
\description{
This is the underlying computing engine called by \code{\link{robslope}} used
  to fit robust slopes. It wraps around the individual functions  \code{\link{TheilSen}}, \code{\link{RepeatedMedian}} or \code{\link{PassingBablok}}. These should usually \emph{not} be used
  directly unless by experienced users.
}
\usage{
robslope.fit(x, y, weights, type, alpha = NULL, beta = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{design matrix of dimension \code{n * p}.}
  \item{y}{vector of observations of length \code{n}, or a matrix with
    \code{n} rows.}
    
    \item{type}{the type of robust slope estimator. Should be one of \code{"TheilSen"} (default), \code{"RepeatedMedian"} or \code{"PassingBablok"}.}
    \item{weights}{vector of weights. Currently not in use.}
    \item{alpha}{Determines the order statistic of the target slope. Defaults to the upper median. See below for details.}
    \item{beta}{Determines the inner order statistic. Only used when \code{type = "RepeatedMedian"}. See below for details.}
    \item{verbose}{
Whether or not to print out the progress of the algorithm. Defaults to \code{TRUE}.
}
}

\details{
This function provides a wrapper around the individual functions \code{\link{TheilSen}}, \code{\link{RepeatedMedian}} or \code{\link{PassingBablok}}. The details on changing the parameters \code{alpha} and \code{beta} can be found in the documentation of those respective functions.
}

\value{
  \code{\link{list}} with components
  \item{coefficients}{\code{p} vector}
  \item{residuals}{\code{n} vector or matrix}
  \item{fitted.values}{\code{n} vector or matrix}
  
}
\references{
Theil, H. (1950), A rank-invariant method of linear and polynomial regression
analysis (Parts 1-3),  \emph{Ned. Akad. Wetensch. Proc. Ser. A}, \bold{53}, 386-392,
521-525, 1397-1412.

Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau.  \emph{Journal of the American statistical association}, \bold{63}(324), 1379-1389.

Dillencourt, M. B., Mount, D. M., & Netanyahu, N. S. (1992). A randomized algorithm for slope selection. \emph{International Journal of Computational Geometry & Applications}, \bold{2}(01), 1-27.

Siegel, A. F. (1982). Robust regression using repeated medians. \emph{Biometrika}, \bold{69}(1), 242-244.  

Matousek, J., Mount, D. M., & Netanyahu, N. S. (1998). Efficient randomized algorithms for the repeated median line estimator. \emph{Algorithmica}, \bold{20}(2), 136-150.


Passing, H., Bablok, W. (1983). A new biometrical procedure for testing the equality of measurements from two different analytical methods. Application of linear regression procedures for method comparison studies in clinical chemistry, Part I,  \emph{Journal of clinical chemistry and clinical biochemistry}, \bold{21},709-720.

Bablok, W., Passing, H., Bender, R., Schneider, B. (1988). A general regression procedure for method transformation. Application of linear
regression procedures for method comparison studies in clinical chemistry, Part III.  \emph{Journal of clinical chemistry and clinical biochemistry}, \bold{26},783-790.

Raymaekers J., Dufey F. (2022). Equivariant Passing-Bablok regression in quasilinear time. \href{https://arxiv.org/abs/2202.08060}{(link to open access pdf)}

Raymaekers (2023). "The R Journal: robslopes: Efficient Computation of the (Repeated) Median Slope", The R Journal. \href{https://journal.r-project.org/articles/RJ-2023-012/}{(link to open access pdf)}
}
\author{
Jakob Raymaekers
}

\seealso{
\code{\link{robslope}}
\code{\link{TheilSen}}
\code{\link{RepeatedMedian}}
\code{\link{PassingBablok}}
}

\examples{
set.seed(123)
x <- rnorm(20)
y <- rnorm(20)

robslope.out <- robslope.fit(x, y, type = "RepeatedMedian", verbose = TRUE)

coef(robslope.out)
plot(fitted.values(robslope.out))

robslope.out <- robslope.fit(x, y, type = "TheilSen", verbose = TRUE)

plot(residuals(robslope.out))
}
