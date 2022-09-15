\name{robslope}
\alias{robslope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust slope estimator
}
\description{
Computes the Theil-Sen median slope, Siegel's repeated median slope or te equivariant Passing-Bablok slope. The algorithms run in an expected linearithmic time while requiring \eqn{O(n)} storage. They are based on Dillencourt et. al (1992), Matousek et. al (1998) and Raymaekers and Dufey (2022).
}
\usage{
robslope(formula, data, subset, weights, na.action,
         type = c("TheilSen", "RepeatedMedian","PassingBablok"),
        alpha = NULL, beta = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    under \sQuote{Details}.}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{robslope} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.  Currently not supported.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s. The default \code{\link{na.exclude}} is applied and an informative message is given in case NAs were removed.}
    \item{type}{the type of robust slope estimator. Should be one of \code{"TheilSen"} (default), \code{"RepeatedMedian"} or \code{"PassingBablok"}.}
    \item{alpha}{Determines the order statistic of the target slope. Defaults to the upper median. See below for details.}
    \item{beta}{Determines the inner order statistic. Only used when \code{type = "RepeatedMedian"}. See below for details.}
    \item{verbose}{
Whether or not to print out the progress of the algorithm. Defaults to \code{TRUE}.
}
}

\details{
This function provides a wrapper around \code{\link{robslope.fit}}, which in turn calls the individual functions \code{\link{TheilSen}}, \code{\link{RepeatedMedian}} or \code{\link{PassingBablok}}. The details on changing the parameters \code{alpha} and \code{beta} can be found in the documentation of those respective functions.
}

\value{
 \code{robslope} returns an object of \code{\link{class}} \code{"lm"}.
 
  The generic accessor functions \code{coefficients},
 \code{fitted.values} and \code{residuals} extract
  various useful features of the value returned by \code{lm}.
  
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
}
\author{
Jakob Raymaekers
}

\seealso{
\code{\link{robslope.fit}}
\code{\link{TheilSen}}
\code{\link{RepeatedMedian}}
\code{\link{PassingBablok}}
}

\examples{
set.seed(123)
df <- data.frame(cbind(rnorm(20), rnorm(20)))
colnames(df) <- c("x", "y")

robslope.out <- robslope(y~x, data = df,
type = "RepeatedMedian", verbose = TRUE)

coef(robslope.out)
plot(fitted.values(robslope.out))

robslope.out <- robslope(y~x, data = df,
type = "TheilSen", verbose = TRUE)

plot(residuals(robslope.out))
}
