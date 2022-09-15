robslope <- function(formula, data, subset, weights, na.action,
                     type = c("TheilSen", "RepeatedMedian", "PassingBablok"),
                     alpha = NULL, beta = NULL,
                     verbose = TRUE) {
  cl        <- match.call()
  type      <- match.arg(type)
  
  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "subset", 
                "weights", "na.action", "offset"), 
              names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if (attr(attr(mf, "terms"), "intercept") == 0) {
    stop("robust slope without intercept is not supported")
  }
  y <- model.response(mf, "numeric")
  if (is.matrix(y)) 
    stop("response should be a vector not a matrix")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w)) 
    stop("'weights' must be a numeric vector")
  
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = numeric(),
              residuals = y,
              fitted.values = 0 *  y,
              weights = weights)
  }
  else {
    x <- model.matrix(mt, mf)
    z <- robslope.fit(x, y, weights, type,
                      alpha, beta, verbose)
  }
  class(z) <- c("lm")
  z$na.action <- attr(mf, "na.omit")
  z$offset <- NULL
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z
}

robslope.fit <- function(x, y, weights, type,
                         alpha = NULL, beta = NULL,
                         verbose = TRUE) {
  x <- as.matrix(x)
  if (is.null(n <- nrow(x))) 
    stop("'x' must be a matrix")
  if (n == 0L) 
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L) {
    return(list(coefficients = numeric(), residuals = y, 
                fitted.values = 0 * y))
  }
  ny <- NCOL(y)
  if (ny > 1)
    stop("y should be univariate")
  if (is.matrix(y) && ny == 1) 
    y <- drop(y)
  if (NROW(y) != n) 
    stop("incompatible dimensions")
  
  if (type == "TheilSen") {
    robslope.out <- TheilSen(x[, ncol(x)], y, alpha, verbose)
  } else if (type == "RepeatedMedian") {
    robslope.out <- RepeatedMedian(x[, ncol(x)], y, alpha, beta, verbose)
  } else if (type == "PassingBablok") {
    robslope.out <- PassingBablok(x[, ncol(x)], y, alpha, verbose)
  } else {
    stop("invalid 'type' argument")
  }
  
  z <- list()
  
  coef  <- c(robslope.out$intercept, robslope.out$slope)
  dn <- colnames(x)[ncol(x)]
  if (is.null(dn)) {dn <- "x"}
  names(coef) <- c("(intercept)", dn)
  
  
  z$coefficients <- coef
  z$residuals <- y - z$coefficients[1] - z$coefficients[2] * x[, ncol(x)]
  z$fitted.values <- y - z$residuals
 
  z
}
