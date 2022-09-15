PassingBablok <- function(x, y, alpha = NULL, verbose = TRUE) {
  
  # Check inputs
  n <- length(x)
  if (length(y) != n) {
    stop("x and y should have the same length")
  }
  naX <- which(is.na(x))
  naY <- which(is.na(y))
  if (length(naX) > 0 | length(naY) > 0) {
    naInds <- sort(unique(c(naX,naY)))
    if (verbose) {
      cat("The data contains observations with NA values. The observations ",
                  naInds, " were removed before computing the Passing-Bablok estimator.")
    }
    x <- x[-naInds]
    y <- y[-naInds]
    n <- length(x)
  }
  
  # Compute correct order statistics
  medind0 <- floor((n + 2) / 2) # upper median (for intercept calculation)
  if (is.null(alpha)) {
    medind1  <- floor(((n * (n - 1) / 2) + 2) / 2) # upper median
  } else {
    medind1  <- max(1, min((n * (n - 1) / 2),
                           round((n * (n - 1) / 2) * alpha)))
  }
  
  PB.out <- rcpp_PassingBablok(x, y,verbose, medind0, medind1)
  
  return(list(intercept = PB.out[1], slope = PB.out[2]))
}