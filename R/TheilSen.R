TheilSen <- function(x, y, alpha = NULL, verbose = TRUE) {
  
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
      warning(cat("The data contains observations with NA values. The observations ",
                  naInds, " were removed before computing the Theil-Sen estimator."))
    }
    x <- x[-naInds]
    y <- y[-naInds]
    n <- length(x)
  }
  
  # Count the number of duplicate values in the predictor variable
  dupnb  <- rle(x)$lengths
  nbdups <- sum(choose(dupnb[dupnb > 1], 2))
  
  # Compute correct order statistics
  if (is.null(alpha)) {
    medind  <- floor(((n * (n - 1) / 2) - nbdups + 1) / 2) # upper median
    medind2 <- floor((n + 1) / 2) # upper median (for intercept calculation)
  } else {
    medind  <- max(1, min(n, round((n * (n - 1) / 2  - nbdups) * alpha)))
    medind2 <- max(1, min(n, round(n * alpha)))
  }
  
  # Compute the Theil-Sen estimator
  TS.out  <- rcpp_TheilSen(x, y, verbose, medind, medind2)
  
  return(list(intercept = TS.out[1], slope = TS.out[2]))
}