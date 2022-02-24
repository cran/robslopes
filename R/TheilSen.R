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
  x.order     <- order(x)
  xs          <- x[x.order]
  dupnb       <- rle(xs)$lengths
  nbdups      <- sum(choose(dupnb[dupnb > 1], 2))
  
  # Compute correct order statistics
  medind0 <- floor((n + 2) / 2) # upper median (for intercept calculation)
 
   if (is.null(alpha)) {
    medind1  <- floor(((n * (n - 1) / 2) - nbdups + 2) / 2) # upper median
  } else {
    medind1  <- max(1, min((n * (n - 1) / 2  - nbdups),
                          round((n * (n - 1) / 2  - nbdups) * alpha)))
  }
  
  # Compute the Theil-Sen estimator
  TS.out  <- rcpp_TheilSen(x, y, verbose, medind0, medind1)
  
  return(list(intercept = TS.out[1], slope = TS.out[2]))
}
