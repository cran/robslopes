RepeatedMedian <- function(x, y, alpha = NULL, beta = NULL, verbose = TRUE) {
  
  # Check inputs
  n <- length(x)
  if (length(y) != n) {
    stop("x and y should have the same length")
  }
  naX <- which(is.na(x))
  naY <- which(is.na(y))
  if (length(naX) > 0 | length(naY) > 0) {
    naInds <- sort(unique(c(naX,naY)))
    warning(paste0("The data contains observations with NA values. The observations ",
                   naInds, " were removed before computing the Theil-Sen estimator."))
    x <- x[-naInds]
    y <- y[-naInds]
    n <- length(x)
  }
  
  # Create a vector with number of duplicate predictor values for each observation
  x.order     <- order(x)
  xs          <- x[x.order]
  nbdups      <- rle(xs)$lengths
  dupnbs      <- which(nbdups > 1)
  nbdups_temp <- rep(1, length(x))
  counter     <- 0
  dup_cumsum  <- cumsum(nbdups)
  if (length(dupnbs) > 0) {
    for (i in 1:length(dupnbs)) {
      dupnbstemp <- dupnbs[i]
      counter    <- max(0, dup_cumsum[dupnbstemp - 1]) + 1
      nbdups_temp[(counter):(counter + nbdups[dupnbs[i]] - 1)] <- nbdups[dupnbs[i]]
    }
    nbdups[x.order] <- nbdups_temp
  }
  
  # Compute correct order statistics
  if (is.null(alpha)) {
    medind1 <- floor((n + 1) / 2) # upper median
  } else {
    medind1 <- pmax(1, pmin(n, round(n * alpha))) # upper median
  }
  if (is.null(beta)) {
    medind2 <- floor((n - nbdups + 1) / 2.0) # vector of upper medians
  } else {
    medind2 <- pmax(1, pmin(n, round((n - nbdups) * beta)))
  }
  
  # Run the RM algorithm
  RM.out  <- rcpp_RepeatedMedian(x, y, verbose, medind1, medind2)
  
  return(list(intercept = RM.out[1], slope = RM.out[2]))
}