#' Soft-thresholding Function
#' 
#' Simple vectorized soft-thresholding utility function.
#' 
#' @details This function is not exported, as it is an internal utility function. Use
#' \code{spca:::soft_thresh} to access it.
#' 
#' @param x Vector to be soft-thresholded
#' @param lambda Threshold value
#' 
#' @examples 
#' set.seed(1)
#' (z <- rnorm(10))
#' spca:::soft_thresh(z, 0.5)
#' 
#' y <- seq(from = -2, to = 2, by = 0.1)
#' plot(y, spca:::soft_thresh(y, 0.4), type = "l")
#' lines(y, y, lty = 3)
#' 
soft_thresh <- function(x, lambda) {
  sign(x) * (pmax(abs(x) - lambda, 0))
}
#' Soft-threshold and normalize
#' 
#' Soft-threshold, and then normalize the result to L2-norm of 1. Care is taken with
#' vectors of L2-norm empirically equal to zero.
#' 
#' @details This function is not exported, as it is an internal utility function. Use
#' \code{spca:::soft_l2norm} to access it.
#' 
#' @param x Vector to be soft-thresholded
#' @param lambda Threshold value
#' 
#' @examples 
#' set.seed(1)
#' (z <- rnorm(10))
#' spca:::soft_l2norm(z, 0.5)
#' sqrt(sum(spca:::soft_l2norm(z, 0.5)^2))
#' 
#' (z <- rnorm(100))
#' spca:::soft_l2norm(z, 0.75)
#' sqrt(sum(spca:::soft_l2norm(z, 0.5)^2))
#' 
soft_l2norm <- function(x, lambda) {
  x <- soft_thresh(x = x, lambda = lambda)
  x_l2 <- sqrt(sum(x^2))
  if (abs(x_l2) > .Machine$double.eps^0.5) {
    x <- x / x_l2
  }
  x
}