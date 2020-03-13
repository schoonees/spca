#' Binary Search for Penalization Parameter
#' 
#' Performs binary search to find a value of the shrinkage parameter
#' lambda which obtains a required L1-norm value after soft-thresholding and
#' normalization to length one. The value of lambda found is the smallest value 
#' that attains the required constraint value.
#' 
#' @details Note that in the \code{\link{PMA}} package, the binary search will find the largest 
#' value of lambda that imposes the required limit on the L1-norm after soft-thresholding in case
#' this norm has regions of lambda for which it is constant. This leads to sparser solutions
#' than what is necessary.
#' 
#' @param x Vector to be soft-thresholded and normalized
#' @param c Desired value of the L1-norm after soft-thresholding and normalization
#' @param maxit Maximum number of iterations for the binary search
#' 
#' @examples 
#' 
#' ## Random vector to soft-threshold and normalize
#' set.seed(1)
#' (z <- rnorm(25))
#' 
#' ## Binary search for L1-norm of 1 (returning lambda)
#' (lambda <- binary_search(z, 1))
#' 
#' ## Graphical illustration
#' pts <- seq(from = 0, to = max(abs(z)), length.out = 100)
#' plot(pts, sapply(pts, function(a) sum(abs(soft_l2norm(z, a)))),
#'      ylab = "L1-norm", xlab = expression(lambda))
#' abline(v = binary_search(z, 1, maxit = 1), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 2), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 3), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 4), col = 2, lty = 2)
#' abline(v = binary_search(z, 1, maxit = 15), col = 2, lty = 2)
#' 
#' ## More tests
#' sum(abs(soft_l2norm(z, lambda)))
#' sum(abs(soft_l2norm(z, binary_search(z, 2))))
#' sum(abs(soft_l2norm(z, binary_search(z, 3))))
#' 
#' @export
binary_search <- function(x, c, maxit = 100) {
  
  ## Check if x already meets the condition, or whether c is zero, and return 0 in either case
  x_l1 <- sum(abs(x))
  if (x_l1 <= c || c == 0) {
    return(0)
  }
  
  ## Initialize search boundaries
  lo <- 0
  hi <- max(abs(x))
  i <- 0
  
  ## Iterate
  while(i < maxit) {
    
    ## Increase iteration counts
    i <- i + 1
    
    ## Get mean of current bounds, and evaluate L1 norm there
    m <- lo + (hi - lo)/2
    m_l1 <- sum(abs(soft_l2norm(x, m)))
    
    ## Update either lower or upper threshold
    if (m_l1 <= c) { ## Note: must be <=
      hi <- m
    } else {
      lo <- m
    }
    
    ## Break if lo and hi have converged
    if (abs(lo - hi) <= sqrt(.Machine$double.eps)) {
      break
    }
  }
  
  return(lo + (hi - lo)/2)
}
