#' Multifactor Penalized Matrix Decomposition
#' 
#' Multifactor version of the rank-1 L1-penalized matrix decomposition in 
#' \code{\link{pmd}}. Components are found by applying \code{\link{pmd}} 
#' to the original matrix after deflation by deducting the components already 
#' found. 
#' 
#' @param Z Matrix to be decomposed
#' @param k Required rank of the result
#' @param c1 L1-norm bound for U (greater than or equal to 1), either length-1, or 
#' with k entries (one for each component). 
#' @param c2 L1-norm bound for V (greater than or equal to 1), either length-1, or 
#' with k entries (one for each component)
#' @param maxit Maximum number of iterations
#' @param eps Stopping criterion, an absolute error tolerance on the mean squared reconstruction error
#' @param centre Logical indicating whether to centre the matrix Z using the overall mean 
#' before analysis
#' 
#' @return A list with the penalized singular value decomposition (\code{d}, \code{U}, \code{V}).
#' 
#' @examples 
#' 
#' ## Simple random test matrix
#' set.seed(1)
#' Z <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' 
#' ## Ordinary SVD (equivalent up to changes in sign)
#' mpmd(Z, c1 = 5, c2 = 3, k = 5)
#' svd(Z)
#' 
#' ## Test with constant c1 and c2
#' mpmd(Z, c1 = 2, c2 = 1.25, k = 5)
#' 
#' ## Test with different c1 and c2 for different components
#' mpmd(Z, c1 = c(3, 1.25), c2 = c(2, 1.25), k = 2)
#' 
#' @export
#' 
mpmd <- function(Z, k = 2, c1 = 1, c2 = 1, maxit = 20, 
                 eps = sqrt(.Machine$double.eps), centre = FALSE) {
  
  ## Check
  if (k > min(dim(Z))) {
    stop(paste("Argument 'k' should be at most", min(dim(Z))))
  }
  
  ## Construct c1 and c2 as vectors, and check
  if (length(c1) == 1) {
    c1 <- rep(c1, k)
  }
  if (length(c2) == 1) {
    c2 <- rep(c2, k)
  }
  if (length(c1) != k) {
    stop("Argument 'c1' should be either length 1 or length 'k'.")
  }
  if (length(c2) != k) {
    stop("Argument 'c2' should be either length 1 or length 'k'.")
  }
  
  ## Initialize
  R <- Z
  i <- 0
  d <- double(k)
  umat <- matrix(NA, nrow = nrow(Z), ncol = k)
  vmat <- matrix(NA, nrow = ncol(Z), ncol = k)
  
  ## Iterate
  for (i in seq_len(k)) {
    
    ## Get rank 1 PMD of R, using potentially different c1 and/or c2
    R_pmd <- pmd(Z = R, c1 = c1[i], c2 = c2[i], maxit = maxit, eps = eps, centre = centre)
    
    ## Update d, u and v
    d[i] <- R_pmd$d
    umat[, i] <- R_pmd$u
    vmat[, i] <- R_pmd$v
    
    ## Deflate R using new component
    R <- R - R_pmd$d[1] * tcrossprod(R_pmd$u, R_pmd$v)
  }
  
  ## Return
  list(d = d, u = umat, v = vmat)
}