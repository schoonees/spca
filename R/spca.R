#' Rank-1 Sparse Principal Components Analysis
#' 
#' Uses the L1-penalized matrix decomposition in \code{\link{pmd}} to obtain sparse principal 
#' components analysis of rank-1. This function is called by \code{\link{mspca}}, which should
#' be used in applications.
#' 
#' @details Note that in the \code{\link[PMA]{SPC}} implementation, the overall mean of the input
#' matrix is set equal to zero before analysis when the (default) \code{center = TRUE}. This 
#' is not done here, or even given as an option.
#' 
#' @param Z Matrix to be decomposed
#' @param c  L1-norm bound for u, the left singular vector. Feasible solutions are available 
#' when values greater than or equal to 1. For values larger than \code{sqrt(nrow(Z)}, it has no effect.
#' @param vstart Optional starting value for the v vector
#' @param maxit Maximum number of iterations
#' @param eps Stopping criterion, and absolute error tolerance on the mean squared reconstruction error
#' @param center Logical indicating whether the columns of \code{Z} should be mean-centered
#' @param scale Logical indicating whether the standard deviations of the column of \code{Z} should be
#' set equal to 1
#' 
#' @return A list with the penalized singular value decomposition(\code{d}, \code{u}, \code{v}), 
#' the vector of errors (\code{error}), and the number of iterations (\code{iteration}). The output from
#' \code{\link{mspca}} is more appropriate for PCA, and that function should be used.
#' 
#' @references 
#' Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized matrix decomposition, with applications to 
#' sparse principal components and canonical correlation analysis. \emph{Biostatistics}, 10(3), 515-534.
#' 
#' @examples 
#' set.seed(1)
#' Z <- matrix(rnorm(600), nrow = 20, ncol = 30)
#' 
#' spca(Z, c = 1.6)
#' spca(Z, c = 1.5)
#' 
#' @export
#' 
spca <- function(Z, c = 1, vstart = NULL, maxit = 100, eps = sqrt(.Machine$double.eps), 
                 center = TRUE, scale = FALSE) {
  
  ## Mean-center and/or scale columns to standard deviation of 1
  if (center || scale) {
    Z <- scale(Z, center = center, scale = scale)
  }
  
  ## Number of rows and columns
  m <- nrow(Z)
  n <- ncol(Z)
  
  ## Check c
  if (c < 1) {
    stop("Argument 'c' is outside the allowed range.")
  }
  
  ## Starting value for v (matrix with 1 column): Use ordinary SVD if not supplied
  if (is.null(vstart)) {
    Z_svd <- svd(Z, nu = 0, nv = 1)
    v <- Z_svd$v
  } else {
    v <- matrix(vstart, ncol = 1, nrow = ncol(Z))
  }
  
  ## Initialize vector for storing reconstruction error, and iteration counter
  err_vec <- double(maxit)
  i <- 0
  
  ## Iteration
  while (i < maxit) {
    
    ## Increase counter
    i <- i + 1
    
    ## Update v (with zero penalty)
    Zv <- Z %*% v
    u <- soft_l2norm(x = Zv, lambda = 0)
    
    ## Find lambda_2
    Zu <- t(Z) %*% u
    lambda2 <- binary_search(x = Zu, c = c)
    
    ## Update v subject to penalty
    v <- soft_l2norm(x = Zu, lambda = lambda2)
    
    ## Update d
    d <- t(u) %*% Z %*% v
    
    ## Mean reconstruction error using Frobenius norm
    err_vec[i] <- mean((Z - drop(d) * tcrossprod(u, v))^2)
    
    ## Break if reduction is not >= eps
    if (i > 1 && err_vec[i - 1] - err_vec[i] < eps) {
      break
    }
  }
  
  ## Warning if it was not monotonically decreasing
  if (any(diff(err_vec) > sqrt(.Machine$double.eps))) {
    warning("Error did not decrease monotonically")
  }
  
  ## Return result
  list(d = d, u = u, v = v, error = err_vec[seq_len(i)], iterations = i)
}