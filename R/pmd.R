#' Rank-1 Penalized Matrix Decomposition
#' 
#' Implementation of the rank-1 penalized matrix decomposition of Witten et al. (2009). This
#' function applies lasso penalties to the left and right singular vectors. This function
#' is called by \code{\link{mpmd}}, which should be used in applications.
#' 
#' @param Z Matrix to be decomposed
#' @param c1 L1-norm bound for u, the left singular vector. Feasible solutions are available 
#' when values greater than or equal to 1. For values larger than \code{sqrt(nrow(Z)}, it has no effect.
#' @param c2 L1-norm bound for v, the right singular vector. Feasible solutions are available 
#' when values greater than or equal to 1. For values larger than \code{sqrt(ncol(Z)}, it has no effect.
#' @param maxit Maximum number of iterations
#' @param eps Stopping criterion, an absolute error tolerance on the mean squared reconstruction error
#' @param centre Logical indicating whether to centre the matrix Z using the overall mean 
#' before analysis
#' 
#' @return A list with the penalized singular value decomposition(\code{d}, \code{u}, \code{v}), 
#' the vector of errors (\code{error}), and the number of iterations (\code{iteration}).
#' 
#' @references 
#' Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized matrix decomposition, with applications to 
#' sparse principal components and canonical correlation analysis. \emph{Biostatistics}, 10(3), 515-534.
#' 
#' @examples 
#' 
#' ## Simple random matrix
#' set.seed(1)
#' Z <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' 
#' ## Default result
#' pmd(Z)
#' 
#' ## Equivalent to rank-1 SVD (almost)
#' pmd_Z <- pmd(Z, c1 = sqrt(nrow(Z)), c2 = sqrt(ncol(Z)))
#' all.equal(pmd_Z[-1], svd(Z, nu = 1, nv = 1)[-1], 
#'           tolerance = (.Machine$double.eps)^0.25)
#'          
#' ## Increasing c2 incrementally
#' pts <- seq(from = 1, to = sqrt(ncol(Z)), length.out = 10)
#' vmat <- matrix(NA, nrow = ncol(Z), ncol = length(pts))
#' for (i in seq_along(pts)) {
#'   vmat[, i] <- pmd(Z, c1 = sqrt(nrow(Z)), c2 = pts[i])$v
#' }
#' vmat
#' 
#' ## Increasing c1 incrementally
#' pts <- seq(from = 1, to = sqrt(nrow(Z)), length.out = 10)
#' umat <- matrix(NA, nrow = nrow(Z), ncol = length(pts))
#' for (i in seq_along(pts)) {
#'   umat[, i] <- pmd(Z, c1 = pts[i], c2 = sqrt(ncol(Z)))$u
#' }
#' umat
#' 
#' ## Testing reconstruction error
#' res <- pmd(Z, c1 = 3, c2 = 1.75)
#' plot(res$error)
#' 
#' @export
#' 
pmd <- function(Z, c1 = 1, c2 = 1, maxit = 100, 
                eps = sqrt(.Machine$double.eps), centre = FALSE) {
  
  ## As in PMA, deduct the overall mean of Z if requested
  if (centre) {
    Z <- Z - mean(Z)
  }
  
  ## Number of rows and columns
  m <- nrow(Z)
  n <- ncol(Z)
  
  ## Check c1 and c2
  if (c1 < 1) { # || c1 > sqrt(m)
    stop("Argument 'c1' is outside the allowed range.")
  }
  if (c2 < 1) { #  || c2 > sqrt(n)
    stop("Argument 'c2' is outside the allowed range.")
  }
  
  ## SVD for starting v (matrix with 1 column)
  Z_svd <- svd(Z, nu = 0, nv = 1)
  v <- Z_svd$v
  
  ## Vector for storing reconstruction error 
  err_vec <- double(maxit)
  
  ## Iteration
  i <- 0
  while (i < maxit) {
    
    ## Increase counter
    i <- i + 1
    
    ## Find lambda_1
    Zv <- Z %*% v
    lambda1 <- binary_search(x = Zv, c = c1)
    
    ## Update u
    u <- soft_l2norm(x = Zv, lambda = lambda1)
    
    ## Find lambda_2
    Zu <- t(Z) %*% u
    lambda2 <- binary_search(x = Zu, c = c2)
    
    ## Update v
    v <- soft_l2norm(x = Zu, lambda = lambda2)
    
    ## Calculate d
    d <- t(u) %*% Z %*% v
    
    ## Mean reconstruction error using Frobenius norm
    err_vec[i] <- mean((Z - drop(d) * tcrossprod(u, v))^2)
    
    ## Break if reduction is not >= eps
    if (i > 1 && err_vec[i - 1] - err_vec[i] < eps) {
      break
    }
  }
  
  ## Return
  list(d = d, u = u, v = v, error = err_vec[1:i], iterations = i)
}