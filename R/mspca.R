#' Rank-k Sparse Principal Components Analysis
#' 
#' Multifactor version of the rank-1 SPCA of \code{\link{spca}}.Cmponents are found by 
#' applying \code{\link{spca}} to the original matrix after deflation by deducting 
#' the components already found. It is possible to apply different penalties to 
#' different components.
#' 
#' @details This is similar to \code{\link[PMA]{SPC}}, but differs especially with respect to
#' the centring of the matrix.
#' 
#' @param Z Matrix to be decomposed
#' @param k Required rank of the result
#' @param c  L1-norm bound for V (greater than or equal to 1), either length-1, or 
#' with k entries (one for each component). Feasible solutions are available 
#' for values greater than or equal to 1. For values larger than \code{sqrt(ncow(Z)}, it has no effect.
#' @param maxit Maximum number of iterations
#' @param eps Stopping criterion, and absolute error tolerance on the mean squared reconstruction error
#' @param center Logical indicating whether to column-centre the matrix Z 
#' @param scale Logical indicating whether to set the standard deviations of the columns of Z equal
#' to one before analysis
#' 
#' @return A list with the following components:
#' \describe{
#'  \item{scores}{The data matrix after projection to the principal component space}
#'  \item{loadings}{The matrix of component loadings}
#'  \item{sdev}{The standard deviations of each of the principal components}
#'  \item{pve}{Data frame giving the proportion of variance explained by the obtained components}
#' }
#' 
#' @examples 
#' 
#' ## Random matrix example, with n > p
#' set.seed(1)
#' Z <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' mspca(Z, c = 1.25, k = 5)
#' 
#' ## Random matrix example with n < p
#' Z2 <- matrix(rnorm(100), nrow = 5, ncol = 20)
#' mspca(Z2, c = 2.5, k = 5)
#' 
#' ## Example with different c for components 1 and 2
#' mspca(Z, k = 2, c = c(1.5, 1.1))
#' 
#' ## Comparison to PCA
#' summary(prcomp(Z2))
#' mspca(Z2, k = 5, c = max(dim(Z2)))$pve
#' 
#' @export
#' 
mspca <- function(Z, k = 2, c = 1, maxit = 20, eps = sqrt(.Machine$double.eps), 
                  center = TRUE, scale = FALSE) {
  
  ## Check
  if (k > min(dim(Z))) {
    stop(paste("Argument 'k' should be at most", min(dim(Z))))
  }
  
  ## Mean-center and/or scale columns to standard devaition of 1
  if (center || scale) {
    Z <- scale(Z, center = center, scale = scale)
  }
  
  ## Construct c as a vector, and check
  if (length(c) == 1) {
    c <- rep(c, k)
  }
  if (length(c) != k) {
    stop("Argument 'c' should be either length 1 or length 'k'.")
  }

  ## Initialization of deflated matrix, iteration counter, SVD components
  R <- Z
  i <- 0
  d <- double(k)
  umat <- matrix(NA, nrow = nrow(Z), ncol = k)
  vmat <- matrix(NA, nrow = ncol(Z), ncol = k)
  
  ## Initialization of variance explained
  cum_var_expl <- double(k)
  
  ## Iterate
  for (i in seq_len(k)) {
    
    ## Get rank 1 SPCA of R (never scale here)
    R_pmd <- spca(Z = R, c = c[i], maxit = maxit, eps = eps, center = FALSE, scale = FALSE)
    
    ## Update d, u, v
    d[i] <- R_pmd$d
    umat[, i] <- R_pmd$u
    vmat[, i] <- R_pmd$v
    
    # Deflate the current residual matrix
    R <- R - R_pmd$d[1] * tcrossprod(R_pmd$u, R_pmd$v)
    
    # Calculate the approximation Z_k and store variance explained
    vcurr <- vmat[, seq_len(i), drop = FALSE]
    Z_k <- Z %*% vcurr %*% solve(t(vcurr) %*% vcurr) %*% t(vcurr)
    cum_var_expl[i] <- sum(Z_k^2)
  }
  
  ## Calculate PC scores, standard deviations of PCs
  colnames(vmat) <- paste0("PC", seq_len(k))
  scores <- as.data.frame(Z %*% vmat %*% solve(t(vmat) %*% vmat) %*% t(vmat))
  colnames(scores) <- paste0("PC", seq_len(k))
  sdev <- d / sqrt((nrow(Z) - 1))
  
  ## Total variance to be explained, and proportion explained
  tot_var <- sum(svd(Z, nu = 0, nv = 0)$d^2)
  pve <- data.frame(PC = seq_len(k), 
                    std_deviation = sdev, 
                    prop_variance = diff(c(0, cum_var_expl)) / tot_var,
                    cum_prop = cum_var_expl / tot_var)
  
  ## Return
  rownames(vmat) <- colnames(Z)
  rownames(scores) <- rownames(Z)
  list(scores = scores, loadings = vmat, sdev = sdev, pve = pve, svd = list(d = d, u = umat, v = vmat))
}