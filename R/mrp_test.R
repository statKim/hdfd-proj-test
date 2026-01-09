library(fda)
Rcpp::sourceCpp("src/mrp_test_cpp.cpp")

#' MRP Test for High-Dimensional Functional Data
#'
#' @param X Array of dimensions c(n1, m, p). Sample 1.
#' @param Y Array of dimensions c(n2, m, p). Sample 2.
#' @param n_cores Integer. Number of cores for parallel computing.
#' @param num_grid Integer. Number of grid points (m). 
#'        If NULL, inferred from X dimensions.
#' @param a_ou Numeric. The parameter 'a' for the Ornstein-Uhlenbeck covariance process.
#'        Default is 1 as per simulation setups in common literature.
#'
#' @return A list containing the test statistic, p-value, and estimated variance.
mrp_test <- function(X, Y, n_cores = 1, num_grid = 51, a_ou = 1) {
  
  # --- Input Validation and Dimensions ---
  dim_x <- dim(X)
  dim_y <- dim(Y)
  data_name <- paste(deparse(substitute(X)), "and", deparse(substitute(Y)))
  
  if (length(dim_x) != 3 || length(dim_y) != 3) {
    stop("Input X and Y must be 3-dimensional arrays (n x m x p).")
  }
  
  n1 <- dim_x[1]
  m  <- dim_x[2]
  p  <- dim_x[3]
  n2 <- dim_y[1]
  gr <- seq(0, 1, length.out = m)
  
  if (dim_y[2] != m || dim_y[3] != p) {
    stop("Dimensions m (timepoints) and p (covariates) must match between X and Y.")
  }
  
  if (n1 < 3 || n2 < 3) {
    stop("Sample size must be at least 3 to calculate variance estimators.")
  }
  
  # Smoothing using B-splines
  n_basis <- min(round(m/5), 20)
  basis_order <- 4  # cubic B-spline
  bspline_basis <- create.bspline.basis(rangeval = c(0, 1),
                                        nbasis = n_basis,
                                        norder = basis_order)
  gr_recon <- seq(0, 1, length.out = num_grid)
  X_recon <- array(NA, dim = c(n1, num_grid, p))
  Y_recon <- array(NA, dim = c(n2, num_grid, p))
  for (j in 1:p) {
    fit <- smooth.basis(gr, t(X[, , j]), bspline_basis)
    X_recon[, , j] <- t(eval.fd(gr_recon, fit$fd))
    
    fit <- smooth.basis(gr, t(Y[, , j]), bspline_basis)
    Y_recon[, , j] <- t(eval.fd(gr_recon, fit$fd))
  }
  m <- num_grid
  gr <- gr_recon
  X <- X_recon
  Y <- Y_recon
  
  
  # --- Covariance Process Construction (V Matrix) ---
  # Stationary Ornstein-Uhlenbeck process: v(s,t) = a^-1 * exp(-a|s-t|)
  # We ignore the scalar constant a^-1 for the test statistic scale invariance, 
  # but strictly following paper def: v(s,t) propto exp(-a|s-t|)
  # We construct the matrix V where V_ij = v(t_i, t_j) * delta
  # Integration approximation: Trapezoidal rule weights
  # Delta = 1/(m-1)
  dist_mat <- abs(outer(gr, gr, "-"))
  Cov_mat <- (1/a_ou) * exp(-a_ou * dist_mat)
  
  # Apply integration weights
  # Weights w = [h/2, h, h, ..., h, h/2]
  h <- 1 / (m - 1)
  weights <- rep(h, m)
  weights[1] <- h / 2
  weights[m] <- h / 2
  
  # To approximate double integral X^T(s) Y(t) v(s,t) ds dt
  # We calculate sum_{i,j} X_i^T Y_j V_ij * w_i * w_j
  # This is equivalent to X^T * (diag(w) * Cov * diag(w)) * Y ? No.
  # Integral is X^T * W * C * W * Y?
  # Matrix form: Integral approx = X_vec^T * Omega * Y_vec
  # Let V_discrete = diag(sqrt(w)) * Cov_mat * diag(sqrt(w))?
  # Simplest form: V_discrete[i,j] = Cov_mat[i,j] * weights[i] * weights[j]
  # Then integral = sum(X[i]*Y[j]*V_discrete[i,j]) = X^T * V_discrete * Y
  W_diag <- diag(weights)
  V_discrete <- W_diag %*% Cov_mat %*% W_diag
  
  # --- Data Preprocessing ---
  # R arrays are (n, m, p). C++ Armadillo cubes are (rows, cols, slices).
  # We want each sample to be a contiguous slice for efficiency.
  # Permute to (m, p, n) so that X[,,i] in R becomes slice i in C++.
  X_perm <- aperm(X, c(2, 3, 1)) # (m, p, n1)
  Y_perm <- aperm(Y, c(2, 3, 1)) # (m, p, n2)
  
  # Pre-compute centered data for Variance calculation efficiency
  # Compute global means
  X_mean <- apply(X_perm, c(1, 2), mean)
  Y_mean <- apply(Y_perm, c(1, 2), mean)
  
  # Center the data
  X_c <- sweep(X_perm, 1:2, X_mean, "-")
  Y_c <- sweep(Y_perm, 1:2, Y_mean, "-")
  
  # --- Run C++ Worker ---
  res <- mrp_cpp_worker(X_perm, Y_perm, X_c, Y_c, V_discrete, n_cores)
  
  mrp_val <- res$mrp_stat
  sigma_sq <- res$sigma_sq
  
  # --- Result Calculation ---
  # According to Corollary 1: Q = MRP / sigma -> N(0, 1) under H0
  # Reject H0 if MRP is large.
  # Since MRP measures distance, it is expected to be positive under H1.
  # The test is one-sided upper tail.
  sigma <- sqrt(max(0, sigma_sq)) # Ensure non-negative before sqrt
  
  # test_stat <- 0
  # p_value <- 0
  
  if (sigma > 1e-8) {
    test_stat <- mrp_val / sigma
    p_value <- 1 - pnorm(test_stat)
  } else {
    warning("Estimated variance is effectively zero.")
    # test_stat <- mrp_val
    # # If variance is 0, and stat > 0, p-value is 0.
    # p_value <- ifelse(mrp_val > 1e-10, 0, 1)
    test_stat <- Inf
    p_value <- 1
  }
  
  # Output
  out <- list(
    statistic = test_stat,
    p.value = p_value,
    method = "Multi-resolution Projection Test for High-Dimensional Functional Data",
    data.name = data_name,
    mrp_value = mrp_val,
    sigma_sq = sigma_sq
  )
  class(out) <- "htest"
  
  # out <- list(
  #   test_stat = test_stat,
  #   p_value = p_value,
  #   mrp_value = mrp_val,
  #   sigma_sq = sigma_sq
  # )
  return(out)
}
