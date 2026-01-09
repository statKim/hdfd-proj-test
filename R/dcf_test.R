# Distribution/correlation-free test
dcf_test <- function(X1,
                     X2,
                     B = 10000,
                     grid = NULL,
                     basis = "bspline",
                     n_basis = NULL,
                     n_basis_list = 4:10,
                     K = 5) {
  n1 <- dim(X1)[1]
  n2 <- dim(X2)[1]
  m <- dim(X1)[2]
  p <- dim(X1)[3]
  data_name <- paste(deparse(substitute(X1)), "and", deparse(substitute(X2)))
  
  # Find the number of basis functions
  if (is.null(n_basis)) {
    # # K-fold CV
    # fold_1 <- sample(1:K, n1, replace = T)
    # fold_2 <- sample(1:K, n2, replace = T)
    # 
    # recon_error <- rep(NA, length(n_basis_list))
    # for (i in 1:length(n_basis_list)) {
    #   n_basis <- n_basis_list[i]
    #   
    #   cv_err <- 0
    #   for (k in 1:K) {
    #     idx_1 <- which(fold_1 == k)
    #     idx_2 <- which(fold_2 == k)
    #     
    #     basis_obj1 <- basis_mfd(X1[-idx_1, , , drop=FALSE],
    #                             grid = grid,
    #                             basis = basis,
    #                             n_basis = n_basis,
    #                             gram = TRUE)
    #     X1_coef_test <- predict(basis_obj1, X1[idx_1, , , drop=FALSE])
    #     
    #     basis_obj2 <- basis_mfd(X2[-idx_2, , , drop=FALSE],
    #                             grid = grid,
    #                             basis = basis,
    #                             n_basis = n_basis,
    #                             gram = TRUE)
    #     X2_coef_test <- predict(basis_obj2, X2[idx_2, , , drop=FALSE])
    #     
    #     basis_mat <- basis_obj1$phi
    #     
    #     for (j in 1:p) {
    #       pred_1 <- X1_coef_test[, ((j-1)*n_basis+1):(j*n_basis)] %*% t(basis_mat)
    #       pred_2 <- X2_coef_test[, ((j-1)*n_basis+1):(j*n_basis)] %*% t(basis_mat)
    #       
    #       cv_err <- cv_err + sum((X1[idx_1, , j] - pred_1)^2) + 
    #         sum((X2[idx_2, , j] - pred_2)^2)
    #     }
    #   }
    #   cv_err <- cv_err / ((n1+n2)*m*p)
    #   recon_error[i] <- cv_err
    # }
    # n_basis <- n_basis_list[which.min(recon_error)]
    
    recon_error <- rep(NA, length(n_basis_list))
    for (i in 1:length(n_basis_list)) {
      n_basis <- n_basis_list[i]
      basis_obj1 <- basis_mfd(X1,
                              grid = grid,
                              basis = basis,
                              n_basis = n_basis,
                              gram = TRUE)
      # basis_mat <- fda::eval.basis(seq(0, 1, length.out = m), basis_obj1$basis_ftn)
      basis_mat <- basis_obj1$phi
      X1_coef <- basis_obj1$X_coef

      basis_obj2 <- basis_mfd(X2,
                              grid = grid,
                              basis = basis,
                              n_basis = n_basis,
                              gram = TRUE)
      X2_coef <- basis_obj2$X_coef

      err <- 0
      for (j in 1:p) {
        pred_1 <- X1_coef[, ((j-1)*n_basis+1):(j*n_basis)] %*% t(basis_mat)
        pred_2 <- X2_coef[, ((j-1)*n_basis+1):(j*n_basis)] %*% t(basis_mat)

        err <- err + sum((X1[, , j] - pred_1)^2) + sum((X2[, , j] - pred_2)^2)
      }
      err <- err / ((n1+n2)*m*p)
      recon_error[i] <- err
      # plot(X1[5, , 2], type = "l")
      # lines(pred[5, ], col = 2)
      # rowMeans((X1[, , j] - pred)^2)
    }
    n_basis <- n_basis_list[which.min(recon_error)]
  }
  
  # Basis representation for each functional covariate
  basis_obj1 <- basis_mfd(X1,
                          grid = grid,
                          basis = basis,
                          n_basis = n_basis,
                          gram = TRUE)
  X1_coef <- basis_obj1$X_coef
  
  basis_obj2 <- basis_mfd(X2,
                          grid = grid,
                          basis = basis,
                          n_basis = n_basis,
                          gram = TRUE)
  X2_coef <- basis_obj2$X_coef
  
  # Compute test statistic
  s_X1 <- colSums(X1_coef) / sqrt(n1)
  s_X2 <- colSums(X2_coef) / sqrt(n2)
  test_stat <- max(abs(s_X1 - sqrt(n1/n2)*s_X2))
  
  # Multiplier bootstrap
  test_stat_boot <- rep(NA, B)
  for (b in 1:B) {
    e1 <- rnorm(n1)
    e2 <- rnorm(n2)
    
    # se_X1 <- (matrix(e1, nrow = 1) %*% scale(X1_coef, center = T, scale = F)) / sqrt(n1)
    # se_X2 <- (matrix(e2, nrow = 1) %*% scale(X2_coef, center = T, scale = F)) / sqrt(n2)
    
    se_X1 <- (matrix(e1, nrow = 1) %*% X1_coef)/sqrt(n1) - sum(e1)/n1*s_X1
    se_X2 <- (matrix(e2, nrow = 1) %*% X2_coef)/sqrt(n2) - sum(e2)/n2*s_X2
    
    test_stat_boot[b] <- max(abs(se_X1 - sqrt(n1/n2)*se_X2))
  }
  
  # p-value for one-sided test (test statisic only takes positive values)
  p_value <- mean(test_stat_boot >= test_stat)
  
  # Output
  out <- list(
    statistic = test_stat,
    p.value = p_value,
    method = "Distribution and Correlation-Free Two-Sample Test for High-Dimensional Functional Means",
    data.name = data_name,
    test_stat_boot = test_stat_boot,
    n_basis = n_basis,
    recon_error = recon_error
  )
  class(out) <- "htest"
  
  return(out)
}


