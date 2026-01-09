##################################################
### Simulation
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
library(progress)
source("R/mrp_test.R")
source("R/dcf_test.R")

### Generate simulated data
generate_mfd <- function(n1 = 50, n2 = 50, m = 101, p = 200,
                         dist = "gaussian",
                         sparsity = 0.1,
                         delta = 0.1,
                         c_k = NULL,
                         df = 3) {
  # 1. Setup constants
  n <- n1 + n2
  rho <- 0.3   # Correlation parameter
  K <- 50   # Number of basis functions
  
  if (is.null(c_k)) {
    c_k <- rep(1, K)
  }
  
  # Define time points T = [0, 1]
  gr <- seq(0, 1, length.out = m)
  
  # 2. Construct the Fourier Basis Matrix (Phi)
  # Dimension: m x K
  # phi_1(t) = 1
  # phi_2l(t) = sqrt(2) * cos(l * pi * (2t - 1))
  # phi_2l_minus_1(t) = sqrt(2) * sin((l - 1) * pi * (2t - 1))
  Phi <- matrix(0, nrow = m, ncol = K)
  Phi[, 1] <- 1  # phi_1
  for (k in 2:K) {
    if (k %% 2 == 0) {
      # Even index: Cosine terms (phi_2, phi_4, ...)
      # If k = 2l, then l = k / 2
      l <- k / 2
      Phi[, k] <- sqrt(2) * cos(l * pi * (2 * gr - 1))
    } else {
      # Odd index: Sine terms (phi_3, phi_5, ...)
      # If k = 2l - 1, then 2l = k + 1 => l = (k + 1) / 2
      # The formula uses (l-1) inside sin
      l <- (k + 1) / 2
      Phi[, k] <- sqrt(2) * sin((l - 1) * pi * (2 * gr - 1))
    }
  }
  
  # 3. Generate Raw Coefficients (tilde_theta)
  # V_ij(t) = sum(tilde_theta_ijk * phi_k(t))
  # tilde_theta_ijk ~ N(0, k^-2) independently
  # Variance decreases as k increases.
  sd_vec <- sqrt(c_k) * 1/(1:K)  # Standard deviation = sqrt(1/k^2) = 1/k
  
  # Storage for raw coefficients: n x p x K
  tilde_theta <- array(0, dim = c(n, p, K))
  for (k in 1:K) {
    # Generate sub-Gaussian random samples
    if (dist == "gaussian") {
      sg_samples <- rnorm(n * p, mean = 0, sd = sd_vec[k])
    } else if (dist == "t") {
      sg_samples <- rt(n * p, df = df) * sqrt((df-2)/df) * sd_vec[k]
    }
    
    # Generate n*p random normal variables for the k-th basis coefficient
    tilde_theta[, , k] <- matrix(sg_samples, n, p)
  }
  
  # 4. Apply Autoregressive Correlation to Predictors
  # X_ij(t) = sum_{j'=1}^p rho^|j-j'| V_ij'(t)
  # This implies the coefficients theta_ijk are linear combinations of tilde_theta
  
  # Create Correlation Matrix Sigma (p x p)
  # Entry (u, v) is rho^|u-v|
  row_idx <- matrix(rep(1:p, p), nrow = p, ncol = p)
  col_idx <- t(row_idx)
  Sigma_rho <- rho^abs(row_idx - col_idx)
  
  # Calculate final coefficients theta (n x p x K)
  theta <- array(0, dim = c(n, p, K))
  
  # Apply correlation for each subject i and basis k
  # We are mixing across the 'p' dimension
  for (i in 1:n) {
    for (k in 1:K) {
      # tilde_theta[i, , k] is a vector of length p
      theta[i, , k] <- Sigma_rho %*% tilde_theta[i, , k]
    }
  }
  
  # 5. Construct the final Functional Data n-m-p X and Y
  # For each subject i and predictor j: X_ij(t) = Phi(t) * theta_ij
  data_arr <- array(0, dim = c(n, m, p))
  for (i in 1:n) {
    # theta[i, , ] is a p x K matrix
    # Phi is an m x K matrix
    # We want Result (m x p) = Phi * t(theta)
    # because X_ij(t) = sum_k theta_ijk * phi_k(t)
    coefs_i <- theta[i, , ] # p x K
    data_arr[i, , ] <- Phi %*% t(coefs_i)
  }
  X <- data_arr[1:n1, , ]
  Y <- data_arr[(n1+1):n, , ]
  
  
  # 6. Mean functions for Y
  if (sparsity > 0) {
    # Number of signal set
    q <- ceiling(p * sparsity)
    if (q < 2) {
      q <- 2
    }
    # Generate sub-Gaussian random samples and generate mean function
    theta_mean <- matrix(0, K, q)
    for (j in 1:q) {
      theta_mean[, j] <- runif(K, -delta, delta)
    }
    mean_ftn <- Phi %*% theta_mean  # m x q matrix
    
    # Add the mean function for Y
    for (i in 1:n2) {
      Y[i, , 1:q] <- Y[i, , 1:q] + mean_ftn
    }
  }
  
  out <- list(
    X1 = X,
    X2 = Y
  )
}




#-----------------------------
# Simulation parameters
#-----------------------------

m <- 101
p <- 200

num_sim <- 200

### Parameters for generating data
set.seed(100)
c_k <- runif(50, 1, 10)

### Combination of simulation parameters
samp_list <- c("balanced","unbalanced")
dist_list <- c("gaussian","t")
delta_list <- c(0.3, 0.45)
sparsity_list <- c(0, 0.05, 0.1, 0.2)  # sparsity
sim_comb <- expand_grid(samp = samp_list,
                        dist = dist_list,
                        delta = delta_list,
                        sparsity = sparsity_list)
# Remove type-I error setting which is duplicated for different delta
sim_comb <- sim_comb %>% 
  filter(!(delta != delta_list[1] & sparsity == 0))


### Parameters for methods
n_basis_list <- 4:8
lambda_list <- 10^seq(-3, -1, length.out = 10)
n_cores <- 20   # number of cores


#-----------------------------
# Empirical size and power
# - Size: sparsity = 0
# - Power: sparsity > 0
#-----------------------------
result <- list()
for (idx in 1:nrow(sim_comb)) {
  # Simulation parameters
  samp <- sim_comb$samp[idx]   # sampling case
  dist <- sim_comb$dist[idx]   # distribution
  delta <- sim_comb$delta[idx]   # signal strength
  sparsity <- sim_comb$sparsity[idx]   # sparsity parameter
  
  if (samp == "balanced") {
    # Balanced samples
    n1 <- 30
    n2 <- 30
  } else if (samp == "unbalanced") {
    # Un-balanced samples
    n1 <- 25
    n2 <- 35
  }
  
  # Progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent [Execution time :elapsedfull || Estimated time remaining: :eta]",
    total = num_sim,   # total number of ticks to complete (default 100)
    clear = FALSE,   # whether to clear the progress bar on completion (default TRUE)
    width = 80   # width of the progress bar
  )
  progress <- function(n){
    pb$tick()
  } 
  
  p_values <- data.frame(
    dspt_hfd = rep(NA, num_sim),
    dcf = rep(NA, num_sim),
    mrp = rep(NA, num_sim)
  )
  test_stat <- data.frame(
    dspt_hfd = rep(NA, num_sim),
    dcf = rep(NA, num_sim),
    mrp = rep(NA, num_sim)
  )
  active_set_len <- data.frame(
    dspt_hfd = rep(NA, num_sim)
  )
  beta_l1_norm <- data.frame(
    dspt_hfd = rep(NA, num_sim)
  )
  for (sim in 1:num_sim) {
    set.seed(sim)
    
    # Show the progress bar
    progress(sim)
    
    # Generate simulation data
    obj <- generate_mfd(n1, n2, m, p,
                        dist = dist,
                        sparsity = sparsity,
                        c_k = c_k,
                        delta = delta)
    X1 <- obj$X1
    X2 <- obj$X2
    
    # DSPT-HFD; 5-fold CV
    set.seed(sim)
    obj_dspt_hfd <- dspt_hfd(X1,
                             X2,
                             penalty = "scad",
                             tuning = TRUE,
                             tune_method = "cv",
                             n_basis_list = n_basis_list,
                             lambda_list = lambda_list,
                             n_cores = n_cores)
    p_values$dspt_hfd[sim] <- obj_dspt_hfd$p.value
    test_stat$dspt_hfd[sim] <- obj_dspt_hfd$obj_test$statistic
    active_set_len$dspt_hfd[sim] <- length(obj_dspt_hfd$obj_proj$active_set)
    beta_l1_norm$dspt_hfd[sim] <- sum(abs(obj_dspt_hfd$obj_proj$beta_hat))
    
    # DCF test
    set.seed(sim)
    obj_dcf <- dcf_test(X1,
                        X2,
                        n_basis_list = n_basis_list)
    p_values$dcf[sim] <- obj_dcf$p.value
    test_stat$dcf[sim] <- obj_dcf$statistic
    
    # MRP test
    set.seed(sim)
    obj_mrp <- mrp_test(X1, X2, n_cores = n_cores)
    p_values$mrp[sim] <- obj_mrp$p.value
    test_stat$mrp[sim] <- obj_mrp$statistic
  }
  print(colMeans(p_values < 0.05))   # empirical power
  
  result[[idx]] <- list(
    p_values = p_values,
    test_stat = test_stat,
    active_set_len = active_set_len,
    beta_l1_norm = beta_l1_norm
  )
  
  save(result, sim_comb,
       file = paste0("RData/sim_num_sim=", num_sim, ".RData"))  
}


cbind(
  sim_comb,
  sapply(result, function(res){
    colMeans(res$p_values < 0.05)
  }) %>%
    t() %>% 
    data.frame()
) %>% 
  print()

sapply(result, function(res){
  colMeans(res$active_set_len)
}) %>%
  t() %>% 
  print()

sapply(result, function(res){
  colMeans(res$beta_l1_norm)
}) %>%
  t() %>% 
  print()


cbind(
  sim_comb,
  sapply(result, function(res){
    colMeans(res$p_values[, -2] < 0.05)
  }) %>%
    t()
) %>% 
  xtable::xtable(digits = 3)

