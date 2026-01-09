##################################################
### Simulation
### - DCF test paper setting
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
library(progress)
source("R/mrp_test.R")
source("R/dcf_test.R")

### Generate simulated data
generate_mfd <- function(n1 = 100, n2 = 150, m = 51, p = 300,
                         tau, delta_k, psi_x, psi_y,
                         n_cores = 1) {
  gr <- seq(0, 1, length.out = m)
  k <- nrow(tau)
  
  # Fourier basis
  phi <- matrix(1, m, k)
  for (i in 1:(k/2)) {
    phi[, 2*i] <- sqrt(2) * cos(i*pi*(2*gr-1))
    if (i > 1) {
      phi[, (2*i)-1] <- sqrt(2) * sin((i-1)*pi*(2*gr-1))
    }
  }
  
  # # Mean functions
  # mu_x <- matrix(0, m, p)
  # mu_y <- matrix(0, m, p)
  # mu_y[, 1:beta_p] <- phi %*% tau
  
  # Generate scores, theta_tilde
  theta_tilde_x <- array(NA, c(n1, p, k))
  theta_tilde_y <- array(NA, c(n2, p, k))
  for (j in 1:p) {
    if (j <= p/2) {
      theta_tilde_x[, j, ] <- t( (8*(delta_k^2))^(-1/2) * matrix(rchisq(n1*k, df = 4) - 4, k, n1) )
      theta_tilde_y[, j, ] <- t( (8*(delta_k^2))^(-1/2) * matrix(rchisq(n2*k, df = 4) - 4, k, n2) )
    } else {
      theta_tilde_x[, j, ] <- t( (5*(delta_k^2)/3)^(-1/2) * matrix(rt(n1*k, df = 5), k, n1) )
      theta_tilde_y[, j, ] <- t( (5*(delta_k^2)/3)^(-1/2) * matrix(rt(n2*k, df = 5), k, n2) )
    }
  }
  
  # AR correlation among the functional predictors
  rho <- 0.3
  
  # Generate X and Y
  # Parallel computing
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  data_obj <- foreach::foreach(i = 1:max(n1, n2)) %dopar% {
    out <- list()
    if (i <= n1) {
      # Generate X
      # k x p
      theta_x_i <- sapply(1:p, function(j) {
        matrix(sqrt(psi_x[i, j] * psi_x[i, ]) * rho^(abs(j - 1:p)), nrow = 1) %*%
          theta_tilde_x[i, , ]
      })
      
      out$X <- phi %*% theta_x_i
    }
    
    if (i <= n2) {
      # Generate Y
      # k x p
      theta_y_i <- sapply(1:p, function(j) {
        matrix(sqrt(psi_y[i, j] * psi_y[i, ]) * rho^(abs(j - 1:p)), nrow = 1) %*%
          theta_tilde_y[i, , ]
      })
      
      out$Y <- phi %*% (tau + theta_y_i)
    }
    
    return(out)
  }
  stopCluster(cl)
  
  X <- array(NA, c(n1, m, p))
  Y <- array(NA, c(n2, m, p))
  for (i in 1:max(n1, n2)) {
    if (i <= n1) {
      X[i, , ] <- data_obj[[i]]$X
    }
    
    if (i <= n2) {
      Y[i, , ] <- data_obj[[i]]$Y
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
p <- 300

num_sim <- 200

### Combination of simulation parameters
samp_list <- c("balanced","unbalanced")
delta_list <- c(0.05, 0.1)
sparsity_list <- c(0, 0.05, 0.1, 0.2)  # sparsity
sim_comb <- expand_grid(samp = samp_list,
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
  delta <- sim_comb$delta[idx]   # signal strength
  beta <- sim_comb$sparsity[idx]   # sparsity parameter
  
  if (samp == "balanced") {
    # Balanced samples
    n1 <- 30
    n2 <- 30
  } else if (samp == "unbalanced") {
    # Un-balanced samples
    n1 <- 25
    n2 <- 35
  }
  
  # Number of different means for Y
  beta_p <- floor(beta * p)
  
  # Mean scores for Y
  set.seed(500)
  k <- 50
  tau <- matrix(0, k, p)
  if (beta_p > 0) {
    tau[, 1:beta_p] <- matrix(runif(k*beta_p, -delta, delta), k, beta_p)
  }
  
  # Functional predictors
  delta_k <- sample(1:k, k)
  psi_x <- matrix(runif(n1*p, 1, 2), n1, p)
  psi_y <- matrix(runif(n2*p, 1, 3), n2, p)
  
  
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
                        tau, delta_k, psi_x, psi_y,
                        n_cores)
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
       file = paste0("RData/sim_dcf_num_sim=", num_sim, ".RData"))
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

