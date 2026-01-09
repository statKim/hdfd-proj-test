##################################################
### Simulation
### - MRP test paper setting
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
library(progress)
source("R/mrp_test.R")
source("R/dcf_test.R")

### Generate simulated data
generate_mfd <- function(n1 = 25, n2 = 25, m = 51, p = 20, 
                         setting = 3, case = 1,
                         prop_same_mean = 0,
                         eps = 0.25, c = 0.45) {
  gr <- seq(0, 1, length.out = m)
  
  # Indices of different mean function
  if (setting == 1) {
    # Mean function for group 1
    if (case == 1) {
      mu <- sapply(gr, function(t){ 
        t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2 + cos(2*pi*t + (1:p)/p)
      })
      mu <- t(mu)
    } else if (case == 2) {
      mu <- sapply(gr, function(t){ 
        (t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2 + cos(2*pi*t + (1:p)/p)) * log(p^2/2)
      })
      mu <- t(mu)
    }
    
    num_diff_mean <- ceiling((1-prop_same_mean) * p)
    idx_diff_mean <- sample(1:p, num_diff_mean)
    
  } else if (setting == 2) {
    # Mean function for group 1
    if (case == 1) {
      mu <- sapply(gr, function(t){ 
        t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2
      })
      mu <- t(mu)
    } else if (case == 2) {
      mu <- sapply(gr, function(t){ 
        (t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2) * log(p^2/2)
      })
      mu <- t(mu)
    }
    
    num_diff_mean <- ceiling((1-prop_same_mean) * p)
    idx_diff_mean <- sample(1:p, num_diff_mean)
    
  } else if (setting == 3) {
    if (prop_same_mean == 1) {
      num_diff_mean <- 0
    } else {
      q <- floor(p^c)
      mu <- matrix(0, m, p)
      # mu[, 1:q] <- eps * sqrt(2*log(p) * gr)
      mu[, 1:q] <- eps * sqrt(2*log(p)) * gr
      
      num_diff_mean <- q
      idx_diff_mean <- 1:q
    }
  }
  
  
  # Generate p-dimensional functional data (1 sample)
  generate_mfd_X_i <- function(group = 1, setting = 1, case = 1) {
    # Standard Brownian motion
    Z_i <- sapply(1:p, function(j){ cumsum(rnorm(m, sd = sqrt(1/m))) })
    # Z_i <- sapply(1:m, function(j){ cumsum(rnorm(p)) })
    # dim(Z_i)
    
    if (setting == 1) {
      # # Mean function for group 1
      # if (case == 1) {
      #   mu <- sapply(gr, function(t){ 
      #     t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2 + cos(2*pi*t + (1:p)/p)
      #   })
      #   mu <- t(mu)
      # } else if (case == 2) {
      #   mu <- sapply(gr, function(t){ 
      #     (t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2 + cos(2*pi*t + (1:p)/p)) * log(p^2/2)
      #   })
      #   mu <- t(mu)
      # }
      # 
      # # Randomly selected functional variables with same mean function
      # num_diff_mean <- ceiling((1-prop_same_mean) * p)
      # idx_diff_mean <- sample(1:p, num_diff_mean)
      
      if (case == 1) {
        c_k <- rep(0, p)
        c_k[1] <- 0.5
        c_k[2] <- 0.3
      } else if (case == 2) {
        c_k <- runif(p, 0.1, 0.6)
      }
      phi <- sapply(c_k, function(c_kk){ c_kk * exp(-gr^2/2) / sqrt(0.7468) })
      # dim(phi)
      
      X_i <- matrix(0, m, p)
      for (j in 1:m) {
        for (k in 1:p) {
          # k=1 => only remain \mu
          z_i_sub <- rep(0, p)
          idx_ma <- which(k - (0:(p-1)) > 0)
          z_i_sub[idx_ma] <- Z_i[j, sort(idx_ma, decreasing = T)]
          X_i[j, k] <- sum(phi[j, ] * z_i_sub)
        }
      }
      
      
    } else if (setting == 2) {
      # # Mean function for group 1
      # if (case == 1) {
      #   mu <- sapply(gr, function(t){ 
      #     t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2
      #   })
      #   mu <- t(mu)
      # } else if (case == 2) {
      #   mu <- sapply(gr, function(t){ 
      #     (t*log((1:p)/p + 1) + (sin(2*pi*t + (1:p)/p))^2) * log(p^2/2)
      #   })
      #   mu <- t(mu)
      # }
      # 
      # # Randomly selected functional variables with same mean function
      # num_diff_mean <- ceiling((1-prop_same_mean) * p)
      # idx_diff_mean <- sample(1:p, num_diff_mean)
      
      if (case == 1) {
        c_k <- rep(0, p)
        c_k[1] <- 0.5
        c_k[2] <- 0.3
      } else if (case == 2) {
        c_k <- runif(p, 0.1, 0.6)
      }
      
      generate_phi <- function(u, t) {
        sapply(c_k, function(c_kk){ c_kk * exp(-(u^2 + t^2)/2) / sqrt(0.7468) })
      }
      
      X_i <- matrix(NA, m, p)
      for (j in 1:m) {
        phi <- generate_phi(gr, gr[j])  # m x p
        
        for (k in 1:p) {
          if (k == 1) {
            # Z + mu
            X_i[j, k] <- Z_i[j, k]
          } else {
            # Z + \sum\int~~ + mu
            idx_ma <- which(k - (1:p) > 0)
            X_i[j, k] <- Z_i[j, k] + 
              sum( colMeans(phi[, idx_ma, drop = FALSE] * Z_i[, sort(idx_ma, decreasing = T), drop = FALSE]) )
          }
        }
      }
      
    } else if (setting == 3) {
      # if (prop_same_mean == 1) {
      #   num_diff_mean <- 0
      # } else {
      #   q <- floor(p^c)
      #   mu <- matrix(0, m, p)
      #   # mu[, 1:q] <- eps * sqrt(2*log(p) * gr)
      #   mu[, 1:q] <- eps * sqrt(2*log(p)) * gr
      #   
      #   num_diff_mean <- q
      #   idx_diff_mean <- 1:q
      # }
      
      X_i <- Z_i
    }
    
    # Additional mean function for group 1
    if ((group == 1) & (num_diff_mean > 0)) {
      X_i[, idx_diff_mean] <- X_i[, idx_diff_mean] + mu[, idx_diff_mean]
    }
    
    return(X_i)
  }
  
  X1 <- array(NA, c(n1, m, p))
  for (i in 1:n1) {
    X1[i, , ] <- generate_mfd_X_i(group = 1, setting = setting, case = case)
  }
  X2 <- array(NA, c(n2, m, p))
  for (i in 1:n2) {
    X2[i, , ] <- generate_mfd_X_i(group = 2, setting = setting, case = case)
  }
  
  out <- list(
    X1 = X1,
    X2 = X2
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
eps_param_list <- c(0.15, 0.2, 0.25)
sparsity_list <- c(0, 0.35, 0.45, 0.55)  # c_param list
sim_comb <- expand_grid(samp = samp_list,
                        eps_param = eps_param_list,
                        sparsity = sparsity_list)
# Remove type-I error setting which is duplicated for different delta
sim_comb <- sim_comb %>% 
  filter(!(eps_param != eps_param_list[1] & sparsity == 0))


### Parameters for methods
n_basis_list <- 4:8
lambda_list <- seq(0.001, 0.2, length.out = 10)
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
  eps_param <- sim_comb$eps_param[idx]   # signal strength
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
    
    # Sampling procedure
    if (sparsity == 0) {
      # Type-1 error
      # Generate simulation data
      obj <- generate_mfd(n1, n2, m, p, 
                          prop_same_mean = 1)
      X1 <- obj$X1
      X2 <- obj$X2
    } else if (sparsity > 0) {
      # Power
      # Generate simulation data
      obj <- generate_mfd(n1, n2, m, p, 
                          c = sparsity,
                          eps = eps_param)
      X1 <- obj$X1
      X2 <- obj$X2
    }
    
    
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
       file = paste0("RData/sim_mrp_num_sim=", num_sim, ".RData"))
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

