##################################################
### EEG data
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
source("R/mrp_test.R")
source("R/dcf_test.R")
library(progress)

# Load EEG dataset
X1 <- data.table::fread("../real_data/eeg/alcoholic_data.txt")
X2 <- data.table::fread("../real_data/eeg/control_data.txt")
# dim(X1)
# dim(X2)

# Combine 2 datasets
X <- rbind(X1, X2)
X <- as.matrix(X)

# Transform to 3D array
n <- 122
m <- 256
p <- 64
X <- array(X, c(n, m, p))
dim(X)

# Class labels
y <- c(rep(1, 77), rep(0, 45))



##################################################
### Test for overall dataset
##################################################
n_cores <- 5
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)
seed <- 123

# DSPT-HFD; 5-fold CV (accuracy)
set.seed(seed)
start_time <- Sys.time()
obj_dspt_hfd <- dspt_hfd(X[y == 0, , ], 
                         X[y == 1, , ], 
                         penalty = "scad",
                         tuning = TRUE,
                         tune_method = "cv",
                         n_basis_list = n_basis_list,
                         lambda_list = lambda_list,
                         n_cores = n_cores)
end_time <- Sys.time()
print(end_time - start_time)
obj_dspt_hfd

# DCF test
set.seed(seed)
start_time <- Sys.time()
obj_dcf <- dcf_test(X[y == 0, , ], 
                    X[y == 1, , ], 
                    n_basis_list = n_basis_list) 
end_time <- Sys.time()
print(end_time - start_time)
obj_dcf

# MRP test
set.seed(seed)
start_time <- Sys.time()
obj_mrp <- mrp_test(X[y == 0, , ],
                    X[y == 1, , ],
                    n_cores = n_cores)
end_time <- Sys.time()
print(end_time - start_time)
obj_mrp




### Mean function surfaces
par(mfrow = c(1, 2))
GA::persp3D(1:256, 1:64,
            apply(X[y == 0, , ], c(2, 3), mean),
            # theta = -70, phi = 30, expand = 1,
            xlab = 't', ylab = 's', zlab = '', main = 'Non-alcoholic Group')
GA::persp3D(1:256, 1:64,
            apply(X[y == 1, , ], c(2, 3), mean),
            # theta = -70, phi = 30, expand = 1,
            xlab = 't', ylab = 's', zlab = '', main = 'Alcoholic Group')



##################################################
### EEG data - simulation
##################################################

num_sim <- 400
sparsity_list <- c(0, 1)  # sparsity (0 is size, 1 is power)

### Parameters for methods
n_basis_list <- 4:8
lambda_list <- seq(0.1, 0.5, length.out = 10)
n_cores <- 5

#-----------------------------
# Empirical size and power
# - Size: sparsity = 0
# - Power: sparsity = 1
#-----------------------------
result <- list()   # containing results of size and power
for (idx in 1:length(sparsity_list)) {
  # Sparsity parameter (0 is size, 1 is power)
  sparsity <- sparsity_list[idx]
  
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
      # Randomly split 77 alcoholic samples (X1: 38 samples, X2: 39 samples)
      idx_control <- which(y == 1)
      idx_sample <- sample(idx_control, round(length(idx_control) / 2))
      X1 <- X[idx_sample, , ]
      X2 <- X[setdiff(idx_control, idx_sample), , ]
    } else if (sparsity == 1) {
      # Power
      # Randomly select EEG samples from each group
      idx_control <- which(y == 1)
      idx_sample <- sample(idx_control, round(length(idx_control) / 2))
      X1 <- X[idx_sample, , ]
      X2 <- X[sample(which(y == 0), 30), , ]
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
  
  save(result,
       file = "RData/eeg_result.RData")
}



sapply(result, function(res){
  colMeans(res$p_values < 0.05)
}) %>%
  t() %>% 
  data.frame()

sapply(result, function(res){
  colMeans(res$active_set_len)
}) %>%
  t()

sapply(result, function(res){
  colMeans(res$beta_l1_norm)
}) %>%
  t()


sapply(result, function(res){
  colMeans(res$p_values < 0.05)
}) %>%
  t() %>% 
  xtable::xtable(digits = 3)

