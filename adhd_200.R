##################################################
### ADHD-200 data
##################################################
library(tidyverse)
library(hdfda)
library(doParallel)
library(progress)
source("R/mrp_test.R")
source("R/dcf_test.R")


# Load preprocessed ADHD-200 (PekingU)
# - 114 control, 77 ADHD
# - X: n-m-p array, 
# - y: class label (0 or 1)
load("../real_data/adhd_200/preprocessed_PU/ADHD-200_Schaefer17_400regions.RData")
dim(X)
table(y)

n <- dim(X)[1]
m <- dim(X)[2]
p <- dim(X)[3]


# Fast Fourier Transform with smoothing splines
# Guo, X., Li, Y., & Hsing, T. (2023). Variable Selection and Minimax Prediction in High-dimensional Functional Linear Models. arXiv preprint arXiv:2310.14419.
X_fft_sm <- X
for (i in 1:p) {
  if (i %% 20 == 0) {
    cat(paste(i, " > "))
  }
  
  X_i_fft_sm <- apply(X[, , i], 1, function(x) {
    smooth.spline(Mod(fft(x)) * (2/m))$y
  })
  X_fft_sm[, , i] <- t(X_i_fft_sm)
}



##################################################
### Test for overall dataset
##################################################
n_cores <- 5
n_basis_list <- 4:8
lambda_list <- 10^seq(-3, -1/2, length.out = 10)
seed <- 123

#------------------------
# Raw data
#------------------------

# DSPT-HFD; 5-fold CV
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


#------------------------
# FFT-smooth
#------------------------

# DSPT-HFD; 5-fold CV
set.seed(seed)
start_time <- Sys.time()
obj_dspt_hfd <- dspt_hfd(X_fft_sm[y == 0, , ], 
                         X_fft_sm[y == 1, , ], 
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
obj_dcf <- dcf_test(X_fft_sm[y == 0, , ], 
                    X_fft_sm[y == 1, , ], 
                    n_basis_list = n_basis_list) 
end_time <- Sys.time()
print(end_time - start_time)
obj_dcf

# MRP test
set.seed(seed)
start_time <- Sys.time()
obj_mrp <- mrp_test(X_fft_sm[y == 0, , ],
                    X_fft_sm[y == 1, , ],
                    n_cores = n_cores)
end_time <- Sys.time()
print(end_time - start_time)
obj_mrp

