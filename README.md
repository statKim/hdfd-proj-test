# Two-Sample Projection Test for High-Dimensional Functional Data

This repository contains the R code to reproduce the simulations and real data analysis presented in the following paper:

> Hyunsung Kim and Junyong Park (2026+). Two-Sample Projection Test for High-Dimensional Functional Data, *submitted*.

## Required R Package

Before you run these reproducible codes in the repository, you need a R package `hdfda` which is in [Here](https://github.com/statKim/hdfda), and it can be installed by the following code:

``` r
# install.packages("devtools")
devtools::install_github("statKim/hdfda")
```

## Description

- **R/dcf_test.R**: R functions to implement *Distribution/correlation-free test*.
- **R/mrp_test.R**: R functions to implement *Multi-resolution projection test*.
- **src/mrp_test_cpp.cpp**: cpp functions to use `mrp_test.R` for fast implementation with parallel computing.
- **sim.R**: R code of the Simulation 1 in the manuscript.
- **sim_dcf.R**: R code of the Simulation 2 in the manuscript.
- **sim_mrp.R**: R code of the Simulation 3 in the manuscript.
- **eeg.R**: R code of the real data analysis (EEG data)
- **adhd_200.R**: R code of the real data analysis (ADHD-200 data)

## Data source

- **EEG**: Download from [UCI Machine Learning Repository](https://archive.ics.uci.edu/dataset/121/eeg+database).
- **ADHD-200**
  - Download from NITRC website: https://www.nitrc.org/frs/?group_id=383
  - Descriptions of preprocessing: https://www.nitrc.org/plugins/mwiki/index.php?title=neurobureau:AthenaPipeline
  - We downloaded the following preprocessed data.
    - Peking_1_preproc_filtfix.tar
    - Peking_2_preproc_filtfix.tar
    - Peking_3_preproc_filtfix.tar
