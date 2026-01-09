#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// Helper function to compute trace(A^T * V * B * B^T * V * A)
// This corresponds to the trace term in the variance calculation:
// tr{ A * B^T * C * D^T } = tr{ (B^T * V * C) * (D^T * V * A) }
// where the matrices passed are already the raw functional data slices.
//
// The integral with v(s,t) corresponds to the quadratic form X^T * V * Y
double compute_trace_integral(const mat& A, const mat& B, const mat& C, const mat& D, const mat& V) {
  // Term 1: integral over (s,t) pairing A and C -> (A^T * V * C)
  mat term1 = A.t() * V * C;
  
  // Term 2: integral over (s1, t1) pairing B and D -> (B^T * V * D)
  mat term2 = B.t() * V * D;
  
  // We need tr(Term1 * Term2)
  // For two square matrices X and Y, tr(X*Y) = sum(sum(X % Y.t()))
  // &: element-wise multiplication
  return accu(term1 % term2.t());
}

// Helper to compute simple quadratic form integral for MRP statistic
// integral X(s)^T Y(t) v(s,t) ds dt  =>  tr(X^T * V * Y)
double compute_mrp_term(const mat& X, const mat& Y, const mat& V) {
  return trace(X.t() * V * Y);
}

// Main function to compute MRP statistics and Variance components
// X and Y are expected to be (m x p x n) cubes (permuted in R)
// X_c and Y_c are centered versions of X and Y
// V is the (m x m) covariance matrix
// [[Rcpp::export]]
Rcpp::List mrp_cpp_worker(const arma::cube& X, const arma::cube& Y, 
                          const arma::cube& X_c, const arma::cube& Y_c,
                          const arma::mat& V, int n_cores) {
  
  omp_set_num_threads(n_cores);
  
  int n = X.n_slices;
  int m_samples = Y.n_slices; // 'm' in paper, but avoiding confusion with timepoints m
  
  // --- 1. Compute MRP Statistic Terms ---
  
  double term1_sum = 0.0;
  double term2_sum = 0.0;
  double term3_sum = 0.0;
  
  // Term 1: Sum over X_i, X_j (i != j)
#pragma omp parallel for reduction(+:term1_sum) collapse(2)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      term1_sum += compute_mrp_term(X.slice(i), X.slice(j), V);
    }
  }
  
  // Term 2: Sum over Y_i, Y_j (i != j)
#pragma omp parallel for reduction(+:term2_sum) collapse(2)
  for (int i = 0; i < m_samples; ++i) {
    for (int j = 0; j < m_samples; ++j) {
      if (i == j) continue;
      term2_sum += compute_mrp_term(Y.slice(i), Y.slice(j), V);
    }
  }
  
  // Term 3: Sum over X_i, Y_j
#pragma omp parallel for reduction(+:term3_sum) collapse(2)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m_samples; ++j) {
      term3_sum += compute_mrp_term(X.slice(i), Y.slice(j), V);
    }
  }
  
  double mrp_stat = term1_sum / (n * (n - 1)) + 
    term2_sum / (m_samples * (m_samples - 1)) - 
    2.0 * term3_sum / (n * m_samples);
  
  // --- 2. Compute Variance Components (Sigma^2) ---
  
  // We need ITR(G1, G1), ITR(G2, G2), and ITR(G1, G2)
  // Using the estimator from the paper involving leave-two-out means.
  // To compute efficiently, we use the pre-calculated centered data (X_c, Y_c)
  // Relations:
  // X_j - X_bar_(j,k) = c1 * Z_j + c2 * Z_k  (where Z is globally centered X)
  // c1 = (n-1)/(n-2), c2 = 1/(n-2)
  
  double tr_11_sum = 0.0;
  double tr_22_sum = 0.0;
  double tr_12_sum = 0.0;
  
  double c1_n = (double)(n - 1) / (n - 2);
  double c2_n = 1.0 / (n - 2);
  
  double c1_m = (double)(m_samples - 1) / (m_samples - 2);
  double c2_m = 1.0 / (m_samples - 2);
  
  // ITR_11: Sum over j != k in Group 1
#pragma omp parallel for reduction(+:tr_11_sum) collapse(2)
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      if (j == k) continue;
      
      // Construct leave-two-out centered vectors using linear combo of global centered
      // A = X_j - mean_(j,k)
      // C = X_k - mean_(j,k)
      // Note: We operate directly on matrices to avoid allocation if possible, 
      // but here temporary mats are needed for the linear combination.
      
      mat A = c1_n * X_c.slice(j) + c2_n * X_c.slice(k);
      mat C = c1_n * X_c.slice(k) + c2_n * X_c.slice(j);
      
      // Compute trace term: tr(A B^T C D^T) with pairing
      // B = X_j (raw), D = X_k (raw)
      tr_11_sum += compute_trace_integral(A, X.slice(j), C, X.slice(k), V);
    }
  }
  
  // ITR_22: Sum over j != k in Group 2
#pragma omp parallel for reduction(+:tr_22_sum) collapse(2)
  for (int j = 0; j < m_samples; ++j) {
    for (int k = 0; k < m_samples; ++k) {
      if (j == k) continue;
      
      mat A = c1_m * Y_c.slice(j) + c2_m * Y_c.slice(k);
      mat C = c1_m * Y_c.slice(k) + c2_m * Y_c.slice(j);
      
      tr_22_sum += compute_trace_integral(A, Y.slice(j), C, Y.slice(k), V);
    }
  }
  
  // ITR_12: Sum over j in G1, k in G2
  // For cross terms, the mean is leave-one-out
  // X_j - X_bar_(j) = (n / (n-1)) * Z_j
  // Y_k - Y_bar_(k) = (m / (m-1)) * Z_k
  double scale_n = (double)n / (n - 1);
  double scale_m = (double)m_samples / (m_samples - 1);
  
#pragma omp parallel for reduction(+:tr_12_sum) collapse(2)
  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < m_samples; ++k) {
      // A = X_j - X_bar_(j)
      // C = Y_k - Y_bar_(k)
      // B = X_j, D = Y_k
      
      // Use pointers to slices to avoid copy where possible, 
      // but scaling needs temp or BLAS.
      mat A = scale_n * X_c.slice(j);
      mat C = scale_m * Y_c.slice(k);
      
      tr_12_sum += compute_trace_integral(A, X.slice(j), C, Y.slice(k), V);
    }
  }
  
  double itr_11 = tr_11_sum / (n * (n - 1));
  double itr_22 = tr_22_sum / (m_samples * (m_samples - 1));
  double itr_12 = tr_12_sum / (n * m_samples);
  
  // Combine for Sigma^2 (Eq 7)
  double sigma_sq = (2.0 / (n * (n - 1))) * itr_11 + 
    (2.0 / (m_samples * (m_samples - 1))) * itr_22 + 
    (4.0 / (n * m_samples)) * itr_12;
  
  return Rcpp::List::create(
    Rcpp::Named("mrp_stat") = mrp_stat,
    Rcpp::Named("sigma_sq") = sigma_sq
  );
}