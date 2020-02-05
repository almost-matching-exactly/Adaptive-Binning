#include <Rcpp.h>
#include <algorithm>
#include <chrono>
using namespace Rcpp;

Function expansion_variance("expansion_variance");

double get_greedy_CATE(IntegerVector MG, LogicalVector test_treatments, 
                       NumericVector test_outcomes) {
  // Computes the CATE, given a MG
  double CATE;
  double control_outcomes = 0;
  double treated_outcomes = 0;
  int n_treated_in_MG = 0;
  int n_control_in_MG = 0;
  for (int i = 0; i < MG.size(); i++) {
    if (test_treatments[MG[i]]) {
      treated_outcomes += test_outcomes[MG[i]];
      n_treated_in_MG += 1;
    }
    else {
      control_outcomes += test_outcomes[MG[i]];
      n_control_in_MG += 1;
    }
  }
  CATE = treated_outcomes / n_treated_in_MG - control_outcomes / n_control_in_MG;
  return(CATE);
}

// [[Rcpp::export]]
NumericVector greedy_cpp(NumericMatrix test_treated_covs, IntegerVector test_control, IntegerVector test_treated,
                         NumericMatrix test_covs, LogicalVector test_treatments, NumericVector test_outcomes,
                         int variation, int n_req_matches, SEXP bart_fit) {
  
  auto start = std::chrono::steady_clock::now();
  
  int n_test_treated = test_treated_covs.nrow(); 
  int p = test_covs.ncol();
  int n_test = test_outcomes.size();
  
  NumericVector CATE(n_test_treated);
  
  std::cout << "Running Greedy" << std::endl;
  // For each test-treated unit
  for (int i = 0; i < n_test_treated; i++) {
    std::cout << "Matching unit " << i + 1 << " of " << n_test_treated << "\r" << std::flush;
    
    // Initialize the unit to be trivially in its own MG
    IntegerVector MG(1, test_treated[i]);
    
    // Lower and upper bounds initialized to be unit's covariate values 
    NumericVector A = test_treated_covs(i, _);
    NumericVector B = test_treated_covs(i, _);
    
    // While we haven't matched enough units 
    while (MG.size() < n_req_matches) {
      // std::cout << "The bins are: \n";
      // for (int l = 0; l < A.size(); l++) {
      //   std::cout << A[l] << " " << B[l] << "\n";
      // }
      NumericVector potential_matches;
      // Find units closest along each axis
      if (variation != 2) { // Only expand to control units
        // Don't consider expanding to any units already in the MG
        potential_matches = setdiff(test_control, MG);
      }
      else {
        // To do
        exit(-1);
      }
      NumericVector bin_var(p); // Variance of expansion along each cov
      NumericVector proposed_bin(p); // Proposed bin endpoint for each cov
      for (int j = 0; j < p; j++) { // Find variance of expanding each cov
        NumericVector test_df_col = test_covs(_, j); 
        double dist_to_bin; // Distance from covariate value to current bin
        double min_dist_to_bin = R_PosInf; 
        int closest_unit; // Closest unit to current bin
        double closest_val_to_bin; // Value of j'th cov for closest_unit
        
        ///////// Test for potential_matches being empty /////////
        
        // 1. Find the closest unit we are allowed to match to
        for (int k = 0; k < n_test; k++) {
          bool contains = false;
          for (int m = 0; m < potential_matches.size(); m++) {
            if (potential_matches[m] == k) {
              contains = true;
              break;
            }
          }
          if (!contains) {
            continue;
          }
          // Consider expanding upper bound
          dist_to_bin = test_df_col[k] - B[j];
          if (dist_to_bin < min_dist_to_bin &&  // Closest option so far
              dist_to_bin > 0) {// Located outside the bin 						
            min_dist_to_bin = dist_to_bin;
            closest_unit = k;
            closest_val_to_bin = test_df_col[k];
          }
          
          // Consider expanding lower bound
          dist_to_bin = A[j] - test_df_col[k];
          if (dist_to_bin < min_dist_to_bin &&  // Closest option so far
              dist_to_bin > 0) { // Located outside the bin 
            min_dist_to_bin = dist_to_bin;
            closest_unit = k;
            closest_val_to_bin = test_df_col[k];
          }
        }
        
        if (traits::is_infinite<REALSXP>(min_dist_to_bin)) {
          bin_var[j] = R_PosInf;
          continue; 
        }
        // 2. Change the bin to the closest unit's covariate value 
        proposed_bin[j] = closest_val_to_bin;
        // 3. Test this new bin 
        // + 1 because C++ --> R indexing
        bin_var[j] = Rcpp::as<double>(expansion_variance(j + 1, A, B, proposed_bin[j], bart_fit));
      }
      // Best covariate to expand along
      int expand_along = which_min(bin_var); // what if all equal
      // Update bin 
      if (proposed_bin[expand_along] < A[expand_along]) {// Expanded downwards
        A[expand_along] = proposed_bin[expand_along];
      }
      else {
        B[expand_along] = proposed_bin[expand_along];
      }
      
      // Find units matched, given the unit's new bin
      LogicalVector in_MG(n_test, true);
      
      for (int k = 0; k < n_test; k++) {
        NumericVector test_unit = test_covs(k, _);
        for (int j = 0; j < p; j++) {
          if (test_unit[j] > B[j] || test_unit[j] < A[j]) {
            in_MG[k] = false;
            break;
          }
        }
        if (in_MG[k]) {
          MG.push_back(k);
        }
      }
      
      MG = unique(MG); // Can also get CATE in running fashion
    CATE[i] = get_greedy_CATE(MG, test_treatments, test_outcomes);
  }
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "Time to match " << n_test_treated << " units: "
       << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
       << " seconds" << std::endl;
  return(CATE);
}
  
  