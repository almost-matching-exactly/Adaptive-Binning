#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <map>
using namespace Rcpp;
using std::find;

Function expansion_variance("expansion_variance");
Function how_curvy("how_curvy");

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
List greedy_cpp(NumericMatrix test_treated_covs, IntegerVector test_control, IntegerVector test_treated,
                NumericMatrix test_covs, LogicalVector test_treatments, NumericVector test_outcomes,
                int variation, int n_req_matches, double multiplier, SEXP bart_fit) {
  
  auto start = std::chrono::steady_clock::now();
  
  int n_test_treated = test_treated_covs.nrow(); 
  int p = test_covs.ncol();
  int n_test = test_outcomes.size();
  
  NumericVector CATE(n_test_treated);
  
  std::cout << "Running Greedy" << std::endl;
  List all_A = List::create();
  List all_B = List::create();
  
  // NumericVector variances; 
  
  double prev_var; 
  int n_matched_controls;
  IntegerVector all_units = seq(0, n_test - 1);
  // For each test-treated unit
  for (int i = 0; i < n_test_treated; i++) {
    std::cout << "Matching unit " << i + 1 << " of " << n_test_treated << "\r" << std::flush;
    // min_ever_var = 10000.0;
    // Initialize the unit to be trivially in its own MG
    IntegerVector MG(1, test_treated[i]);
    
    n_matched_controls = 0; 
    
    // Lower and upper bounds initialized to be unit's covariate values 
    NumericVector A = test_treated_covs(i, _);
    NumericVector B = test_treated_covs(i, _);
    NumericVector bin_var(p, 10000.0); // Variance of expansion along each cov
    // While we haven't matched enough units 
    do {
      prev_var = min(bin_var) + 0.01;
      // if (min(bin_var) < min_ever_var) {
      //   min_ever_var = min(bin_var) + 0.02;
      // }

      NumericVector potential_matches;
      // Find units closest along each axis
      if (variation != 2) { 
        // Don't consider expanding to any units already in the MG
        // Don't take union every time 
        potential_matches = setdiff(all_units, MG);
      }
      else {
        // To do
        exit(-1);
      }
      
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
        // bin_var[j] = Rcpp::as<double>(how_curvy(j + 1, A, B, proposed_bin[j], bart_fit));
        bin_var[j] = Rcpp::as<double>(expansion_variance(j + 1, A, B, proposed_bin[j], bart_fit));
      }
      // Best covariate to expand along
      int expand_along = which_min(bin_var); // what if all equal
      // Update bin 
      // std::cout << "i: bin_var prev_var relerror " << i << " " << min(bin_var) << " " << prev_var << " " <<
      //   abs(min(bin_var) - prev_var) / prev_var << std::endl;
      // if (abs(min(bin_var) - prev_var) / prev_var < multiplier) {
      if (min(bin_var) < multiplier * prev_var || n_matched_controls < 1) {
        // std::cout << "min_ever_var: " << min_ever_var << " and minbinvar: " << min(bin_var) << std::endl;
        // std::cout << i << "Got all up in here" << std::endl;
        if (proposed_bin[expand_along] < A[expand_along]) {// Expanded downwards
          A[expand_along] = proposed_bin[expand_along];
        }
        else {
          B[expand_along] = proposed_bin[expand_along];
        }
      }
      else {
        break;
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
          if (std::binary_search(test_control.begin(), test_control.end(), k)) {
            n_matched_controls += 1; 
          }
          MG.push_back(k);
        }
      }
      
      MG = unique(MG); // Can also get CATE in running fashion
      CATE[i] = get_greedy_CATE(MG, test_treatments, test_outcomes);
      // std::cout << "Unit " << i << "has min_var = " << min(bin_var) << std::endl;
    }
    // while (abs(min(bin_var) - prev_var) / prev_var < multiplier);
    while (min(bin_var) < multiplier * prev_var || n_matched_controls < 1);
    all_A.push_back(A);
    all_B.push_back(B);
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "Time to match " << n_test_treated << " units: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
            << " seconds" << std::endl;
  List ret = List::create(CATE, all_A, all_B);
  return(ret);
}

