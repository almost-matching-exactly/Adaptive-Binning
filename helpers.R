expit <- function(a, x) {
  exp(a * x) / (1 + exp(a * x))
}

expansion_variance <- function(cov, current_bin, expanded_bin, df, bart_fit, n_grid_pts = 8) {
  p <- ncol(df) - 2
  
  expanded <- which(current_bin[cov, ] != expanded_bin[cov, ])
  grid_pts <- seq(current_bin[cov, expanded], expanded_bin[cov, expanded], 
                  length.out = n_grid_pts)
  
  bin_centers <- rowMeans(current_bin)
  
  pred_data <- 
    sapply(grid_pts, function(x) {
      bin_centers[cov] <- x
      bin_centers
    }) %>%
    t() %>% 
    cbind(1) # Treatment = TRUE
  
  return(var(predict(bart_fit, pred_data))) ## For XGBoost
  return(var(colMeans(predict(bart_fit, pred_data)))) ## For BART
}

get_greedy_CATE <- function(n_test_treated, test_covs, bins, test_df) {
  CATE <- vector('numeric', n_test_treated)
  for (i in 1:n_test_treated) {
    in_MG <- which(apply(test_covs, 1, function(x) all(x >= bins[i, , 1]) & all(x <= bins[i, , 2])))
    treated <- in_MG[which(test_df$treated[in_MG])]
    control <- in_MG[which(!test_df$treated[in_MG])]
    CATE[i] <- mean(test_df$Y[treated]) - mean(test_df$Y[control])
  }
  return(CATE)
}

summarize_CATEs <- function(all_CATEs) {
  all_CATEs %>%
    group_by(estimator) %>%
    summarize(MAE = mean(abs(actual - predicted), na.rm = TRUE),
              percent_missing = 100 * mean(is.na(predicted))) %>%
    arrange(MAE) %>%
    return()
}

format_CATEs <- function(CATE_obj, true_CATE, estimators) {
  n_estimators <- ncol(CATE_obj)
  CATE_out <- NULL
  for (i in 1:n_estimators) {
    CATE_out %<>%
      rbind(cbind(true_CATE, CATE_obj[, i]))
  }
  CATE_out %<>%
  as.data.frame() %>%
    `colnames<-`(c('actual', 'predicted')) %>%
    mutate(estimator = rep(estimators, each = nrow(.) / n_estimators))
  return(CATE_out)
}

make_mg <- function(X, lbs, ubs){
  return(which(colMeans((t(X) <= ubs)*(t(X) >= lbs))==1))
}
