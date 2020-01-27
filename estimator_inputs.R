estimator_inputs <- function(df, n_train, n_units) {
  p <- ncol(df) - 2
  f <- formula(paste('treated ~', paste(colnames(df)[1:p], collapse = ' + '))) 
  
  train_df <-
    df %>%
    dplyr::slice(1:n_train)
  
  train_covs <- 
    train_df %>%
    dplyr::select(1:p)
  
  train_control <- which(!train_df$treated)
  train_treated <- which(train_df$treated)
  
  test_df <-
    df %>%
    dplyr::slice((n_train + 1):n_units)
  
  test_covs <- 
    test_df %>%
    dplyr::select(1:p)
  
  test_control <- which(!test_df$treated)
  test_treated <- which(test_df$treated)
  
  n_test_control <- length(test_control)
  n_test_treated <- length(test_treated)
  
  # bart_fit <- bart(x.train = dplyr::select(train_df, -Y),
  #                  y.train = train_df$Y,
  #                  x.test = mutate(test_df[test_df$treated, 1:p], treated = 0), # Prognostic score on test units 
  #                  keeptrees = TRUE,
  #                  verbose = FALSE)
  # counterfactuals <- bart_fit$yhat.test.mean
  
  # Not actually a bart fit
  bart_fit <- xgboost(data = as.matrix(dplyr::select(train_df, -Y)),
                      label = train_df$Y,
                      nround = 50,
                      verbose = 0)
  
  counterfactuals <- 
    bart_fit %>%
    predict(newdata = as.matrix(mutate(test_df[test_df$treated, 1:p], treated = 0)))
  
  # browser()
  
  
  return(list(df = df,
              f = f,
              n = n_units,
              n_train = n_train,
              p = p,
              train_df = train_df, 
              train_covs = train_covs, 
              train_control = train_control,
              train_treated = train_treated, 
              test_df = test_df, 
              test_covs = test_covs, 
              test_control = test_control,
              test_treated = test_treated,
              n_test_control = n_test_control,
              n_test_treated = n_test_treated,
              bart_fit = bart_fit,
              counterfactuals = counterfactuals))
}