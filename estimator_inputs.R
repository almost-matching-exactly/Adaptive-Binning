estimator_inputs <- function(df, n_train, n_units, black_box) {
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
    
  if (black_box == 'BART') {
    black_box_fit <- bart(x.train = dplyr::select(train_df, -Y),
                     y.train = train_df$Y,
                     x.test = mutate(test_df[test_df$treated, 1:p], treated = 0), # Prognostic score on test units
                     keeptrees = TRUE,
                     verbose = FALSE)
    counterfactuals <- black_box_fit$yhat.test.mean
  }
  else if (black_box == 'xgb') {
    eta <- c(.01, .05, .1, .2, .3, .5)
    max_depth <- c(2, 3, 4, 6, 8)
    alpha <- c(.01, .1, .5, 1, 5)
    nrounds <- c(50, 100, 200)
    subsample <- c(0.1, 0.3, 0.5, 0.75, 1)
    param_combs <- 
      expand.grid(eta, max_depth, alpha, nrounds, subsample) %>%
      as.data.frame() %>%
      `colnames<-`(c('eta', 'max_depth', 'alpha', 'nrounds', 'subsample'))
    
    tmp <- as.matrix(dplyr::select(train_df, -Y))
    
    RMSE <- vector(mode = 'numeric', length = length(param_combs))
    for (i in 1:length(param_combs)) {
      params <- list(objective = 'reg:squarederror', 
                     eta = param_combs$eta[i],
                     max_depth = param_combs$max_depth[i],
                     alpha = param_combs$alpha[i],
                     subsample = param_combs$subsample[i])
      cv <- xgb.cv(data = tmp, 
                   label = train_df$Y,
                   params = params, metrics = list('rmse'), 
                   nrounds = param_combs$nrounds[i], nfold = 5, verbose = 0)
      RMSE[i] <- cv$evaluation_log$test_rmse_mean[param_combs$nrounds[i]]
    }
    
    best_params <- param_combs[which.min(RMSE), ]
    
    black_box_fit <- xgboost(data = tmp,
                        label = train_df$Y,
                        params = list(objective = 'reg:squarederror', 
                                      eta = best_params$eta,
                                      max_depth = best_params$max_depth,
                                      alpha = best_params$alpha,
                                      subsample = best_params$subsample),
                        nround = best_params$nrounds,
                        verbose = 0)
    
    counterfactuals <-
      black_box_fit %>%
      predict(newdata = as.matrix(mutate(test_df[test_df$treated, 1:p], treated = 0)))  
  }
  else if (black_box == 'LASSO') {
    
  }
  else {
    stop("black_box must be one of: BART, xgb, LASSO")
  }
  
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
              bart_fit = black_box_fit,
              counterfactuals = counterfactuals))
}