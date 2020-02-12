estimator_inputs <- function(df, n_train, n_units, black_box, cv) {
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
    # chisq <- function(df, quant) {
    #   new('dbartsChiSqPrior', df = df, quantile = quant)
    # }
    k <- 2
    n.trees <- 200
    if (cv) {
      message('Cross-validating BART; should take a minute')
      flush.console()
      n.trees <- c(2, 5, 10, 20, 50, 200)
      nu_q <- list(c(3, 0.9), c(3, 0.99), c(10, 0.75))
      k <- c(1, 2, 3, 5, 7)
      alpha <- c(0.5, 0.95)[2]
      beta <- c(.5, 2)[2]
      
      cv <- xbart(formula = dplyr::select(train_df, -Y),
                  data = train_df$Y,
                  verbose = FALSE,
                  method = 'k-fold',
                  n.test = 10,
                  n.reps = 10,
                  loss = 'rmse',
                  n.trees = n.trees,
                  k = k,
                  power = beta,
                  base = alpha) %>% 
        apply(c(2, 3), mean)
      
      best_params <- arrayInd(which.min(cv), dim(cv))
      n.trees <- n.trees[best_params[1]]
      k <- k[best_params[2]]
    }
    
    black_box_fit <- bart(x.train = dplyr::select(train_df, -Y),
                     y.train = train_df$Y,
                     x.test = mutate(test_df[test_df$treated, 1:p], treated = 0), # Prognostic score on test units
                     keeptrees = TRUE,
                     keepevery = 10,
                     verbose = FALSE,
                     k = k,
                     ntree = n.trees)
    
    # pred_points <- 
    #   expand.grid(seq(-5, 5, length.out = 50), 
    #               seq(-5, 5, length.out = 50)) %>% 
    #   `colnames<-`(c('X1', 'X2')) %>% 
    #   mutate(treated = 1)
    # 
    # preds <- colMeans(predict(black_box_fit, pred_points))
    # pred_points %<>% mutate(preds = preds)
    # # browser()
    # ggplot(data = pred_points, aes(x = X1, y = X2, fill = preds)) +
    #   geom_raster() +
    #   geom_rect(xmin = -1, xmax = 2, ymin = 2, ymax= 4, color = 'red', fill = NA)
    counterfactuals <- black_box_fit$yhat.test.mean
  }
  else if (black_box == 'xgb') {
    eta <- c(.01, .05, .1, .2, .3, .5)
    max_depth <- c(2, 3, 4, 6, 8)
    alpha <- c(.01, .1, .5, 1, 5)
    nrounds <- c(5, 10, 50, 100, 200)
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