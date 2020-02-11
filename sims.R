 # Libraries ---------------------------------------------------------------
require(beepr)
require(zeallot)
require(tidyverse)
source('estimator_inputs.R')
source('helpers.R')
source('estimators.R')
source('plot_results.R')
 
simulate_data <- function(n_units=100, p=3, n_train=floor(n_units/2), 
                          X_dgp=NULL, e_dgp=NULL, y_dgp=NULL){
   n_test <- n_units - n_train
  # n_estimators <- 7 # Full Matching, Prognostic, CEM, Mahalanobis, Propensity Score, Greedy, MIP

  beta0 <- 2 # Baseline response
  beta_tilde <- 5 # treatment effect
  beta <- runif(p, -1, 1) # To go from X to propensity score 
  
  # Effect of covariates on outcome; not always used 
  gamma <- matrix(runif(p, -3, 3), nrow = p)

  ## For generating propensity scores and assigning treatment
  if (is.null(X_dgp))
    X <- matrix(runif(p * n_units, -5, 5), nrow = n_units)
  else
    X <- X_dgp(n_units, p)
  #X <- matrix(rexp(p * n_units, 0.5), nrow = n_units)
  if(is.null(e_dgp)) {
    e <- expit(.01, X %*% beta)
  }
  else {
    e <- e_dgp(X)
  }
  Z <- rbinom(n_units, 1, e)
  
  ## Generate outcome 
  if(is.null(y_dgp)){
    eps <- rnorm(n_units, 0, 1)
    Y1 <- beta0 + (X[, 1] > 1.5) * beta_tilde + eps
    Y0 <- beta0 + eps
  }else{
    eps <- rnorm(n_units, 0, .1^0.5)
    Y1 <- y_dgp(X, rep(1, n_units), eps)
    Y0 <- y_dgp(X, rep(0, n_units), eps)
  }
  HTE <- Y1 - Y0
  Y = Y1 * Z + Y0 * (1-Z)
  
  df <- cbind(data.frame(X), data.frame(Y = Y, treated = as.logical(Z)))
  return(list(df, HTE))
}

matching_sim <- function(n_sims = 10, n_units = 100, p = 3, n_train = floor(n_units / 2),
                         estimators = c('Full Matching', 'Prognostic', 'CEM',
                                        'Mahalanobis', 'Nearest Neighbor', 'Greedy',
                                        'MIP-Explain', 'MIP-Predict', 'MIQP-Variance',
                                        'MIQP-Fhat', 'MIQP-Grid'),
                         black_box = 'xgb',
                         X_dgp = NULL, e_dgp=NULL, y_dgp = NULL,
                         lambda = 1, alpha=0, beta=2, gamma=1, 
                         lambda0=0, lambda1=0, gamma0=2, gamma1=2, m=1, M=1e5) {
  all_CATEs <- NULL
  all_bins <- vector('list', length = n_sims)
  # all_train_dfs <- vector('list', length = n_sims)
  all_test_dfs <- vector('list', length = n_sims)
  
  for (sim in 1:n_sims) {
    c(df, HTE) %<-% simulate_data(n_units, p, n_train, X_dgp, e_dgp, y_dgp)  
    
    inputs <- estimator_inputs(df, n_train, n_units, black_box)
    hyperparameters <- list(lambda, alpha, beta, gamma, lambda0, lambda1, gamma0, gamma1, m, M)
      
    est_out <- 
      inputs %>%
      get_CATEs(estimators, hyperparameters)
    
    HTE_test_treated = HTE[(n_train + 1):n_units][test_df$treated]
    
    this_sim_CATEs <- 
      est_out %>%
      extract2('CATEs') %>%
      format_CATEs(HTE_test_treated, estimators)
    
    this_sim_bins <- 
      est_out %>%
      extract2('bins')
    
    all_CATEs <- rbind(all_CATEs, this_sim_CATEs)
    all_bins[[sim]] <- this_sim_bins
    # all_train_dfs[[sim]] <- inputs$train_df
    all_test_dfs[[sim]] <- inputs$test_df
    
    print(sprintf('%d of %d simulations completed', sim, n_sims))
  }
  # beep()
  
  return(list(CATEs = all_CATEs, 
              bins = all_bins,
              test_df = all_test_dfs))
}

