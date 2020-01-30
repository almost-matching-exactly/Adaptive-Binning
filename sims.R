# Libraries ---------------------------------------------------------------
require(beepr)
require(zeallot)
require(tidyverse)
source('estimator_inputs.R')
source('helpers.R')
source('estimators.R')
source('plot_results.R')

matching_sim <- function(n_sims = 10, n_units = 100, p = 3, n_train = floor(n_units / 2),
                         estimators = c('Full Matching', 'Prognostic', 'CEM',
                                        'Mahalanobis', 'Nearest Neighbor', 'Greedy',
                                        'MIP-Explain', 'MIP-Predict', 'MIQP-Variance'),
                         X_dgp = NULL, e_dgp=NULL, HTE_dgp = NULL, y_dgp = NULL, ... ) {
  
  n_test <- n_units - n_train
  # n_estimators <- 7 # Full Matching, Prognostic, CEM, Mahalanobis, Propensity Score, Greedy, MIP

  alpha <- 2 # Baseline response
  beta_tilde <- 5 # treatment effect
  beta <- runif(p, -1, 1) # To go from X to propensity score 
  
  # Effect of covariates on outcome; not always used 
  gamma <- matrix(runif(p, -3, 3), nrow = p)
  all_CATEs <- NULL
  all_bins <- vector('list', length = n_sims)
  for (sim in 1:n_sims) {
    ## For generating propensity scores and assigning treatment
    if (is.null(X_dgp))
      X <- matrix(runif(p * n_units, 0, 5), nrow = n_units)
    else
      X <- X_dgp(n, p, ...)
    #X <- matrix(rexp(p * n_units, 0.5), nrow = n_units)
    if(is.null(e_dgp))
      e <- expit(.01, X %*% beta)
    else
      e <- e_dgp(X, ...)
    
    Z <- rbinom(n_units, 1, e)
    
    ## Generate outcome 
    if(is.null(HTE_dgp))
      HTE <- (X[, 1] > 1.5) * beta_tilde * Z
    else
      HTE <- HTE_dgp(X, Z, ...)

    # HTE <- beta_tilde * Z
    HTE_true <- HTE[intersect((n_train + 1):n_units, which(Z == 1))]

    if(is.null(y_dgp))
      Y <- alpha + HTE + rnorm(n_units, 0, 1)
    else
      Y <- y_dgp(X, Z, HTE, ...)

    df <- cbind(data.frame(X), data.frame(Y = Y, treated = as.logical(Z)))
    
    inputs <- estimator_inputs(df, n_train, n_units)
    
    est_out <- 
      inputs %>%
      get_CATEs(estimators)
    
    this_sim_CATEs <- 
      est_out %>%
      extract2('CATEs') %>%
      format_CATEs(HTE_true, estimators)
    
    this_sim_bins <- 
      est_out %>%
      extract2('bins')
    # this_sim_CATEs <-
    #   inputs %>%
    #   get_CATEs(estimators) %>%
    #   format_CATEs(HTE_true, estimators)
    
    all_CATEs <- rbind(all_CATEs, this_sim_CATEs)
    all_bins[[sim]] <- this_sim_bins
    
    print(sprintf('%d of %d simulations completed', sim, n_sims))
  }
  beep()
  
  all_CATEs 
}

