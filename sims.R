
rm(list = ls())

# Libraries ---------------------------------------------------------------
require(beepr)
require(tidyverse)
source('estimator_inputs.R')
source('helpers.R')
source('estimators.R')

matching_sim <- function(n_sims = 10, n_units = 100, p = 3, n_train = floor(n_units / 2)) {
  
  n_test <- n_units - n_train
  # n_estimators <- 7 # Full Matching, Prognostic, CEM, Mahalanobis, Propensity Score, Greedy, MIP

  alpha <- 2 # Baseline response
  beta_tilde <- 5 # treatment effect
  beta <- runif(p, -1, 1) # To go from X to propensity score 
  
  # Effect of covariates on outcome; not always used 
  gamma <- matrix(runif(p, -3, 3), nrow = p)
  store_all_CATEs <- NULL
  for (sim in 1:n_sims) {
    ## For generating propensity scores and assigning treatment
    X <- matrix(runif(p * n_units, 0, 5), nrow = n_units)
    e <- expit(.01, X %*% beta)
    Z <- rbinom(n_units, 1, e)
    
    ## Generate outcome 
    HTE <- (X[, 1] > 1.5) * beta_tilde * Z
    # HTE <- beta_tilde * Z
    HTE_true <- HTE[intersect((n_train + 1):n_units, which(Z == 1))]
    Y <- alpha + HTE + rnorm(n_units, 0, 1)
    
    df <- cbind(data.frame(X), data.frame(Y = Y, treated = as.logical(Z)))
    
    inputs <- estimator_inputs(df, n_train, n_units)

    this_sim_CATEs <- 
      inputs %>%
      get_CATEs() %>%
      format_CATEs(HTE_true)
    
    store_all_CATEs = rbind(store_all_CATEs, this_sim_CATEs)
    
    print(sprintf('%d of %d simulations completed', sim, n_sims))
  }
  beep()
  
  store_all_CATEs %>%
    group_by(estimator) %>%
    summarize(MSE = mean((actual - predicted) ^ 2, na.rm = TRUE),
              percent_missing = 100 * mean(is.na(predicted))) %>%
    arrange(MSE) %>%
    return()
}

# Analysis ----------------------------------------------------------------
res = matching_sim(n_sims = 3, n_units = 100, p = 2)

# unique_HTEs <- unique(HTE$actual)
# if (length(unique_HTEs) < nrow(HTE) / n_estimators) { # Constant treatment effect 
#   gg <- 
#     ggplot(HTE, aes(x = as.factor(actual), y = predicted, color = estimator)) + 
#     geom_boxplot()
#   for (i in 1:length(unique_HTEs)) {
#     gg <- gg + 
#       geom_hline(yintercept = unique_HTEs[i])
#   }
#   gg <- gg + 
#     labs(x = 'Actual', 
#          y = 'Predicted',
#          title = '(Piecewise-)Constant Treatment Effect')
#   print(gg)
# } else {
#   ggplot(HTE, aes(x = actual, y = predicted)) + 
#     geom_point(aes(color = estimator)) + 
#     geom_abline(intercept = 0, slope = 1, color = 'black') + 
#     labs(title = 'Heterogeneous Treatment Effect')
# }
# 
# partitions <- 
#   bins %>%
#   c() %>%
#   unique()
# 
# g <- ggplot(data = df, aes(x = X, y = Y)) + 
#   geom_point(aes(color = treated))
# for (partition in partitions) {
#   g <- g + geom_vline(xintercept = partition)
# }
# plot(g)