# Libraries ---------------------------------------------------------------
require(dbarts)
require(MatchIt)
require(beepr)
require(cem)
require(tidyverse)
source("matching_estimators.R")
source("AB_MIP.R")

# Functions ---------------------------------------------------------------
expit <- function(a, x) {
  exp(a * x) / (1 + exp(a * x))
}

expansion_variance <- function(cov, current_bin, expanded_bin, df, bart_fit) {
  p <- ncol(df) - 2
  
  expanded <- which(current_bin[cov, ] != expanded_bin[cov, ])
  grid_pts <- seq(current_bin[cov, expanded], expanded_bin[cov, expanded], 
                  length.out = 8)
  bin_centers <- rowMeans(current_bin)
  
  pred_data <- 
    sapply(grid_pts, function(x) {
      bin_centers[cov] <- x
      bin_centers
    }) %>%
    t() %>% 
    as.data.frame() %>%
    mutate(treatment = 1) %>%
    `colnames<-`(colnames(dplyr::select(df, -Y)))
  
  return(var(colMeans(predict(bart_fit, pred_data))))
}
# Algorithm ---------------------------------------------------------------


matching_sim <- function(n_sims = 10, n_units = 100, p = 3, n_train = floor(n_units / 2)) {
  
  n_test <- n_units - n_train
  n_estimators <- 7 # CEM, Dynamic Binning, Full Matching, 
                    # Mahalanobis, Nearest Neighbor, Prognostic, MIP DB

  alpha <- 2 # Baseline response
  beta_tilde <- 5 # treatment effect
  beta <- runif(p, -1, 1) # To go from X to propensity score 
  
  # Effect of covariates on outcome; not always used 
  gamma <- matrix(runif(p, -3, 3), nrow = p)
  
  for (sim in 1:n_sims) {
    ## For generating propensity scores and assigning treatment
    X <- matrix(runif(p * n_units, 0, 5), nrow = n_units)
    e <- expit(.01, X %*% beta)
    Z <- rbinom(n_units, 1, e)
    
    ## Generate outcome 
    HTE <- (X[, 1] > 1.5) * beta_tilde * Z
    HTE_true <- HTE[intersect((n_train + 1):n_units, which(Z == 1))]
    Y <- alpha + HTE + rnorm(n_units, 0, 1)
    
    df <- cbind(data.frame(X), data.frame(Y = Y, treated = as.logical(Z)))
    
    # Formula for matching methods
    f <- formula(paste('treated ~', paste(colnames(df)[1:p], collapse = ' + '))) 
    
    train_df <-
      df %>%
      slice(1:n_train)
    
    train_covs <- 
      train_df %>%
      dplyr::select(1:p)
    
    train_control <- which(!train_df$treated)
    train_treated <- which(train_df$treated)
    
    test_df <-
      df %>%
      slice((n_train + 1):n_units)
    
    test_covs <- 
      test_df %>%
      dplyr::select(1:p)
    
    test_control <- which(!test_df$treated)
    test_treated <- which(test_df$treated)
    
    n_test_control <- length(test_control)
    n_test_treated <- length(test_treated)
    
    # sorted_test_control_covs <- 
    #   test_covs %>% 
    #   slice(test_control) %>% 
    #   lapply(sort, index.return = TRUE)
    
    store_all_HTEs <- NULL
    
    bart_fit <- bart(x.train = dplyr::select(train_df, -Y),
                     y.train = train_df$Y,
                     x.test = mutate(test_df[test_df$treated, 1:p], treated = 0), # Prognostic score on test units 
                     keeptrees = TRUE,
                     verbose = FALSE)
    
    message("Running alternative estimators...")
    ################### Alternative Matching Estimators ################### 
    HTE_fullmatch <- est_fullmatch(test_df, f)
    HTE_prog <- est_prog(test_df, bart_fit$yhat.test.mean)
    HTE_cem <- est_cem(test_df)
    HTE_mahal <- est_mahal(f, test_df)
    HTE_nn <- est_nn(f, test_df, ratio = 1)
    # HTE_caliper <- est_caliper(f, df, ratio = 1) # Doesn't work...? 
    ######################################################################
    
    
    # Initialize bins to be centered right on treated unit values 
    bins <- array(dim = c(n_test_treated, p, 2))
    bins[, , 1] <- as.matrix(test_df[test_treated, 1:p])
    bins[, , 2] <- as.matrix(test_df[test_treated, 1:p])
    
    # Indices of test units already matched to, for each test, treated unit
    # Every unit will trivially be in their own MG 
    already_matched <- as.list(test_treated)
    
    message("Running Greedy AB...")
    for (i in 1:n_test_treated) {
      message(paste("Matching unit", i, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
      bin_copy <- matrix(bins[i, ,], ncol = 2) # For case of 1 covariate
      counter <- 0
      while (length(already_matched[[i]]) < 5) { # Variance is not too big
        min_var <- Inf
        min_size_increase <- Inf
        
        # Find units closest along each axis
        potential_matches <- setdiff(test_control, already_matched[[i]])
        bin_var <- vector('numeric', length = p)
        proposed_bin <- vector('list', length = p)
        for (j in 1:p) {
          expand_up <- intersect(potential_matches, which(test_df[, j] > bin_copy[j, 2]))
          if (length(expand_up) == 0) {
            distance_up <- Inf
          }
          else {
            distance_up <- min(abs(bin_copy[j, 2] - test_df[expand_up, j])) ###### Consider sorting covariates so we don't have to do this   
          }
          
          expand_down <- intersect(potential_matches, which(test_df[, j] < bin_copy[j, 1]))
          
          if (length(expand_down) == 0) { # Can't go any lower
            distance_down <- Inf
          }
          else {
            distance_down <- min(abs(bin_copy[j, 1] - test_df[expand_down, j]))  
          }
          
          if (is.infinite(distance_up) & is.infinite(distance_down)) { 
            bin_var[j] <- Inf
            next
          }
          
          if (distance_up < distance_down) {
            proposed_bin[[j]] <- bin_copy
            proposed_bin[[j]][j, 2] <- test_df[expand_up[which.min(abs(bin_copy[j, 2] - test_df[expand_up, j]))], j]
          }
          else {
            proposed_bin[[j]] <- bin_copy
            proposed_bin[[j]][j, 1] <- test_df[expand_down[which.min(abs(bin_copy[j, 1] - test_df[expand_down, j]))], j]
          }
          
          bin_var[j] <- expansion_variance(j, bin_copy, proposed_bin[[j]], train_df, bart_fit)
        }
        
        expand_along <- which.min(bin_var)
        bin_copy <- proposed_bin[[expand_along]]
        in_MG <- apply(test_covs, 1, function(x) all(x >= bin_copy[, 1]) & all(x <= bin_copy[, 2]))
        already_matched[[i]] <- unique(c(already_matched[[i]], which(in_MG)))
        
        counter <- counter + 1
      }
      
      bins[i, , ] <- bin_copy
    }
    message("\n")
    
    CATE <- vector('numeric', n_test_treated)
    size <- vector('numeric', n_test_treated)
    for (i in 1:n_test_treated) {
      in_MG <- which(apply(test_covs, 1, function(x) all(x >= bins[i, , 1]) & all(x <= bins[i, , 2])))
      size[i] <- length(in_MG)
      treated <- in_MG[which(test_df$treated[in_MG])]
      control <- in_MG[which(!test_df$treated[in_MG])]
      CATE[i] <- mean(test_df$Y[treated]) - mean(test_df$Y[control])
    }
    ATE <- sum(CATE * size) / sum(size)
    
    ## MIP DB
    MIP_cates = vector('numeric', n_test_treated)
    message("Running MIP AB...")
    for (l in 1:n_test_treated){
      i = test_treated[l]
      message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
      
      mip_pars =  create_unit_mip(xi = as.numeric(test_covs[i, ]), zi =  1, y_train = train_df$Y, 
                                  x_train = as.matrix(train_covs), z_train = train_df$treated, 
                                  x_test = as.matrix(test_covs), z_test = test_df$treated, 
                                  alpha=5, lambda=20, m=1, M=1e5)
      sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
      mip_out = recover_pars(sol, n_train, n_test, p)
    
      MIP_cates[l] = test_df$Y[i] - mean(test_df$Y[mip_out$w>=0.1 & test_df$treated==0])
    }
    message("\n")
    this_sim <- 
      rbind(cbind(HTE_true, CATE), # dynamic binning
            cbind(HTE_true, MIP_cates),
            cbind(HTE_true, HTE_fullmatch), 
            cbind(HTE_true, HTE_prog),
            cbind(HTE_true, HTE_cem), 
            cbind(HTE_true, HTE_mahal), 
            cbind(HTE_true, HTE_nn)) %>%
      as.data.frame() %>%
      `colnames<-`(c('actual', 'predicted')) %>%
      mutate(estimator = rep(c('Dynamic Binning','MIP DB',
                               'Full Matching', 'Prognostic',
                               'CEM', 'Mahalanobis', 'Nearest Neighbor'),
                             each = nrow(.) / n_estimators))
    store_all_HTEs = rbind(store_all_HTEs, this_sim)
    
    print(sprintf('%d of %d simulations completed', sim, n_sims))
  }
  beep()
  
  store_all_HTEs %>%
    group_by(estimator) %>%
    summarize(MSE = mean((actual - predicted) ^ 2, na.rm = TRUE),
              percent_missing = 100 * mean(is.na(predicted))) %>%
    arrange(MSE) %>%
    return()
}

# Analysis ----------------------------------------------------------------

res = matching_sim(n_sims = 10, n_units = 100, p = 5)

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