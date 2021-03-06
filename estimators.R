source('helpers.R')
source("AB_MIPs.R")
Rcpp::sourceCpp('greedy2.cpp')

require(MatchIt)
require(cem)
require(dbarts)
require(magrittr)
require(xgboost)
require(Matching)

est_true <- function(df, Y0true){
  return(est_prog(df, Y0true))
}

est_fullmatch <- function(df, n, n_train, f) {
  m.out <- suppressWarnings(matchit(f, data = df, method = 'full'))
  HTE <- vector('numeric', length = n)
  for (i in 1:n) {
    in_MG <- which(m.out$subclass == m.out$subclass[i])
    treated <- intersect(in_MG, which(df$treated))
    control <- intersect(in_MG, which(!df$treated))
    HTE[i] <- mean(df$Y[treated]) - mean(df$Y[control])
  }
  HTE <- HTE[intersect((n_train + 1):n, which(df$treated))]
  return(HTE)
}

est_gen_match <- function(df, n_train, f, ratio) {
  capture.output(suppressWarnings(gmatch_out <- 
                   Matching::GenMatch(df$treated, df[, ncol(df) - 2], M = ratio)), file = 'NUL')
  
  n <- nrow(df)
  matches <- gmatch_out$matches
  test_treated_inds <- intersect(which(df$treated), (n_train + 1):n)
  HTE <- vector('numeric', length = length(test_treated_inds))
  for (i in 1:length(test_treated_inds)) {
    matching_units <- matches[which(matches[, 1] == test_treated_inds[i]), 2]
    HTE[i] <- df$Y[test_treated_inds[i]] - mean(df$Y[matching_units])
  }
  return(HTE)
}

est_prog <- function(df, counterfactuals, ratio, replace = TRUE) {
  treated_outcomes <- df$Y[df$treated]
  control_outcomes <- df$Y[!df$treated]
  HTE <- rep(NA, length(treated_outcomes))
  for (i in 1:length(treated_outcomes)) {
    if (length(control_outcomes) == 0) {
      break
    }
    match_outcomes <- vector('numeric', length = ratio)
    tmp <- control_outcomes
    for (j in 1:ratio) {
      closest <- which.min(abs(counterfactuals[i] - tmp))
      # matches[j] <- closest
      match_outcomes[j] <- tmp[closest]
      tmp <- tmp[-closest]
      # if (!replace) {
      #   control_outcomes <- control_outcomes[-match]
      # }
    }
    HTE[i] <- treated_outcomes[i] - mean(match_outcomes)
  }
  return(HTE)
}

est_cem <- function(df, n_train) {
  n <- nrow(df)
  capture.output(mat <- tryCatch({cem(treatment = 'treated', drop = 'Y', data = df)},
           error = function(e){return(rep(NA, sum(df[(n_train + 1):n, "treated"])))}), file = 'NUL')
  if (all(is.na(mat))) {
    # browser()
    return(mat)
  }
    
  HTE <- vector('numeric', length = n)
  for (i in 1:n) {
    in_MG <- which(mat$strata == mat$strata[i])
    treated <- intersect(in_MG, which(df$treated))
    control <- intersect(in_MG, which(!df$treated))
    HTE[i] <- mean(df$Y[treated]) - mean(df$Y[control])
  }
  HTE <- HTE[intersect(which(df$treated), (n_train + 1):n)]
  # browser()
  return(HTE)
}

est_mahal <- function(df, n_train, f, ratio) {
  m.out <- matchit(f, data = df, 
                   distance = 'mahalanobis',
                   method = "nearest", 
                   ratio = 3,
                   replace = TRUE)
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - apply(matches, 1, function(x) mean(df$Y[as.numeric(x)]))
  HTE <- HTE[which(as.numeric(rownames(matches)) > n_train)]
  return(HTE)
}

est_nn <- function(df, n_train, f, ratio) {
  m.out <-suppressWarnings(matchit(f, data = df, 
                                   method = "nearest", 
                                   ratio = 3, 
                                   replace = TRUE))
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - apply(matches, 1, function(x) mean(df$Y[as.numeric(x)]))
  HTE <- HTE[which(as.numeric(rownames(matches)) > n_train)]
  return(HTE)
}

est_greedy <- function(train_df, 
                       test_df, 
                       test_covs,
                       test_control, 
                       test_treated,
                       n_test_treated,
                       p,
                       bart_fit, 
                       variation = 1,
                       n_req_matches = 5) {
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
    bin_copy <- matrix(bins[i, ,], ncol = 2) # Conversion to matrix handles case of 1 covariate
    # print(sprintf('i = %d', i))
    while (length(already_matched[[i]]) < n_req_matches) { 
      # Find units closest along each axis
      if (variation != 2) { ## Only expand to control units
        potential_matches <- setdiff(test_control, already_matched[[i]])
      }
      else { ## Expand to any units 
        potential_matches <- setdiff(1:nrow(test_df), already_matched[[i]])
      }
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
        
        bin_var[j] <- R_expansion_variance(j, bin_copy, proposed_bin[[j]], train_df, bart_fit)
      }
      
      expand_along <- which.min(bin_var)
      bin_copy <- proposed_bin[[expand_along]]
      in_MG <- apply(test_covs, 1, function(x) all(x >= bin_copy[, 1]) & all(x <= bin_copy[, 2]))
      already_matched[[i]] <- unique(c(already_matched[[i]], which(in_MG)))
    }
    bins[i, , ] <- bin_copy
  }
  message("\n")
  CATE <- get_greedy_CATE(n_test_treated, test_covs, bins, test_df)
  return(list(CATE = CATE, 
              bins = bins))
}

est_MIP_predict <- function(train_df, test_df, 
                    test_treated, n_test_treated, 
                    train_covs, test_covs,
                    n_train, p, lambda=10, alpha=1, m=1, M=1e05) {
  MIP_cates <- vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  message("Running MIP-Predict...")
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    mip_pars =  setup_mip_predict(xi = as.numeric(test_covs[i, ]), yi = test_df$Y[i], 
                                  zi =  1,  y_train = train_df$Y, 
                                  x_train = as.matrix(train_covs), z_train = train_df$treated, 
                                  x_test = as.matrix(test_covs[test_df$treated==0, ]),  
                                  alpha=alpha, lambda=lambda, m=m, M=M)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, sum(test_df$treated==0), p)
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b

    MIP_cates[l] = test_df$Y[i] - mean(test_df$Y[test_df$treated==0][mip_out$w>=0.1])
  }
  message("\n")
  return(list(CATE = MIP_cates, bins = mip_bins))
}

est_MIP_explain <- function(train_df, test_df, 
                            test_treated, n_test_treated, 
                            train_covs, test_covs,
                            n_train, p, lambda=10, alpha=0, m=1, M=1e05) {
  
  mip_cates = vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  message("Running MIP-Explain...")
  bt = bart(data.frame(train_covs, treated=train_df$treated), 
            train_df$Y, keeptrees = T, verbose = FALSE)
  fhat = predict(bt, data.frame(test_covs, treated=test_df$treated))
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    mip_pars =  setup_mip_explain(fhat = fhat[i], xi = as.numeric(test_covs[i, ]), 
                                  y_test = test_df$Y[test_df$treated==0],
                                  x_test = as.matrix(test_covs[test_df$treated==0, ]),  
                                  lambda=lambda, alpha=alpha, m=m, M=M)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, sum(test_df$treated==0), p)
    
    mip_cates[l] = test_df$Y[i] - mean(test_df$Y[test_df$treated==0][mip_out$w>=0.1])
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b
  }
  message("\n")
  return(list(CATE = mip_cates, bins = mip_bins))
}

est_MIQP_variance <- function(train_df, test_df, 
                              test_treated, n_test_treated, 
                              train_covs, test_covs,
                              n_train, p, lambda=10, alpha=1, beta=1, gamma=1, m=1, M=1e05) {
  
  mip_cates = vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  message("Running MIQP-Variance")
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    mip_pars =  setup_miqp_variance(xi = as.numeric(test_covs[i, ]),
                                    y_train = train_df$Y,
                                    x_train = train_covs,
                                    z_train = train_df$treated,
                                    x_test = as.matrix(test_covs[test_df$treated==0, ]),  
                                    alpha=alpha, lambda=lambda, beta=beta, gamma=gamma, m=m, M=M)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, sum(test_df$treated==0), p)
    
    mip_cates[l] = test_df$Y[i] - mean(test_df$Y[test_df$treated==0][mip_out$w>=0.1])
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b
  }
  message("\n")
  return(list(CATE = mip_cates, bins = mip_bins))
}

est_MIQP_fhat <- function(test_df, test_treated, n_test_treated,test_covs, bart_fit0, bart_fit1,
                          n_train, p, lambda0=10, lambda1=10, alpha=1, 
                          beta=1, gamma0=1, gamma1=1, m=1, M=1e05) {
  mip_cates = vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  # fhat1 = predict(bart_fit, newx=as.matrix(cbind(test_covs, treated=1)))
  # fhat1 = fhat1[ ,ncol(fhat1)]
  # fhat0 = predict(bart_fit, newx=as.matrix(cbind(test_covs, treated=0)))
  # fhat0 = fhat0[, ncol(fhat0)]
  fhat0 <- predict(bart_fit0, as.matrix(test_covs))
  fhat1 <- predict(bart_fit1, as.matrix(test_covs))
  # fhat1 = colMeans(predict(bart_fit1, test = as.matrix(test_covs)))
  # fhat0 = colMeans(predict(bart_fit0, test = as.matrix(test_covs)))
  message("Running MIQP-Fhat")
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    mip_pars =  setup_miqp_fhat(xi = as.numeric(test_covs[i, ]),
                                x_test = as.matrix(test_covs[-i, ]),
                                z_test = test_df$treated[-i],
                                fhati1=fhat1[i],
                                fhati0=fhat0[i],
                                fhat1=fhat1[-i],
                                fhat0=fhat0[-i], 
                                alpha=alpha, lambda0=lambda0, lambda1=lambda1,
                                gamma0=gamma0/sd(fhat0), gamma1=gamma1/sd(fhat1), beta=beta, m=m, M=M)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, nrow(test_covs), p)
    
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b
    mg = make_mg(test_covs, mip_out$a, mip_out$b)
    mip_cates[l] = mean(test_df$Y[mg][test_df$treated[mg]]) - mean(test_df$Y[mg][!test_df$treated[mg]])
  }
  message("\n")
  return(list(CATE = mip_cates, bins = mip_bins))
}

est_MIQP_grid <- function(test_df, test_treated, n_test_treated,test_covs, bart_fit,
                          n_train, p, lambda0 = 10, lambda1 = 10, alpha = 0, 
                          beta = 1, gamma0 = 1, gamma1 = 1, m = 1, M = 1e05, grid_size = 3) {
  mip_cates = vector('numeric', n_test_treated)
  mip_bins = array(NA, c(n_test_treated, p, 2))
  
  message("Running MIQP-Grid")
  for (l in 1:n_test_treated){
    i = test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
    
    grid = generate_grid(test_covs[i, ], test_covs[-i, ], grid_size)
    fhat1 = predict(bart_fit, newdata=as.matrix(cbind(grid, treated=1)))
    fhat0 = predict(bart_fit, newdata=as.matrix(cbind(grid, treated=0)))
    fhati1 = predict(bart_fit, newdata=as.matrix(cbind(test_covs, treated=1)))[i]
    fhati0 = predict(bart_fit, newdata=as.matrix(cbind(test_covs, treated=0)))[i]
    mip_pars =  setup_miqp_fhat(xi = as.numeric(test_covs[i, ]),
                                x_test = grid,
                                z_test = rep(1, nrow(grid)),
                                fhati1=fhati1,
                                fhati0=fhati0,
                                fhat1=fhat1,
                                fhat0=fhat0, 
                                alpha=alpha, lambda0=lambda0, lambda1=lambda1,
                                gamma0=gamma0, gamma1=gamma1, beta=beta, m=m, M=M)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense="max", control=list(trace=0))))
    mip_out = recover_pars(sol, n_train, nrow(test_covs), p)
    
    mip_bins[l, ,1] = mip_out$a
    mip_bins[l, ,2] = mip_out$b
    
    mg = make_mg(test_covs, mip_out$a, mip_out$b)
    mip_cates[l] = test_df$Y[i] - mean(test_df$Y[mg][!test_df$treated[mg]])
  }
  message("\n")
  return(list(CATE = mip_cates, bins = mip_bins))
}



get_CATEs <- function(inputs, ratio, estimators, hyperparameters) {
  n_estimators <- length(estimators)
  
  bins <- vector(mode = 'list', length = 3)
  names(bins) <- c('Greedy', 'MIQP-Fhat', 'MIQP-Grid')
  
  c(df, f, n, n_train, p,
    train_df, train_covs, train_control, train_treated,
    test_df, test_covs, test_control, test_treated, 
    n_test_control, n_test_treated, 
    bart_fit0, bart_fit1, counterfactuals) %<-% inputs
  
  c(lambda, alpha, beta, gamma, lambda0, lambda1, gamma0, gamma1, m, M) %<-% hyperparameters
  
  CATEs <- matrix(nrow = n_test_treated, ncol = n_estimators)
  for (i in 1:n_estimators) {
    if (estimators[i] == 'Full Matching') {
      CATEs[, i] <- est_fullmatch(df, n, n_train, f)    
    }
    else if (estimators[i] == 'Prognostic') {
      CATEs[, i] <- est_prog(test_df, counterfactuals, ratio = ratio)
    }
    else if (estimators[i] == 'CEM') {
      CATEs[, i] <- est_cem(df, n_train)
    }
    else if (estimators[i] == 'Mahalanobis') {
      CATEs[, i] <- est_mahal(df, n_train, f, ratio = ratio)
    }
    else if (estimators[i] == 'Nearest Neighbor') {
      CATEs[, i] <- est_nn(df, n_train, f, ratio = ratio)
    }
    else if (estimators[i] == 'GenMatch') {
      CATEs[, i] <- est_gen_match(df, n_train, f, ratio = ratio)
    }
    else if (estimators[i] == 'Greedy') {
      greedy_out <- greedy_cpp(as.matrix(test_covs[test_treated, ]), 
                               test_control - 1, test_treated - 1,
                               as.matrix(test_covs), as.logical(test_df$treated), test_df$Y,
                               1, 15, 1.1, bart_fit0, bart_fit1)
      # (mean(bart_fit0$sigma) + mean(bart_fit1$sigma))
      CATEs[, i] <- greedy_out[[1]]
      lower_bounds <- do.call(rbind, greedy_out[[2]])
      upper_bounds <- do.call(rbind, greedy_out[[3]])
      bins[['Greedy']] <- array(c(lower_bounds, upper_bounds), dim = c(dim(lower_bounds), 2))
    } else if (estimators[i] == 'True'){
      CATEs[, i] <- est_true(df, )
    }
    else if (estimators[i] == 'MIP-Predict') {
      mip_predict_out <- 
        est_MIP_predict(train_df, test_df, test_treated, 
                        n_test_treated, train_covs, test_covs, n_train, p, 
                        lambda=lambda, alpha=alpha, m=m, M=M)
      CATEs[, i] <- mip_predict_out$CATE
      bins[['MIP-Predict']] <- mip_predict_out$bins
    }
    else if (estimators[i] == 'MIP-Explain') {
      mip_explain_out <- 
        est_MIP_explain(train_df, test_df, test_treated, 
                        n_test_treated, train_covs, test_covs, n_train, p,
                        lambda = lambda, alpha=alpha, m=m, M=M)
      CATEs[, i] <- mip_explain_out$CATE
      bins[['MIP-Explain']] <- mip_explain_out$bins
    }
    else if (estimators[i] == 'MIQP-Variance') {
      miqp_variance_out <- 
        est_MIQP_variance(train_df, test_df, test_treated, 
                          n_test_treated, train_covs, test_covs, n_train, p,
                          lambda=lambda, alpha=alpha, beta=beta, gamma=gamma, m=m, M=M)
      CATEs[, i] <- miqp_variance_out$CATE
      bins[['MIQP-Variance']] <- miqp_variance_out$bins
    }    
    else if (estimators[i] == 'MIQP-Fhat') {
      # lasso_fit = glmnet(as.matrix(dplyr::select(train_df, -Y)), 
      #                    train_df$Y, family="gaussian", alpha=1)
      miqp_fhat_out <- 
        est_MIQP_fhat(test_df, test_treated, n_test_treated,test_covs, bart_fit0, bart_fit1,
                      n_train, p, lambda0=lambda0, lambda1=lambda1, 
                      alpha=alpha, beta=beta, gamma0=gamma0, gamma1=gamma1, m=m, M=M)
      CATEs[, i] <- miqp_fhat_out$CATE
      bins[['MIQP-Fhat']] <- miqp_fhat_out$bins
    }
    else if (estimators[i] == 'MIQP-Grid') {
      miqp_grid_out <- 
        est_MIQP_grid(test_df, test_treated, n_test_treated ,test_covs, bart_fit,
      n_train, p, lambda0 = lambda0, lambda1 = lambda1, alpha=alpha, 
      beta = beta, gamma0 = gamma0, gamma1 = gamma1, m = m, M = M, grid_size = 3)
      CATEs[, i] <- miqp_grid_out$CATE
      bins[['MIQP-Grid']] <- miqp_grid_out$bins
    }
    else {
      stop('Unrecognized Estimator')
    }
  }
  return(list(CATEs = CATEs,
              bins = bins))
}

# est_caliper <- function(f, df, ratio = 1) {
#   p <- ncol(df) - 2
#   m.out <- matchit(f, data = df, 
#                    method = 'nearest', 
#                    ratio = ratio, 
#                    mahvars = c('X1' = 'X1'))
#   matches <- m.out$match.matrix
#   HTE <- df$Y[as.numeric(rownames(matches))] - df$Y[as.numeric(matches)]
#   HTE <- HTE[which(as.numeric(rownames(matches)) > n_train)]
#   return(HTE)
# }

