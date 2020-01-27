
source('helpers.R')
# source("AB_MIP.R")

require(MatchIt)
require(cem)
require(dbarts)

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

est_prog <- function(df, counterfactuals) {
  treated_outcomes <- df$Y[df$treated]
  control_outcomes <- df$Y[!df$treated]
  HTE <- rep(NA, length(treated_outcomes))
  for (i in 1:length(treated_outcomes)) {
    if (length(control_outcomes) == 0) {
      break
    }
    match <- which.min(abs(counterfactuals[i] - control_outcomes))
    HTE[i] <- treated_outcomes[i] - control_outcomes[match]
    control_outcomes <- control_outcomes[-match]
  }
  return(HTE)
}

est_cem <- function(df, n_train) {
  n <- nrow(df)
  capture.output(mat <- cem(treatment = 'treated', drop = 'Y', data = df), file = 'NUL')
  HTE <- vector('numeric', length = n)
  for (i in 1:n) {
    in_MG <- which(mat$strata == mat$strata[i])
    treated <- intersect(in_MG, which(df$treated))
    control <- intersect(in_MG, which(!df$treated))
    HTE[i] <- mean(df$Y[treated]) - mean(df$Y[control])
  }
  HTE <- HTE[intersect(which(df$treated), (n_train + 1):n)]
  
  return(HTE)
}

est_mahal <- function(df, n_train, f, ratio = 1) {
  m.out <- matchit(f, data = df, 
                   distance = 'mahalanobis',
                   method = "nearest", 
                   ratio = ratio)
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - df$Y[as.numeric(matches)]
  HTE <- HTE[which(as.numeric(rownames(matches)) > n_train)]
  return(HTE)
}

est_nn <- function(df, n_train, f, ratio = 1) {
  m.out <- matchit(f, data = df, 
                   method = "nearest", 
                   ratio = ratio) %>%
    suppressWarnings()
  
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - df$Y[as.numeric(matches)]
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
                       bart_fit) {
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
    while (length(already_matched[[i]]) < 5) { 
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
    }
    
    bins[i, , ] <- bin_copy
  }
  message("\n")
  CATE <- get_greedy_CATE(n_test_treated, test_covs, bins, test_df)
  return(CATE = CATE, 
         bins = bins)
}

est_MIP <- function(train_df, test_df, 
                    test_treated, n_test_treated, 
                    train_covs, test_covs,
                    n_train) {
  MIP_cates <- vector('numeric', n_test_treated)
  message("Running MIP AB...")
  for (l in 1:n_test_treated){
    i <- test_treated[l]
    message(paste("Matching unit", l, "of", n_test_treated), "\r", appendLF = FALSE); flush.console()
  
    mip_pars <- create_unit_mip(xi = as.numeric(test_covs[i, ]), zi =  1, y_train = train_df$Y,
                                x_train = as.matrix(train_covs), z_train = train_df$treated,
                                x_test = as.matrix(test_covs), z_test = test_df$treated,
                                alpha = 5, lambda = 20, m = 1, M = 1e5)
    sol <- do.call(Rcplex, c(mip_pars, list(objsense = "max", control = list(trace=0))))
    mip_out <- recover_pars(sol, n_train, n_test, p)
  
    MIP_cates[l] <- test_df$Y[i] - mean(test_df$Y[mip_out$w >= 0.1 & test_df$treated == 0])
  }
  message("\n")
  return(MIP_cates)
}

get_CATEs <- function(inputs) {
  n_estimators <- 7 - 1
  n_MIPs <- 0
  bins <- vector(mode = 'list', length = n_MIPs + 1)
  ab_names <- c('Greedy')
  if (n_MIPs > 0) {
    ab_names <- c(ab_names, paste0(rep('MIP ', n_MIPs), 1:n_MIPs))
  }
  names(bins) <- ab_names
  
  c(df, f, n, n_train, p,
    train_df, train_covs, train_control, train_treated,
    test_df, test_covs, test_control, test_treated, 
    n_test_control, n_test_treated, 
    bart_fit, counterfactuals) %<-% inputs
  
  CATEs <- matrix(nrow = n_test_treated, ncol = n_estimators)
  CATEs[, 1] <- est_fullmatch(df, n, n_train, f)
  CATEs[, 2] <- est_prog(test_df, counterfactuals)
  CATEs[, 3] <- est_cem(df, n_train)
  CATEs[, 4] <- est_mahal(df, n_train, f, ratio = 1)
  CATEs[, 5] <- est_nn(df, n_train, f, ratio = 1)
  greedy_out <- est_greedy(train_df, 
                           test_df, 
                           test_covs, 
                           test_control, 
                           test_treated, 
                           n_test_treated, 
                           p, 
                           bart_fit)
  CATEs[, 6] <- greedy_out$CATE
  bins[['Greedy']] <- greedy_out$bins
  
  # CATEs[, 7] <- est_MIP(df)$CATEs
  # bins[['MIP 1']] <- est_MIP(df)$bins
  # CATEs[, 8] <- est_caliper(f, df) # Doesn't work...? 
  return(CATEs)
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

