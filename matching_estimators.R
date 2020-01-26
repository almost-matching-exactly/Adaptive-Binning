est_fullmatch <- function(df, f) {
  m.out <- suppressWarnings(matchit(f, data = df, method = 'full'))
  n = nrow(df)
  HTE <- vector('numeric', length = n)
  for (i in 1:n) {
    in_MG <- which(m.out$subclass == m.out$subclass[i])
    treated <- intersect(in_MG, which(df$treated))
    control <- intersect(in_MG, which(!df$treated))
    HTE[i] <- mean(df$Y[treated]) - mean(df$Y[control])
  }
  HTE <- HTE[which(df$treated)]
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

est_cem <- function(df) {
  n = nrow(df)
  
  HTE <- tryCatch({capture.output(mat <- cem(treatment = 'treated', drop = 'Y', data = df), file = 'NUL')
    HTE <- vector('numeric', length = n)
    size <- vector('numeric', length = n)
    for (i in 1:n) {
      in_MG <- which(mat$strata == mat$strata[i])
      treated <- intersect(in_MG, which(df$treated))
      control <- intersect(in_MG, which(!df$treated))
      size[i] <- length(in_MG)
      HTE[i] <- mean(df$Y[treated]) - mean(df$Y[control])
    }
    HTE[which(df$treated)]},
  error=function(e){
    NA
  })
  
  return(HTE)
}

est_mahal <- function(f, df, ratio = 1) {
  m.out <- matchit(f, data = df, 
                   distance = 'mahalanobis',
                   method = "nearest", 
                   ratio = ratio)
  n= nrow(df)
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - df$Y[as.numeric(matches)]
  HTE <- HTE[which(df$treated)]
  return(HTE)
}

est_nn <- function(f, df, ratio = 1) {
  m.out <- matchit(f, data = df, 
                   method = "nearest", 
                   ratio = ratio)
  
  matches <- m.out$match.matrix
  HTE <- df$Y[as.numeric(rownames(matches))] - df$Y[as.numeric(matches)]
  HTE <- HTE[which(df$treated)]
  return(HTE)
}
