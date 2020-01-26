library(Rcplex)
library(MatchIt)
library(data.table)
library(tidyverse)

create_unit_mip = function(xi, zi, y_train, x_train, z_train, x_test, z_test, 
                           lambda=1, alpha=1,  m=1, M=1e10){
  n_test = nrow(x_test)
  n_train = nrow(x_train)
  p = length(xi)
  
  #Unit level mip:
  # Objective: max c'X
  cvec = c(alpha * rep(1, p),                               # aij
           alpha * rep(-1, p),                               # bij
           rep(1, n_test),                          # wij
           rep(0, n_test * p),                      # qij
           rep(0, n_test * p),                      # rij
           rep(0, n_train * p),                      # uij
           rep(0, n_train * p),                      # vij
           rep(0, n_train),                          # sk
           rep(0, n_train),                          # dk
           rep(0, n_train),                          # ek
           -lambda, lambda
)
  
  a_start = 0
  b_start = p
  w_start = b_start + p
  q_start = w_start + n_test
  r_start = q_start + n_test * p
  u_start = r_start + n_test * p
  v_start = u_start + n_train * p
  s_start = v_start + n_train * p
  d_start = s_start + n_train
  e_start = d_start + n_train
  ymax_start = e_start + n_train
  ymin_start = ymax_start + 1
  n_vars = length(cvec)
  
  # Constraint 2 a_j < x_j 
  a2 = matrix(0, p, n_vars)
  for (j in 1:p){
    a2[j, a_start + j] = 1
  }
  rownames(a2) = rep("C2", nrow(a2))
  b2 = xi
  names(b2) = rep("C2", length(b2))
  s2 = rep("L", p)
  
  # Constraint 3 -b_j < -x_j 
  a3 =  matrix(0, p, n_vars)
  for (j in 1:p){
    a3[j, b_start + j] = 1
  }
  rownames(a3) = rep("C3", nrow(a3))
  b3 = xi
  names(b3) = rep("C3", length(b3))
  s3 = rep("G", p)
  
  ### Training set constraints ###
  # Constraint 4: uikj = 0 if xkj < aij, either 0 or 1 otherwise
  a4 = matrix(0, n_train * p, n_vars)
  b4 = rep(NA, n_train * p)
  l = 1
  for (k in 1:n_train){
    for (j in 1:p){
      a4[l, u_start + l] = M
      a4[l, a_start + j] = 1 
      b4[l] = x_train[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a4) = rep("C4", nrow(a4))
  names(b4) = rep("C4", length(b4))
  s4 = rep("L", n_train * p)
  
  # Constraint 5: uikj = 1 if xkj > aij, either 0 or 1 otherwise
  a5 = matrix(0, n_train * p, n_vars)
  b5 = rep(NA, n_train * p)
  l = 1
  for (k in 1:n_train){
    for (j in 1:p){
      a5[l, u_start + l] = -M
      a5[l, a_start + j] = -1 
      b5[l] = - x_train[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a5) = rep("C5", nrow(a5))
  names(b5) = rep("C5", length(b5))
  s5 = rep("L", n_train * p)
  
  # Constraint 6: vikj = 0 if xkj > bij, either 0 or 1 otherwise
  a6 = matrix(0, n_train * p, n_vars)
  b6 = rep(NA, n_train * p)
  l = 1
  for (k in 1:n_train){
    for (j in 1:p){
      a6[l, v_start + l] = M
      a6[l, b_start + j] = -1 
      b6[l] = -x_train[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a6) = rep("C6", nrow(a6))
  names(b6) = rep("C6", length(b6))
  s6 = rep("L", n_train * p)
  
  # Constraint 7: vikj = 1if xkj < bij, either 0 or 1 otherwise
  a7 = matrix(0, n_train * p, n_vars)
  b7 = rep(NA, n_train * p)
  l = 1
  for (k in 1:n_train){
    for (j in 1:p){
      a7[l, v_start + l] = -M
      a7[l, b_start + j] = 1 
      b7[l] = x_train[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a7) = rep("C7", nrow(a7))
  names(b7) = rep("C7", length(b7))
  s7 = rep("L", n_train * p)
  
  # Constraint 8: sik = 1 if sum over j uijk + sum over j vijk = 4P, 
  # either 0 or 1 otherwise
  a8 = matrix(0, n_train, n_vars)
  b8 = rep(NA, n_train)
  for (k in 1:n_train){
    a8[k, s_start + k] = -M
    a8[k, (v_start + 1 + p*(k-1)):(v_start + p*k)] = 1 
    a8[k, (u_start + 1 + p*(k-1)):(u_start + p*k)] = 1 
    b8[k] = 2 * p - 1
  }
  rownames(a8) = rep("C8", nrow(a8))
  names(b8) = rep("C8", length(b8))
  s8 = rep("L", n_train)
  
  # Constraint 9: wik = 0 if sum over j uijk + sum over j vijk < 2P, 
  # either 0 or 1 otherwise
  a9 = matrix(0, n_train, n_vars)
  b9 = rep(NA, n_train)
  for (k in 1:n_train){
    a9[k, s_start + k] = M
    a9[k, (v_start + 1 + p*(k-1)):(v_start + p*k)] = -1 
    a9[k, (u_start + 1 + p*(k-1)):(u_start + p*k)] = -1 
    b9[k] = -2 * p + M
  }
  rownames(a9) = rep("C9", nrow(a9))
  names(b9) = rep("C9", length(b9))
  s9 = rep("L", n_train)
  
  ### Test set constraints ###
  # Constraint 10: qikj = 0 if xkj < aij, either 0 or 1 otherwise
  a10 = matrix(0, n_test * p, n_vars)
  b10 = rep(NA, n_test * p)
  l = 1
  for (k in 1:n_test){
    for (j in 1:p){
      a10[l, q_start + l] = M
      a10[l, a_start + j] = 1 
      b10[l] = x_test[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a10) = rep("C10", nrow(a10))
  names(b10) = rep("C10", length(b10))
  s10 = rep("L", n_test * p)
  
  # Constraint 11: qikj = 1 if xkj > aij, either 0 or 1 otherwise
  a11 = matrix(0, n_test * p, n_vars)
  b11 = rep(NA, n_test * p)
  l = 1
  for (k in 1:n_test){
    for (j in 1:p){
      a11[l, q_start + l] = -M
      a11[l, a_start + j] = -1 
      b11[l] = - x_test[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a11) = rep("C11", nrow(a11))
  names(b11) = rep("C11", length(b11))
  s11 = rep("L", n_test * p)
  
  # Constraint 12: rikj = 0 if xkj > bij, either 0 or 1 otherwise
  a12 = matrix(0, n_test * p, n_vars)
  b12 = rep(NA, n_test * p)
  l = 1
  for (k in 1:n_test){
    for (j in 1:p){
      a12[l, r_start + l] = M
      a12[l, b_start + j] = -1 
      b12[l] = -x_test[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a12) = rep("C12", nrow(a12))
  names(b12) = rep("C12", length(b12))
  s12 = rep("L", n_test * p)
  
  # Constraint 13: rikj = 1if xkj < bij, either 0 or 1 otherwise
  a13 = matrix(0, n_test * p, n_vars)
  b13 = rep(NA, n_test * p)
  l = 1
  for (k in 1:n_test){
    for (j in 1:p){
      a13[l, r_start + l] = -M
      a13[l, b_start + j] = 1 
      b13[l] = x_test[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a13) = rep("C13", nrow(a13))
  names(b13) = rep("C13", length(b13))
  s13 = rep("L", n_test * p)
  
  # Constraint 14: wik = 1 if sum over j uijk + sum over j vijk = 4P, 
  # either 0 or 1 otherwise
  a14 = matrix(0, n_test,  n_vars)
  b14 = rep(NA, n_test)
  for (k in 1:n_test){
    a14[k, w_start + k] = -M
    a14[k, (r_start + 1 + p*(k-1)):(r_start + p*k)] = 1 
    a14[k, (q_start + 1 + p*(k-1)):(q_start + p*k)] = 1 
    b14[k] = 2 * p - 1
  }
  rownames(a14) = rep("C14", nrow(a14))
  names(b14) = rep("C14", length(b14))
  s14 = rep("L", n_test)
  
  # Constraint 15: wik = 0 if sum over j uijk + sum over j vijk < 2P, 
  # either 0 or 1 otherwise
  a15 = matrix(0, n_test, n_vars)
  b15 = rep(NA, n_test)
  for (k in 1:n_test){
    a15[k, w_start + k] = M
    a15[k, (r_start + 1 + p*(k-1)):(r_start + p*k)] = -1 
    a15[k, (q_start + 1 + p*(k-1)):(q_start + p*k)] = -1 
    b15[k] = -2 * p + M
  }
  rownames(a15) = rep("C15", nrow(a15))
  names(b15) = rep("C15", length(b15))
  s15 = rep("L", n_test)

  # Constraint 16 sum(wik) >= m
  a16 = matrix(0, 1, n_vars)
  a16[1, (w_start+1):(w_start + n_test)] = abs(zi - z_test)
  rownames(a16) = "C16"
  b16 = m
  names(b16) = "C16"
  s16 = "G"

  # Constraint 17 sum(sikc) >= m
  a17 = matrix(0, 1,  n_vars)
  a17[1, (s_start+1):(s_start + n_train)] = 1
  rownames(a17) = "C17"
  b17 = m
  names(b17) = "C17"
  s17 = "G"
  
  # Constraint 18 yksk <= ymax
  a18 = matrix(0, n_train, n_vars)
  b18 = rep(NA, n_train)
  for (k in 1:n_train){
    a18[k, s_start + k] = y_train[k] + M
    a18[k, ymax_start + 1] = -1
    b18[k] = M
  }
  rownames(a18) = rep("C18", nrow(a18))
  names(b18) = rep("C18", length(b18))
  s18 = rep("L", n_train)
  
  
  # Constraint 19 ymax <= yksk + M(1-dk)
  a19 = matrix(0, n_train, n_vars)
  b19 = rep(NA, n_train)
  for(k in 1:n_train){
    a19[k, s_start + k] = -y_train[k] 
    a19[k, ymax_start + 1] = 1
    a19[k, d_start + k] = M
    b19[k] = M
  }
  rownames(a19) = rep("C19", nrow(a19))
  names(b19) = rep("C19", length(b19))
  s19 = rep("L", n_train)
  
  # Constraint 20 sum dk = 1
  a20 = matrix(0, 1, n_vars)
  a20[1, (d_start+1):(d_start + n_train)] = 1
  b20 = 1
  s20 = "E"
  rownames(a20) = rep("C20", nrow(a20))
  names(b20) = rep("C20", length(b20))
  
  # Constraint 21 dk <= sk
  a21 = matrix(0, n_train, n_vars)
  b21 = rep(0, n_train)
  for (k in 1:n_train){
    a21[k, s_start + k] = -1
    a21[k, d_start + k] = 1
  }
  s21 = rep("L", n_train)
  rownames(a21) = rep("C21", nrow(a21))
  names(b21) = rep("C21", length(b21))

  #Constraint 22 yksk >= ymin
  a22 = matrix(0, n_train, n_vars)
  b22 = rep(0, n_train)
  for(k in 1:n_train){
    a22[k, s_start + k] = y_train[k] - M
    a22[k, ymin_start + 1] = -1
    b22[k] = -M
  }
  s22 = rep("G", n_train)
  rownames(a22) = rep("C22", nrow(a22))
  names(b22) = rep("C22", length(b22))
  
  # Constraint 23 ymin > ykwk if ek = 1
  a23 = matrix(0, n_train, n_vars)
  b23 = rep(0, n_train)
  for(k in 1:n_train){
    a23[k, s_start + k] = - y_train[k]
    a23[k, ymin_start + 1] = 1
    a23[k, e_start + k] = -M
    b23[k] = -M
  }
  s23 = rep("G", n_train)
  rownames(a23) = rep("C23", nrow(a23))
  names(b23) = rep("C23", length(b23))
  
  # Constraint 24 sum ek = 1
  a24 = matrix(0, 1, n_vars)
  a24[1, (e_start+1):(e_start+n_train)] = 1
  rownames(a24) = "C24"
  b24 = 1
  names(b24) = "C24"
  s24 = "E"
  
  # Constraint 25 ek <= sk
  a25 = matrix(0, n_train, n_vars)
  b25 = rep(0, n_train)
  for(k in 1:n_train){
    a25[k, s_start + k] = -1
    a25[k, e_start + k] = 1
  }
  s25 = rep("L", n_train)
  rownames(a25) = rep("C25", nrow(a25))
  names(b25) = rep("C25", length(b25))
  
  Amat = rbind(a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15,
               a16, a17, a18, a19, a20, a21, a22, a23, a24, a25)
  bvec = c(b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, 
           b17, b18, b19, b20, b21, b22, b23, b24, b25)
  svec = c(s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, 
           s17, s18, s19, s20, s21, s22, s23, s24, s25)
  
  lbs = c(sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          rep(0, n_test),  rep(0, n_test * p), rep(0, n_test * p),
          rep(0, n_train * p), rep(0, n_train * p), rep(0, n_train),  
          rep(0, n_train),  rep(0, n_train), -M, -M)
  ubs = c(sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          rep(1, n_test),  rep(1, n_test * p), rep(1, n_test * p),
          rep(1, n_train * p), rep(1, n_train * p), rep(1, n_train),
          rep(1, n_train),  rep(1, n_train), M, M)
  
  vtype = c(rep("C", p), rep("C", p), rep("B", n_test), 
            rep("B", n_test * p), rep("B", n_test * p),
            rep("B", n_train * p), rep("B", n_train * p), rep("B", n_train),
            rep("B", n_train), rep("B", n_train), "C", "C")
  list(Amat=Amat, bvec=bvec, cvec=cvec, lb=lbs, ub=ubs, vtype=vtype, sense=svec)
}

recover_pars = function(sol, n_train, n_test, p){
   a_start = 0
    b_start = p
    w_start = b_start + p
    q_start = w_start + n_test
    r_start = q_start + n_test * p
    u_start = r_start + n_test * p
    v_start = u_start + n_train * p
    s_start = v_start + n_train * p
    d_start = s_start + n_train
    e_start = d_start + n_train
    ymax_start = e_start + n_train
    ymin_start = ymax_start + 1
    
    list(
        a = sol$xopt[(a_start+1):b_start],
        b = sol$xopt[(b_start+1):w_start],
        w = sol$xopt[(w_start+1):(q_start)],
        q = matrix(sol$xopt[(q_start+1):(r_start)], n_test, p, byrow=T),
        r = matrix(sol$xopt[(r_start+1):(u_start)], n_test, p, byrow=T),
        u = matrix(sol$xopt[(u_start+1):(v_start)], n_train, p, byrow=T),
        v = matrix(sol$xopt[(v_start+1):s_start], n_train, p, byrow=T),
        s = sol$xopt[(s_start + 1):(d_start)],
        d = sol$xopt[(d_start + 1):(e_start)],
        e = sol$xopt[(e_start + 1):(ymax_start)],
        ymax = sol$xopt[ymax_start+1],
        ymin = sol$xopt[ymin_start + 1]
    )
 
}


make_mgs = function(X, bins){
  MGs = list()
  for(i in 1:length(bins)){
    b = bins[[i]]
    MGs[[i]] = which(apply(X, 1, function(x) all(x >= b$lbs) && all(x <= b$ubs)))
    if (!i %in% MGs[[i]])
      MGs[[i]] = c(i, MGs[[i]])
  }
  return(MGs)
}

get_ate = function(Y, Z, MGs){
  cates = c()
  for (mg in MGs){
    if(length(mg)==0)
      next
    if(sum(Z[mg]==1)==0 || sum(Z[mg]==0)==0)
      next
    cates = c(cates, mean(Y[mg][Z[mg]==1]) - mean(Y[mg][Z[mg]==0]))
  }
  return(mean(cates))
}

get_training_cates = function(Y, Z, bins){
  cates = rep(NA, length(bins))
  for (i in seq_along(bins)){
    if(length(bins[[i]]$mg) == 0)
      next
    if(Z[i]==1)
      cates[i] = mean(Y[i] - Y[Z==0][bins[[i]]$mg])
    else
      cates[i] = mean(Y[Z==1][bins[[i]]$mg] - Y[i])
  }
  cates
}

make_mmgs_matching_on_Y = function(y_train, z_train, y_test, z_test, X_test, 
                            lbs, ubs, eps=1e-10){
  all_mmgs = make_all_mmgs(X_test, z_test, lbs, ubs, eps)
  mgs = list()
  for(i in seq_along(all_mmgs)){
    candidates = as.integer(names(all_mmgs[[i]]))
    candidates = candidates[z_train[candidates] == z_test[i]]
    win = as.character(candidates[which.min(abs(y_test[i] - y_train[candidates]))])
    mgs[[i]] = all_mmgs[[i]][win]
  }
  mgs
}

make_mmgs_best_mq = function(z_train, y_test, z_test, X_test, lbs, ubs, mq, eps=1e-10){
  all_mmgs = make_all_mmgs(X_test, z_test, lbs, ubs, eps)
  mgs = list()
  for(i in seq_along(all_mmgs)){
    candidates = as.integer(names(all_mmgs[[i]]))
    candidates = candidates[z_train[candidates] == z_test[i]]
    win = as.character(candidates[which.min(mq[candidates])])
    mgs[[i]] = all_mmgs[[i]][win]
  }
  mgs
}


make_all_mmgs = function(X_test, z_test, lbs, ubs, eps=1e-10){
  mmgs = list()
  for (i in seq_along(z_test)){
    mmgs[[i]] = list()
    candidates = which(colSums((X_test[i, ] >= t(lbs) - eps) * 
                                 (X_test[i, ] <= t(ubs) + eps)) == ncol(X_test))
    Xopp = X_test[z_test != z_test[i], ]
    all_mmgs = sapply(candidates, function(c) {
      which(colSums((t(Xopp) >= lbs[c, ] - eps) * 
                      (t(Xopp) <= ubs[c, ] + eps)) == ncol(X_test))}, simplify = F)
    mmgs[[i]] = all_mmgs[unlist(lapply(all_mmgs, length) > 0)]
    names(mmgs[[i]]) = candidates[unlist(lapply(all_mmgs, length) > 0)]
  }
  mmgs
}

get_all_cates_from_mmgs = function(Y, Z, mmgs){
  cates = list()
  for(i in seq_along(Y)){
    cates[[i]] = rep(NA, length(mmgs[[i]]))
    for(k in seq_along(mmgs[[i]])){
      mg = mmgs[[i]][[k]]
      if(Z[i]==1)
        cates[[i]][k] = mean(Y[i] - Y[Z==0][mg])
      else
        cates[[i]][k] = mean(Y[Z==1][mg] - Y[i])
    }
  }
  cates
}

get_cates_from_mmgs = function(Y, Z, mgs){
  cates = rep(NA, length(mgs))
  for (i in seq_along(mgs)){
    if(length(mgs[[i]]) == 0)
      next
    if(Z[i]==1)
      cates[i] = mean(Y[i] - Y[Z==0][unlist(mgs[[i]])])
    else
      cates[i] = mean(Y[Z==1][unlist(mgs[[i]])] - Y[i])
  }
  cates
}

make_mg_by_interpolation = function(xi, X, lbs, ubs, eps=1e-10){
  if(all(lbs < xi)){
    closest_lbs = max(X)
  }else{
    c_lbs = matrix(lbs[lbs >= xi, ])
    closest_lbs = c_lbs[which.min(sqrt(colSums((xi - t(c_lbs))^2))), ]
  }
  if(all(ubs > xi)){
    closest_ubs = min(X)
  }else{
    c_ubs = matrix(ubs[ubs <= xi, ])
    closest_ubs = ubs[which.min(sqrt(colSums((xi - t(c_ubs))^2))), ]
  }
  mg = which(colSums((t(X) >= closest_ubs - eps) * (t(X) <= closest_lbs + eps)) == length(xi))
  mg
}

bin_covariate = function(lbs, ubs, max=3, interpolate=T, eps=1e-3){
  glbs = c()
  gubs = c()
  while(length(lbs) > 0 && length(ubs) > 0){
    lb = min(lbs)
    ubs = ubs[ubs > lb + eps]
    ub = ubs[which.min(abs(ubs - lb))]
    lbs = lbs[lbs > ub + eps]
    ub2 = lbs[which.min(abs(lbs - ub))]
    lbs = lbs[lbs > ub2 + eps]
    if(interpolate)
      lbs = c(lbs, ub2)
    ubs = ubs[ubs > ub2 + eps]
    glbs = c(glbs, lb, ub)
    gubs = c(gubs, ub, ub2)
  }
  list(glbs, c(gubs, max))
}

make_mg_with_bounds = function(X, lbs, ubs, eps=1e-10){
  which(colSums((t(X) >= lbs - eps) * (t(X) <= ubs + eps)) == ncol(X))
}

find_bounds = function(x, lbs, ubs){
  lbsi = rep(NA, length(x))
  ubsi = rep(NA, length(x))
  for (j in seq_along(x)){
    lbj = lbs[[j]]
    ubj = ubs[[j]]
    ubsi[j] = min(ubj[ubj >= x[j]])
    lbsi[j] = max(lbj[lbj <= x[j]])
  }
  cbind(lbsi, ubsi)
}

bin_all_covariates = function(lbs, ubs, max=3, interpolate=T, eps=1e-3){
  bin_lbs = list()
  bin_ubs = list()
  for(j in 1:ncol(lbs)){
    bins = bin_covariate(lbs[!is.na(lbs[, j]), j], ubs[!is.na(ubs[, j]), j], max, interpolate, eps)
    bin_lbs[[j]] = bins[[1]]
    bin_ubs[[j]] = bins[[2]]
  }
  list(bin_lbs, bin_ubs)
}



