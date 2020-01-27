library(Rcplex)

setup_mip_explain = function(fhat, xi, y_test, x_test, lambda=1, m=1, M=1e10){
  Nop = length(y_test)
  p = length(xi)
  
  #Unit level mip:
  # Objective: max c'X
  cvec = c(rep(0, p),                            # aij
           rep(0, p),                            # bij
           1-lambda*abs(fhat - y_test),          # wij
           rep(0, Nop * p),                      # uij
           rep(0, Nop * p)                       # vij
  )
  
  a_start = 0
  b_start = p
  w_start = p + p
  u_start = p + p + Nop
  v_start = p + p + Nop + Nop * p
  
  # Constraint 2 
  #a_j < x_j 
  a2 = matrix(0, p, p + p + Nop + Nop * p + Nop * p)
  for (j in 1:p){
    a2[j, a_start + j] = 1
  }
  rownames(a2) = rep("C2", nrow(a2))
  b2 = xi
  names(b2) = rep("C2", length(b2))
  s2 = rep("L", p)
  
  # Constraint 3 -b_j < -x_j 
  a3 =  matrix(0, p, p + p + Nop + Nop * p + Nop * p)
  for (j in 1:p){
    a3[j, b_start + j] = 1
  }
  rownames(a3) = rep("C3", nrow(a3))
  b3 = xi
  names(b3) = rep("C3", length(b3))
  s3 = rep("G", p)
  
  # Constraint 4: uikj = 0 if xkj < aij, either 0 or 1 otherwise
  a4 = matrix(0, Nop * p, p + p + Nop + Nop * p + Nop * p)
  b4 = rep(NA, Nop * p)
  l = 1
  for (k in 1:Nop){
    for (j in 1:p){
      a4[l, u_start + l] = M
      a4[l, a_start + j] = 1 
      b4[l] = x_test[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a4) = rep("C4", nrow(a4))
  names(b4) = rep("C4", length(b4))
  s4 = rep("L", Nop * p)
  
  # Constraint 5: uikj = 1 if xkj > aij, either 0 or 1 otherwise
  a5 = matrix(0, Nop * p, p + p + Nop + Nop * p + Nop * p)
  b5 = rep(NA, Nop * p)
  l = 1
  for (k in 1:Nop){
    for (j in 1:p){
      a5[l, u_start + l] = -M
      a5[l, a_start + j] = -1 
      b5[l] = - x_test[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a5) = rep("C5", nrow(a5))
  names(b5) = rep("C5", length(b5))
  s5 = rep("L", Nop * p)
  
  # Constraint 6: vikj = 0 if xkj > bij, either 0 or 1 otherwise
  a6 = matrix(0, Nop * p, p + p + Nop + Nop * p + Nop * p)
  b6 = rep(NA, Nop * p)
  l = 1
  for (k in 1:Nop){
    for (j in 1:p){
      a6[l, v_start + l] = M
      a6[l, b_start + j] = -1 
      b6[l] = -x_test[k, j] + M - 1/M
      l = l + 1
    }
  }
  rownames(a6) = rep("C6", nrow(a6))
  names(b6) = rep("C6", length(b6))
  s6 = rep("L", Nop * p)
  
  # Constraint 7: vikj = 1if xkj < bij, either 0 or 1 otherwise
  a7 = matrix(0, Nop * p, p + p + Nop + Nop * p + Nop * p)
  b7 = rep(NA, Nop * p)
  l = 1
  for (k in 1:Nop){
    for (j in 1:p){
      a7[l, v_start + l] = -M
      a7[l, b_start + j] = 1 
      b7[l] = x_test[k, j] - 1/M
      l = l + 1
    }
  }
  rownames(a7) = rep("C7", nrow(a7))
  names(b7) = rep("C7", length(b7))
  s7 = rep("L", Nop * p)
  
  # Constraint 8: wik = 1 if sum over j uijk + sum over j vijk = 2P, 
  # either 0 or 1 otherwise
  a8 = matrix(0, Nop, p + p + Nop + Nop * p + Nop * p)
  b8 = rep(NA, Nop)
  for (k in 1:Nop){
    a8[k, w_start + k] = -M
    a8[k, (v_start + 1 + p*(k-1)):(v_start + p*k)] = 1 
    a8[k, (u_start + 1 + p*(k-1)):(u_start + p*k)] = 1 
    b8[k] = 2 * p - 1
  }
  rownames(a8) = rep("C8", nrow(a8))
  names(b8) = rep("C8", length(b8))
  s8 = rep("L", Nop)
  
  # Constraint 9: wik = 0 if sum over j uijk + sum over j vijk < 2P, 
  # either 0 or 1 otherwise
  a9 = matrix(0, Nop, p + p + Nop + Nop * p + Nop * p)
  b9 = rep(NA, Nop)
  for (k in 1:Nop){
    a9[k, w_start + k] = M
    a9[k, (v_start + 1 + p*(k-1)):(v_start + p*k)] = -1 
    a9[k, (u_start + 1 + p*(k-1)):(u_start + p*k)] = -1 
    b9[k] = -2 * p + M
  }
  rownames(a9) = rep("C9", nrow(a9))
  names(b9) = rep("C9", length(b9))
  s9 = rep("L", Nop)
  
  # Constraint 10 sum(wik) >= m
  a10 = c(rep(0, p), rep(0, p), rep(1, Nop), rep(0, Nop * p), rep(0, Nop * p))
  b10 = m
  s10 = "G"
  
  Amat = rbind(a2, a3, a4, a5, a6, a7, a8, a9, a10)
  bvec = c(b2, b3, b4, b5, b6, b7, b8, b9, b10)
  svec = c(s2, s3, s4, s5, s6, s7, s8, s9, s10)
  
  lbs = c(apply(x_test, 2, min)-1/M, apply(x_test, 2, min)-1/M, rep(0, Nop),  rep(0, Nop * p), rep(0, Nop * p))
  ubs = c(apply(x_test, 2, max)+1/M, apply(x_test, 2, max)+1/M, rep(1, Nop),  rep(1, Nop * p), rep(1, Nop * p))
  
  vtype = c(rep("C", p), rep("C", p), rep("B", Nop), rep("B", Nop * p), rep("B", Nop * p))
  list(Amat=Amat, bvec=bvec, cvec=cvec, lb=lbs, ub=ubs, vtype=vtype, sense=svec)
}

setup_miqp_variance = function(xi, y_train, x_train, x_test, lambda=1, alpha=0,  m=1, M=1e10){
  n_test = nrow(x_test)
  n_train = nrow(x_train)
  p = length(xi)
  
  #Unit level mip:
  # Objective: max c'X
  cvec = c(alpha * rep(1, p),                       # aij
           alpha * rep(-1, p),                      # bij
           rep(1, n_test),                          # wij
           rep(0, n_test * p),                      # qij
           rep(0, n_test * p),                      # rij
           rep(0, n_train * p),                     # uij
           rep(0, n_train * p),                     # vij
           rep(0, n_train)                          # sk
          )
  
  a_start = 0
  b_start = p
  w_start = b_start + p
  q_start = w_start + n_test
  r_start = q_start + n_test * p
  u_start = r_start + n_test * p
  v_start = u_start + n_train * p
  s_start = v_start + n_train * p
  n_vars = length(cvec)
  
  Qmat = matrix(0, n_vars, n_vars)
  Qmat[(s_start+1):n_vars, (s_start+1):n_vars] = outer(y_train, y_train, 
                                                       FUN=function(x1, x2) - lambda * sqrt((x1 - x2)^2))

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
  a16[1, (w_start+1):(w_start + n_test)] = 1
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
  
  Amat = rbind(a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15,
               a16, a17)
  bvec = c(b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, 
           b17)
  svec = c(s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, 
           s17)
  
  lbs = c(sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          rep(0, n_test),  rep(0, n_test * p), rep(0, n_test * p),
          rep(0, n_train * p), rep(0, n_train * p), rep(0, n_train))
  ubs = c(sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          rep(1, n_test),  rep(1, n_test * p), rep(1, n_test * p),
          rep(1, n_train * p), rep(1, n_train * p), rep(1, n_train))
  
  vtype = c(rep("C", p), rep("C", p), rep("B", n_test), 
            rep("B", n_test * p), rep("B", n_test * p),
            rep("B", n_train * p), rep("B", n_train * p), rep("B", n_train))
  list(Amat=Amat, bvec=bvec, cvec=cvec, Qmat=Qmat, lb=lbs, ub=ubs, vtype=vtype, sense=svec)
}

setup_mip_predict = function(xi, yi, zi, y_train, x_train, z_train, x_test, lambda1=1, 
                             lambda2=1, alpha=0,  m=1, M=1e10){
  n_test = nrow(x_test)
  n_train = nrow(x_train)
  p = length(xi)
  
  #Unit level mip:
  # Objective: max c'X
  cvec = c(alpha * rep(1, p),                       # aij
           alpha * rep(-1, p),                      # bij
           rep(1, n_test),                          # wij
           rep(0, n_test * p),                      # qij
           rep(0, n_test * p),                      # rij
           rep(0, n_train * p),                     # uij
           rep(0, n_train * p),                     # vij
           {- lambda1 * abs(y_train - yi) * abs(zi - z_train)
            - lambda2 * abs(y_train - yi) * (1-abs(zi - z_train))}             # sk
          )
  
  a_start = 0
  b_start = p
  w_start = b_start + p
  q_start = w_start + n_test
  r_start = q_start + n_test * p
  u_start = r_start + n_test * p
  v_start = u_start + n_train * p
  s_start = v_start + n_train * p
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
  a16[1, (w_start+1):(w_start + n_test)] = 1
  rownames(a16) = "C16"
  b16 = m
  names(b16) = "C16"
  s16 = "G"

  # Constraint 17 sum(sikc * t) >= m
  a17 = matrix(0, 1,  n_vars)
  a17[1, (s_start+1):(s_start + n_train)] = z_train
  rownames(a17) = "C17"
  b17 = m
  names(b17) = "C17"
  s17 = "G"

  # Constraint 18 sum(sikc * 1-t) >= m
  a18 = matrix(0, 1,  n_vars)
  a18[1, (s_start+1):(s_start + n_train)] = 1- z_train
  rownames(a18) = "C18"
  b18 = m
  names(b18) = "C18"
  s18 = "G"

  Amat = rbind(a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15,
               a16, a17, a18)
  bvec = c(b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, 
           b17, b18)
  svec = c(s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, 
           s17, s18)
  
  lbs = c(sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          sapply(1:p, function(j) min(c(x_train[, j], x_test[, j])))-1/M, 
          rep(0, n_test),  rep(0, n_test * p), rep(0, n_test * p),
          rep(0, n_train * p), rep(0, n_train * p), rep(0, n_train))
  ubs = c(sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          sapply(1:p, function(j) max(c(x_train[, j], x_test[, j])))+1/M, 
          rep(1, n_test),  rep(1, n_test * p), rep(1, n_test * p),
          rep(1, n_train * p), rep(1, n_train * p), rep(1, n_train))
  
  vtype = c(rep("C", p), rep("C", p), rep("B", n_test), 
            rep("B", n_test * p), rep("B", n_test * p),
            rep("B", n_train * p), rep("B", n_train * p), rep("B", n_train))
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
    
    list(
        a = sol$xopt[(a_start+1):b_start],
        b = sol$xopt[(b_start+1):w_start],
        w = sol$xopt[(w_start+1):(q_start)],
        q = matrix(sol$xopt[(q_start+1):(r_start)], n_test, p, byrow=T),
        r = matrix(sol$xopt[(r_start+1):(u_start)], n_test, p, byrow=T),
        u = matrix(sol$xopt[(u_start+1):(v_start)], n_train, p, byrow=T),
        v = matrix(sol$xopt[(v_start+1):s_start], n_train, p, byrow=T),
        s = sol$xopt[(s_start + 1):(n_train + s_start)]
     )
}


