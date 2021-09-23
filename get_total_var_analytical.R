require(mvtnorm)

get_tot_var = function(B = NULL, var_Ek = 0.2, var_cj = 0.2, var_Bi = NULL, S = NULL, D = NULL) {
  # B: a SxD matrix, or NULL (B is generated automatically)
  # S: number of species
  # D: dimensionality of resources
  
  if(!is.null(B)) {
    S = dim(B)[1]
    D = dim(B)[2]
  } else {
    if(is.null(S) | is.null(D) | is.null(var_Bi)) {
      return("Error: if B is NULL, then S, D, and var_Bi must be specified")
    } else {
      B = rmvnorm(S, sigma = diag(D)*var_Bi)
    }
  }
  
  varTot = mean(diag(B %*% (diag(D)*var_Ek) %*% t(B)) + var_cj) +
    mean((rowSums(B)/D-mean(rowSums(B)/D))^2)
  
  return(varTot)
}

iterate_totvar = function(niter = 1000, B = NULL, var_Ek = 0.2, var_cj = 0.2, var_Bi = NULL, S = NULL, D = NULL) {
  totvar = numeric(niter)
  for(i in 1:niter) {
    totvar[i] = get_tot_var(B, var_Ek, var_cj, var_Bi, S, D)
  }
  return(totvar)
}

find_opt_par = function(tot_var = NULL, nmcmc = 200, niter = 1000, B = NULL, var_Ek = NULL, var_cj = NULL, var_Bi = NULL, S = NULL, D = NULL, minv = 0, maxv = 2, findB = FALSE) {
  # tot_var: desired total variance
  # nmcmc: number of mcmc iterations
  # minv: minimum value for paramter to be optimised
  # maxv: maximum value for paramter to be optimised
  # note: if no B mat is given, then function simulates niter differnet B matrices and uses the average attribute variance
  # if findB is TRUE then a B will be identified that attempts to satisfy the tot_var requirement
  
  if(is.null(tot_var)) {
    if(is.null(B)) {
      totvar = mean(iterate_totvar(niter = niter, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D))
    } else {
      totvar = mean(iterate_totvar(niter = 1, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D))
    }
    return(totvar)
  } else {
    estout = numeric(nmcmc)
    diffout = numeric(nmcmc)
    
    # find missing value
    which_missing = which(c(is.null(var_Ek), is.null(var_cj), is.null(var_Bi), is.null(B)))
    
    # first step
    estout[1] = runif(1, minv, maxv)
    if(1 %in% which_missing) {
      var_Ek = estout[1]
    } else if(2 %in% which_missing) {
      var_cj = estout[1]
    } else if(3 %in% which_missing) {
      var_Bi = estout[1]
    }
    
    if(is.null(B) & !findB) {
      diffout[1] = abs(tot_var-mean(iterate_totvar(niter = niter, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D)))
    } else {
      if(findB){
        B = rmvnorm(S, sigma = diag(D)*var_Bi)
      }
      diffout[1] = abs(tot_var-mean(iterate_totvar(niter = 1, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D)))
    }
    
    for(i in 2:nmcmc) {
      estnew = runif(1, minv, maxv)
      if(1 %in% which_missing) {
        var_Ek = estnew
      } else if(2 %in% which_missing) {
        var_cj = estnew
      } else if(3 %in% which_missing) {
        var_Bi = estnew
      }
      
      if(is.null(B)) {
        diffoutnew = abs(tot_var-mean(iterate_totvar(niter = niter, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D)))
      } else {
        if(findB){
          Bold = B
          B = rmvnorm(S, sigma = diag(D)*var_Bi)
        }
        diffoutnew = abs(tot_var-mean(iterate_totvar(niter = 1, B = B, var_Ek = var_Ek, var_cj = var_cj, var_Bi = var_Bi, S = S, D = D)))
      }
        
      if(diffoutnew<diffout[i-1]) {
        estout[i] = estnew
        diffout[i] = diffoutnew
      } else {
        estout[i] = estout[i-1]
        diffout[i] = diffout[i-1]
        B = Bold
      }
    }
    if(findB){
      return(list(opt = data.frame(estout=estout, diffout=diffout), B = B))
    } else {
      return(data.frame(estout=estout, diffout=diffout))
    }
  }
}


if(FALSE) {
  # demonstrate optimiser function
  S = 14; D = 3; var_Bi = 0.15
  B = rmvnorm(S, sigma = diag(D)*var_Bi)
  
  ## with fixed B
  # find var_Ek
  var_Ek_est = find_opt_par(tot_var = 0.5, B = B, var_Ek = NULL, var_cj = 0.2, maxv = 2)
  var_Ek_est[nrow(var_Ek_est),]
  # estimate, and absolute difference between achieved and desired total variance
  plot(var_Ek_est$estout, type = "l")  # estimate over iterations
  plot(var_Ek_est$diffout, type = "l") # distance over iterations
  
  # find var_cj
  var_cj_est = find_opt_par(tot_var = 0.5, B = B, var_Ek = 0.2, var_cj = NULL, maxv = 2)
  var_cj_est[nrow(var_cj_est),]
  
  # find B
  var_Bi_est = find_opt_par(tot_var = 0.5, B = NULL, var_Ek = 0.2, var_cj = 0.2, S=S, D=D, maxv = 2, findB = TRUE)
  var_Bi_est$opt[nrow(var_Bi_est$opt),]
  Bnew = var_Bi_est$B
  tot_var = find_opt_par(B = Bnew, var_Ek = 0.2, var_cj = 0.2, maxv = 2)
  
  #get total variance
  tot_var = find_opt_par(B = B, var_Ek = 0.2, var_cj = 0.38, maxv = 2)
  
  ## without fixed B
  # find var_Ek
  var_Ek_est = find_opt_par(tot_var = 0.5, var_Ek = NULL, var_cj = 0.2, var_Bi = 0.15, S = S, D = D, maxv = 2)
  var_Ek_est[nrow(var_Ek_est),]
  
  # find var_cj
  var_cj_est = find_opt_par(tot_var = 0.5, var_Ek = 0.2, var_cj = NULL, var_Bi = 0.15, S = S, D = D, maxv = 2)
  var_cj_est[nrow(var_cj_est),]
  
  # find var_Bi
  var_Bi_est = find_opt_par(tot_var = 0.5, var_Ek = 0.2, var_cj = 0.2, var_Bi = NULL, S = S, D = D, maxv = 2)
  var_Bi_est[nrow(var_Bi_est),]
  
  # get total variance
  var_Bi_est = find_opt_par(var_Ek = 0.2, var_cj = 0.2, var_Bi = 0.33, B = NULL, S = S, D = D, maxv = 2)
}




if(FALSE) {
  #get distribution
  niter = 1e3
  totvar = iterate_totvar(niter = niter, var_Ek = 0.2, var_cj = 0.2, var_Bi = 0.2, S = 14, D = 3) 
  mean(totvar)
  
  hist(totvar, breaks = 20)
}




if(FALSE) {
  #check correctness
  
  d = 3   # dimensionality of E
  S = 14  # total number of species
  N0 = 200 # initial individuals per species
  M = N0*S # total number of initial occupied sites
  niter = 1000
  
  var_Bi = 0.13
  var_cj = 0.08
  var_Ek = 0.09
  
  var_saved = numeric(niter)
  var_est = numeric(S)
  var_est2 = numeric(S)
  var_saved_S = matrix(nrow = niter, ncol = S)
  var_est_S = matrix(nrow = niter, ncol = S)
  mean_saved_S = matrix(nrow = niter, ncol = S)
  
  var_in_spAijk_saved = numeric(niter)
  var_in_spAijk_est = numeric(niter)
    
  for(i in 1:niter) {
    # make environment
    E = rmvnorm(M, mean = rep(1/d, d), sigma = diag(d)*var_Ek)
    
    # make B
    B = rmvnorm(S, sigma = diag(d)*var_Bi)
    
    # make individuals
    cM = matrix(nrow = S, ncol = N0, rnorm(N0*S, sd = sqrt(var_cj)))
    
    # place individuals
    L = rep(1:S, each=N0)
    
    # calculate attributes
    Aijk = rowSums(B[L,]*E) + c(cM)
    
    # calculate variance
    var_saved[i] = var(Aijk)
    
    Z = diag(d)*var_Ek
    
    for(j in 1:S) {
      var_saved_S[i,j] = var(Aijk[L==j])
      var_est_S[i,j] = t(B[j,]) %*% Z %*% B[j,] + var_cj
      mean_saved_S[i,j] = mean(Aijk[L==j])
    }
    
    tmp = mean_saved_S[i,]
    var_in_spAijk_saved[i] = mean((mean(tmp)-tmp)^2)
    
    tmp = rowSums(B)/d
    var_in_spAijk_est[i] =  (mean((mean(tmp)-tmp)^2))
    
    var_est[i] = mean(var_est_S[i,]) + var_in_spAijk_est[i]
    
    
    #var_est2[i] = mean(diag(B %*% (diag(d)*var_Ek) %*% t(B)) + var_cj) +
    #  mean((rowSums(B)/d-mean(rowSums(B)/d))^2)
    
    
    var_est2[i] = get_tot_var(B, var_Ek, var_cj)
    
    if(i/100 == floor(i/100)) {
      print(i/niter)
    }
  }
  
  
  est_var_tot_mom = rowSums(1/S*(var_saved_S+mean_saved_S^2))-
    rowSums(1/S*mean_saved_S)^2
  plot(est_var_tot_mom, var_saved); abline(a=0, b=1)
  
  est_var_tot_mom = rowMeans(var_est_S)+
    var_in_spAijk_est
  plot(est_var_tot_mom, var_saved); abline(a=0, b=1)
  
  plot(var_est, var_saved); abline(a=0, b=1)
  
  plot(var_est2, var_saved); abline(a=0, b=1)
}
