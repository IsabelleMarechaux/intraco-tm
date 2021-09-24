if(tot_variance == "variable") {
  # matrix for storing vars
  var_par_mat = matrix(nrow = 6, ncol = 6)
  row.names(var_par_mat) = c("Scen1a", "Scen1b", "Scen1c", "Scen2", "Scen3", "Scen4")
  colnames(var_par_mat) = c("var_beta", "var_intra", "var_env", "var_it", "TotVar", "TotVar_achieved")
  var_par_mat[,"TotVar"] = NA
  
  var_par_mat[,"var_beta"] = var_beta
  
  if(is.null(sp_beta_save)) {
    sp_beta <- matrix(rnorm(n_species*env_dim, mean=0, sd=sqrt(var_beta)), ncol=env_dim)
  }
  var_par_mat["Scen1a","TotVar_achieved"] = find_opt_par(B = sp_beta, var_Ek = 0, var_cj = 0)
  
  var_par_mat[c("Scen1a", "Scen1b", "Scen1c", "Scen3"),"var_intra"] = 0
  var_par_mat[c("Scen1a", "Scen1b", "Scen1c", "Scen2"),"var_env"] = 0
  var_par_mat[c("Scen1a", "Scen1c", "Scen2", "Scen3", "Scen4"),"var_it"] = 0
  var_par_mat[c("Scen1b"),"var_it"] = var_par_mat["Scen1a","TotVar_achieved"]/2
  
  # find B that meets twice the achieved tot_var
  var_Bi_est <- find_opt_par(tot_var = var_par_mat["Scen1a","TotVar_achieved"]*2, B = NULL, var_Ek = 0, var_cj = 0, S=n_species, D=env_dim, maxv = 2*var_par_mat["Scen1a","TotVar_achieved"]*2, findB = TRUE)
  sp_beta_1c <- var_Bi_est$B
  var_beta_1c <- tail(var_Bi_est$opt$estout,1)
  var_par_mat[c("Scen1c"),"var_beta"] = var_beta_1c
  
  # scenario 2
  var_par_mat["Scen2","var_intra"] = var_intra
  
  # scenario 3
  var_par_mat["Scen3","var_env"] = var_env
  
  # scenario 4
  var_par_mat["Scen4","var_intra"] = var_intra
  var_par_mat["Scen4","var_env"] = var_env
  
  #test all
  for(i in 1:nrow(var_par_mat)) {
    if(row.names(var_par_mat)[i] == "Scen1c") {
      Buse = sp_beta_1c
    } else {
      Buse = sp_beta
    }
    
    var_par_mat[i, "TotVar_achieved"] = find_opt_par(B = Buse, var_Ek = var_par_mat[i,"var_env"], var_cj = var_par_mat[i,"var_intra"])
    
    if(row.names(var_par_mat)[i] == "Scen1b") {
      var_par_mat[i, "TotVar_achieved"] = var_par_mat[i, "TotVar_achieved"] + var_par_mat[c("Scen1b"),"var_it"]
    }
  }
} else {
  # matrix for storing vars
  var_par_mat = matrix(nrow = 6, ncol = 6)
  row.names(var_par_mat) = c("Scen1a", "Scen1b", "Scen1c", "Scen2", "Scen3", "Scen4")
  colnames(var_par_mat) = c("var_beta", "var_intra", "var_env", "var_it", "TotVar", "TotVar_achieved")
  var_par_mat[,"TotVar"] = c(tot_var/2, rep(tot_var, 5))
  
  var_par_mat[c("Scen1a", "Scen1b", "Scen1c", "Scen3"),"var_intra"] = 0
  var_par_mat[c("Scen1a", "Scen1b", "Scen1c", "Scen2"),"var_env"] = 0
  var_par_mat[c("Scen1a", "Scen1c", "Scen2", "Scen3", "Scen4"),"var_it"] = 0
  var_par_mat[c("Scen1b"),"var_it"] = tot_var/2
  
  # find B that meets desired tot_var/2
  if(is.null(sp_beta_save)) {
    var_Bi_est <- find_opt_par(tot_var = tot_var/2, B = NULL, var_Ek = 0, var_cj = 0, S=n_species, D=env_dim, maxv = 2*tot_var, findB = TRUE)
    sp_beta <- var_Bi_est$B
  } else {
    sp_beta = sp_beta_save
  }
  var_beta = var(c(sp_beta))
  var_par_mat[c("Scen1a", "Scen1b", "Scen2", "Scen3", "Scen4"),"var_beta"] = var_beta
  
    
  # find B that meets desired tot_var
  var_Bi_est <- find_opt_par(tot_var = tot_var, B = NULL, var_Ek = 0, var_cj = 0, S=n_species, D=env_dim, maxv = 2*tot_var, findB = TRUE)
  sp_beta_1c <- var_Bi_est$B
  var_beta_1c <- tail(var_Bi_est$opt$estout,1)
  var_par_mat[c("Scen1c"),"var_beta"] = var_beta_1c
  
  # scenario 2
  var_par_mat["Scen2","var_intra"] = tot_var/2
  
  # scenario 3
  var_Ek_est = find_opt_par(tot_var = tot_var, B = sp_beta, var_Ek = NULL, var_cj = 0, maxv = 2*tot_var)
  var_par_mat["Scen3","var_env"] = tail(var_Ek_est$estout, 1)
  
  # scenario 4
  var_par_mat["Scen4","var_intra"] = tot_var/4
  
  var_Ek_est = find_opt_par(tot_var = tot_var, B = sp_beta, var_Ek = NULL, var_cj = var_par_mat["Scen4","var_intra"], maxv = 2*tot_var)
  var_par_mat["Scen4","var_env"] = tail(var_Ek_est$estout, 1)
  
  #test all
  for(i in 1:nrow(var_par_mat)) {
    if(row.names(var_par_mat)[i] == "Scen1c") {
      Buse = sp_beta_1c
    } else {
      Buse = sp_beta
    }
    
    var_par_mat[i, "TotVar_achieved"] = find_opt_par(B = Buse, var_Ek = var_par_mat[i,"var_env"], var_cj = var_par_mat[i,"var_intra"])
    
    if(row.names(var_par_mat)[i] == "Scen1b") {
      var_par_mat[i, "TotVar_achieved"] = var_par_mat[i, "TotVar_achieved"] + var_par_mat[c("Scen1b"),"var_it"]
    }
  }
}
if(is.null(sp_beta_save)) {
  sp_beta_save = sp_beta
}

