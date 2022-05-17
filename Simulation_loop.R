# INTRACO theoretical model: simulation loop

rm(list = ls())
setwd("/Users/marechaux/OneDrive/INTRACO/Theoretical_model/workshop2")

source("get_total_var_analytical.R")

option_attribute <- c("Scen1a", "Scen1b", "Scen1c", "Scen2", "Scen3", "Scen4")
option_perfvar <- c("fecundity", "mortality", "competitive_hierarchy_lin", "competitive_hierarchy_step")
option_nonlin <- c( 5, -5) # default is one, would have been already explored
option_totvariance <- c( "high", "variable") #default is "low", would have been already explored

default_nonlin <- 1
default_perfvar <- "competitive_hierarchy_lin"
default_totvariance <- "low"

n_beta <-10
n_repet <- 5
env_dim <- 2
n_species <- 20

outputs <- data.frame()


####### loop of options for demographic rates that vary with attributes
## for each option, we test all 6 options of attribute variation, while keeping the nonlin and totalvariance options to default

for (p in 1:length(option_perfvar)) {
  
  perf_var <- option_perfvar[p]
  nonlin <- default_nonlin
  tot_variance <- default_totvariance
  
  compet_fun <- "-"
  if (perf_var=="competitive_hierarchy_lin") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "lin"
  } 
  if (perf_var=="competitive_hierarchy_step") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "step"
  } 
  
  for (br in 1:n_beta) {
  sp_beta_save = NULL
  
  for (a in 1:length(option_attribute)) {
    
    variance_scenario <- option_attribute[a]
  
    if(tot_variance == "variable") {
      var_beta = 0.5
      var_intra = 0.5
      var_env = 0.5
      
      source("get_vars_for_scenarios.R")
      
      if(variance_scenario=="Scen1c") {
        sp_beta = sp_beta_1c
      } else {
        sp_beta = sp_beta_save
      }
      var_it = var_par_mat[variance_scenario,"var_it"]
      var_beta = var_par_mat[variance_scenario,"var_beta"]
      var_intra = var_par_mat[variance_scenario,"var_intra"]
      var_env = var_par_mat[variance_scenario,"var_env"]
      totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
      
    } else {
      if(tot_variance == "low") {
        tot_var = 1
      } else {
        tot_var = 2
      }
      source("get_vars_for_scenarios.R")
      
      var_beta = var_par_mat[variance_scenario,"var_beta"]
      var_intra = var_par_mat[variance_scenario,"var_intra"]
      var_env = var_par_mat[variance_scenario,"var_env"]
      var_it = var_par_mat[variance_scenario,"var_it"]
      totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
      
      #sp_beta and sp_beta_1c have been automatically updated
      if(variance_scenario=="Scen1c") {
        sp_beta = sp_beta_1c
      } else {
        sp_beta = sp_beta_save
      }
    }
    
    if (variance_scenario=="Scen1a"){
      sp_beta_save <- sp_beta
    } ## so that for the other scenarios (except Scen1b), the same betas are used
    
    for (r in 1:n_repet) {
      
      output <- c(variance_scenario, perf_var, compet_fun, nonlin, tot_variance, var_beta, var_intra, var_env, var_it,totvar, br, r)
      for (sb in 1:20) {
        output <- c(output, sp_beta[sb, 1], sp_beta[sb, 2] )
      }
      outputs <- rbind(outputs, 
                       t(data.frame(output)))
      
    }
 
    
  }
  }
  

}

## this leads to 4 options_perfvar * 10 betas * 6 option_attribute * 5 repetitions each = 1200 simulations



####### loop over options for the shapes of the function that relates the attribute to performance
## for each option, we test all 6 options of attribute variation, while keeping the perfvar and totalvariance options to default

for (n in 1:length(option_nonlin)) {
  
  perf_var <- default_perfvar
  nonlin <- option_nonlin[n]
  tot_variance <- default_totvariance
  
  compet_fun <- "-"
  if (perf_var=="competitive_hierarchy_lin") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "lin"
  }
  if (perf_var=="competitive_hierarchy_step") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "step"
  }
  
  for (br in 1:n_beta) {
    sp_beta_save = NULL
    
    for (a in 1:length(option_attribute)) {
      
      variance_scenario <- option_attribute[a]
      
      if(tot_variance == "variable") {
        var_beta = 0.5
        var_intra = 0.5
        var_env = 0.5
        
        source("get_vars_for_scenarios.R")
        
        if(variance_scenario=="Scen1c") {
          sp_beta = sp_beta_1c
        } else {
          sp_beta = sp_beta_save
        }
        var_it = var_par_mat[variance_scenario,"var_it"]
        var_beta = var_par_mat[variance_scenario,"var_beta"]
        var_intra = var_par_mat[variance_scenario,"var_intra"]
        var_env = var_par_mat[variance_scenario,"var_env"]
        totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
        
      } else {
        if(tot_variance == "low") {
          tot_var = 1
        } else {
          tot_var = 2
        }
        source("get_vars_for_scenarios.R")
        
        var_beta = var_par_mat[variance_scenario,"var_beta"]
        var_intra = var_par_mat[variance_scenario,"var_intra"]
        var_env = var_par_mat[variance_scenario,"var_env"]
        var_it = var_par_mat[variance_scenario,"var_it"]
        totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
        
        #sp_beta and sp_beta_1c have been automatically updated
        if(variance_scenario=="Scen1c") {
          sp_beta = sp_beta_1c
        } else {
          sp_beta = sp_beta_save
        }
      }
      
      if (variance_scenario=="Scen1a"){
        sp_beta_save <- sp_beta
      }
      
      
      for (r in 1:n_repet) {
        
        output <- c(variance_scenario, perf_var, compet_fun, nonlin, tot_variance, var_beta, var_intra, var_env, var_it,totvar, br, r)
        for (sb in 1:20) {
          output <- c(output, sp_beta[sb, 1], sp_beta[sb, 2] )
        }
        outputs <- rbind(outputs, 
                         t(data.frame(output)))
        
      }
      
      
    }
  }
}

## this leads to 2 options_nonlin * 10 betas * 6 option_attribute * 5 repetitions each = 600 simulations


####### loop over options on the constraints on total variance
## for each option, we test all 6 options of attribute variation, while keeping the perfvar and nonlin options to default

for (v in 1:length(option_totvariance)) {
  
  perf_var <- default_perfvar
  nonlin <- default_nonlin
  tot_variance <- option_totvariance[v]
  
  compet_fun <- "-"
  if (perf_var=="competitive_hierarchy_lin") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "lin"
  }
  if (perf_var=="competitive_hierarchy_step") {
    perf_var <- "competitive_hierarchy"
    compet_fun <- "step"
  }
  
  for (br in 1:n_beta) {
    sp_beta_save = NULL
    
    for (a in 1:length(option_attribute)) {
      
      variance_scenario <- option_attribute[a]
      
      if(tot_variance == "variable") {
        var_beta = 0.5
        var_intra = 0.5
        var_env = 0.5
        
        source("get_vars_for_scenarios.R")
        
        if(variance_scenario=="Scen1c") {
          sp_beta = sp_beta_1c
        } else {
          sp_beta = sp_beta_save
        }
        var_it = var_par_mat[variance_scenario,"var_it"]
        var_beta = var_par_mat[variance_scenario,"var_beta"]
        var_intra = var_par_mat[variance_scenario,"var_intra"]
        var_env = var_par_mat[variance_scenario,"var_env"]
        totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
        
      } else {
        if(tot_variance == "low") {
          tot_var = 1
        } else {
          tot_var = 2
        }
        source("get_vars_for_scenarios.R")
        
        var_beta = var_par_mat[variance_scenario,"var_beta"]
        var_intra = var_par_mat[variance_scenario,"var_intra"]
        var_env = var_par_mat[variance_scenario,"var_env"]
        var_it = var_par_mat[variance_scenario,"var_it"]
        totvar = var_par_mat[variance_scenario, "TotVar_achieved"]
        
        #sp_beta and sp_beta_1c have been automatically updated
        if(variance_scenario=="Scen1c") {
          sp_beta = sp_beta_1c
        } else {
          sp_beta = sp_beta_save
        }
      }
      
      if (variance_scenario=="Scen1a"){
        sp_beta_save <- sp_beta
      }
      
      
      
      for (r in 1:n_repet) {
        
        output <- c(variance_scenario, perf_var, compet_fun, nonlin, tot_variance, var_beta, var_intra, var_env, var_it,totvar, br, r)
        for (sb in 1:20) {
          output <- c(output, sp_beta[sb, 1], sp_beta[sb, 2] )
        }
        outputs <- rbind(outputs, 
                         t(data.frame(output)))
        
      }
      
      
    }
  }
}

## this leads to 2 options_totvariance * 10 betas * 6 option_attribute * 5 repetitions each = 600 simulations
##### In total, this leads to 1200+ 600 + 600 = 2400 simulations


colnames(outputs) <- c("variance_scenario", "perf_var", "compet_fun", "nonlin", "tot_variance", "var_beta", "var_intra", "var_env", "var_it","totvar", "br", "r", 
                       "sp_beta1_1", "sp_beta1_2",
                       "sp_beta2_1", "sp_beta2_2",
                       "sp_beta3_1", "sp_beta3_2",
                       "sp_beta4_1", "sp_beta4_2",
                       "sp_beta5_1", "sp_beta5_2",
                       "sp_beta6_1", "sp_beta6_2",
                       "sp_beta7_1", "sp_beta7_2",
                       "sp_beta8_1", "sp_beta8_2",
                       "sp_beta9_1", "sp_beta9_2",
                       "sp_beta10_1", "sp_beta10_2",
                       "sp_beta11_1", "sp_beta11_2",
                       "sp_beta12_1", "sp_beta12_2",
                       "sp_beta13_1", "sp_beta13_2",
                       "sp_beta14_1", "sp_beta14_2",
                       "sp_beta15_1", "sp_beta15_2",
                       "sp_beta16_1", "sp_beta16_2",
                       "sp_beta17_1", "sp_beta18_2",
                       "sp_beta18_1", "sp_beta18_2",
                       "sp_beta19_1", "sp_beta19_2",
                       "sp_beta20_1", "sp_beta20_2")

summary(outputs)
nrow(outputs)

write.table(outputs, file="simulation_scenarios.txt", sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)


