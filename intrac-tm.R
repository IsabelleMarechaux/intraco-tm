# INTRACO theoretical model
source("att_to_perf_function.R")
source("compet_kernel.R")
source("get_total_var_analytical.R")

# =========================
# Scenario options 
# =========================

## which demographic rates ('performance') vary with the attributes A?
perf_var <- "fecundity"  # to be chosen among "mortality", "competitive_hierarchy", "fecundity"
## how the attribute vary across individuals ?
attribute_scenario <- 1 # to be chosen among 1, 2, 3 or 4 (eg. 1: species differences, no IV)
attriute_letter <- "a"
# options for scenario 1: Scen1a, Scen1b, Scen1c, Scen2, Scen3, Scen4
# see var_par_mat below for details

## how does the performance vary with the attribute ?
nonlin <- 1 # changed to new function of Adams att_to_perf_function 1 linear, -5 sublinear convex 5 suplinear concave (MUST BE EITHER BELOW 0 or ABOVE 1)
## is the total variance in all individuals' attribute constant across scenarios or not ? 
tot_variance <- "low" # to be chosen among: "low", "high", "variable"
## if perf_var="competitive_hierarchy", how competitive hierarchy translates into probability of recruitement?
compet_fun <- "lin"
compet_gen <- compet_lin # to be chosen among "compet_lin", "compet_logistic", "compet_step"

# save name of total scenario type
if(attribute_scenario!=1) {
  variance_senario <- paste("Scen", attribute_scenario, sep = "")
} else {
  variance_senario <- paste("Scen", attribute_scenario, attriute_letter, sep = "")
}


# =========================
# Global parameters
# =========================

gridsize <- 50
env_dim <- 2
n_species <- 20
n_ind_per_species <- 10
n_ind <- n_species * n_ind_per_species
n_cell <- gridsize * gridsize
n_iter <- 1500

# =========================
# Demographic parameters
# =========================

mortality_rate <- 0.05
recruitment_rate <- 0.5
recruitment_rate_max <- 5

# =========================
# Variability
# =========================
sp_beta_save = NULL # input a specific beta matrix here if desired

if(tot_variance == "variable") {
  var_beta = 0.5
  var_intra = 0.5
  var_env = 0.5
  
  source("get_vars_for_scenarios.R")
  
  if(variance_senario=="Scen1c") {
    sp_beta = sp_beta_1c
  } else {
    sp_beta = sp_beta_save
  }
  var_it = var_par_mat[variance_senario,"var_it"]
  var_beta = var_par_mat[variance_senario,"var_beta"]
  var_intra = var_par_mat[variance_senario,"var_intra"]
  var_env = var_par_mat[variance_senario,"var_env"]
  
} else {
  if(tot_variance == "low") {
    tot_var = 1
  } else {
    tot_var = 2
  }
  source("get_vars_for_scenarios.R")
  
  var_beta = var_par_mat[variance_senario,"var_beta"]
  var_intra = var_par_mat[variance_senario,"var_intra"]
  var_env = var_par_mat[variance_senario,"var_env"]
  var_it = var_par_mat[variance_senario,"var_it"]
  
  #sp_beta and sp_beta_1c have been automatically updated
  if(variance_senario=="Scen1c") {
    sp_beta = sp_beta_1c
  } else {
    sp_beta = sp_beta_save
  }
}

# =========================
# Outputs
# =========================

output_cells <- matrix(NA, n_cell, n_iter+1)
species_abundance <- matrix(0, n_species, n_iter+1)

# =========================
# Landscape
# =========================

# Set up the spatialized environement (which is now of dimension env_dim -- IM)
env_xpos <- env_ypos <- seq(0, 1, length=gridsize+1)
env_var <- array(rnorm(gridsize^2*env_dim, mean=0, sd=sqrt(var_env)), dim=c(gridsize, gridsize, env_dim)) # normal distribution (cf. doc)
# !! The id number of each cell in the lanscape are given by column (R style)
env_var_df <- matrix(c(env_var), ncol=env_dim)

# Vector of cells
cell_state <- rep(0, n_cell)

# Data-frame with cells and individuals characteristics
cell_char <- data.frame(cell_state=cell_state, env_var_df, attr_ind=0, perf_ind=-1)

# Positioning the individuals in the grid
id_cell_start <- sample(n_cell, n_ind, replace = FALSE)

# Fill the cells with individuals of a given species
cell_char$cell_state[id_cell_start] <- sample(rep(1:n_species, each=n_ind_per_species))


# Assigning attribute to individuals
cell_char$attr_ind[cell_char$cell_state!=0] <- apply((as.matrix(cell_char[cell_char$cell_state!=0, 2:(1+env_dim)]) + 1/env_dim) * sp_beta[cell_char$cell_state[cell_char$cell_state!=0], ], 1, sum) +
  rnorm(n_ind, mean=0, sd=sqrt(var_intra))
if(variance_senario == "Scen1b") {
  cell_char$attr_ind[cell_char$cell_state!=0] = cell_char$attr_ind[cell_char$cell_state!=0]+rnorm(sum(cell_char$cell_state!=0),0,sqrt(var_it))
}

# ALREADY DONE IN att_to_perf_function
# Computing performance from attribute
cell_char$perf_ind[cell_char$cell_state!=0] <-  att_to_perf_function(cell_char$attr_ind[cell_char$cell_state!=0], nonlin = nonlin)

# Outputs
# Species richness
output_cells[, 1] <- cell_char$cell_state
# Abundance
abundance <- as.data.frame(table(cell_char$cell_state))
names(abundance) <- c("sp", "ab")
abundance$sp <- as.numeric(as.character(abundance$sp))
if(sum(abundance$sp==0)>0) {
  abundance <- abundance[-which(abundance$sp==0), ]
}
species_abundance[abundance$sp, 1] <- abundance$ab 

# Iterations
for (iter in 1:n_iter) {
  # =========================
  # Epsilon (IT) variability: only works for scenario 1b
  # =========================
  if(variance_senario == "Scen1b") {
    cell_char$attr_ind[cell_char$cell_state!=0] = cell_char$attr_ind[cell_char$cell_state!=0]+rnorm(sum(cell_char$cell_state!=0),0,sqrt(var_it))
    cell_char$perf_ind[cell_char$cell_state!=0] <-  att_to_perf_function(cell_char$attr_ind[cell_char$cell_state!=0], nonlin = nonlin)
  }
  
  # =========================
  # Mortality
  # =========================
  
  # Which cell has an individual
  occupied_cells <- which(cell_char$cell_state != 0)
  n_occupied_cells <- length(occupied_cells)
  # Mortality envents
  if (perf_var=="mortality") {
    mort_events <- rbinom(n_occupied_cells, size=1, prob=mortality_rate*cell_char$perf_ind[occupied_cells])
  }else{
  mort_events <- rbinom(n_occupied_cells, size=1, prob=mortality_rate)
  }
  cell_char$cell_state[occupied_cells[mort_events==1]] <- 0
  cell_char$perf_ind[occupied_cells[mort_events==1]] <- -1
  
  # =========================
  # Recruitment
  # =========================
  if (perf_var=="fecundity") {
    
    # Number of new individuals per species
    n_recruit_sp <- as.data.frame(tapply(cell_char$perf_ind[which(cell_char$cell_state!=0)]*recruitment_rate_max, INDEX = as.factor(cell_char$cell_state[which(cell_char$cell_state!=0)]), FUN = sum))
    n_recruit_sp <- cbind(as.factor(levels(as.factor(cell_char$cell_state[which(cell_char$cell_state!=0)]))), n_recruit_sp)
    names(n_recruit_sp) <- c("sp", "fecundity")
    n_recruit_sp$n_recruit <- floor(n_recruit_sp$fecundity)
    
  } else {
    #Number of new individuals per species
    n_recruit_sp <- as.data.frame(table(cell_char$cell_state))
    names(n_recruit_sp) <- c("sp", "abundance")
    if(sum(n_recruit_sp$sp==0)>0) {
      n_recruit_sp = n_recruit_sp[n_recruit_sp$sp!=0,]
    }

    n_recruit_sp$n_recruit <- floor(n_recruit_sp$abundance*recruitment_rate)
  }
 
   
  # Data-frame of recruits
  n_tot_recruit <- sum(n_recruit_sp$n_recruit)
  recruit_char <- data.frame(id_recruit=1:n_tot_recruit, sp=rep(n_recruit_sp$sp, times=n_recruit_sp$n_recruit))
  recruit_char$pot_cell <- sample(1:n_cell, n_tot_recruit, replace=TRUE)
  recruit_sp <- as.numeric(as.character(recruit_char$sp))
  
  recruit_char$attr <- apply(as.matrix(cell_char[recruit_char$pot_cell, 2:(1+env_dim)]) * sp_beta[recruit_sp, ], 1, sum) + rnorm(n_tot_recruit, mean=0, sd=sqrt(var_intra))
  
  if(variance_senario == "Scen1b") {
    recruit_char$attr = recruit_char$attr+rnorm(length(recruit_char$attr),0,sqrt(var_it))
  }
  
  # Computing performance from attribute
  recruit_char$perf <-att_to_perf_function(recruit_char$attr, nonlin = nonlin)  
  
  # Looping on recruits
  
  if (perf_var =="competitive_hierarchy") {
    for (i in 1:n_tot_recruit) {
      if (cell_char$cell_state[recruit_char$pot_cell[i]]==0) {
        cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
        cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
        cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
      }
      if ((cell_char$cell_state[recruit_char$pot_cell[i]]!=0) & 
          runif(1) < compet_gen(recruit_char$perf[i] , cell_char$perf_ind[recruit_char$pot_cell[i]]))  {
       cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
       cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
       cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
      }
    }
  } else { # not competitive hierarchy
    tbltmp = table(recruit_char$pot_cell)
    if (length(tbltmp)==n_tot_recruit){ # no multiple new recruits per cell
      for (i in 1:n_tot_recruit) {
        if (cell_char$cell_state[recruit_char$pot_cell[i]]==0 | (cell_char$cell_state[recruit_char$pot_cell[i]]!=0 & runif(1)<0.5)) {
          cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
          cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
          cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
        }
      }
    } else {
      for (c in 1:length(which(tbltmp!=1))) {
        cell=as.integer(row.names(tbltmp))[which(tbltmp!=1)[c]]
        
        if (cell_char$cell_state[cell]==0 | (cell_char$cell_state[cell]!=0 & runif(1)>1/(length(which(recruit_char$pot_cell==cell))+1))) {
          winner=sample(which(recruit_char$pot_cell==cell), 1)
          cell_char$cell_state[cell] <- as.numeric(as.character(recruit_char$sp[winner]))
          cell_char$attr_ind[cell] <- recruit_char$attr[winner]
          cell_char$perf_ind[cell] <- recruit_char$perf[winner]
        }
      }
      for (c in 1:length(which(tbltmp==1))){
        cell=as.integer(row.names(tbltmp))[which(tbltmp==1)[c]]
        if (cell_char$cell_state[cell]==0 | (cell_char$cell_state[cell]!=0 & runif(1)<0.5)) {
          r=which(recruit_char$pot_cell==cell)
          cell_char$cell_state[cell] <- as.numeric(as.character(recruit_char$sp[r]))
          cell_char$attr_ind[cell] <- recruit_char$attr[r]
          cell_char$perf_ind[cell] <- recruit_char$perf[r]
        }
      }
   }
  }
  
  # Output
  # Species richness
  output_cells[, iter+1] <- cell_char$cell_state
  # Abundance
  abundance <- as.data.frame(table(cell_char$cell_state))
  names(abundance) <- c("sp", "ab")
  abundance$sp <- as.numeric(as.character(abundance$sp))
  if(sum(abundance$sp==0)>0) {
    abundance <- abundance[-which(abundance$sp==0), ]
  }
  species_abundance[abundance$sp, iter+1] <- abundance$ab 
  
  message(paste("progress:", round(iter/n_iter,3)),"\r",appendLF=FALSE)
}

# ===========================
# Graphics
# ===========================

# Species richness
species_richness <- apply(output_cells, 2, FUN=function(x){length(unique(x[x!=0]))})

#summary(species_abundance)

tot_number_ind=apply(species_abundance, 2, sum)

# Shannon index
prop <- species_abundance/tot_number_ind
shannon_index <- - apply(prop*log(prop), 2, sum, na.rm=TRUE)

par(mfrow=c(1,3))
plot(1:(n_iter+1), species_richness, pch=20, xlab="Iterations", ylab="Species richness", ylim=c(0,22), cex=0.5, las=1)
plot(1:(n_iter+1), shannon_index, pch=20, xlab="Iterations", ylab="shannon_index", ylim=c(0,8), cex=0.5, las=1)
plot(1:(n_iter+1), species_abundance[1,]/tot_number_ind, type="l", xlab="Iterations", ylab="Species relative abundance", ylim=c(0,1), las=1)
for (i in (2:n_species)) {
  par(new=T)
  plot(1:(n_iter+1), species_abundance[i,]/tot_number_ind, type="l", xlab="", ylab="", ylim=c(0,1), col=i, xaxt="n", yaxt="n")
}

sp_rank <- order(species_abundance[,(n_iter+1)], decreasing=T)
sp_rank[(species_richness[(n_iter+1)]+1):n_species] <- NA
rank<-rep(NA, n_species)
rank[sp_rank[1:species_richness[n_iter+1]]] <- 1:species_richness[n_iter+1]


# Species richness by scenario
results <- data.frame(matrix(NA, nrow=1, ncol=36))
results[1,] <- c(perf_var, attribute_scenario, nonlin, tot_variance, compet_fun,env_dim, 
             species_richness[floor(n_iter/5)], species_richness[floor(2*n_iter/5)], species_richness[floor(3*n_iter/5)], species_richness[floor(4*n_iter/5)], species_richness[floor(n_iter)],
             shannon_index[floor(n_iter/5)], shannon_index[floor(2*n_iter/5)], shannon_index[floor(3*n_iter/5)], shannon_index[floor(4*n_iter/5)], shannon_index[floor(n_iter)],
             rank)
colnames(results)[1:16]=c("dem_rate", "attr_scen", "nonlin", "tot_var", "compet_fun", "env_dim", 
                    "sp_rich1", "sp_rich2", "sp_rich3", "sp_rich4", "sp_rich5",
                    "shan1", "shan2", "shan3", "shan4", "shan5")

#=========================================================================================================================

# Recruits arriving on empty cells
list_empty_cell <- which(cell_char$cell_state==0)
id_recruit_by_default <- which(recruit_char$pot_cell %in% list_empty_cell)
cell_char$cell_state[recruit_char$pot_cell] 

