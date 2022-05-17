# =============================#
# ==== Community dynamics =====#
# =============================#

rm(list = ls())

#directory_reading <- "/Users/marechaux/OneDrive/INTRACO/Theoretical_model/workshop2"
#directory_writing <- "/Users/marechaux/OneDrive/INTRACO/Theoretical_model/workshop2"

directory_reading <- "/home/marechauxi/INTRACO"
directory_writing <- "/lustre/marechauxi" # for the cluster

setwd(directory_reading) 

source("att_to_perf_function.R")
source("compet_kernel.R")

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

simulation_scenarios <- read.table(file=paste(directory_reading, "simulation_scenarios.txt", sep="/"), header=TRUE, dec=".", sep="\t")
#s <- args[1] # for cluster
#s <- 1

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# coerce the value to an integer
s <- as.integer(slurm_arrayid)


#message(paste("simulation:", s),"\r",appendLF=FALSE)

perf_var <- simulation_scenarios$perf_var[s]
compet_fun <- simulation_scenarios$compet_fun[s]
variance_scenario <- simulation_scenarios$variance_scenario[s]
var_intra <- simulation_scenarios$var_intra[s]
var_env <- simulation_scenarios$var_env[s]
var_it <- simulation_scenarios$var_it[s]
var_it <- simulation_scenarios$var_it[s]
sp_beta <- matrix(nrow=n_species, ncol=env_dim)
for (ss in 1:n_species) {
  for (d in 1:env_dim) {
    sp_beta[ss,d] <- simulation_scenarios[s,13+env_dim*(ss-1)+d-1]
  }
}
nonlin <- simulation_scenarios$nonlin[s]




# =========================
# Demographic parameters
# =========================

mortality_rate <- 0.05
recruitment_rate <- 0.5
recruitment_rate_max <- 5

#output_cells <- matrix(NA, n_cell, n_iter+1)
species_abundance <- matrix(0, n_species, n_iter+1)
species_richness=NULL
shannon_index=NULL

if (perf_var =="competitive_hierarchy") {
  if (compet_fun=="lin") {
    compet_gen <- compet_lin
  }
  if (compet_fun=="step") {
    compet_gen <- compet_step
  }
}

# =========================
# Landscape
# =========================

# Set up the spatialized environement (which is now of dimension env_dim -- IM)
env_xpos <- env_ypos <- seq(0, 1, length=gridsize+1)
env_var <- array(rnorm(gridsize^2*env_dim, mean=0, sd=sqrt(var_env)), dim=c(gridsize, gridsize, env_dim)) # normal distribution (cf. doc)
# !! The id number of each cell in the landscape are given by column (R style)
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
if(variance_scenario == "Scen1b") {
  cell_char$attr_ind[cell_char$cell_state!=0] = cell_char$attr_ind[cell_char$cell_state!=0]+rnorm(n_ind,0,sqrt(var_it))
}

# Computing performance from attribute
cell_char$perf_ind[cell_char$cell_state!=0] <-  att_to_perf_function(cell_char$attr_ind[cell_char$cell_state!=0], nonlin = nonlin)

# Outputs
# Species richness
#output_cells[, 1] <- cell_char$cell_state
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
  //if(variance_scenario == "Scen1b") {
    //cell_char$attr_ind[cell_char$cell_state!=0] = cell_char$attr_ind[cell_char$cell_state!=0]+rnorm(sum(cell_char$cell_state!=0),0,sqrt(var_it))
    //cell_char$perf_ind[cell_char$cell_state!=0] <-  att_to_perf_function(cell_char$attr_ind[cell_char$cell_state!=0], nonlin = nonlin)
  //}
  
  # =========================
  # Mortality
  # =========================
  
  # Which cell has an individual
  occupied_cells <- which(cell_char$cell_state != 0)
  n_occupied_cells <- length(occupied_cells)
  # Mortality events
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
  
  if(variance_scenario == "Scen1b") {
    recruit_char$attr = recruit_char$attr+rnorm(length(recruit_char$attr),0,sqrt(var_it))
  }
  
  # Computing performance from attribute
  recruit_char$perf <-att_to_perf_function(recruit_char$attr, nonlin = nonlin)  
  
  # Looping on recruits
  
  if (perf_var =="competitive_hierarchy") {
    for (i in 1:n_tot_recruit) {
      if (cell_char$cell_state[recruit_char$pot_cell[i]]==0) { # if the potential cell is empty, then the individual is recruited
        cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
        cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
        cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
      }
      if ((cell_char$cell_state[recruit_char$pot_cell[i]]!=0) &  # if the potential cell is occupied, then use the competition function to decide if it wins the site or not
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
      for (c in 1:length(which(tbltmp!=1))) { # loop over the cells are in competition for recruitement among new recruits
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
  
  #### Outputs
  
  # Species richness
  species_richness=c(species_richness, length(unique(cell_char$cell_state[cell_char$cell_state!=0])))
  
  # Abundance
  abundance <- as.data.frame(table(cell_char$cell_state))
  names(abundance) <- c("sp", "ab")
  abundance$sp <- as.numeric(as.character(abundance$sp))
  if(sum(abundance$sp==0)>0) {
    abundance <- abundance[-which(abundance$sp==0), ]
  }
  species_abundance[abundance$sp, iter+1] <- abundance$ab 
  
  #Shanon_index
  tot_ind=sum(abundance$ab)
  prop=abundance$ab/tot_ind
  shannon <- - sum(prop*log(prop))
  shannon_index=c(shannon_index, shannon)
  
  #message(paste("progress:", round(iter/n_iter,3)),"\r",appendLF=FALSE)
}


write.table(species_abundance, file=paste( directory_writing, paste(paste("species_abundance", as.character(s), sep="_"), ".txt", sep=""), sep="/"), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(species_richness, file=paste( directory_writing, paste(paste("species_richness", as.character(s), sep="_"), ".txt", sep=""), sep="/"), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(shannon_index, file=paste(directory_writing, paste(paste("shannon_index", as.character(s), sep="_"), ".txt", sep=""), sep="/"), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)




