# INTRACO theoretical model

# Parameters
gridsize <- 100
env_dim <- 2
n_species <- 20
n_ind_per_species <- 10
n_ind <- n_species * n_ind_per_species
n_cell <- gridsize * gridsize
n_iter <- 100

# Demographic rates
mortality_rate <- 0.05
recruitment_rate <- 0.5

# Scenarios
attribute_scenario <- 1 # Scenario 1, 2, 3 or 4 (eg. 1: species differences, no IV)
nonlin <- 0 # Performance = (1-nonlin)*A + nonlin*A^2

# Interspecific variability
var_inter <- 1

# =========================
# Variables
# =========================

output_cells <- matrix(NA, n_cell, n_iter+1)
species_abundance <- matrix(0, n_species, n_iter+1)

# =========================
# Landscape
# =========================

# Set up the spatialized environement (which is now of dimension env_dim -- IM)
env_xpos <- env_ypos <- seq(0, 1, length=gridsize+1)
env_var <- array(rnorm(gridsize^2*env_dim, mean=0, sd=1), dim=c(gridsize, gridsize, env_dim)) # normal distribution (cf. doc)
# !! The id number of each cell in the lanscape are given by column (R style)
env_var_df <- matrix(c(env_var), ncol=env_dim)

# Vector of cells
cell_state <- rep(0, n_cell)

# Data-frame with cells and individuals characteristics
cell_char <- data.frame(cell_state=cell_state, env_var_df, attr_ind=0, perf_ind=0)

# Positioning the individuals in the grid
id_cell_start <- sample(n_cell, n_ind, replace = FALSE)

# Fill the cells with individuals of a given species
cell_char$cell_state[id_cell_start] <- sample(rep(1:n_species, each=n_ind_per_species))

# Assigning attribute to species
sp_attr <- rnorm(n_species, mean=0, sd=sqrt(var_inter))

# Assigning attribute to individuals
if (attribute_scenario==1) {
  cell_char$attr_ind[cell_char$cell_state!=0] <- sp_attr[cell_char$cell_state[cell_char$cell_state!=0]]
}

# Computing performance form attribute
cell_char$perf_ind <- (1-nonlin)*cell_char$attr_ind + nonlin*cell_char$attr_ind^2

# Output
# Species richness
output_cells[, 1] <- cell_char$cell_state
# Abundance
abundance <- as.data.frame(table(cell_char$cell_state))
names(abundance) <- c("sp", "ab")
abundance$sp <- as.numeric(as.character(abundance$sp))
abundance <- abundance[-which(abundance$sp==0), ]
species_abundance[abundance$sp, 1] <- abundance$ab 

# Iterations
for (iter in 1:n_iter) {
  
  # =========================
  # Mortality
  # =========================
  
  # Which cell has an individual
  occupied_cells <- which(cell_char$cell_state != 0)
  n_occupied_cells <- length(occupied_cells)
  # Mortality envents
  mort_events <- rbinom(n_occupied_cells, size=1, prob=mortality_rate)
  cell_char$cell_state[occupied_cells[mort_events==1]] <- 0
  # Setting zero for attribute and performance to dead individuals
  cell_char$attr_ind[occupied_cells[mort_events==1]] <- 0
  cell_char$perf_ind[occupied_cells[mort_events==1]] <- 0
  
  # =========================
  # Recruitment
  # =========================
  
  # Number of new individuals per species
  n_recruit_sp <- as.data.frame(table(cell_char$cell_state))[-1, ]
  names(n_recruit_sp) <- c("sp", "abundance")
  n_recruit_sp$n_recruit <- floor(n_recruit_sp$abundance*recruitment_rate)
  
  # Data-frame of recruits
  n_tot_recruit <- sum(n_recruit_sp$n_recruit)
  recruit_char <- data.frame(id_recruit=1:n_tot_recruit, sp=rep(n_recruit_sp$sp, times=n_recruit_sp$n_recruit))
  recruit_char$pot_cell <- sample(1:n_cell, n_tot_recruit, replace=TRUE) 
  if (attribute_scenario==1) {
    recruit_char$attr <- sp_attr[as.numeric(as.character(recruit_char$sp))]
  }
  # Computing performance form attribute
  recruit_char$perf <- (1-nonlin)*recruit_char$attr + nonlin*recruit_char$attr^2
  
  # Looping on recruits
  for (i in 1:n_tot_recruit) {
    if (cell_char$cell_state[recruit_char$pot_cell[i]]==0) {
      cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
      cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
      cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
    }
    if ( (cell_char$cell_state[recruit_char$pot_cell[i]]!=0) & (recruit_char$perf[i] > cell_char$perf_ind[recruit_char$pot_cell[i]]) ) {
      cell_char$cell_state[recruit_char$pot_cell[i]] <- as.numeric(as.character(recruit_char$sp[i]))
      cell_char$attr_ind[recruit_char$pot_cell[i]] <- recruit_char$attr[i]
      cell_char$perf_ind[recruit_char$pot_cell[i]] <- recruit_char$perf[i]
    }
  }
  
  # Output
  # Species richness
  output_cells[, iter+1] <- cell_char$cell_state
  # Abundance
  abundance <- as.data.frame(table(cell_char$cell_state))
  names(abundance) <- c("sp", "ab")
  abundance$sp <- as.numeric(as.character(abundance$sp))
  abundance <- abundance[-which(abundance$sp==0), ]
  species_abundance[abundance$sp, iter+1] <- abundance$ab 

}

# ===========================
# Graphics
# ===========================

# Species richness
species_richness <- apply(output_cells, 2, FUN=function(x){length(unique(x[x!=0]))})



#=========================================================================================================================

# Recruits arriving on empty cells
list_empty_cell <- which(cell_char$cell_state==0)
id_recruit_by_default <- which(recruit_char$pot_cell %in% list_empty_cell)
cell_char$cell_state[recruit_char$pot_cell] 

