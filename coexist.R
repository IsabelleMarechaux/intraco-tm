## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr
## web             :https://ecology.ghislainv.fr
## license         :GPLv3
## ==============================================================================

# Libraries
library(raster)
library(here)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(glue)
library(geoR) # variog()
library(plot3D) # scatter3D()
library(Rcpp)
library(RcppArmadillo)

# Create output directories
dir.create(here("outputs/m0"), recursive=TRUE)

# Seed for reproducibility
seed <- 1234

# Figure width
fig_width <- 16.6 # in cm

# logit/inv_logit functions
logit <- function(x, min=0, max=1) {
  p <- (x-min)/(max-min)
  return(log(p/(1-p)))
}

inv_logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse( is.na(p) & !is.na(x), 1, p ) # fix problems with +Inf
  return(p * (max-min) + min)
}

# Cpp function to compute distance between Sites and Species
Rcpp::sourceCpp(here("_src", "dist_Site_Sp.cpp"))

# =========================
# Landscape and environment
# =========================

# Landscape
nsite_side <- 25
mat <- matrix(0, nrow=nsite_side, ncol=nsite_side)
r <- raster(mat, crs="+proj=utm +zone=1")
nsite <- ncell(r)
coords <- coordinates(r)

# Environment on each site
Sites <- data.frame(site=1:nsite, V1_env=NA, V2_env=NA, V3_env=NA)

# Neighbourhood matrix
neighbors.mat <- adjacent(r, cells=c(1:nsite), directions=8,
                          pairs=TRUE, sorted=TRUE)
# Number of neighbours by site
n.neighbors <- as.data.frame(table(as.factor(neighbors.mat[,1])))[,2]
# Adjacent sites
adj <- neighbors.mat[,2]
# Generate symmetric adjacency matrix, A
A <- matrix(0,nsite,nsite)
index.start <- 1
for (i in 1:nsite) {
  index.end <- index.start+n.neighbors[i]-1
  A[i,adj[c(index.start:index.end)]] <- 1
  index.start <- index.end+1
}

# Function to draw in a multivariate normal
rmvn <- function(n, mu=0, V=matrix(1), seed=1234) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) {
    stop("Dimension problem!")
  }
  D <- chol(V)
  set.seed(seed)
  t(matrix(rnorm(n*p),ncol=p)%*%D+rep(mu,rep(n,p)))
}

# Generate spatial random effects
Vrho.target <- 1 # Variance of spatial random effects
d <- 1 # Spatial dependence parameter = 1 for intrinsic CAR
Q <- diag(n.neighbors)-d*A + diag(.0001,nsite) # Add small constant to make Q non-singular
covrho <- Vrho.target*solve(Q) # Covariance of rhos

# Number of axis for the niche
n_axis <- 3

# Environment on each site
sites <- data.frame(V1_env=rep(NA, nsite), V2_env=NA, V3_env=NA)
env <- list()
for (i in 1:n_axis) {
  seed <- 1234 + i - 1
  rho <- c(rmvn(1, mu=rep(0, nsite), V=covrho, seed=seed)) # Spatial Random Effects
  rho <- rho-mean(rho) # Centering rhos on zero
  rho <- scales::rescale(rho, to=c(0, 1))
  sites[,i] <- rho
  env[[i]] <- matrix(rho, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
}
env_stack <- stack(raster(env[[1]]*255), raster(env[[2]]*255), raster(env[[3]]*255))
crs(env_stack) <- "+proj=utm +zone=1"

# Plot
png(file=here("outputs", "m0", "environment.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2))
plot(raster(env[[1]]), main="Environment var1", col=topo.colors(255))
plot(raster(env[[2]]), main="Environment var2", col=topo.colors(255))
plot(raster(env[[3]]), main="Environment var3", col=topo.colors(255))
plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
dev.off()

# Habitat frequency
png(file=here("outputs", "m0", "hab_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(env[[1]], main="", xlab="Environment var1")    
dev.off()

# =========================================
# Species niche
# =========================================

# Niche width
niche_width <- 0.25
# Niches per axis
n_niche <- 1/niche_width
# Number of species
nsp <- n_niche^n_axis 
# Species coordinates on the three niche axis (x, y, z)
base_coord <- seq(0, 1, length.out=n_niche+1)[-(n_niche+1)]+niche_width/2
sp_x <- rep(rep(base_coord, n_niche), n_niche)
sp_y <- rep(rep(base_coord, each=n_niche), n_niche)
sp_z <- rep(base_coord, each=n_niche^2)
niche_optimum <- as.data.frame(cbind(sp_x, sp_y, sp_z))

# Random optimum for species
randomOptSp <- TRUE
if (randomOptSp) {
  set.seed(seed)
  niche_optimum <- data.frame(sp_x=runif(nsp), sp_y=runif(nsp), sp_z=runif(nsp))
  niche_optimum <- niche_optimum[order(niche_optimum$sp_z, niche_optimum$sp_y, niche_optimum$sp_x), ]
}

# Plot the species niche
png(file=here("outputs", "m0", "species_niche.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
par(mar=c(1,1,2,2))
scatter3D(niche_optimum$sp_x, niche_optimum$sp_y, niche_optimum$sp_z,
          pch=16, 
          colvar=1:64, col=rainbow(nsp),
          bty = "f", main ="Three-dimensional niche", phi=0,
          xlim=c(0,1), ylim=c(0,1), zlim=c(0,1))
dev.off()

# Matrix of species performance on each site (distances)
# Sites in rows, Species in columns
dist_E_Sp <- dist_Site_Sp(as.matrix(sites), as.matrix(niche_optimum))
dprim_E_Sp <- (dist_E_Sp-mean(dist_E_Sp))/sd(dist_E_Sp)
perf_E_Sp <- -dprim_E_Sp

# Function to identify the species with the highest performance
high_perf_sp <- function(dist, sp_pres) {
  dist_pres <- dist[sp_pres]
  min_dist <- min(dist_pres)
  sp_high_perf <- sp_pres[which(dist_pres==min_dist)]
  # If more than one species, selection at random
  if (length(sp_high_perf)>1) {
    # Random permutation
    sp_high_perf <- sample(sp_high_perf)
    sp_high_perf <- sp_high_perf[1]
  }
  return(sp_high_perf)
}

# Probability of dying of each species on each site
# Strength of unsuitability
b <- -0.5
mortality_E_Sp <- inv_logit(logit(0.1) + b * perf_E_Sp)
# Mortality rate distribution
png(file=here("outputs", "m0", "hist_mortality.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
hist(mortality_E_Sp)
dev.off()

# Habitat frequency for each species
rank_dist_E <- t(apply(dist_E_Sp, 1, rank, ties.method="min"))
sp_hab_freq <- apply(rank_dist_E, 2, function(x){sum(x==1)})
sp_hab_freq <- as.table(sp_hab_freq)
names(sp_hab_freq) <- 1:nsp
png(file=here("outputs", "m0", "species_habitat_freq.png"),
    width=fig_width, height=fig_width, units="cm", res=300)
plot(sp_hab_freq, xlab="Species", ylab="Habitat frequency")
dev.off()

# Species with no habitat
sp_no_habitat <- as.vector(which(sp_hab_freq==0))
nsp_no_habitat <- length(sp_no_habitat)

# =========================================
# Repetitions
# =========================================

# Number of repetitions
nrep <- 10
# Number of generations
ngen <- 100

# Species richness
sp_rich <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Species rank at the end of the generations
rank_sp <- matrix(NA, nrow=nrep, ncol=nsp)
# Environmental filtering
env_filt <- matrix(NA, nrow=ngen+1, ncol=nrep)
# Mean mortality rate in the community
theta_comm <- matrix(NA, nrow=ngen+1, ncol=nrep)

# Loop on repetitions
for (r in 1:nrep) {

  # -----------------------------------------
  # Initial conditions
  # -----------------------------------------
  
  # Draw species at random in the landscape (one individual per site)
  sp <- sample(1:nsp, size=nsite, replace=TRUE)
  #hist(sp)
  community_start <- matrix(sp, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
  if (r==1) {
    png(file=here("outputs", "m0", "community_start.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    plot(raster(community_start), main="Species - Start", zlim=c(0, nsp),
         col=c("black", rainbow(nsp)))
    dev.off()
  }
  community <- community_start
  
  # Species richness
  sp_rich[1, r] <- length(unique(c(community)))
  
  # Abundances
  abund <- matrix(NA, ncol=nsp, nrow=ngen+1)
  abund[1,] <- table(factor(c(community), levels=1:nsp))
  
  # Environmental filtering
  dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
  env_filt[1, r] <- mean(dist_site)
  
  # Mean mortality rate
  theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
  theta_comm[1, r] <- mean(theta_site) 

  # -----------------------------------------
  # Dynamics
  # -----------------------------------------
  
  # Simulating generation
  for (g in 1:ngen) {
    
    # ******************
    # Mortality
    # ******************
    
    # Mortality rate on each site
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    
    # Mortality events
    mort_ev <- rbinom(nsite, size=1, prob=theta_site)
    mortality <- matrix(mort_ev, nrow=nsite_side, ncol=nsite_side, byrow=TRUE)
    
    # Update community
    community[mortality==1] <- 0
    # Plot once
    if (r==1 & g==1) {
      png(file=here("outputs", "m0", "mortality_events.png"),
          width=fig_width, height=fig_width, units="cm", res=300)
      plot(raster(community), main="Species - with vacant sites", zlim=c(0, nsp),
           col=c("black", rainbow(nsp)))
      dev.off()
    }

    # *********************
    # Fecundity/Recruitment
    # *********************
    
    # Species present in the community
    sp_present <- sort(unique(community[community!=0]))
    nsp_present <- length(sp_present)
    
    # Vacant sites
    community_rast <- raster(community)
    sites_vacant <- which(values(community_rast)==0)
    
    # Performance of species on vacant sites
    dist_E_Sp_vacant <- dist_E_Sp[sites_vacant, ]
    
    # Identify the species with the highest performance
    new_ind <- apply(dist_E_Sp_vacant, 1, high_perf_sp, sp_pres=sp_present)
    
    # Recruitment
    community_rast[sites_vacant] <- new_ind
    community <- as.matrix(community_rast)
    
    # *********************
    # Diversity
    # *********************
    
    # Species richness
    sp_rich[g+1, r] <- length(unique(as.vector(community)))
    abund[g+1, ] <- table(factor(as.vector(community), levels=1:nsp))
    
    # Environmental filtering
    dist_site <- diag(dist_E_Sp[, as.vector(t(community))])
    env_filt[g+1, r] <- mean(dist_site)

    # Mean mortality rate in the community
    theta_site <- diag(mortality_E_Sp[, as.vector(t(community))])
    theta_comm[g+1, r] <- mean(theta_site)
    
  } # End ngen
  
  # Plot final community once
  if (r==1) {
    png(file=here("outputs", "m0", "community_end.png"),
        width=fig_width, height=fig_width, units="cm", res=300)
    plot(raster(community), main=glue("Species - End (ngen={ngen})"),
         zlim=c(0, nsp), col=c("black", rainbow(nsp)))
    dev.off()
  }
  
  # Species rank
  rank_sp[r, ] <- rank(-abund[ngen, ], ties.method="min")
  
} # End nrep

# =========================
# Diversity analysis
# =========================

sp_rich
rank_sp

# ---------------------------------------------
# Plot species richness
# ---------------------------------------------

sp_rich <- data.frame(sp_rich)
sp_rich_long <- sp_rich %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="sp_rich")
p <- ggplot(data=sp_rich_long, aes(x=gen, y=sp_rich, col=rep)) +
  geom_line() + 
  xlab("Generations") + 
  ylab("Species richness")
ggsave(p, filename=here("outputs", "m0", "species_richness_with_time.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ---------------------------------------------
# Link between final rank and habitat frequency
# ---------------------------------------------

# Mean final rank
sp_mean_rank <- apply(rank_sp, 2, mean)
# Plot
df <- data.frame(cbind(sp_mean_rank, sp_hab_freq))
p <- ggplot(data=df, aes(x=sp_hab_freq, y=sp_mean_rank)) +
  geom_point() +
  geom_smooth(method="gam", formula=y~s(x, bs = "cs"), color="red", fill="#69b3a2", se=TRUE) +
  xlab("Species habitat frequency") +
  ylab("Species mean rank (higher rank = lower abundance)") +
  theme(axis.title=element_text(size=16))
ggsave(p, filename=here("outputs", "m0", "mean_rank-habitat_freq.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# ---------------------------------------------
# Environmental filtering
# ---------------------------------------------

env_filt <- data.frame(env_filt)
env_filt_long <- env_filt %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="env_filt")
p <- ggplot(data=env_filt_long, aes(x=gen, y=env_filt, col=rep)) +
  geom_line() +
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean env-species perf difference")
ggsave(p, filename=here("outputs", "m0", "environmental_filtering.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# Plot
png(file=here("outputs", "m0", "spatial_comp_env_sp.png"), 
    width=fig_width, height=fig_width, units="cm", res=300)
par(mfrow=c(2,2))
plot(raster(community_start), main="Species - Start", zlim=c(0, nsp),
     col=c("black", rainbow(nsp)))
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", rainbow(nsp)))
plotRGB(env_stack, main="Environment RGB", axes=TRUE, margins=TRUE)
plot(raster(community), main="Species - End", zlim=c(0, nsp),
     col=c("black", rainbow(nsp)))
dev.off()

# ---------------------------------------------
# Theta community
# ---------------------------------------------

theta_comm <- data.frame(theta_comm)
theta_comm_long <- theta_comm %>%
  mutate(gen=1:(ngen+1)) %>%
  pivot_longer(cols=X1:X10, names_to="rep",
               names_prefix="X", values_to="theta_comm")
p <- ggplot(data=theta_comm_long, aes(x=gen, y=theta_comm, col=rep)) +
  geom_line() +
  labs(title="Environmental filtering") +
  xlab("Generations") + 
  ylab("Mean mortality rate in the community")
ggsave(p, filename=here("outputs", "m0", "mortality_rate_community.png"),
       width=fig_width, height=fig_width/2, units="cm", dpi=300)

# ----------------------------------
# Spatial autocorrelation of species
# ----------------------------------

# Species autocorrelation
sp_XY <- data.frame(rasterToPoints(raster(community)))
names(sp_XY) <- c("x", "y", "sp")
vario_sp <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=sp_XY$sp)
# Environment autocorrelation
# 3D voxel for each site
x_site <- pmin(floor(sites$V1_env/niche_width)+1, 4)
y_site <- pmin(floor(sites$V2_env/niche_width)+1, 4)
z_site <- pmin(floor(sites$V3_env/niche_width)+1, 4)
class_site <- (z_site-1)*n_niche^2+(y_site-1)*n_niche+(x_site-1)+1
vario_env <- variog(coords=cbind(sp_XY$x, sp_XY$y), data=class_site)
# Plot with correlation
png(file=here("outputs", "m0", "sp_autocorrelation.png"),
    width=fig_width, height=fig_width*0.8, units="cm", res=300)
par(mfrow=c(2,2))
plot(vario_sp, main="Species - End")
plot(vario_env, main="Environment")
plot(vario_env$v, vario_sp$v,
     xlab="Semivariance for environment",
     ylab="Semivariance for species")
m <- lm(vario_sp$v ~ vario_env$v-1)
abline(a=0, b=coef(m), col="red")
dev.off()

## Conclusions

# 1. Because mortality rate is equal for all species, species with the lowest habitat frequency have a higher probability to disappear
# 2. Stable coexistence (cf. species rank from one repetition to the other)
# 3. Species abundance at the end correlated with species habitat frequency
# 4. Species autocorrelation with correspondence between species and environment
# 5. For a given habitat distribution, the number of species at the equilibrium depends on the mortality rate
#    If the mortality rate is higher, a higher number of low abundance species disappear.   

# ===========================
# "Intraspecific variability"
# ===========================

# Data-set
df <- data.frame(perf_E_Sp)
names(df) <- sprintf("Sp_%03d", 1:nsp)
df_perf <- tibble(df) %>%
  mutate(Env=values(raster(env[[1]]))) %>%
  mutate(Env2=Env^2) %>%
  pivot_longer(cols=Sp_001:glue("Sp_0{nsp}"), names_to="Species", values_to="Perf")

# Observed niche
# Select 8 species at random
set.seed(seed)
sp_sel <- sample(unique(df_perf$Species), 9, replace=FALSE)
df_sp_sel <- df_perf %>% filter(Species %in% sp_sel)
p <- ggplot(data=df_sp_sel, aes(x=Env, y=Perf)) +
  geom_point() +
  geom_smooth(method="lm", formula=y~poly(x,2), se=TRUE) +
  facet_wrap(vars(Species), nrow=3) +
  xlab("Environment (first axis)") +
  ylab("Performance")
ggsave(p, filename=here("outputs", "m0", "infering_species_niche.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)
  
# Observed intraspecific variability
lm_fit <- lm(Perf~Species+Species*Env+Species*Env2, data=df_perf)
V_intra <- df_perf %>%
  mutate(res=lm_fit$residuals) %>%
  group_by(Species) %>%
  summarise(V=var(res))
p <- ggplot(data=V_intra, aes(x=Species, y=V)) +
  geom_col() +
  theme(axis.text.x=element_text(angle=90, size=6)) +
  ylab("Intraspecific variance")
ggsave(p, filename=here("outputs", "m0", "intraspecific_variance.png"),
       width=fig_width, height=fig_width, units="cm", dpi=300)

# =========================
# End of file
# =========================