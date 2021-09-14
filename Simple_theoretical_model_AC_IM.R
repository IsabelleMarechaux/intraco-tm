#####################################################################
################ INTRACO workshop  ################################## 
####### A first simple theoretical model to start with ##############
#####################################################################


###### Adam Clark - sept 2020
###### Isabelle Maréchaux - dec 2020


error # prevent auto-run
rm(list = ls())

# =========================================
# tests (functions can be found below)
# =========================================

# n: number of individuals (zero-sum model)
# N: number of species
# niter: number of death events and replacement # CHANGE MADE IM.
# gridsize: the space is discretized in m*m square # CHANGE TO BE MADE (spatial autocorrelation).
# fractintr: the fraction of initial total trait variance (here always=1) that is intraspecific
# dcrowding: control the strength of crowding (competition) effect on mortality
# denv: control the strength of the effect of environemntal filtering
# r_crowding: interaction distance (beyond which individuals do not compete with each other)
# env_dim: number of independent environmental dimensions

#default values:
n=1000
N=100
niter=60000
gridsize=10
fracintr=0.01
dcrowding=0.0001
r_crowding=0.2 
env_dim=1
denv=1
useenv=TRUE


par(mfrow=c(3, 4))

simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.3, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.7, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.99, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)

simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=5, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=10, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=15, denv=1, useenv=TRUE)

simulation(n=1000, N=100, niter=60000, gridsize=2,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=10,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=20,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)
simulation(n=1000, N=100, niter=60000, gridsize=100,
           fracintr=0.01, dcrowding=0.0001, r_crowding=0.2, 
           env_dim=1, denv=1, useenv=TRUE)


# =========================================
# FUNCTIONS
# =========================================


### jaccard: compute the jaccard index (similarity of two ensembles: J=1 --> perfect similarity, intersection=union; J=0 --> intesection=0 ) 
### (used to assess the stability of coexistence between replicates)

jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

### tdist: compute the distances on a taurus
### (used to avoid edge effects due to finite landscape)

tdist <- function(xydata, xynew = NULL) { 
    if(is.null(xynew)) {
        xd = as.matrix(dist(xydata[,1]))
        yd = as.matrix(dist(xydata[,2]))
    } else {
        xd = abs(xynew[1] - xydata[,1])
        yd = abs(xynew[2] - xydata[,2])
    }
    xd[xd > 0.5] = 1 - xd[xd > 0.5]
    yd[yd > 0.5] = 1 - yd[yd > 0.5]
    
    distout = sqrt(xd^2 + yd^2)
    return(distout)
}


tdist_mnew <- function(xydata, xynew = NULL) { 
  if(is.null(xynew)) {
    xd = as.matrix(dist(xydata[,1]))
    yd = as.matrix(dist(xydata[,2]))
  } else {
    l=nrow(xynew)
    xd = abs(xynew[,1] - t(replicate(l, xydata[,1])))
    yd = abs(xynew[,2] - t(replicate(l, xydata[,2])))
  }
  xd[xd > 0.5] = 1 - xd[xd > 0.5]
  yd[yd > 0.5] = 1 - yd[yd > 0.5]
  
  distout = sqrt(xd^2 + yd^2)
  return(distout)
}

envdist <- function(trait, env) {
  diff=abs(trait-env)
  diff[diff>0.5]=1-diff[diff>0.5]
  dist=sqrt(sum(diff^2))
  return(dist)
}


### logit and ilogit
### (used to compute the probability of mortility)

logit<-function(x, ...) qlogis(x, ...)
ilogit<-function(x, ...) plogis(x, ...)


### simulation: the model itself, performs one simulation given a certain number of parameters

simulation <- function(n, N, niter, gridsize,fracintr, dcrowding, r_crowding, env_dim, denv=1, useenv=TRUE) {
  # n: number of individuals (zero-sum model)
  # N: number of species
  # niter: number of death events and replacement # CHANGE MADE IM.
  # gridsize: the space is discretized in m*m square # CHANGE TO BE MADE (spatial autocorrelation).
  # fractintr: the fraction of initial total trait variance (here always=1) that is intraspecific
  # dcrowding: control the strength of crowding (competition) effect on mortality
  # denv: control the strength of the effect of environemntal filtering
  # r_crowding: interaction distance (beyond which individuals do not compete with each other)
  # env_dim: number of independent environmental dimensions
    
    #Simulation_setup<-rbind(Simulation_setup, c(n, N, niter, gridsize,fracintr, dcrowding, r_crowding, env_dim, denv, useenv))
    
    tottrvar<-1 # total trait variance
    
    # set up the spatialized environement (which is now of dimension env_dim -- IM)
    env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
    #env_var<-array(rnorm(gridsize^2*env_dim, mean=0,sd=sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim)) # normal distributions to be discussed
    env_var=array(runif(gridsize^2*env_dim, min=0,max=1), dim=c(gridsize,gridsize,env_dim)) # uniform distributions to be tested
    
    # assign species and trait values to individuals
    #Ntrait <- matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, mean=0, sd=sqrt(tottrvar*(1-fracintr)))) # species traits
    Ntrait<- matrix(nrow=env_dim, ncol=N, data=runif(N*env_dim, min=0, max=1)) # uniform distributions to be tested
    spid <- sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
    trid <- matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid], sqrt(tottrvar*fracintr))) # individual-level traits 
    trid[which(trid<0)]<- trid[which(trid<0)]-floor(trid[which(trid<0)]) # if used with uniform distribution for the environment, use a taurus 
    trid[which(trid>1)]<- trid[which(trid>1)]-floor(trid[which(trid>1)]) # if used with uniform distribution for the environment, use a taurus 
    
    # assign locations to individuals
    xymat <- cbind(x = runif(n), y = runif(n))
    
    # get distances among individuals
    distmat <- tdist(xymat)
    is_neighborhood<- ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath <- (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath <- (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath<-(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath<-(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
        gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
        )
        env_eff=NULL
        for (j in 1:n) { # would be great to replace this for loop
            #env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
            env_eff=c(env_eff, envdist(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])*denv)
        }
        ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
        ptot=pdeath
    }
    
    #outputs
    nsmp <- 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    Nsp<-NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
        
        # save status ==> output
        if(ncount%%nsmp==0) {
            #datout[ncount%/%nsmp+1,2] = length(unique(spid))
            Nsp=c(Nsp, length(unique(spid)))
        }
        
        # which individual dies?
        nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE -- IM
        dead_neighborhood= is_neighborhood[nextdeath,]
        
        # replace the individual with a new one
        spid[nextdeath] = sample(spid[-nextdeath], 1)
        xymat[nextdeath,] = runif(2)
        
        newdist <- tdist(xymat, xymat[nextdeath,])
        distmat[nextdeath,] <- newdist
        distmat[,nextdeath] <- newdist
        new_neighborhood<-ifelse(newdist<r_crowding  & newdist>0, 1, NA)
        is_neighborhood[nextdeath,]<- new_neighborhood
        is_neighborhood[,nextdeath]<- new_neighborhood
        
        trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
        trid[which(trid[,nextdeath]<0),nextdeath]=trid[which(trid[,nextdeath]<0),nextdeath]-floor(trid[which(trid[,nextdeath]<0),nextdeath]) # if used with uniform distribution for the environment, use a taurus 
        trid[which(trid[,nextdeath]>1),nextdeath]=trid[which(trid[,nextdeath]>1),nextdeath]-floor(trid[which(trid[,nextdeath]>1),nextdeath])
        
        pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))

        # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
        updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
        if (length(which(updated_neighborhood==1))>0) {
            if (length(which(updated_neighborhood==1))==1) {
                pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
            } else {
                pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
            }
        }
        
        # Compute the environmental effect on the probability of dying of the newly recruited individual
        if(useenv) {
            gridloc[nextdeath,] = c(
            x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
            y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
            )
            #env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
            env_eff[nextdeath]= envdist(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])*denv
            
            # update the probabilities of dying
            ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
            if (length(which(updated_neighborhood==1))>0) {
                ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
            }
        } else {
            ptot[nextdeath] = pdeath[nextdeath]
            if (length(which(updated_neighborhood==1))>0) {
                ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
            }
        }
        
        ncount = ncount + 1
    }
    
    plot(Nsp, type = "l", ylim=c(0, N),
    main = paste("UNIF",paste(paste("fracintr =", fracintr), paste("gridsize =", gridsize), paste("envdim =", env_dim))))
    abline(h = 0, lty = 3)
    #Diversity_dynamics<<-rbind(Diversity_dynamics, Nsp)
    
}



###Test 24/02/2021: several trees can die at each iteration (independent Bernouilli draw for each tree)

simulation <- function(n, N, niter, gridsize,fracintr, dcrowding, r_crowding, env_dim, denv=1, useenv=TRUE) {
  # n: number of individuals (zero-sum model)
  # N: number of species
  # niter: number of death events and replacement # CHANGE MADE IM.
  # gridsize: the space is discretized in m*m square # CHANGE TO BE MADE (spatial autocorrelation).
  # fractintr: the fraction of initial total trait variance (here always=1) that is intraspecific
  # dcrowding: control the strength of crowding (competition) effect on mortality
  # denv: control the strength of the effect of environemntal filtering
  # r_crowding: interaction distance (beyond which individuals do not compete with each other)
  # env_dim: number of independent environmental dimensions
  
  #Simulation_setup<-rbind(Simulation_setup, c(n, N, niter, gridsize,fracintr, dcrowding, r_crowding, env_dim, denv, useenv))
  
  tottrvar<-1 # total trait variance
  
  # set up the spatialized environement (which is now of dimension env_dim -- IM)
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  #env_var<-array(rnorm(gridsize^2*env_dim, mean=0,sd=sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim)) # normal distributions to be discussed
  env_var=array(runif(gridsize^2*env_dim, min=0,max=1), dim=c(gridsize,gridsize,env_dim)) # uniform distributions to be tested
  
  # assign species and trait values to individuals
  #Ntrait <- matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, mean=0, sd=sqrt(tottrvar*(1-fracintr)))) # species traits
  Ntrait<- matrix(nrow=env_dim, ncol=N, data=runif(N*env_dim, min=0, max=1)) # uniform distributions to be tested
  spid <- sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid <- matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid], sqrt(tottrvar*fracintr))) # individual-level traits 
  trid[which(trid<0)]<- trid[which(trid<0)]-floor(trid[which(trid<0)]) # if used with uniform distribution for the environment, use a taurus 
  trid[which(trid>1)]<- trid[which(trid>1)]-floor(trid[which(trid>1)]) # if used with uniform distribution for the environment, use a taurus 
  
  # assign locations to individuals
  xymat <- cbind(x = runif(n), y = runif(n))
  
  # get distances among individuals
  distmat <- tdist(xymat)
  is_neighborhood<- ifelse(distmat<r_crowding & distmat>0, 1, NA)
  
  # calculate each individual probability to die
  ## effect of crowding
  #pdeath <- (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
  #pdeath <- (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
  #pdeath<-(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
  pdeath<-(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
  ## effect of environemental filtering
  if(useenv) {
    gridloc = cbind( # in which cell is each individual?
      x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
      y = as.numeric(cut(xymat[,2], breaks = env_ypos))
    )
    env_eff=NULL
    for (j in 1:n) { # would be great to replace this for loop
      #env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
      env_eff=c(env_eff, envdist(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])*denv)
    }
    ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
  } else {
    ptot=pdeath
  }
  
  #outputs
  nsmp <- 100  # sampling interval
  #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
  Nsp<-NULL
  
  #run simulation
  ncount = 0
  while(ncount < niter) {
    
    # save status ==> output
    if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      Nsp=c(Nsp, length(unique(spid)))
    }
    
    # which individuals dies? (to accelerate the dynamics, I here pick 10 individuals at each iteration)
    nextdeath=which(rmultinom(1,10, ptot)==1) ## NEW CHANGE -- IM
    dead_neighborhood <- ifelse(apply(is_neighborhood[nextdeath,], 2, sum, na.rm=TRUE)>0, 1, NA)
    
    # replace the individuals with new ones
    spid[nextdeath] <- sample(spid[-nextdeath], length(nextdeath))
    xymat[nextdeath,] <- cbind(x = runif(length(nextdeath)), y = runif(length(nextdeath)))
    
    newdist <- tdist_mnew(xymat, xymat[nextdeath,])
    distmat[nextdeath,] <- newdist
    distmat[,nextdeath] <- newdist
    new_neighborhood <- ifelse(newdist<r_crowding  & newdist>0, 1, NA)
    is_neighborhood[nextdeath,]<- new_neighborhood
    is_neighborhood[,nextdeath]<- new_neighborhood
    new_neighborhood<-ifelse(apply(new_neighborhood, 2, sum, na.rm=TRUE)>0, 1, NA)
    
    trid[,nextdeath] = rnorm(length(nextdeath)*env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
    trid[,nextdeath][which(trid[,nextdeath]<0)] <- trid[,nextdeath][which(trid[,nextdeath]<0)]-floor(trid[,nextdeath][which(trid[,nextdeath]<0)]) # if used with uniform distribution for the environment, use a taurus 
    trid[,nextdeath][which(trid[,nextdeath]>1)] <- trid[,nextdeath][which(trid[,nextdeath]>1)]-floor(trid[,nextdeath][which(trid[,nextdeath]>1)])
    
    pdeath[nextdeath] <- (1 - exp(-dcrowding*rowSums(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
    
    # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
    updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
    if (length(which(updated_neighborhood==1))>0) {
      if (length(which(updated_neighborhood==1))==1) {
        pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
      } else {
        pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
      }
    }
    
    # Compute the environmental effect on the probability of dying of the newly recruited individuals
    if(useenv) {
      gridloc[nextdeath,] = c(
        x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
      )
      #env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
      for (l in 1:length(nextdeath)) {
      env_eff[nextdeath[l]]= envdist(trid[,nextdeath[l]], env_var[gridloc[nextdeath[l],1], gridloc[nextdeath[l],2],])*denv
      }
      
      # update the probabilities of dying
      ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
      if (length(which(updated_neighborhood==1))>0) {
        ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
      }
    } else {
      ptot[nextdeath] <- pdeath[nextdeath]
      if (length(which(updated_neighborhood==1))>0) {
        ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
      }
    }
    
    ncount = ncount + 1
  }
  
  plot(Nsp, type = "l", ylim=c(0, N),
       main = paste("UNIF",paste(paste("fracintr =", fracintr), paste("gridsize =", gridsize), paste("envdim =", env_dim))))
  abline(h = 0, lty = 3)
  #Diversity_dynamics<<-rbind(Diversity_dynamics, Nsp)
  
}








### Test 9/12/2020 with variation in fractintr, env_dim and gridsize, and 10 replicates for each parameter combination to investigate the stability of the "coexistence"

f=c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99)
ed=c(1,2,5,10,15,20)
g=c(2,5,10,20,30,40)

Jaccard_f=NULL
Jaccard_ed=NULL
Jaccard_g=NULL
nr=10 # number of replicates with same initial conditions
Nsp_f=matrix(nrow = nr, ncol=length(f), data=NA)
Nsp_ed=matrix(nrow = nr, ncol=length(ed), data=NA)
Nsp_g=matrix(nrow = nr, ncol=length(g), data=NA)

n=1000
N=100
niter=60000
dcrowding=0.0001
r_crowding=0.2 
denv=1 
useenv=TRUE
tottrvar=1 # total trait variance

### Variation of fracintr

for (i in 1:length(f)) {
  
  fracintr=f[i]
  gridsize=10
  env_dim=1
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  
  # assign initial species and trait values to individuals
  Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
  
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
        #datout[ncount%/%nsmp+1,2] = length(unique(spid))
        #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_f[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_f=c(Jaccard_f, mean(jac))
  
}


### Variation of environmental dimensions

for (i in 1:length(ed)) {
  
  fracintr=0.01
  gridsize=10
  env_dim=ed[i]
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  
  # assign initial species and trait values to individuals
  Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
    
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_ed[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_ed=c(Jaccard_ed, mean(jac))
  print(mean(jac))
  
}


### Variation of gridsize

for (i in 1:length(g)) {
  
  fracintr=0.01
  gridsize=g[i]
  env_dim=1
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  
  # assign initial species and trait values to individuals
  Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
    
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_g[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_g=c(Jaccard_g, mean(jac))
  print(mean(jac))
  
}


setwd("/Users/marechaux/OneDrive/INTRACO/Theoretical_model")
write.table(Jaccard_f, file = "Jaccard_f.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)
write.table(Jaccard_ed, file = "Jaccard_ed.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)
write.table(Jaccard_g, file = "Jaccard_g.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)
write.table(Nsp_f, file = "Nsp_f.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)
write.table(Nsp_ed, file = "Nsp_ed.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)
write.table(Nsp_g, file = "Nsp_g.txt", quote = FALSE, dec=".", sep = "\t", na = "-", row.names=T)




### Test 18/02/2021 with variation in fractintr, env_dim and gridsize, and 10 replicates for each parameter combination to investigate the stability of the "coexistence"
### with uniform distributions for environmental space.

f=c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99)
ed=c(1,2,5,10,15,20)
g=c(2,5,10,20,30,40)

Jaccard_f=NULL
Jaccard_ed=NULL
Jaccard_g=NULL
nr=10 # number of replicates with same initial conditions
Nsp_f=matrix(nrow = nr, ncol=length(f), data=NA)
Nsp_ed=matrix(nrow = nr, ncol=length(ed), data=NA)
Nsp_g=matrix(nrow = nr, ncol=length(g), data=NA)

n=1000
N=100
niter=60000
dcrowding=0.0001
r_crowding=0.2 
denv=1 
useenv=TRUE
tottrvar=1 # total trait variance

### Variation of fracintr

for (i in 1:length(f)) {
  
  fracintr=f[i]
  gridsize=10
  env_dim=1
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  #env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  env_var=array(runif(gridsize^2*env_dim, min=0,max=1), dim=c(gridsize,gridsize,env_dim)) # uniform distributions to be tested
  
  # assign initial species and trait values to individuals
  #Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  Ntrait<- matrix(nrow=env_dim, ncol=N, data=runif(N*env_dim, min=0, max=1)) # uniform distributions to be tested
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  trid_init[which(trid_init<0)]<- trid_init[which(trid_init<0)]-floor(trid_init[which(trid_init<0)]) # if used with uniform distribution for the environment, use a taurus 
  trid_init[which(trid_init>1)]<- floor(trid_init[which(trid_init>1)]) # if used with uniform distribution for the environment, use a taurus 
  
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
    
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        #env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
        env_eff=c(env_eff, envdist(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      trid[which(trid[,nextdeath]<0),nextdeath]=trid[which(trid[,nextdeath]<0),nextdeath]-floor(trid[which(trid[,nextdeath]<0),nextdeath]) # if used with uniform distribution for the environment, use a taurus 
      trid[which(trid[,nextdeath]<1),nextdeath]=floor(trid[which(trid[,nextdeath]<1),nextdeath])
      
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        #env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        env_eff[nextdeath]= envdist(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_f[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_f=c(Jaccard_f, mean(jac))
  
}


### Variation of environmental dimensions

for (i in 1:length(ed)) {
  
  fracintr=0.01
  gridsize=10
  env_dim=ed[i]
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  #env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  env_var=array(runif(gridsize^2*env_dim, min=0,max=1), dim=c(gridsize,gridsize,env_dim)) # uniform distributions to be tested
  
  # assign initial species and trait values to individuals
  #Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  Ntrait<- matrix(nrow=env_dim, ncol=N, data=runif(N*env_dim, min=0, max=1)) # uniform distributions to be tested
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  trid_init[which(trid_init<0)]<- trid_init[which(trid_init<0)]-floor(trid_init[which(trid_init<0)]) # if used with uniform distribution for the environment, use a taurus 
  trid_init[which(trid_init>1)]<- floor(trid_init[which(trid_init>1)]) # if used with uniform distribution for the environment, use a taurus 
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
    
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        #env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
        env_eff=c(env_eff, envdist(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      trid[which(trid[,nextdeath]<0),nextdeath]=trid[which(trid[,nextdeath]<0),nextdeath]-floor(trid[which(trid[,nextdeath]<0),nextdeath]) # if used with uniform distribution for the environment, use a taurus 
      trid[which(trid[,nextdeath]<1),nextdeath]=floor(trid[which(trid[,nextdeath]<1),nextdeath])
      
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        #env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        env_eff[nextdeath]= envdist(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_ed[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_ed=c(Jaccard_ed, mean(jac))
  print(mean(jac))
  
}


### Variation of gridsize

for (i in 1:length(g)) {
  
  fracintr=0.01
  gridsize=g[i]
  env_dim=1
  
  spid_replicate=matrix(nrow=nr, ncol=n, data=NA)
  
  # set up the spatialized environement
  env_xpos = env_ypos = seq(0, 1, length=gridsize+1)
  env_var=array(rnorm(gridsize^2*env_dim, 0,sqrt(tottrvar)), dim=c(gridsize,gridsize,env_dim))
  
  # assign initial species and trait values to individuals
  Ntrait = matrix(nrow=env_dim, ncol=N, data=rnorm(N*env_dim, 0, sqrt(tottrvar*(1-fracintr)))) # species traits
  spid_init = sample(rep(1:N, n/N), n, replace = FALSE) # assign species to individual (works when n is a multiple of N) # CHANGE MADE so that species are equally abundant at the beginning of the simulation
  trid_init = matrix(nrow=env_dim, ncol= n, data=rnorm(n*env_dim, Ntrait[,spid_init], sqrt(tottrvar*fracintr))) # individual-level traits
  
  # assign initial locations to individuals
  xymat_init = cbind(x = runif(n), y = runif(n))
  
  for (r in 1:nr) {
    
    spid=spid_init
    trid=trid_init
    xymat=xymat_init
    
    # get distances among individuals
    distmat = tdist(xymat)
    is_neighborhood= ifelse(distmat<r_crowding & distmat>0, 1, NA)
    
    # calculate each individual probability to die
    ## effect of crowding
    #pdeath = (1 - exp(-dcrowding*1/rowSums(distmat))) # effect of crowding on individual mortality probability # CHANGE MADE (see line below), test to be made
    #pdeath = (1 - exp(-dcrowding*rowSums(distmat))) # pb in this term, as the probability to die is closer to 1 when the sum of distance increases (and hence crowding decreases...)
    #pdeath=(1 - exp(-dcrowding*rowSums(is_crowding))) # simply counts the number of individuals within r_crowding
    pdeath=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood*distmat), na.rm=TRUE))) # sum of the inverse of distances with individual within the neighborhood
    ## effect of environemental filtering
    if(useenv) {
      gridloc = cbind( # in which cell is each individual?
        x = as.numeric(cut(xymat[,1], breaks = env_xpos)),
        y = as.numeric(cut(xymat[,2], breaks = env_ypos))
      )
      env_eff=NULL
      for (j in 1:n) { # would be great to replace such for loop
        env_eff=c(env_eff, as.numeric(dist(rbind(trid[,j], env_var[gridloc[j,1], gridloc[j,2],])))*denv)
      }
      ptot=ilogit(logit(pdeath) + env_eff) # shape to be discussed
    } else {
      ptot=pdeath
    }
    
    #outputs
    #nsmp = 100  # sampling interval
    #datout = data.frame(iter = seq(0, niter, by = nsmp), Nsp = NA)
    #Nsp=NULL
    
    #run simulation
    ncount = 0
    while(ncount < niter) {
      
      # save status ==> output
      #if(ncount%%nsmp==0) {
      #datout[ncount%/%nsmp+1,2] = length(unique(spid))
      #Nsp=c(Nsp, length(unique(spid)))
      #}
      
      # which individual dies?
      nextdeath = which(rmultinom(1,1, ptot)==1) # CHANGE MADE
      dead_neighborhood= is_neighborhood[nextdeath,]
      
      # replace the individual with a new one
      spid[nextdeath] = sample(spid[-nextdeath], 1)
      xymat[nextdeath,] = runif(2)
      
      newdist = tdist(xymat, xymat[nextdeath,])
      distmat[nextdeath,] = newdist
      distmat[,nextdeath] = newdist
      new_neighborhood=ifelse(newdist<r_crowding  & newdist>0, 1, NA)
      is_neighborhood[nextdeath,]= new_neighborhood
      is_neighborhood[,nextdeath]= new_neighborhood
      
      trid[,nextdeath] = rnorm(env_dim, Ntrait[,spid[nextdeath]], sqrt(tottrvar*fracintr))
      
      pdeath[nextdeath] = (1 - exp(-dcrowding*sum(1/(is_neighborhood[nextdeath,]*distmat[nextdeath,]), na.rm=TRUE)))
      
      # update the crowding effect on the probabilities of dying of individuals within the neighborhood of the dead or newly recruited individuals
      updated_neighborhood=ifelse(dead_neighborhood==1 | new_neighborhood==1, 1, NA)
      if (length(which(updated_neighborhood==1))>0) {
        if (length(which(updated_neighborhood==1))==1) {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*sum(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        } else {
          pdeath[which(updated_neighborhood==1)]=(1 - exp(-dcrowding*rowSums(1/(is_neighborhood[which(updated_neighborhood==1),]*distmat[which(updated_neighborhood==1),]), na.rm=TRUE)))
        }
      }
      
      # Compute the environmental effect on the probability of dying of the newly recruited individual
      if(useenv) {
        gridloc[nextdeath,] = c(
          x = as.numeric(cut(xymat[nextdeath,1], breaks = env_xpos)),
          y = as.numeric(cut(xymat[nextdeath,2], breaks = env_ypos))
        )
        env_eff[nextdeath]= as.numeric(dist(rbind(trid[,nextdeath], env_var[gridloc[nextdeath,1], gridloc[nextdeath,2],])))*denv
        
        # update the probabilities of dying 
        ptot[nextdeath] =  ilogit(logit(pdeath[nextdeath]) + env_eff[nextdeath]) # for the newly recruited individual
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=ilogit(logit(pdeath[which(updated_neighborhood==1)]) + env_eff[which(updated_neighborhood==1)])
        }
      } else {
        ptot[nextdeath] = pdeath[nextdeath]
        if (length(which(updated_neighborhood==1))>0) {
          ptot[which(updated_neighborhood==1)]=pdeath[which(updated_neighborhood==1)]
        }
      }
      
      ncount = ncount + 1
    }
    
    Nsp_g[r,i]=length(unique(spid))
    spid_replicate[r,]=spid
    print(length(unique(spid)))
    
  }
  
  jac=NULL
  for (a in 1:(nr-1)) {
    for (b in (a+1):nr) {
      jac=c(jac, jaccard(unique(spid_replicate[a,]), unique(spid_replicate[b,])))
    }
  }
  
  Jaccard_g=c(Jaccard_g, mean(jac))
  print(mean(jac))
  
}





















