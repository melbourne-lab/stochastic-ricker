# These functions calculate the variance and expected value of some extended
# stochastic Ricker models.

# Simulation is used - we don't yet know the variance formulas for these.

source("source/RickerStBS_EhSdi.r") #Poisson-beta (env)
source("source/RickerStBS_DhSdi.r") #Poisson-beta (dem)
source("source/RickerStBS_EhSdd_gamma.r") #gamma alpha (env)
source("source/RickerStBS_DhSdd_gamma.r") #gamma alpha (dem)


#----Variance and expected values by simulation for the Poisson-beta (env)
# Of course the expected value of N(t+1) is (alternatively) easily obtained, as
# it is just the deterministic Ricker.
# 
Ricker_poisbeta_e.var <- function(Nt,R,m,alpha,var_m,reps=1e6) {
  V <- NA*(1:length(Nt))
  E <- NA*(1:length(Nt))
  for (i in 1:length(Nt)){
    if (Nt[i]==0) { #If Nt=0, then p of any offspring = 0, and V,E=0.
      V[i] <- E[i] <- 0
      next
    }
    Ntp1 <- RickerStBS_EhSdi(rep(Nt[i],reps),R,m,alpha,var_m)
    V[i] <- var(Ntp1)
    E[i] <- mean(Ntp1)
  }
  return(list(var=V,E=E))
}

#----Variance and expected values by simulation for the Poisson-beta (dem)
#
Ricker_poisbeta_d.var <- function(Nt,R,m,alpha,var_m,reps=1e6) {
  V <- NA*(1:length(Nt))
  E <- NA*(1:length(Nt))
  for (i in 1:length(Nt)){
    if (Nt[i]==0) { #If Nt=0, then p of any offspring = 0, and V,E=0.
      V[i] <- E[i] <- 0
      next
    }
    Ntp1 <- RickerStBS_DhSdi(rep(Nt[i],reps),R,m,alpha,var_m)
    V[i] <- var(Ntp1)
    E[i] <- mean(Ntp1)
  }
  return(list(var=V,E=E))
}

#----Variance and expected values by simulation for the gamma alpha (env)
#
Ricker_gammaalpha_e.var <- function(Nt,R,alpha,kalpha,reps=1e6) {
  V <- NA*(1:length(Nt))
  E <- NA*(1:length(Nt))
  for (i in 1:length(Nt)){
    if (Nt[i]==0) { #If Nt=0, then p of any offspring = 0, and V,E=0.
      V[i] <- E[i] <- 0
      next
    }
    Ntp1 <- RickerStBS_EhSdd_gamma(rep(Nt[i],reps),R,alpha,kalpha)
    V[i] <- var(Ntp1)
    E[i] <- mean(Ntp1)
  }
  return(list(var=V,E=E))
}

#----Variance and expected values by simulation for the gamma alpha (dem)
#
Ricker_gammaalpha_d.var <- function(Nt,R,alpha,kalpha,reps=1e6) {
  V <- NA*(1:length(Nt))
  E <- NA*(1:length(Nt))
  for (i in 1:length(Nt)){
    if (Nt[i]==0) { #If Nt=0, then p of any offspring = 0, and V,E=0.
      V[i] <- E[i] <- 0
      next
    }
    Ntp1 <- RickerStBS_DhSdd_gamma(rep(Nt[i],reps),R,alpha,kalpha)
    V[i] <- var(Ntp1)
    E[i] <- mean(Ntp1)
  }
  return(list(var=V,E=E))
}