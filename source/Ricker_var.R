# These functions calculate the variance of the stochastic Ricker models. The
# variances are calculated directly from the additive variance components.

####----Variances-----------------------------------------


#----Variance of the Poisson model
#
Ricker_pois.var <- function(Nt,R,alpha) {
  return(Ricker(Nt,R,alpha))
}

#----Variance of the Negative binomial (demographic)
#
Ricker_nbinom_d.var <- function(Nt,R,alpha,kD) {
  mu <- Ricker(Nt,R,alpha)
  return( mu + (mu*mu)/(Nt*kD) )
}

#----Variance of the Negative binomial (environmental)
#
Ricker_nbinom_e.var <- function(Nt,R,alpha,kE) {
  mu <- Ricker(Nt,R,alpha)
  return( mu + mu*mu/kE )
}

#----Variance of the Negative binomial-gamma
#
Ricker_nbinomgamma.var <- function(Nt,R,alpha,kD,kE) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  dhvar <- Ricker(Nt,R,alpha)^2 / (kD*Nt)   #dem hetero variance
  evar <- Ricker(Nt,R,alpha)^2 / kE         #env variance
  return( poisvar + dhvar + evar )
}

#----Variance of the Poisson-binomial
#
Ricker_poisbinom.var <- function(Nt,R,alpha) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  return( poisvar + sexvar )
}

#----Variance of the NB-binomial (demographic)
#
Ricker_nbinombinom_d.var <- function(Nt,R,alpha,kD) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  dhvar <- Ricker(Nt,R,alpha)^2 / (kD*Nt)   #raw dhvar for no-sex model
  return( poisvar + sexvar + 2*dhvar )
}

#----Variance of the NB-binomial (environmental)
#
Ricker_nbinombinom_e.var <- function(Nt,R,alpha,kE) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  evar <- Ricker(Nt,R,alpha)^2 / kE         #env variance
  return( poisvar+sexvar+evar )
}

#----Variance of the NB-binomial-gamma
#
Ricker_nbinombinomgamma.var <- function(Nt,R,alpha,kD,kE) {
  poisvar <- Ricker(Nt,R,alpha)             #Poisson variance
  sexvar <- Ricker(Nt,R,alpha)^2 / Nt       #sex variance
  dhvar <- Ricker(Nt,R,alpha)^2 / (kD*Nt)   #raw dhvar for no-sex model
  evar <- Ricker(Nt,R,alpha)^2 / kE         #env variance
  return( poisvar + sexvar + 2*dhvar + evar )
}
