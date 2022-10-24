# Model name: Poisson-beta-demographic Ricker
# Model abbreviation: PBd
# Ricker with stochasticity in birth (B), survival (S), and between individual
# heterogeneity in DI survival (S).
# i.e. Density independent demographic heterogeneity
# Beta variation in survival results in a Poisson-beta distribution?
# It is a sum of Poisson-betas - not sure if that's P-beta
# m is the probability of mortality
# var_m is the variance in mortality/survival (will be the same)
#
# Current version is an IBM but if a sum of PBs is also PB, this can be reduced
# (eliminating the for loop).
# 
RickerStBS_DhSdi <- function(Nt, R, m, alpha, var_m) {
  if ( length(c(R,m,alpha,var_m)) > 4 ) stop("R,m,alpha,var_m must be scalar")
  birthmn <- R/(1-m)
  births <- rpois( length(Nt), Nt * birthmn )
# Heterogeneity in survival between individuals
  survivors <- NA*Nt
  const <- m*(1-m)/var_m - 1 #shape1 + shape2 = constant
  for (x in 1:length(Nt)) {
    si <- rbeta(births[x],shape1=(1-m)*const,shape2=m*const) #beta heterogeneity
  # Density independent and density dependent survival
    survivors[x] <- sum( rbinom( births[x], 1, si * exp( -1 * alpha * Nt[x] ) ) )
  }
  return(survivors)
}

