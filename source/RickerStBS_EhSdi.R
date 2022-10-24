# Model name: Poisson-beta-environmental Ricker
# Model abbreviation: PBe
# Ricker with stochasticity in birth (B), survival (S), and between patch or
# time (ie environmental) heterogeneity in DI survival (S).
# i.e. Density independent environmental stochasticity
# Beta variation in survival results in a Poisson-beta distribution.
# m is the probability of mortality
# var_m is the variance in mortality/survival (will be the same)
#
RickerStBS_EhSdi <- function(Nt, R, m, alpha, var_m) {
  birthmn <- R/(1-m)
  births <- rpois( length(Nt), Nt * birthmn )
# Heterogeneity in survival between times or locations
  const <- m*(1-m)/var_m - 1 #shape1 + shape2 = constant
  sx <- rbeta(length(Nt),shape1=(1-m)*const,shape2=m*const)
# Density independent and density dependent survival
  survivors <- rbinom( length(Nt), births, sx * exp( -1 * alpha * Nt ) )
  return(survivors) #Poisson-beta
}




