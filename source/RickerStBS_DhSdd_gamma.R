# Model name: Gamma-alpha-demographic Ricker
# Model abbreviation: GAd
# Ricker with stochasticity in birth (B), survival (S), and between individual
# heterogeneity in DD survival (S).
# i.e. Density dependent demographic heterogeneity
#
RickerStBS_DhSdd_gamma <- function(Nt, R, alpha, kalpha) {
  if ( length(c(R,alpha,kalpha)) > 3 ) stop("R,alpha,kalpha must be scalar")
  births <- rpois( length(Nt), Nt * R ) #This includes DI mortality
# Heterogeneity in survival between individuals
  survivors <- NA*Nt
  for (x in 1:length(Nt)) {
    ai <- rgamma( births[x],shape=kalpha,scale=alpha/kalpha )
    ci <- exp(-ai*Nt[x])
  # Density dependent survival
    survivors[x] <- sum( rbinom( births[x], 1, ci ) )
  }
  return(survivors)
}





