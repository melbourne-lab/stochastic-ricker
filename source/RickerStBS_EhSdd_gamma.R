# Model name: Gamma-alpha-environmental Ricker
# Model abbreviation: GAe
# Ricker with stochasticity in birth (B), survival (S), and between patch or
# time (ie environmental) heterogeneity in DD survival (S).
# i.e. Density dependent environmental stochasticity
# Gamma variation in alpha.
# Final distribution is unknown.
# kalpha: shape parameter of the gamma for alpha
#
RickerStBS_EhSdd_gamma <- function(Nt, R, alpha, kalpha) {
  births <- rpois( length(Nt), Nt * R ) #This includes DI mortality
# Heterogeneity in survival between times or locations
  ax <- rgamma( length(Nt),shape=kalpha,scale=alpha/kalpha )
  cx <- exp(-ax*Nt)
  survivors <- rbinom( length(Nt), births, cx )
  return(survivors)
}

