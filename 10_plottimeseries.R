# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922

# Simulate time series of different stochastic Ricker models

source("source/Ricker.r")             #Deterministic Ricker
source("source/RickerStBS.r")         #Stochastic birth and survival
source("source/RickerStBS_DhB.r")     #Stochastic b, s, + dem het in b.
source("source/RickerStBS_EhB.r")     #Stochastic b, s, + env het in b.
source("source/RickerStBS_DEhB.r")    #Stochastic b, s, + dem & env het in b.
source("source/RickerStSexBS.r")      #Stochastic sex ratio, birth, survival
source("source/RickerStSexBS_DhB.r")  #St sr, b, s, + dem het in b.
source("source/RickerStSexBS_EhB.r")  #St sr, b, s, + env het in b.
source("source/RickerStSexBS_DEhB.r") #St sr, b, s, + dem & env het in b.

#Simulate one time series
R <- 2.6
alpha <- 0.00373
kE <- 17.62
kD <- 1.07

t <- 0:200
N <- t*NA
N[1] <- 20
for (i in 1:max(t)) {
#  N[i+1] <- Ricker(N[i],R,alpha)
#  N[i+1] <- RickerStBS(N[i],R,alpha)
#  N[i+1] <- RickerStBS_DhB(N[i],R,alpha,kD)
#  N[i+1] <- RickerStBS_EhB(N[i],R,alpha,kE)
#  N[i+1] <- RickerStBS_DEhB(N[i],R,alpha,kD,kE) 
#  N[i+1] <- RickerStSexBS(N[i],R,alpha)
#  N[i+1] <- RickerStSexBS_DhB(N[i],R,alpha,kD)
#  N[i+1] <- RickerStSexBS_EhB(N[i],R,alpha,kE)
   N[i+1] <- RickerStSexBS_DEhB(N[i],R,alpha,kD,kE)
}
plot(t,N,type="l")
