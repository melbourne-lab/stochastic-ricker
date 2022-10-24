# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Supplementary Fig. 4 

# This script runs the simulations for Fig. 4 (see separate script for
# plotting). See how the parameter estimates are biased by measurement error.
# This is hypothetical (since no error) but was suggested by reviewer 2. Take
# the true model and add error. Assume the error is lognormal but rounded to
# integer values.

source("source/RickerStSexBS_DEhB.r")
source("source/Ricker_nll.r")      #Likelihood functions
library(stats4) #For mle

mydata <- read.csv("data/ricker_data.csv")
Nt_tr <- round(mydata$At) #The "true" Nt

# For "true" values, we use the fitted values from the data for the NBBg model
R <- 2.613
alpha <- 0.003731
kD <- 1.1475
kE <- 26.6221

sigmas <- seq(0.05,0.4,by=0.05)
reps <- 100 #Number of simulated datasets for each value of sigma

NBBg_R <- NBBg_a <- NBBg_kD <- NBBg_kE <- NA

for (i in 1:reps) {

  for (s in sigmas) {
  # Calc true data and add error
    Ntp1 <- RickerStSexBS_DEhB(Nt_tr,R,alpha,kD,kE) #NBBg model (true data)
    Ntp1 <- round(exp( log(Ntp1)+rnorm(length(Ntp1),sd=s) ))
    Nt <- round(exp( log(Nt_tr)+rnorm(length(Ntp1),sd=s) ))

  # NBBg model
    llfit <-try(mle( Ricker_nbinombinomgamma.nll,
                 start=list(lnR=log(R),lnalpha=log(alpha),lnkD=log(kD),lnkE=log(kE)) ) )
    if(class(llfit)=="try-error") next
    NBBg_R <- exp(coef(llfit)[1])
    NBBg_a <- exp(coef(llfit)[2])
    NBBg_kD <- exp(coef(llfit)[3])
    NBBg_kE <- exp(coef(llfit)[4])

    datm <- cbind(s,NBBg_R,NBBg_a,NBBg_kD,NBBg_kE)
    write.table(datm,"output/m_error.csv",sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
    print(c(i,s))
  }
}
