# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Table 1

# This script fits the stochastic Ricker models to data from a density
# experiment with Tribolium castaneum.

# In addition to the models described in the paper, included here are "zero
# inflated" model fits, which account for the case when all individuals at the
# beginning of a generation are female and do not have access to males to mate
# with. This situation did not occur for the experimental design in the paper
# since females had access to males before establishing the initial generation.

source("source/Ricker.R")
source("source/Ricker_nll.R")
library(stats4) #mle

tribdata <- read.csv("data/ricker_data.csv")
# Round abundance to nearest integer. We'll use "N" as the symbol for adults.
tribdata$Nt <- round(tribdata$At)
tribdata$Ntp1 <- round(tribdata$Atp1)
head(tribdata)
attach(tribdata)

# Basic plot
plot(Nt,Ntp1)

# Plot by batch (nothing here indicates an effect of batch)
clr <- round((batch-3)*10)
plot(Nt,Ntp1,col=clr)
plot(jitter(Nt),Ntp1,col=clr,xlim=c(0,50),ylim=c(0,150))
plot(log(Nt),log(Ntp1),col=clr)

# Try leaving out influential/outliers (case 19, 24, 58)
# tribdata <- tribdata[-58,]
# tribdata <- tribdata[-c(19,24,58),]
# Or try leaving out all the points at the max of prod func (viz reviewer 1)
# tribdata <- tribdata[!(tribdata$Nt %in% c(358,363,360)),]


# Fitting the Ricker models by maximum likelihood

# Likelihood for Poisson
llfit <- mle( Ricker_pois.nll,start=list(lnR=log(3),lnalpha=log(0.003)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
Poiss_AIC <- AIC(llfit)
Poiss_AIC

# Likelihood for Negative binomial (demographic)
llfit <- mle( Ricker_nbinom_d.nll,start=list(lnR=log(3),lnalpha=log(0.003),
              lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NB_AIC <- AIC(llfit)
NB_AIC

# Likelihood for Negative binomial (environmental)
llfit <- mle( Ricker_nbinom_e.nll,start=list(lnR=log(3),lnalpha=log(0.003),
              lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NB_AIC <- AIC(llfit)
NB_AIC

# Likelihood for Negative binomial-gamma
llfit <- mle( Ricker_nbinomgamma.nll,start=list(lnR=log(3),lnalpha=log(0.003),
              lnkD=log(1),lnkE=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBg_AIC <- AIC(llfit)
NBg_AIC

# Likelihood for Poisson-binomial
llfit <- mle( Ricker_poisbinom.nll,start=list(lnR=log(3),lnalpha=log(0.003)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
PB_AIC <- AIC(llfit)
PB_AIC


# Likelihood for Negative binomial-binomial (demographic)
llfit <- mle( Ricker_nbinombinom_d.nll,start=list(lnR=log(2.64),
              lnalpha=log(0.0037),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBBd_AIC <- AIC(llfit)
NBBd_AIC

# Likelihood for Negative binomial-binomial (environmental)
llfit <- mle( Ricker_nbinombinom_e.nll,start=list(lnR=log(2.64),
              lnalpha=log(0.0037),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBBe_AIC <- AIC(llfit)
NBBe_AIC

# Likelihood for Negative binomial-binomial-gamma
llfit <- mle( Ricker_nbinombinomgamma.nll,start=list(lnR=log(2.6),
              lnalpha=log(0.0037),lnkD=log(1),lnkE=log(30)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE) #10-20 mins
confint(llfit) #10-20 mins
NBBg_AIC <- AIC(llfit)
NBBg_AIC

# Likelihood for ZIP-binomial
llfit <- mle( Ricker_zipbinom.nll,start=list(lnR=log(3),lnalpha=log(0.003)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
ZIPB_AIC <- AIC(logLik(llfit))
ZIPB_AIC

# Likelihood for ZINB-binomial (demographic)
llfit <-mle( Ricker_zinbbinom_d.nll,start=list(lnR=log(2.64),
             lnalpha=log(0.0037),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
ZINBd_AIC <- AIC(logLik(llfit))
ZINBd_AIC

# Likelihood for ZINB-binomial (environmental)
llfit <- mle( Ricker_zinbbinom_e.nll,start=list(lnR=log(2.64),
              lnalpha=log(0.0037),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
ZINBe_AIC <- AIC(logLik(llfit))
ZINBe_AIC

# Likelihood for ZINB-binomial-gamma
llfit <-mle( Ricker_zinbbinomgamma.nll,start=list(lnR=log(2.6),
             lnalpha=log(0.0037),lnkD=log(1),lnkE=log(30)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
ZINBBg_AIC <- AIC(logLik(llfit))
ZINBBg_AIC
