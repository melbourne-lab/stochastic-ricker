# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Figure 3 (this script produces the simulation data used for the plot)

# Extinction times for stochastic Ricker models found by direct simulation
# Two standardizations are possible (plus no standardization):
# 1) Same equilibrium abundance (carrying capacity) by adjusting alpha
# 2) Equal total variance of kE and kD (kE adjusted each time to a/c for alpha)
# To do the standardizations, set their flags to TRUE.

source("source/Ricker.r")             #Deterministic Ricker
source("source/RickerStBS.r")         #Stochastic birth and survival
source("source/RickerStBS_DhB.r")     #Stochastic b, s, plus dem het in b.
source("source/RickerStBS_EhB.r")     #Stochastic b, s, plus env het in b.
source("source/RickerStBS_DEhB.r")    #Stochastic b, s, plus env & dem het in b.
source("source/RickerStSexBS.r")      #Stochastic sex ratio, birth, survival
source("source/RickerStSexBS_DhB.r")  #St sr, b, s, plus dem het in b.
source("source/RickerStSexBS_EhB.r")  #St sr, b, s, plus env het in b.
source("source/RickerStSexBS_DEhB.r") #St sr, b, s, plus env & dem het in b.

alpha <- 0.05
Rseq <- c(1.05,1.1,seq(1.2,20,by=0.2)) #R sequence for sims
kdem <- 0.5     #k for Negative binomial (demographic) model
kenv <- kdem/alpha
maxt <- 1000 #round(exp(11+0.05*(11-2))) #Maximum time for 1 simulation
reps <- 2000   #Replicate sims per R value (suggest 2000)
Neq <- 30 #Mean equilibrium population size
Neqadjust <- TRUE #Adjusts alpha to hold Neq constant
kEadjust <- TRUE #Adjusts kE to give same total variance as kD

# Time series of different models
R <- 5
# Adjust alpha to get desired Neq
if (Neqadjust == TRUE) alpha <- log(R)/Neq
# Adjust kE to same total var as kD
if (kEadjust == TRUE) kenv <- kdem/alpha

t <- 0:10000
N <- t*NA
N[1] <- 50
for (i in 1:max(t)) {
# N[i+1] <- Ricker(N[i],R,alpha)
# N[i+1] <- RickerStBS(N[i],R,alpha)
  N[i+1] <- RickerStBS_DhB(N[i],R,alpha,kdem)
#  N[i+1] <- RickerStSexBS(N[i],R,alpha)
#  N[i+1] <- RickerStSexBS_DhB(N[i],R,alpha,kdem)
  if (N[i+1]==0) {
    endi <- i + 1
    break
  }
  endi <- i + 1 #If we get to here on final iteration
}
plot(t[1:endi],N[1:endi],type="l")
abline(h=round(log(R) / alpha),col="red")



#--NB (dem)----------------------
print("NB-dem")
NBd_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStBS_DhB(N,R,alpha,kdem)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBd_Tm <- c(NBd_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--NB (env)----------------------
print("NB-env")
NBe_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    if (kEadjust == TRUE) kenv <- kdem/alpha   #Adjust kE to same total var as kD
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStBS_EhB(N,R,alpha,kenv)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBe_Tm <- c(NBe_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--NB (gamma)----------------------
print("NB-gamma")
NBg_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    if (kEadjust == TRUE) kenv <- kdem/alpha   #Adjust kE to same total var as kD
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStBS_DEhB(N,R,alpha,kdem,kenv)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBg_Tm <- c(NBg_Tm,Tm)
  print(cbind(R,Tm,c1))
}


#--Poisson-binomial----------------------
print("Pois-bin")
PB_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStSexBS(N,R,alpha)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  PB_Tm <- c(PB_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--NB-binomial-dem----------------------
print("NB-bin-dem")
NBBd_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStSexBS_DhB(N,R,alpha,kdem)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBBd_Tm <- c(NBBd_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--NB-binomial-env----------------------
print("NB-bin-env")
NBBe_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    if (kEadjust == TRUE) kenv <- kdem/alpha   #Adjust kE to same total var as kD
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStSexBS_EhB(N,R,alpha,kenv)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBBe_Tm <- c(NBBe_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--NB-binomial (gamma)----------------------
print("NBB-gamma")
NBBg_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    if (kEadjust == TRUE) kenv <- kdem/alpha   #Adjust kE to same total var as kD
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStSexBS_DEhB(N,R,alpha,kdem,kenv)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<maxt] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  NBBg_Tm <- c(NBBg_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#--Poisson----------------------
Rseq <- c(1.05,1.1,1.2,1.4,1.6,seq(8,20,by=0.2)) #Adjust R sequence for sims
print("Poisson")
P_Tm <- NULL #Vector to hold extinction times
for (R in Rseq) {
  et <- NA*1:reps #vector for extinction time
  for (i in 1:reps) {
    if (Neqadjust == TRUE) alpha <- log(R) / Neq   #Adjust alpha to maintain same equil N across R
    N <- ceiling(log(R) / alpha) #To start N at equilibrium
    for (t in 0:(maxt-1)) {
      N <- RickerStBS(N,R,alpha)
    # If extinct, record extinction time and stop here.
      if (N == 0){
        et[i] <- t+1
        break
      }
    if (t==(maxt-1)) et[i] <- maxt
    }
  }
#  et <- et[et<(maxt/2)] #Trim off the cases where we hit the simulation maximum time
  cdf <- ecdf(et)
  t <- unique(et)
  t <- t[t<max(t)] #Trim the highest point, since P0=1, FnP0=Inf
  P0 <- cdf(t)
  FnP0 <- -log(1-P0)
  plot(t,FnP0,main=paste("R=",R,sep=" "),ylab="-ln(1-Po)")
  fit <- lm(FnP0~t)
  abline(fit)
  c1 <- exp(-coef(fit)[1]) #c1 should be close to 1 if there is no effect of initial pop
  Tm <- 1/coef(fit)[2]
  P_Tm <- c(P_Tm,Tm)
  print(cbind(R,Tm,c1))
}

#Save simulation data to file
data <- cbind(Rseq,P_Tm,NBd_Tm,NBe_Tm,NBg_Tm,PB_Tm,NBBd_Tm,NBBe_Tm,NBBg_Tm)
rownames(data) <- NULL
data <- as.data.frame(data)
write.csv(data,file="output/extinctiondata_xxx.csv")
