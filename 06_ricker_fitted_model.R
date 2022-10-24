# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Supplementary Fig. 2

source("source/Ricker.r")
source("source/Ricker_nll.r")
source("source/Ricker_var.r") #Variances of stochastic Ricker models
library(stats4) #mle

tribdata <- read.csv("data/ricker_data.csv")
tribdata$Nt <- round(tribdata$At)
tribdata$Ntp1 <- round(tribdata$Atp1)
attach(tribdata)

# Fit Negative binomial-binomial-gamma

llfit <- mle( Ricker_nbinombinomgamma.nll,start=list(
                                             lnR=log(2.613),
                                             lnalpha=log(0.003731),
                                             lnkD=log(1.1475),
                                             lnkE=log(26.6221) ) )
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat

# Figure

m_Nt <- expression( italic(N) [ italic(t) ] )
m_Ntp1 <- expression( italic(N) [ italic(t)+1 ] )
plot(Nt,Ntp1,xlab=m_Nt,ylab=m_Ntp1,col="black")
lines(1:1100,Ricker(1:1100,phat["R"], phat["alpha"]),col="red") #col="grey75"
x <- c(2,6,10,18,28,46,78,130,215,360,600,1000)
y <- Ricker(x,phat["R"], phat["alpha"])
std <- sqrt(Ricker_nbinombinomgamma.var(x,phat["R"], phat["alpha"],
                                          phat["kD"],phat["kE"]))
segments(x,y-std,x,y+std,col="red") #col="grey75"
