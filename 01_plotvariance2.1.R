# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Figure 2

#This is to plot the variance of the different Ricker models

source("source/Ricker.R")       #Deterministic Ricker
source("source/Ricker_var.R")   #Variances of stochastic Ricker models

#Parameters
R <- 5
alpha <- 0.05
kdem <- 0.5 #Shape parameter of gamma for birth rate (between individuals)
kenv <- kdem/alpha #k for env variation: Var same as NBdem-Ricker @ stat pt
Nt <- c(seq(1,50,by=1),seq(55,105,by=5)) #For final figure

#----Calculate the variance------------

#Poisson model
varp <- Ricker_pois.var(Nt,R,alpha)

#Negative binomial - demographic
varnbd <- Ricker_nbinom_d.var(Nt,R,alpha,kdem)

#Negative binomial - environmental
varnbe <- Ricker_nbinom_e.var(Nt,R,alpha,kenv)

#Negative binomial-gamma
varnbg <- Ricker_nbinomgamma.var(Nt,R,alpha,kdem,kenv)

#Poisson-binomial
varpb <- Ricker_poisbinom.var(Nt,R,alpha)

#Negative binomial-binomial demographic
varnbbd <- Ricker_nbinombinom_d.var(Nt,R,alpha,kdem)

#Negative binomial-binomial environmental
varnbbe <- Ricker_nbinombinom_e.var(Nt,R,alpha,kenv)

#Negative binomial-binomial-gamma
varnbbg <- Ricker_nbinombinomgamma.var(Nt,R,alpha,kdem,kenv)


#----Plot the variance----------------

windows()

xlim <- c(0,105) #x-axis limits
ylim <- c(0,650)

labsz <- 0.8 #cex for label size within the plot
vo <- 15 #vertical offset for curve labels in panel A

#Mathematical expressions for axis labels
m_Nt <- expression( italic(N) [ italic(t) ] )
sigsq <- expression( italic(sigma)[ italic(N) [ italic(t)+1 ] ]^2 )

#Plot axes and labels
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab=m_Nt,ylab="")
mtext(sigsq,2,cex=1.2,line=2)
axis(1, at = seq(from = 0 , to = 90 , by = 30 ))
axis(1, at = seq(from = 0 , to = 100 , by = 10 ),labels=FALSE,tcl=-0.25 )
axis(2, at = seq(from = 0 , to = 650 , by = 200 ))
axis(2, at = seq(from = 0 , to = 650 , by = 50 ),labels=FALSE)
axis(2, at = seq(from = 0 , to = 650 , by = 250 ),labels=FALSE,tcl=-0.25 )
abline(v=1/alpha,col="grey")
box()

#Poisson model
lines(Nt,varp)
text(1/alpha,max(varp)+vo,"P",pos=4,offset=0.1,cex=labsz)

#Negative binomial - demographic
lines(Nt,varnbd)
text(12,max(varnbd)+vo/2,"NBd",pos=4,offset=0,cex=labsz)

#Negative binomial - environmental
lines(Nt,varnbe)
text(1/alpha,max(varnbe)+vo,"NBe",pos=4,offset=0.1,cex=labsz)

#Poisson-binomial
lines(Nt,varpb,lty=2)
text(12,max(varpb)+vo/2,"PB",pos=4,offset=0,cex=labsz)

#Negative binomial-binomial demographic
lines(Nt,varnbbd,lty=2)
text(11,max(varnbbd)+vo,"NBBd",pos=4,offset=0.1,cex=labsz)

#Negative binomial-binomial environmental
lines(Nt,varnbbe,lty=2)
text(15,max(varnbbe)+vo,"NBBe",pos=4,offset=0.1,cex=labsz)

#Negative binomial-gamma
lines(Nt,varnbg)
text(13,max(varnbg)+vo,"NBg",pos=4,offset=0.1,cex=labsz)

#Negative binomial-binomial-gamma
lines(Nt,varnbbg,lty=2)
text(11,max(varnbbg)+vo,"NBBg",pos=4,offset=0.1,cex=labsz)
