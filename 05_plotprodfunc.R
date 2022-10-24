# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Supplementary Fig 1.

# Simulate the and plot the production function for the family of stochastic
# Ricker models

source("source/Ricker.r")             #Deterministic Ricker
source("source/RickerStBS.r")         #Stochastic birth and survival
source("source/RickerStBS_DhB.r")     #Stochastic b, s, + dem het in b.
source("source/RickerStBS_EhB.r")     #Stochastic b, s, + env het in b.
source("source/RickerStBS_DEhB.r")    #Stochastic b, s, + dem & env het in b.
source("source/RickerStSexBS.r")      #Stochastic sex ratio, birth, survival
source("source/RickerStSexBS_DhB.r")  #St sr, b, s, + dem het in b.
source("source/RickerStSexBS_EhB.r")  #St sr, b, s, + env het in b.
source("source/RickerStSexBS_DEhB.r") #St sr, b, s, + dem & env het in b.
source("source/Ricker_var.r")         #Variances of stochastic Ricker models

# Parameters for simulating data
R <- 5
alpha <- 0.05
kD <- 0.5 #Shape parameter of gamma for birth rate
kE <- kD/alpha #k for env variation: Var same as NBdem-Ricker @ stat pt
Nt <- c(2,6,12,seq(20,160,by=10)) 
reps <- 100

# Plot setup
windows(width=7,height=7*2/4)
par(mfrow=c(2,4),mar=c(0,0,0,0),oma=c(4,5,1,1))
xlim <- c(0,160) #x-axis limits
ylim <- c(0,100)
xtks <- seq(from = 0 , to = 160 , by = 50 ) #x ticks
ytks <- seq(from = 0 , to = 100 , by = 20 )
xpl <- 150 #x position of panel label
ypl <- 90

# Poisson model
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks, labels = FALSE )
axis(2, at = ytks, las=1 )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStBS(Nt,R,alpha)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(y)
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"Poisson",pos=2,offset=0)

# Negative binomial - environmental
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks, labels = FALSE )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStBS_EhB(Nt,R,alpha,kE)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(y+(y*y)/kE)
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-environmental",pos=2,offset=0)

# Negative binomial - demographic
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks, labels = FALSE )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStBS_DhB(Nt,R,alpha,kD)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(y+(y*y)/(x*kD))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-demographic",pos=2,offset=0)

# Negative binomial - gamma
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks, labels = FALSE )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStBS_DEhB(Nt,R,alpha,kD,kE)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_nbinomgamma.var(x,R,alpha,kD,kE))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-gamma",pos=2,offset=0)

# Poisson-binomial
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks )
axis(2, at = ytks, las=1 )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS(Nt,R,alpha)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_poisbinom.var(x,R,alpha))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"Poisson-binomial",pos=2,offset=0)

# Negative binomial-binomial environmental
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS_EhB(Nt, R, alpha, kE)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_nbinombinom_e.var(x,R,alpha,kE))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-binomial-env",pos=2,offset=0)

# Negative binomial-binomial demographic
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS_DhB(Nt, R, alpha, kD)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_nbinombinom_d.var(x,R,alpha,kD))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-binomial-dem",pos=2,offset=0)

# Negative binomial-binomial gamma
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE)
axis(1, at = xtks )
axis(2, at = ytks, labels = FALSE )
box()
for (i in 1:reps) {
  Ntp1 <- RickerStSexBS_DEhB(Nt, R, alpha, kD, kE)
  points(jitter(Nt),Ntp1,col="grey75")
}
lines(0:1100,Ricker(0:1100,R,alpha),col="grey50")
x <- Nt
y <- Ricker(x,R,alpha)
std <- sqrt(Ricker_nbinombinomgamma.var(x,R,alpha,kD,kE))
segments(x,y-std,x,y+std,col="black")
text(xpl,ypl,"NB-binomial-gamma",pos=2,offset=0)

# Axis labels
m_Nt <- expression( italic(N) [ italic(t) ] )
m_Ntp1 <- expression( italic(N) [ italic(t)+1 ] )
mtext(m_Nt,side=1,line=2.8,outer=TRUE)
mtext(m_Ntp1,side=2,line=3,outer=TRUE)
