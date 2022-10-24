# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Supplementary Fig. 4

# This program makes the plots from saved simulation output and calculates
# bias-corrected parameter values.

simdata <- read.csv("output/m_error_dat.csv")

# For "true" values, we use the fitted values from the data for the NBBg model
R <- 2.613
alpha <- 0.003731
kD <- 1.1475
kE <- 26.6221

# First check histograms for problems
sub <- subset(simdata,sigma==0.15)
par(mfrow=c(2,3))
attach(sub)
hist(NBBg_R)
abline(v=R,col="red")
hist(NBBg_a)
abline(v=alpha,col="red")
hist(NBBg_kD)
abline(v=kD,col="red")
hist(log(NBBg_kD))
abline(v=log(kD),col="red")
hist(NBBg_kE)
abline(v=kE,col="red")
hist(log(NBBg_kE))
abline(v=log(kE),col="red")
detach(sub)

# Make the graphs
windows(width=7,height=7*0.8)
par(mfrow=c(2,2),mar=c(1,4,0,0),oma=c(3,0.5,0.5,0.5))

R_hat <- expression( hat(italic(R)) )
a_hat <- expression( hat(italic(alpha)) )
kD_hat <- expression( hat(italic(k))[D] )
kE_hat <- expression( hat(italic(k))[E] )
m_err <- expression( Measurement*~~error*~~(italic(sigma)) )

mn <- NULL
for (s in unique(simdata$sigma)) {
  sub <- subset(simdata,sigma==s,select=NBBg_R)
  attach(sub)
  mn <- c(mn,mean(NBBg_R,trim=0))
  detach(sub)
}
plot(unique(simdata$sigma),mn,xlim=c(0,0.4),ylim=c(2.3,2.9),xlab="",ylab="",axes=FALSE)
axis(1, labels = FALSE )
axis(1, at=seq(0,0.4,by=0.05 ),labels=FALSE,tcl=-0.25)
axis(2)
abline(h=R,col="gray")
smooth <- smooth.spline(unique(simdata$sigma),mn)
xx  <- seq(0,0.4,length.out=201)
lines(predict(smooth,xx),col="red")
box()
mtext(R_hat,side=2,las=2,line=2.5)

mn <- NULL
for (s in unique(simdata$sigma)) {
  sub <- subset(simdata,sigma==s,select=NBBg_a)
  attach(sub)
  mn <- c(mn,mean(NBBg_a,trim=0))
  detach(sub)
}
plot(unique(simdata$sigma),mn,xlim=c(0,0.4),ylim=c(0.0031,0.0038),xlab="",ylab="",axes=FALSE)
axis(1, labels = FALSE )
axis(1, at=seq(0,0.4,by=0.05 ),labels=FALSE,tcl=-0.25)
axis(2)
abline(h=alpha,col="gray")
smooth <- smooth.spline(unique(simdata$sigma),mn,df=5)
xx  <- seq(0,0.4,length.out=201)
lines(predict(smooth,xx),col="red")
box()
mtext(a_hat,side=2,las=2,line=2.5)

mn <- NULL
for (s in unique(simdata$sigma)) {
  sub <- subset(simdata,sigma==s,select=NBBg_kD)
  attach(sub)
  mn <- c(mn,mean(NBBg_kD,trim=0.05))
  detach(sub)
}
plot(unique(simdata$sigma),mn,xlim=c(0,0.4),ylim=c(0,25),xlab="sigma",ylab="",axes=FALSE)
axis(1)
axis(1, at=seq(0,0.4,by=0.05 ),labels=FALSE,tcl=-0.25)
axis(2)
abline(h=kD,col="gray")
smooth <- smooth.spline(unique(simdata$sigma),mn,df=7)
xx  <- seq(0,0.4,length.out=201)
lines(predict(smooth,xx),col="red")
box()
mtext(kD_hat,side=2,las=2,line=2.5)

mn <- NULL
for (s in unique(simdata$sigma)) {
  sub <- subset(simdata,sigma==s,select=NBBg_kE)
  attach(sub)
  mn <- c(mn,mean(NBBg_kE,trim=0.05))
  detach(sub)
}
plot(unique(simdata$sigma),mn,xlim=c(0,0.4),ylim=c(0,45),xlab="sigma",ylab="",axes=FALSE)
axis(1)
axis(1, at=seq(0,0.4,by=0.05 ),labels=FALSE,tcl=-0.25)
axis(2)
abline(h=kE,col="gray")
smooth <- smooth.spline(unique(simdata$sigma),mn,df=7)
xx  <- seq(0,0.4,length.out=201)
lines(predict(smooth,xx),col="red")
box()
mtext(kE_hat,side=2,las=2,line=2.5)
mtext(m_err,side=1,line=1.5,outer=TRUE)

# Bias corrected estimates for kD and kE
sub <- subset(simdata,sigma==0,select=NBBg_kE)
mn <- mean(sub$NBBg_kE,trim=0.05)
rd <- (log(mn)-log(kE))/log(mn) #ln scale relative deviation
exp(log(kE) - rd*log(kE))
rd <- (mn-kE)/mn #natural scale relative deviation
kE - rd*kE

sub <- subset(simdata,sigma==0,select=NBBg_kD)
mn <- mean(sub$NBBg_kD,trim=0.05)
rd <- (log(mn)-log(kD))/log(mn) #ln scale relative deviation
exp(log(kD) - rd*log(kD))
rd <- (mn-kD)/mn #natural scale relative deviation
kD - rd*kD

