# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Supplementary Fig. 3

# This is to plot the variance of some stochastic Ricker models that have
# stochasticity in survival instead of in births.

source("source/Ricker.R")           #Deterministic Ricker
source("source/Ricker_var.R")       #Variances of stochastic Ricker models
source("source/Ricker_var_extra.R") #Variances of extra stochastic Ricker models

#----Parameters-------------
R <- 5
m <- 0.9  #Probability of mortality
var_mE <- 0.02
var_mD <- 0.08 #Maximum variance possible in var_mD is < m(1-m)
alpha <- 0.05
sd_alpha <- 0.01
kalpha <- 25 #Shape parameter of gamma for cannibalism parameter
kdem <- 0.5 #Shape parameter of gamma for birth rate (between individuals)
kenv <- kdem/alpha #k for env variation: Total Var same as NBdem-Ricker
var_mE <- (1-m)^2/kenv

#----Calculate the variance------------

NtL <- 1:110 #long version of Nt
NtM <- c(1,seq(5,110,by=5)) #medium length version of Nt
NtS <- c(1,seq(10,110,by=10)) #short version of Nt

# Poisson model
varp <- Ricker_pois.var(NtL,R,alpha)

# Negative binomial - environmental
varnbe <- Ricker_nbinom_e.var(NtL,R,alpha,kenv)

# Poisson beta - environmental
VEpbeta_e <- Ricker_poisbeta_e.var(NtS,R,m,alpha,var_mE,reps=100000)
varpbetae <- VEpbeta_e$var

# Poisson beta - demographic
VEpbeta_d <- Ricker_poisbeta_d.var(NtS,R,m,alpha,var_mD,reps=10000)
varpbetad <- VEpbeta_d$var

# Gamma alpha - demographic
VEga_d <- Ricker_gammaalpha_d.var(NtM,R,alpha,kalpha,reps=10000)
vargad <- VEga_d$var

# Gamma alpha - environmental
VEga_e <- Ricker_gammaalpha_e.var(NtM,R,alpha,kalpha,reps=100000)
vargae <- VEga_e$var


#----Plot the variance----------------

xlim <- c(0,105) #x-axis limits
ylim <- c(0,200)

labsz <- 0.8 #cex for label size within the plot
vo <- 5 #vertical offset for curve labels in panel A

# Mathematical expressions for axis labels
m_Nt <- expression( italic(N) [ italic(t) ] )
sigsq <- expression( italic(sigma)[ italic(N) [ italic(t)+1 ] ]^2 )

# Plot axes and labels
plot(1,1,xlim=xlim,ylim=ylim,type="n",axes=FALSE,xlab=m_Nt,ylab="")
mtext(sigsq,2,cex=1.2,line=2.3)
axis(1, at = seq(from = 0 , to = 90 , by = 30 ))
axis(1, at = seq(from = 0 , to = 100 , by = 10 ),labels=FALSE,tcl=-0.25 )
axis(2, at = seq(from = 0 , to = 200 , by = 50 ))
axis(2, at = seq(from = 0 , to = 200 , by = 25 ),labels=FALSE,tcl=-0.25 )
abline(v=1/alpha,col="grey")
box()

# Poisson model
lines(NtL,varp)
text(1/alpha,max(varp)+vo,"P",pos=4,offset=0.1,cex=labsz)

# Negative binomial - environmental
lines(NtL,varnbe)
text(1/alpha,max(varnbe)+vo,"NBe",pos=4,offset=0.1,cex=labsz)

# Gamma alpha - environmental
#points(NtM,vargae,col="red")
smooth <- smooth.spline(NtM,vargae,df=12)
xx  <- seq(min(NtM),max(NtM),length.out=201)
lines(predict(smooth,xx),col="red")
#ae <- expression( italic(alpha)~(env) )
#text(40,max(vargae)+vo,ae,pos=4,offset=0.1,cex=labsz,col="pink")

# Gamma alpha - demographic
#points(NtM,vargad,col="green")
smooth <- smooth.spline(NtM,vargad,df=12)
xx  <- seq(min(NtM),max(NtM),length.out=201)
preds <- predict(smooth,xx)
lines(c(NtL[1:20],preds$x[36:201]),c(varp[1:20],preds$y[36:201]),col="blue")
#ad <- expression( italic(alpha)~(dem) )
#text(30,max(vargad)+vo,ad,pos=4,offset=0.1,cex=labsz,col="green")

# Poisson beta - environmental
points(NtS,varpbetae,col="red")

# Poisson beta - demographic
points(NtS,varpbetad,col="blue")

legend(80,200, legend=expression(italic(m)~(env),italic(m)~(dem),
                                 italic(alpha)~(env),italic(alpha)~(dem)),
       pch=c(1,1,NA,NA), lty=c(NA,NA,1,1),
       col=c("red","blue","red","blue"))
