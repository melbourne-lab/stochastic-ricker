# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922

# Bifurcation diagrams for stochastic Ricker models

source("source/Ricker.r")             #Deterministic Ricker
source("source/RickerStBS.r")         #Stochastic birth and survival
source("source/RickerStBS_DhB.r")     #Stochastic b, s, + dem het in b.
source("source/RickerStBS_EhB.r")     #Stochastic b, s, + env het in b.
source("source/RickerStBS_DEhB.r")    #Stochastic b, s, + dem & env het in b.
source("source/RickerStSexBS.r")      #Stochastic sex ratio, birth, survival
source("source/RickerStSexBS_DhB.r")  #St sr, b, s, + dem het in b.
source("source/RickerStSexBS_EhB.r")  #St sr, b, s, + env het in b.
source("source/RickerStSexBS_DEhB.r") #St sr, b, s, + dem & env het in b.


windows(width=9.5,height=7*2/3)
par(mfrow=c(2,4),mar=c(0,0,0,0),oma=c(4,5,1,1))

alpha <- 0.005
kdem <- 0.5     #k for Negative binomial (demographic) model
ylim <- c(0,2000)
xpl <- 2        #Panel label position
ypl <- 0.95*max(ylim)
dres <- 0.02    #Resolution in R for deterministic Ricker
sres <- 0.2     #Resolution in R for stochastic Ricker
dcol <- "grey"  #Color for deterministic Ricker
scol <- "red"   #Color for stochastic Ricker



#--Poisson----------------------
k <- kdem
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5),labels=FALSE)
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,las=2)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStBS(N[i],R,alpha)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"Poisson",pos=4,offset=0)


#--Negative binomial (dem)----------------------
k <- kdem
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5),labels=FALSE)
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStBS_DhB(N[i],R,alpha,k)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-demographic",pos=4,offset=0)


#--Negative binomial (env)----------------------
k <- kdem / alpha
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5),labels=FALSE)
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStBS_EhB(N[i],R,alpha,k)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-environmental",pos=4,offset=0)


#--Negative binomial gamma----------------------
kE <- kdem / alpha
kD <- kdem
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5),labels=FALSE)
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStBS_DEhB(N[i],R,alpha,kD,kE)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-gamma",pos=4,offset=0)


#--Poisson-binomial----------------------
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5))
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,las=2)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStSexBS(N[i],R,alpha)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"Poisson-binomial",pos=4,offset=0)


#--NB-binomial (dem)----------------------
k <- kdem
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5))
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStSexBS_DhB(N[i],R,alpha,k)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-binomial-dem",pos=4,offset=0)


#--NB-binomial (env)----------------------
k <- kdem / alpha
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5))
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStSexBS_EhB(N[i],R,alpha,k)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-binomial-env",pos=4,offset=0)


#--Negative binomial binomial gamma----------------------
kE <- kdem / alpha
kD <- kdem
plot(1,1,xlim=c(2,20),ylim=ylim,xlab="R",ylab="N",type="n",axes=FALSE)
axis(1,at=seq(5,20,by=5),labels=FALSE)
axis(1,at=seq(2,20,by=1),tcl=-0.25,labels=FALSE)
axis(2,labels=FALSE)
box()
for (R in seq(2,20,by=dres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  for (i in 1:max(t)) {
    N[i+1] <- Ricker(N[i],R,alpha)
  }
  points(rep(R,201),N[100:300],pch=".",col=dcol)
}
for (R in seq(2,20,by=sres)) {
  t <- 0:300
  N <- t*NA
  N[1] <- 20
  while ( N[300]==0|is.na(N[300]) ) { #Keep trying until we get no extinction
    for (i in 1:max(t)) {
      N[i+1] <- RickerStSexBS_DEhB(N[i],R,alpha,kD,kE)
    }
    points(rep(R,201),N[100:300],pch=".",col=scol)
  }
}
text(xpl,ypl,"NB-binomial-gamma",pos=4,offset=0)

m_R <- expression(italic(R))
m_N <- expression(italic(N))
mtext(m_R,side=1,line=2.5,outer=TRUE)
mtext(m_N,side=2,line=3.5,outer=TRUE,las=2)
