# Melbourne BA, Hastings A (2008). Extinction risk depends strongly on factors
# contributing to stochasticity. Nature 454:100–103.
# https://doi.org/10.1038/nature06922
#
# Figure 3

# Graph extinction times for stochastic Ricker models from previously run and
# saved simulations. This is the graph for the standardised run (standardised to
# the same equilibrium abundance (K) and same total variance.

simdata <- read.csv("output/extinctiondata_stdKVsim.csv")
attach(simdata)

windows()

#----Graph-----------------------------------------------
makeminorticks <- function(majorticks) {
  tks <- majorticks[1]
  for ( i in 1:(length(majorticks)-1) ) {
    l <- majorticks[i]
    u <- majorticks[i+1]
    tks <- c(tks,seq(l,u,((u-l)/10))[-1])
  }
  return(tks)
}

# Ticks for natural numbers on ln scale
majticks <- c(2,7,20,55,150,400,1100,3000,8100,22000,60000)
lnmajticks <- log(majticks)
majtickslb <- c(2,7,20,55,150,400,"1,100","3,000","8,100","22,000","60,000")
minticks <- makeminorticks(majticks)
lnminticks <- log(minticks)

ylim <- c(1,11)
xlim <- c(1,20)
labsz <- 0.8 #Size for curve labels

# Mathematical expressions for axis labels
m_ylab <- expression( T[m]  )
m_ylab2 <- expression( log(T[m]) )
m_xlab <- expression( italic(R) )

plot(1,1,ylim=ylim,xlim=xlim,xlab="",ylab="",type="n",axes=FALSE)
axis(1)
axis(1,at=seq(min(xlim),max(xlim),by=1),tcl=-0.25,labels=FALSE)
axis(2,at=lnmajticks,labels=majtickslb,cex.axis=0.8,las=2)
axis(2,at=lnminticks,tcl=-0.25,labels=FALSE)
axis(2,las=2,tcl=0.5,labels=FALSE) #Interior ticks for ln scale
axis(2,at=ylim[1]:ylim[2],tcl=0.25,labels=FALSE)
axis(4,las=2,tcl=0.5,mgp=c(1, 0.3, 0))
axis(4,at=ylim[1]:ylim[2],tcl=0.25,labels=FALSE)
mtext(m_xlab,1,line=2.7)
mtext(m_ylab,2,line=3)
mtext(m_ylab2,4,line=1)
box()

#points(Rseq,log(P_Tm))
lines(Rseq[1:5],log(P_Tm[1:5]))
smooth <- smooth.spline(Rseq[37:97],log(P_Tm[37:97]),df=5)
lines(smooth)
text(9,11,"Pois",pos=4,offset=0.1,cex=labsz)

#points(Rseq,log(NBd_Tm))
smooth <- smooth.spline(Rseq,log(NBd_Tm),df=20)
bit1 <- smooth$y[1:46]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "NBd",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(NBd_Tm),df=8)
bit2 <- smooth$y[47:97]
lines(Rseq,c(bit1,bit2))

#points(Rseq,log(NBe_Tm))
lines(Rseq[1:7],log(NBe_Tm[1:7]))
smooth <- smooth.spline(Rseq[10:97],log(NBe_Tm[10:97]),df=8)
lines(smooth)
text(4.2,11,"NBe",pos=4,offset=0.1,cex=labsz)

#points(Rseq,log(NBg_Tm))
smooth <- smooth.spline(Rseq,log(NBg_Tm),df=20)
bit1 <- smooth$y[1:46]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "NBg",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(NBg_Tm),df=8)
bit2 <- smooth$y[47:97]
lines(Rseq,c(bit1,bit2))

#points(Rseq,log(PB_Tm),lty=2)
smooth <- smooth.spline(Rseq,log(PB_Tm),df=20)
bit1 <- smooth$y[1:46]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "PB",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(PB_Tm),df=8)
bit2 <- smooth$y[47:97]
lines(Rseq,c(bit1,bit2),lty=2)

#points(Rseq,log(NBBd_Tm),lty=2)
smooth <- smooth.spline(Rseq,log(NBBd_Tm),df=15)
bit1 <- smooth$y[1:46]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "NBBd",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(NBBd_Tm),df=8)
bit2 <- smooth$y[47:97]
lines(Rseq,c(bit1,bit2),lty=2)

#points(Rseq,log(NBBe_Tm))
smooth <- smooth.spline(Rseq,log(NBBe_Tm),df=20)
bit1 <- smooth$y[1:46]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "NBBe",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(NBBe_Tm),df=8)
bit2 <- smooth$y[47:97]
lines(Rseq,c(bit1,bit2),lty=2)

#points(Rseq,log(NBBg_Tm))
smooth <- smooth.spline(Rseq,log(NBBg_Tm),df=15)
bit1 <- smooth$y[1:61]
text(Rseq[which.max(smooth$y)],max(smooth$y)+0.15,
     "NBBg",pos=4,offset=0.1,cex=labsz)
smooth <- smooth.spline(Rseq,log(NBBg_Tm),df=5)
bit2 <- smooth$y[62:97]
lines(Rseq,c(bit1,bit2),lty=2)
