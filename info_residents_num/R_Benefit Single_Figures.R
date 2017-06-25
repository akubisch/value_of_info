# Version for new Alex Data with modiefied mutation (mutating logit transform) and with explicit data on dispersal trait. 
# Estimating of emigration based on distribution statistics (mean, s.d. thus not necessary)

load("processed_data.RDat")
Cap=100
mu=0.1
Lambda=2
sigma=1

################ emergent emigration prob. for invasive TA - blue lines indicate strategy with max. fitness
par(mfrow=c(5,4), ps=18, mar=c(3,3,1,1), oma=c(3,3,3,1))
for (i in 1:length(TAeInfo))
{
  plot(vCT, M_vemip_DD[,i], type="l", bty="l", lwd=2, yaxs='i',
       xlab="", ylab="")
  #abline(h=M_geo0_DI[i], lty=3, col="grey")
  #abline(h=M_geo1_DI[i], lty=3, col="orange3")
  if(M_bestemi_DD[i]>0.01)
  {
    segments(0, M_bestemi_DD[i], M_bestCT[i], M_bestemi_DD[i], col="blue3")
    segments(M_bestCT[i],0,M_bestCT[i], M_bestemi_DD[i], col="blue3")
  }
  mtext(paste("lev. info.=", TAeInfo[i], sep=''), side=3, line=0, adj=1, cex=0.6)
}
mtext("emigration threshold", side=1, outer=T, line=1)
mtext("emigration prob. (TA)", side=2, outer=T, line=1)
M_bestemi_DD; M_bestCT

###### plot emergent fitness over threhold parameter
par(mfrow=c(5,4), ps=18, mar=c(3,3,1,1), oma=c(3,3,3,1))
for (i in 1:length(TAeInfo))
{
  plot(vCT, M_vmfit_DD[,i]-1, type="l", bty="l", lwd=2,
       xlab="", ylab="", ylim=c(-0.06, 0.06))
  #abline(h=M_geo0_DI[i]-1, lty=3, col="grey")
  #abline(h=M_geo1_DI[i]-1, lty=3, col="orange3")
  #abline(h=M_mfit_DI[i]-1, col="green3")
  abline(h=0, lty=3, col="grey70")
 # abline(h=M_fitMax_DD[i]-1, lty=3, col="salmon")
  if(M_bestemi_DD[i]>0.01)
  {
    # indicating optimal CT and maximum fitness 
    segments(M_bestCT[i],0,M_bestCT[i], M_fitMax_DD[i]-1, col="blue3", lwd=1)
  }
  mtext(paste("lev. info.=", TAeInfo[i], sep=''), side=3, line=0, adj=1, cex=0.6)
  mtext(paste("bestCT=", M_bestCT[i], sep=''), side=3, line=-1, adj=1, cex=0.6)
  
}
mtext("emigration threshold", side=1, outer=T, line=1)
mtext("fitness gain invasive TA", side=2, outer=T, line=1)

### plot maximum fitness gain (arithmetic) and opt emigr. prob over information accuracy in resident population
par(mfrow=c(2,1), mar=c(2,4,1,1), oma=c(3,2,1,1))
plot(TAeInfo, M_fitMax_DD-1, pch=16, bty="l", ylim=c(0, 0.06), cex=1.2, yaxs="i", type="b",
     xlab="", ylab="fitness gain")
#lines(smooth.spline(TAeInfo, M_fitMax_DD-1), col="red3", lwd=2)
#abline(h=0, lty=3, lwd=2, col="grey70")

plot(TAeInfo, M_bestemi_DD, pch=16, bty="l", cex=1.2,yaxs="i", ylim=c(0, 0.5), type="b",
     xlab="", ylab="emigration prob.")
#lines(smooth.spline(TAeInfo, M_bestemi_DD, df=10), col="red3", lwd=2)
points(TAeInfo, M_resDIE, pch=16, col="grey70", cex=1) # mean emigration of resident population
mtext("information accuracy", 1, outer=T, line=1)

