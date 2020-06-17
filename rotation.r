#!/usr/bin/env Rscript

data <- read.table("gasslit.txt", head=TRUE)
pts = 0.5*(data$vr+(-data$vl))
if (sum(pts)<0) {
	pts <- -pts
	data$vr <- -data$vr
	data$vl <- -data$vl
}

phi <- read.table("halophi.txt", head=TRUE)
phifunc <- splinefun(phi$r, phi$pot)
realv <- sqrt(phifunc(phi$r, deriv=1)*phi$r)


pole1phi <- read.table("polar1phi.txt", head=TRUE)
phifunc <- splinefun(pole1phi$r, pole1phi$pot)
realv1 <- sqrt(phifunc(phi$r, deriv=1)*phi$r)


pole2phi <- read.table("polar2phi.txt", head=TRUE)
phifunc <- splinefun(pole2phi$r, pole2phi$pot)
realv2 <- sqrt(phifunc(phi$r, deriv=1)*phi$r)

gas <- read.table("gaspoints.txt", head=TRUE)
star <- read.table("starpoints.txt", head=TRUE)

#dv <- numeric(0)
#for(r in data$r) {
#	n <- which(abs(gas$r-r)<data$r[2]-data$r[1])
#	dv <- c(dv, median(sqrt(gas$vzsq[n])))
#}

spherical <- read.table("gmoverr.txt", head=TRUE)
tring <- read.table("trcurve.txt", head=TRUE)
galid <- read.table("galid.txt", head=FALSE)[1,1]

pdf("rc.pdf")
plot(data$r*1e3, pts, ylim=c(0,1.5*max(c(data$vr, -data$vl))), xlim=c(0,max(data$r*1e3)/2), xlab="R[kpc]", ylab=expression(v[c]~"[kpc]"))
segments(data$r*1e3, data$vr, data$r*1e3, -data$vl)
#plot(phi$r*1000, realv, lty=1, col="purple", lwd=2, xlim=c(0,10), type="l", xlab="R[kpc]", ylab=expression(v[c]~"[kpc]"))
lines(phi$r*1000, realv, lty=1, col="purple", lwd=2)
lines(phi$r*1000, realv1, lty=2, col="purple", lwd=2)
lines(phi$r*1000, realv2, lty=3, col="purple", lwd=2)
lines(tring$r*1e3, tring$v,lty=2, col="green", lwd=2)
#points(data$r*1e3, pts+dv, pch=4)
lines(spherical$r*1e3, spherical$v, lty=1, col="red", lwd=2)
#points(gas$r*1e3, sqrt(gas$vzsq), cex=0.1, xlim=c(0,10))
if (galid==4) {
	segments(0, 207.79, 15, 207.79, lty=2)
}
if (galid==1) {
	segments(0, 221.36, 15, 221.36, lty=2)
}

dev.off()
