#!/usr/bin/env Rscript

tring <- read.table("../trcurve.txt", head=TRUE)
mcmc <-read.table("curvelog.txt", head=FALSE)

pdf("curves.pdf")
plot(tring$r,mcmc[1e5,],col="light grey", xlab="R [kpc]", ylab="Vc (km/s)")
for (n in seq(1e5,1.01e5)) {
  lines(tring$r,mcmc[n,],col="light grey")
}

lines(tring$r, tring$v, lty=2, lwd=2, col="green")
dev.off() -> dumpvar
