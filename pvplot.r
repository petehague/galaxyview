#!/usr/bin/env Rscript

cat("Making plots...\n")

tring <- read.table("trcurve.txt", head=TRUE)
gmoverr <- read.table("gmoverr.txt", head=TRUE)

data <- read.table("maxv.txt",head=FALSE)
vrange <- data[1,1]
dv <- 4
data <- read.table("rsize.txt", head=FALSE)
rsize <- data[1,1]

data <- as.matrix(read.table("observation.txt", head=FALSE))
pdf("pvplot.pdf")
image(log10(data), axes=0, xaxs="i", yaxs="i")
axis(1, at=seq(0,1,length.out=9),labels=1e3*c(seq(-rsize, 0, length.out=5),seq(rsize/4, rsize,length.out=4)))
axis(2, at=seq(0,1,length.out=9),labels=c(seq(-vrange,0,length.out=5),seq(vrange/4,vrange,length.out=4)))
segments(-1,0.5,2,0.5,lty=2)
segments(0.5,-1,0.5,2,lty=2)
dev.off() -> dumpvar

distmedian <- function(vels, distribution) {
	for (i in seq(2,length(distribution))) {
		distribution[i] <- distribution[i] + distribution[i-1]
	}
	midpoint <- distribution[length(distribution)]*0.5
	medpoint <- 0
	for (i in seq(1,length(distribution))) {
		if (distribution[i]<midpoint) {
			medpoint <- i
		}
	}
	return (vels[medpoint])
}

weightedmean <- function(vels, distribution) {
	v <-  sum(vels*distribution)/sum(distribution)
	return(v)
}


velseq <- seq((dv/2)-vrange,vrange-(dv/2),by=dv)
radius <- gmoverr$r[seq(2,1+(length(data[,1])-1)/2)]*1e3
leftcurve <- numeric(0)
rightcurve <- numeric(0)
for (i in seq(1,(length(data[,1])-1)/2 ) ) {
	leftdist <- as.numeric(data[(length(data[,1])+1)/2-i,])
	rightdist <- as.numeric(data[(length(data[,1]+1))/2+i,])
	#leftcurve <- c(leftcurve, distmedian(velseq,leftdist))
	#rightcurve <- c(rightcurve, distmedian(velseq,rightdist))
	leftcurve <- c(leftcurve, weightedmean(velseq,leftdist))
	rightcurve <- c(rightcurve, weightedmean(velseq,rightdist))
}

data <- read.table("finalcurve.txt", head=TRUE)
radius <- data$radius
rightcurve <- data$low
leftcurve <- -data$high

rc <- data.frame(r=radius, vr=rightcurve, vl=leftcurve)

pts <- (rightcurve-leftcurve)/2

pdf("tiltedringcurve.pdf")
plot(rc$r, pts, ylim=c(0,max(tring$v)*1.5), xlim=c(0,max(rc$r)), xlab="R[kpc]", ylab=expression(v[c]~"[kpc]"))
segments(rc$r, rc$vr, rc$r, -rc$vl)
lines(tring$r*1e3, tring$v,lty=2, col="green", lwd=2)
lines(gmoverr$r*1e3, gmoverr$v, lty=2, col="red", lwd=2)
#lines(tring$r*1e3, tring$gas,col="blue",lwd=2)
lines(tring$r*1e3, tring$stars,col="purple",lwd=2)
dev.off() -> dumpvar

errs <- numeric(0)
for (i in seq(1,length(rc$r))) {
	errs[i] = abs((rightcurve[i]+leftcurve[i])/2)
	if (is.nan(errs[i])) {
		errs[i] = 4
	} else {
		if (errs[i] == 0) { errs[i]=4 }
	}
}

sink("curvedata.txt")
cat("# Radius  v_obs  err  v_gas  v_disk  v_bulge\n")
cat("#  kpc    km/s   km/s km/s   km/s    km/s\n")
for (i in seq(1,length(rc$r))) {
	cat(rc$r[i])
	cat("  ")
	cat(pts[i])
	cat("  ")
	cat(errs[i])
	cat("  ")
  cat(tring$gas[i])
	cat("  ")
  cat(tring$stars[i])
	cat("  ")
	cat("0\n")
}
sink()

sink("realcurvedata.txt")
cat("# Radius  v_obs  err  v_gas  v_disk  v_bulge\n")
cat("#  kpc    km/s   km/s km/s   km/s    km/s\n")
for (i in seq(1,length(rc$r))) {
	cat(rc$r[i])
	cat("  ")
	cat(tring$v[i])
	cat("  ")
	cat(errs[i])
	cat("  ")
  cat(tring$gas[i])
	cat("  ")
  cat(tring$stars[i])
	cat("  ")
	cat("0\n")
}
sink()
