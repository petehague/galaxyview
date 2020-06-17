#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
	path = "."
} else {
	path = args[1]
}

peakfind <- function(dist) {
	x <- dist
	level <- 0
	for (y in seq(sort(x))) {
		z <- sum(x[which(x>y)])
		if (z/sum(x) > 0.68) { level <- y }
	}
	x[which(x<level)] <- 0
	return (x)	
}

flux <- as.matrix(read.table(paste(path, "/gasimage.txt", sep=""), head=FALSE))
vel <- as.matrix(read.table(paste(path, "/gasvel.txt", sep=""), head=FALSE))
params <- read.table(paste(path,"/imagedata.txt", sep=""), head=FALSE)[[1]]

vel[which(is.nan(vel))] <- 0

vlimits <- c(vel[1,1], vel[1,length(vel[1,])])

vbins <- vel[1,]

nbins <- 25

size <- params[1]
windowSize <- 1.0/params[2]
baser <- 0.02/params[2]

inclination <- params[3]
if (inclination>(pi/2)) { inclination <- pi-inclination }

vc <- function(x,y) {
	return(1+y+(x-1)*size)
}

vplot <- function() {
	dev.new(width=4, height=4)
	vimage <- matrix(0,size,size)
	for (y in seq(size)) {
		for (x in seq(size)) {
			n <- 1+x+(y-1)*size
			dist <- vel[n,]
			medv <- vbins[which.max(dist)]
			if (sum(dist)<100) { medv = 0 }
			vimage[x,y] = medv
		}
	}
	#par(mfrow=c(1,1))
	image(vimage)
}

vplot2 <- function() {
	vimage <- matrix(0,size,size)
	vlimit <- 0

	fluxmin <- 0.10*max(flux)
	print(fluxmin)

	for (y in seq(size)) {
		for (x in seq(size)) {
			n <- 1+x+(y-1)*size
			dist <- peakfind(vel[n,])
			red <- sum(dist[which(vbins>vlimit)])
			blue <- sum(dist[which(vbins<(-vlimit))])
			if (sum(dist)>fluxmin) {
				vimage[y,x]=blue-red
				#if (blue-red==0) {
				#	vimage[y,x]=10
				#}
			} else {
				vimage[y,x]=-999999
			}
		}
	}
	#vimage[which(abs(vimage)>1e10)] <- 0
	#par(mfrow=c(1,1))
	colseq <- rgb(seq(1,0, by=-0.05), 0, seq(0,1,by=0.05))
	colseq <- c(rgb(seq(0,0.5,by=0.05), 0, seq(1,0.5,by=-0.05)), rgb(0,0,0), rgb(seq(0.5,1, by=0.05), 0, seq(0.5,0,by=-0.05)))

	rangesize <- 40
	rangeclip <- 50

	colseq <- c(rgb(0,0,0), rgb(seq(0,0.5,length.out=rangesize), 0, seq(1,0.5,length.out=rangesize)), rgb(seq(0.5,1, length.out=rangesize), 0, seq(0.5,0,length.out=rangesize)))

	vsize <- max(abs(vimage[which(vimage>-999999)]))+1
	breakseq <- c(-1000000, seq(-vsize,vsize,length.out=1+rangesize*2))

	image(vimage, col=colseq, breaks=breakseq, xaxt= "n", yaxt= "n", xlab="x(kpc)", ylab="y(kpc)")
	w <- round(1e3*windowSize/2)
	axis(1, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
	axis(2, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
}

residplot <- function() {
	velmodel <- as.matrix(read.table(paste(path, "/model/gasvel.txt", sep=""), head=FALSE))
	vimage <- matrix(0,size,size)
	for (y in seq(size)) {
		for (x in seq(size)) {
			n <- 1+x+(y-1)*size
			dist1 <- cumsum(vel[n,])
			dist2 <- cumsum(velmodel[n,])
			med1 <- which(dist1>=0.5*max(dist1))[1]
			med2 <- which(dist2>=0.5*max(dist2))[1]
			if (sum(dist1)>0 && sum(dist2)>0) {
				vimage[y,x] <- med2-med1
			} else {
				vimage[y,x] <- -999999
			}
		}
	}
	#vimage[which(abs(vimage)>1e10)] <- 0
	#par(mfrow=c(1,1))
	rangesize <- 40
	rangeclip <- 50

	#colseq <- c(rgb(seq(0,0.5,length.out=rangesize), 0, seq(1,0.5,length.out=rangesize)), rgb(0,0,0), rgb(seq(0.5,1, length.out=rangesize), 0, seq(0.5,0,length.out=rangesize)))

	#breakseq <- seq(-size/rangeclip,size/rangeclip,length.out=2+rangesize*2)

	colseq <- c(rgb(0,0,0), rgb(seq(0,0.5,length.out=rangesize), 0, seq(1,0.5,length.out=rangesize)), rgb(seq(0.5,1, length.out=rangesize), 0, seq(0.5,0,length.out=rangesize)))

	breakseq <- c(-1000000, seq(-size/rangeclip,size/rangeclip,length.out=1+rangesize*2))

	image(vimage, col=colseq, breaks=breakseq, xaxt= "n", yaxt= "n", xlab="x(kpc)", ylab="y(kpc)")
	w <- round(1e3*windowSize/2)
	axis(1, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
	axis(2, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
}

centreplot <- function() {
	n <- 1+(size/2)+((size/2)-1)*size
	plot(vbins, vel[n,]/sum(vel[n,]), type="l", xlab=expression(v[los]~(km/s)), ylab=expression(n/n[tot]))
	segments(0,0,0,1e10, lty=2)
}

plainfluxplot <- function() {
	image(flux, xaxt= "n", yaxt= "n", xlab="x(kpc)", ylab="y(kpc)")
	w <- round(1e3*windowSize/2)
	axis(1, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
	axis(2, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
}

fluxplot <- function(x,y, angle) {
	image(flux, xaxt= "n", yaxt= "n", xlab="x(kpc)", ylab="y(kpc)")
	w <- round(1e3*windowSize/2)
	axis(1, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
	axis(2, at=seq(0,1, length.out=5), labels=c(-w, -w/2, 0, w/2, w))
	xend <- (0.05*cos(angle))/windowSize
	yend <- (0.05*sin(angle))/windowSize
	x <- x/windowSize + 0.5
	y <- y/windowSize + 0.5
	points(x,y)
	segments(x-1000*cos(angle), y-1000*sin(angle), x+1000*cos(angle), y+1000*sin(angle), lty=3)
	points(c(x+xend, x-xend), c(y+yend, y-yend), pch=4)
}

dplot <- function(x,y, angle, dptype) {
	par(mfrow=c(2,2))
	fluxplot(x,y,angle)
	rotgrid <- matrix(0, nbins, length(vel[1,]))
	rotcur <- rep(0,nbins)
	radius <- seq(baser, baser*nbins/2,by=baser)
	radius <- c(sort(-radius), 0, radius)

 	sink(paste(path,"/curve.txt", sep=""))
	cat("r v\n")
	for(i in seq(length(radius))) {
		px <- radius[i]*cos(angle)
		py <- radius[i]*sin(angle)
		gx <- px * size/windowSize + size*0.5
		gy <- py * size/windowSize + size*0.5
		n <- c(vc(floor(gx), floor(gy)), vc(floor(gx), ceiling(gy)), vc(ceiling(gx), floor(gy)), vc(ceiling(gx), ceiling(gy)))
		xdist <- gx-floor(gx)
		ydist <- gy-floor(gy)
		v1 <- vel[n[1],]*xdist + vel[n[3],]*(1-xdist)
		v2 <- vel[n[2],]*xdist + vel[n[4],]*(1-xdist)
		vhist <- v1*(1-ydist) + v2*(ydist)
		peak1 <- vel[1,which.max(vel[n[1],])]*xdist + vel[1,which.max(vel[n[3],])]*(1-xdist)
		peak2 <- vel[1,which.max(vel[n[2],])]*xdist + vel[1,which.max(vel[n[4],])]*(1-xdist)
		peak <- peak1*(1-ydist) + peak2*(ydist)
		peak <- seq(vlimits[1], vlimits[2], length.out=length(vhist))[which.max(vhist)]
		rotgrid[i,] = vhist
		rotcur[i] = peak
		cat(c(radius[i]," ",rotcur[i]/sin(inclination),"\n"))
	}
	sink()

	vplot2()
	image(rotgrid, xaxt= "n", yaxt= "n", xlab="dR(kpc)", ylab=expression(v[rot]~(km/s)))
	axis(1, at=seq(0,1, length.out=5), labels=1000*as.numeric(quantile(radius)))
	axis(2, at=seq(0,1, length.out=5), labels=seq(round(vlimits[1]), round(vlimits[2]),length.out=5))
	segments(-1.0, -1.0/(vlimits[2]/vlimits[1] - 1.0), 2.0, -1.0/(vlimits[2]/vlimits[1] - 1.0), lty=2)
	segments(0.5, -1.0, 0.5, 2.0, lty=2)
	#lines(seq(0,1,length.out=nbins), (rotcur-vlimits[1])/(vlimits[2]-vlimits[1]), lty=1)
	if (dptype==0) {
		centreplot()
	} else {
		residplot()
	}
}

rotation <- function() {
	data <- read.table(paste(path,"/curve.txt", sep=""), head=TRUE)
	n <- which(data$v>0)[1]
	r0 <- data$r[n-1]+(data$r[n]-data$r[n-1])*(-data$v[n-1]/(data$v[n-1]+data$v[n]))
	data$r <- abs(data$r-r0)
	data$v <- abs(data$v)
	r1 <- data$r[(n-1):1]*1e3
	v1 <- data$v[(n-1):1]
	r2 <- data$r[n:length(data$r)]*1e3
	v2 <- data$v[n:length(data$r)]
	smooth1 <- splinefun(r1,v1)
	smooth2 <- splinefun(r2,v2)
	par(mfrow=c(1,1))
	rs <- seq(1,12, by=0.5)
	centre <- (smooth1(rs) + smooth2(rs))/2
	err <- abs(smooth1(rs) - smooth2(rs))/2
	err[which(err<1)] = 1
	plot(rs, centre, ylim=c(0,max(data$v)*1.3))
	segments(rs, smooth1(rs), rs, smooth2(rs))
	sink(paste(path,"/obs.csv", sep=""))
	for(i in seq(length(rs))) {
		cat(c(rs[i],", ",centre[i],"\n"))
	}
	sink()
	sink(paste(path,"/err.csv", sep=""))
	for(i in seq(length(rs))) {
		cat(c(rs[i],", ",err[i],"\n"))
	}
	sink()
}

baryons <- function() {
	data <- read.table(paste(path, "/gasphi.txt", sep=""), head=TRUE)
	curve <- splinefun(data$r*1e3, data$phi)
	rs <- seq(1,12, by=0.5)
	vc <- sqrt(curve(rs, deriv=1)*rs)
	sink(paste(path, "/gas.csv", sep=""))
	for(i in seq(length(rs))) {
		cat(c(rs[i], ", ", vc[i], "\n"))
	}
	sink()

	data <- read.table(paste(path,"/starphi.txt", sep=""), head=TRUE)
	curve <- splinefun(data$r*1e3, data$phi)
	rs <- seq(1,12, by=0.5)
	vc <- sqrt(curve(rs, deriv=1)*rs)
	sink(paste(path, "/star.csv", sep=""))
	for(i in seq(length(rs))) {
		cat(c(rs[i], ", ", vc[i], "\n"))
	}
	sink()

	sink(paste(path,"/star2.csv", sep=""))
	for(i in seq(length(rs))) {
		cat(c(rs[i], ", 0\n"))
	}
	sink()
}

curveplot <- function() {
	data <- read.csv(paste(path, "/obs.csv", sep=""), head=FALSE)
	radius <- data[[1]]
	obs <- data[[2]]
	data <- read.csv(paste(path, "/err.csv", sep=""), head=FALSE)
	err <- data[[2]]
	data <- read.csv(paste(path, "/gas.csv", sep=""), head=FALSE)
	gas <- data[[2]]
	data <- read.csv(paste(path, "/star.csv", sep=""), head=FALSE)
	stars <- data[[2]]
	plot(radius, obs, ylim=c(0,max(obs)*1.2), ylab=expression(v[c]), xlab="R(kpc)")
	segments(radius, obs-err, radius, obs+err)
	lines(radius, stars, col="red")
	lines(radius, gas, col="blue")
}

pdf("output.pdf")
dplot(0,0,(pi/2.0)-params[5], 0)
dev.off() -> dumpvar

pdf("flux.pdf")
plainfluxplot()
dev.off() -> dumpvar

pdf("output2.pdf")
dplot(0,0,(pi/2.0)-params[5], 1)
dev.off() -> dumpvar
