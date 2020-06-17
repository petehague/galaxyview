library("rgl")

crossdot <- function(a, b, c) {
	return (c[1]*(a[2]*b[3]-a[3]*b[2]) + c[2]*(a[3]*b[1]-a[1]*b[3]) + c[3]*(a[1]*b[2]-a[2]*b[1]) )
}

cross <- function(a, b) {
	return (c( a[2]*b[3]-a[3]*b[2] , a[3]*b[1]-a[1]*b[3] , a[1]*b[2]-a[2]*b[1] ))
}

arrow3d <- function(start, end, cl) {
	x <- (cross(end-start, c(0,0,1))/10)
	lines3d(c(start[1], end[1]), c(start[2], end[2]), c(start[3], end[3]), col=cl, lwd=3)
	head <- start+(end-start)*0.8
	lines3d(c(end[1], head[1]+x[1]), c(end[2], head[2]+x[2]), c(end[3], head[3]+x[3]), col=cl, lwd=3)
	lines3d(c(end[1], head[1]-x[1]), c(end[2], head[2]-x[2]), c(end[3], head[3]-x[3]), col=cl, lwd=3)
	lines3d(c(head[1]-x[1], head[1]+x[1]), c(head[2]-x[2], head[2]+x[2]), c(head[3]-x[3], head[3]+x[3]), col=cl, lwd=3)
}

if (rerun==0) {
	gas <- read.table("halo_gas.txt", head=TRUE)
	#n <- intersect(which(gas$m>2e5), which(sqrt(gas$x^2 + gas$y^2 + gas$z^2)<0.01))
	n <- which(sqrt(gas$x^2 + gas$y^2 + gas$z^2)<0.04)
	gas <- gas[n,]

	star <- read.table("halo_star.txt", head=TRUE)
	n <- which(sqrt(star$x^2 + star$y^2 + star$z^2)<0.04)
	star <- star[n,]
	n <- which(star$m > quantile(star$m, probs=0.8))
	#n <- sample(seq(length(star[,1])), 10000)
	star <- star[n,]

	phi <- read.table("halophi.txt", head=TRUE)
	phi$pot[which(is.nan(phi$pot))] <- 0
	phifunc <- splinefun(phi$r, phi$pot)
	realv <- sqrt(phifunc(phi$r, deriv=1)*phi$r)

	j0 <- c(0,0,0)
	for (i in seq(length(gas$x))) {
	  	j0 <- j0 + cross(c(gas$x[i], gas$y[i], gas$z[i]), c(gas$vx[i], gas$vy[i], gas$vz[i]))
	}
	j0 <- j0/sqrt(sum(j0^2))

	ratio <- numeric(length(gas$x))
	for (i in seq(length(gas$x))) {
		radius <- sqrt(gas$x[i]^2 + gas$y[i]^2 + gas$z[i]^2)
		ratio[i] <- crossdot(c(gas$x[i], gas$y[i], gas$z[i]), c(gas$vx[i], gas$vy[i], gas$vz[i]), j0)/radius
		radbin <- which.min(abs(phi$r - radius))
		ratio[i] <- ratio[i]/realv[radbin]
	}

	pcols <- ratio
	pcols[which(is.infinite(pcols))] = 0
	pcols <- pcols/max(pcols)
	pcols[which(pcols<0)] = 0

	sratio <- numeric(length(star$x))
	for (i in seq(length(star$x))) {
		radius <- sqrt(star$x[i]^2 + star$y[i]^2 + star$z[i]^2)
		sratio[i] <- crossdot(c(star$x[i], star$y[i], star$z[i]), c(star$vx[i], star$vy[i], star$vz[i]), j0)/radius
		radbin <- which.min(abs(phi$r - radius))
		sratio[i] <- sratio[i]/realv[radbin]
	}

	spcols <- sratio
	spcols[which(is.infinite(spcols))] = 0
	spcols <- spcols/max(spcols)
	spcols[which(spcols<0)] = 0

	force <- read.table("forcemap.txt", head=TRUE)
}

if (type==0) {
	n <- which(ratio>=1.0)
	plot3d(gas$x[n], gas$y[n], gas$z[n], col=rgb(pcols[n], 0, 0), xlim=c(-0.01, 0.01), ylim=c(-0.01, 0.01), zlim=c(-0.01, 0.01))
	n <- which(ratio<1.0)
	points3d(gas$x[n], gas$y[n], gas$z[n], col=rgb(0, 0, pcols[n]))
}

if (type==1) {
	n <- which(ratio>=1.0)
	plot3d(gas$x[n], gas$y[n], gas$z[n], col=rgb(pcols[n], 0, 0), xlim=c(-0.01, 0.01), ylim=c(-0.01, 0.01), zlim=c(-0.01, 0.01))
	n <- which(ratio<1.0)
	points3d(gas$x[n], gas$y[n], gas$z[n], col=rgb(0, 0, pcols[n]))
	arrow3d(-j0/500,j0/500, "green")
	for (i in seq(length(force$x))) {
		fscale <- 4e10
		arrow3d(c(force$x[i], force$y[i], force$z[i]), c(force$x[i], force$y[i], force$z[i])+c(force$Fx[i], force$Fy[i], force$Fz[i])/fscale, "black")
	}
}

if (type==2) {
	n <- which(sratio>=1.0)
	plot3d(star$x[n], star$y[n], star$z[n], col=rgb(spcols[n], 0, 0), xlim=c(-0.01, 0.01), ylim=c(-0.01, 0.01), zlim=c(-0.01, 0.01))
	n <- which(sratio<1.0)
	points3d(star$x[n], star$y[n], star$z[n], col=rgb(0, 0, spcols[n]))
	arrow3d(-j0/500,j0/500, "green")
}

if (type==3) {
	n <- which(sratio>=1.0)
	plot3d(star$x[n], star$y[n], star$z[n], col=rgb(spcols[n], 0, 0), xlim=c(-0.01, 0.01), ylim=c(-0.01, 0.01), zlim=c(-0.01, 0.01))
	n <- which(sratio<1.0)
	points3d(star$x[n], star$y[n], star$z[n], col=rgb(0, 0, spcols[n]))
	tr <- read.table("trcurve.txt", head=TRUE)
	theta <- seq(0,2*pi,by=pi/100)
	for (i in seq(length(tr$r)/2)) {
		xhat <- cross(c(tr$Lx[i], tr$Ly[i], tr$Lz[i]), c(1,0,0))
		yhat <- cross(c(tr$Lx[i], tr$Ly[i], tr$Lz[i]), xhat)
		xhat <- (xhat*tr$r[i])/sqrt(xhat[1]^2 + xhat[2]^2 + xhat[3]^2)
		yhat <- (yhat*tr$r[i])/sqrt(yhat[1]^2 + yhat[2]^2 + yhat[3]^2)
		lines3d(xhat[1]*cos(theta)+yhat[1]*sin(theta), xhat[2]*cos(theta)+yhat[2]*sin(theta), xhat[3]*cos(theta)+yhat[3]*sin(theta), col="green", lwd=2)
	}
}
