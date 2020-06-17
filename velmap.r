#!/usr/bin/env Rscript

data <- read.table("gasdisk.txt", head=TRUE)
trmodel <- read.table("trcurve.txt", head=TRUE)

rc <- splinefun(trmodel$r, trmodel$v)
theta <- seq(0,2*pi, length.out=10000)

cross <- function(a, b) {
	return (c( a[2]*b[3]-a[3]*b[2] , a[3]*b[1]-a[1]*b[3] , a[1]*b[2]-a[2]*b[1] ))
}

dotprod <- function(a,b) {
  return ( a[1]*b[1] + a[2]*b[2] + a[3]*b[3] )
}

lhat <- c(-0.137364,0.332812,0.932935)
slit0 <- c(-0.000159857,0.00147301,0.000698394)
slit90 <- c(0.00197305,-0.000506444,-0.00141966)
i <- cross(lhat, c(1,0,0))
j <- cross(lhat, i)
slitpos0 <- c(dotprod(slit0, i), dotprod(slit0, j))
slitpos90 <- c(dotprod(slit90, i), dotprod(slit90, j))

slitpos0 <- 10 * slitpos0/sqrt(sum(slitpos0^2))
slitpos90 <- 10 * slitpos90/sqrt(sum(slitpos90^2))

pdf("velmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vtot <- sqrt(data$vx^2 + data$vy^2)
vscale <- 1e-3 * max(vtot)/50
radius <- sqrt(data$x^2 + data$y^2)*1e-3
n <- which(vtot>rc(radius))
m <- which(vtot<rc(radius))
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
segments(data$x[m], data$y[m], data$x[m]+data$vx[m]*vscale, data$y[m]+data$vy[m]*vscale, col="blue")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
segments(-slitpos0[1], -slitpos0[2], slitpos0[1], slitpos0[2], lty=2)
segments(-slitpos90[1], -slitpos90[2], slitpos90[1], slitpos90[2], lty=2)
dev.off()

pdf("radmotion.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vr <- abs(data$x*data$vx + data$y*data$vy)/sqrt(data$x^2 + data$y^2)
for (ratio in seq(0,0.9,by=0.1)) {
	n <- intersect(which(vr/vtot > ratio), which(vr/vtot < ratio+0.1))
	segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(ratio,0,1-ratio^0.5))
}
n <- which(vr/vtot > 1)
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("lzmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
Lz <- data$x*data$vy - data$y*data$vx
for (i in seq(10,1000, by=10)) {
	n <- intersect(which(Lz<i), which(Lz>i-100))
	segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(i/1e3,0,1-i/1e3))
}
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("philz.pdf")
phi <- atan2(data$x, data$y)
m <- intersect(which(radius<2.6e-3),which(radius>2.5e-3))
m <- m[order(phi[m])]
plot(phi[m], Lz[m], cex=0.1, ylim=c(0,4e-3), type="l")
dL <- 500
for (L in seq(dL,2000,by=dL)) {
	m <- intersect(which(Lz<L+dL),which(Lz>L))
	m <- m[order(phi[m])]
	if (length(m)>0) { lines(phi[m], radius[m]) }
}

n <- which(radius<4e-3)
coldmass <- sum(data$m[n])
data <- read.table("hotdisk.txt", head=TRUE)
radius <- sqrt(data$x^2 + data$y^2)*1e-3
n <- which(radius<4e-3)
hotmass <- sum(data$m[n])
dev.off()

pdf("fullvelmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vtot <- sqrt(data$vx^2 + data$vy^2)
n <- which(vtot>rc(radius))
m <- which(vtot<rc(radius))
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
segments(data$x[m], data$y[m], data$x[m]+data$vx[m]*vscale, data$y[m]+data$vy[m]*vscale, col="blue")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("fullradmotion.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vr <- abs(data$x*data$vx + data$y*data$vy)/sqrt(data$x^2 + data$y^2)
for (ratio in seq(0,0.9,by=0.1)) {
	n <- intersect(which(vr/vtot > ratio), which(vr/vtot < ratio+0.1))
	segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(ratio,0,1-ratio^0.5))
}
n <- which(vr/vtot > 1)
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("fulllzmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
Lz <- data$x*data$vy - data$y*data$vx
for (i in seq(10,1000, by=10)) {
	n <- intersect(which(Lz<i), which(Lz>i-100))
	segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(i/1e3,0,1-i/1e3))
}
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

data <- read.table("stellardisk.txt", head=TRUE)
radius <- sqrt(data$x^2 + data$y^2)*1e-3
n <- which(radius<4e-3)
starmass <- sum(data$m[n])
n <- sample(seq(length(data$x)), min(c(50000,length(data$m))))
data <- data[n,]

pdf("starvelmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vtot <- sqrt(data$vx^2 + data$vy^2)
radius <- sqrt(data$x^2 + data$y^2)*1e-3
n <- which(vtot>rc(radius))
m <- which(vtot<rc(radius))
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
segments(data$x[m], data$y[m], data$x[m]+data$vx[m]*vscale, data$y[m]+data$vy[m]*vscale, col="blue")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("starradmotion.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
vr <- abs(data$x*data$vx + data$y*data$vy)/sqrt(data$x^2 + data$y^2)
for (ratio in seq(0,0.9,by=0.1)) {
	n <- intersect(which(vr/vtot > ratio), which(vr/vtot < ratio+0.1))
	segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(ratio,0,1-ratio^0.5))
}
n <- which(vr/vtot > 1)
segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col="red")
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()

pdf("starlzmap.pdf")
plot(data$x, data$y, cex=0.1, xlim=c(-5,5), ylim=c(-5,5), xlab="x[kpc]", ylab="y[kpc]")
Lz <- data$x*data$vy - data$y*data$vx
for (i in seq(10,1000, by=10)) {
	n <- intersect(which(Lz<i), which(Lz>i-100))
	#segments(data$x[n], data$y[n], data$x[n]+data$vx[n]*vscale, data$y[n]+data$vy[n]*vscale, col=rgb(i/1e3,0,1-i/1e3))
	points(data$x[n], data$y[n], col=rgb(i/1e3,0,1-i/1e3), cex=0.2, pch=19)
}
lines(3*cos(theta), 3*sin(theta), lty=3, col="green", lwd=3)
points(0,0, pch=16)
dev.off()


hotmass <- hotmass-coldmass
totmass <- coldmass+hotmass+starmass
print(coldmass/totmass)
print(hotmass/totmass)
print(starmass/totmass)
