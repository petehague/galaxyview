hirespath <- "survey_Mar312358"
medrespath <- "surveymrhi"

#hires <- c("v4-7-0-0.3-0", "v4-7-0-0.2-0", "v4-7-0-0.4-0")
#hiresreal <- c("v4-7-0-0.3-0")

hires <- c("v4-7-0-0.3-0", "v4-7-0-0.3-0.5", "v4-8-0-0.3-0", "v4-8-0-0.3-0.5" )
hiresreal <- c("v4-7-0-0.3-0", "v4-8-0-0.3-0")

medres <- c("v1-4-0-0.3-0", "v1-7-0-0.3-0", "v1-4-0-0.3-0.5", "v1-7-0-0.3-0.5", "v1-8-0-0.3-0", "v1-8-0-0.3-0.5", "v1-9-0-0.3-0", "v1-9-0-0.3-0.5", "v1-11-0-0.3-0", "v1-11-0-0.3-0.5", "v1-12-0-0.3-0", "v1-12-0-0.3-0.5")
medresreal <- c("v1-4-0-0.3-0", "v1-7-0-0.3-0", "v1-8-0-0.3-0","v1-9-0-0.3-0","v1-11-0-0.3-0","v1-12-0-0.3-0")

gamma <- numeric(0)
gammaerr1 <- numeric(0)
gammaerr2 <- numeric(0)

vmax <- numeric(0)
vmaxerr1 <- numeric(0)
vmaxerr2 <- numeric(0)

r2 <- numeric(0)
r2err1 <- numeric(0)
r2err2 <- numeric(0)

for (file in hires) {
	path <- paste(hirespath, "/", file, "/curve/mcmc.txt", sep="")
	data <- read.table(path, head=FALSE)
	gamma <- c(gamma, median(abs(data[[4]])))
	vmax <- c(vmax, median(data[[6]]))
	r2 <- c(r2, median(data[[5]]))
	
	x <- quantile(abs(data[[4]]), c(0.16, 0.86))
	gammaerr1 <- c(gammaerr1, x[[1]])
	gammaerr2 <- c(gammaerr2, x[[2]])
	
	x <- quantile(data[[6]], c(0.16, 0.86))
	vmaxerr1 <- c(vmaxerr1, x[[1]])
	vmaxerr2 <- c(vmaxerr2, x[[2]])
	
	x <- quantile(data[[5]], c(0.16, 0.86))
	r2err1 <- c(r2err1, x[[1]])
	r2err2 <- c(r2err2, x[[2]])
}

print(r2)
print(gamma)

plot(gamma, r2, pch=5, col="red", xlim=c(0,2), ylim=c(0,3), xlab=expression(gamma), ylab=expression(r[2]))
segments(gammaerr1, r2, gammaerr2, r2, col="red")
segments(gamma, r2err1, gamma, r2err2, col="red")

gamma <- numeric(0)
gammaerr1 <- numeric(0)
gammaerr2 <- numeric(0)

vmax <- numeric(0)
vmaxerr1 <- numeric(0)
vmaxerr2 <- numeric(0)

r2 <- numeric(0)
r2err1 <- numeric(0)
r2err2 <- numeric(0)

for (file in hiresreal) {	
	path <- paste(hirespath, "/", file, "/realcurve/mcmc.txt", sep="")
	data <- read.table(path, head=FALSE)
	gamma <- c(gamma, median(abs(data[[4]])))
	vmax <- c(vmax, median(data[[6]]))
	r2 <- c(r2, median(data[[5]]))
	
	x <- quantile(abs(data[[4]]), c(0.16, 0.86))
	gammaerr1 <- c(gammaerr1, x[[1]])
	gammaerr2 <- c(gammaerr2, x[[2]])
	
	x <- quantile(data[[6]], c(0.16, 0.86))
	vmaxerr1 <- c(vmaxerr1, x[[1]])
	vmaxerr2 <- c(vmaxerr2, x[[2]])
	
	x <- quantile(data[[5]], c(0.16, 0.86))
	r2err1 <- c(r2err1, x[[1]])
	r2err2 <- c(r2err2, x[[2]])
	
}

points(gamma, r2, pch=5, col="black")
segments(gammaerr1, r2, gammaerr2, r2, col="black")
segments(gamma, r2err1, gamma, r2err2, col="black")

gamma <- numeric(0)
gammaerr1 <- numeric(0)
gammaerr2 <- numeric(0)

vmax <- numeric(0)
vmaxerr1 <- numeric(0)
vmaxerr2 <- numeric(0)

r2 <- numeric(0)
r2err1 <- numeric(0)
r2err2 <- numeric(0)

for (file in medres) {
	path <- paste(medrespath, "/", file, "/curve/mcmc.txt", sep="")
	data <- read.table(path, head=FALSE)
	gamma <- c(gamma, median(abs(data[[4]])))
	vmax <- c(vmax, median(data[[6]]))
	r2 <- c(r2, median(data[[5]]))
	
	x <- quantile(abs(data[[4]]), c(0.16, 0.86))
	gammaerr1 <- c(gammaerr1, x[[1]])
	gammaerr2 <- c(gammaerr2, x[[2]])
	
	x <- quantile(data[[6]], c(0.16, 0.86))
	vmaxerr1 <- c(vmaxerr1, x[[1]])
	vmaxerr2 <- c(vmaxerr2, x[[2]])
	
	x <- quantile(data[[5]], c(0.16, 0.86))
	r2err1 <- c(r2err1, x[[1]])
	r2err2 <- c(r2err2, x[[2]])
}

points(gamma, r2, pch=5, col="blue")
segments(gammaerr1, r2, gammaerr2, r2, col="blue")
segments(gamma, r2err1, gamma, r2err2, col="blue")

gamma <- numeric(0)
gammaerr1 <- numeric(0)
gammaerr2 <- numeric(0)

vmax <- numeric(0)
vmaxerr1 <- numeric(0)
vmaxerr2 <- numeric(0)

r2 <- numeric(0)
r2err1 <- numeric(0)
r2err2 <- numeric(0)

for (file in medresreal) {	
	path <- paste(medrespath, "/", file, "/realcurve/mcmc.txt", sep="")
	data <- read.table(path, head=FALSE)
	gamma <- c(gamma, median(abs(data[[4]])))
	vmax <- c(vmax, median(data[[6]]))
	r2 <- c(r2, median(data[[5]]))
	
	x <- quantile(abs(data[[4]]), c(0.16, 0.86))
	gammaerr1 <- c(gammaerr1, x[[1]])
	gammaerr2 <- c(gammaerr2, x[[2]])
	
	x <- quantile(data[[6]], c(0.16, 0.86))
	vmaxerr1 <- c(vmaxerr1, x[[1]])
	vmaxerr2 <- c(vmaxerr2, x[[2]])
	
	x <- quantile(data[[5]], c(0.16, 0.86))
	r2err1 <- c(r2err1, x[[1]])
	r2err2 <- c(r2err2, x[[2]])
	
}

points(gamma, r2, pch=5, col="grey")
segments(gammaerr1, r2, gammaerr2, r2, col="grey")
segments(gamma, r2err1, gamma, r2err2, col="grey")