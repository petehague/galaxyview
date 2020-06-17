#!/usr/bin/env Rscript

args <- commandArgs(trailingonly = TRUE)
if (length(args)==0) {
  stop ("Path not supplied\n")
}
datapath <- args[1]

halolist <- c(1,2,4,7,8,9,11,12)

getdata <- function(id) {
  fileroot0 <- paste(datapath,"/v1-",id,"-0-0.3-0/", sep="")
  fileroot1 <- paste(datapath,"/v1-",id,"-0-0.3-0.5/", sep="")

  data <- read.table(paste(fileroot0,"obscurve/r2-histogram.txt", sep=""), head=FALSE)
  r2obs1 <- data[[1]][which.max(data[[2]])]

  data <- read.table(paste(fileroot0,"obscurve/gamma-histogram.txt", sep=""), head=FALSE)
  gammaobs1 <- abs(data[[1]][which.max(data[[2]])])

  data <- read.table(paste(fileroot0,"obscurve/vmax-histogram.txt", sep=""), head=FALSE)
  vmax1 <- data[[1]][which.max(data[[2]])]

  data <- read.table(paste(fileroot1,"obscurve/r2-histogram.txt", sep=""), head=FALSE)
  r2obs2 <- data[[1]][which.max(data[[2]])]

  data <- read.table(paste(fileroot1,"obscurve/gamma-histogram.txt", sep=""), head=FALSE)
  gammaobs2 <- abs(data[[1]][which.max(data[[2]])])

  data <- read.table(paste(fileroot1,"obscurve/vmax-histogram.txt", sep=""), head=FALSE)
  vmax2 <- data[[1]][which.max(data[[2]])]

  data <- read.table(paste(fileroot0,"realcurve/r2-histogram.txt", sep=""), head=FALSE)
  r2real <- data[[1]][which.max(data[[2]])]

  data <- read.table(paste(fileroot0,"realcurve/gamma-histogram.txt", sep=""), head=FALSE)
  gammareal <- abs(data[[1]][which.max(data[[2]])])

  data <- read.table(paste(fileroot0,"trcurve.txt", sep=""), head=TRUE)
  Rd <- 1000*data$r[which.max(data$stars)]
  vmax <- max(data$v)

  return (c(Rd,vmax,r2real,gammareal,r2obs1,gammaobs1,vmax1,r2obs2,gammaobs2,vmax2))
}

sink("dwarfdata.txt")
cat("Rd vmax r2 gamma r2_a gamma_a vmax_a r2_b gamma_b vmax_b\n")
for (id in halolist) {
  cat(getdata(id))
  cat("\n")
}
sink()

data <- read.table("dwarfdata.txt", head=TRUE)

pdf("dwplot1.pdf")
plot(data$r2/data$Rd, data$gamma, pch=19, cex=0.5, ylim=c(0,2), xlab=expression(r[2]/R[D]), ylab=expression(gamma))
points(data$r2_a/data$Rd, data$gamma_a, pch=19, cex=0.5,col="red")
segments(data$r2/data$Rd,data$gamma, data$r2_a/data$Rd, data$gamma_a, col="red")
dev.off()

pdf("realgammahist.pdf")
hist(data$gamma, xlim=c(0,2), breaks=10, xlab=expression(gamma), ylab="n", main=NULL)
dev.off()

pdf("obsgammahist.pdf")
hist(c(data$gamma_a, data$gamma_b), xlim=c(0,2), breaks=10, xlab=expression(gamma), ylab="n", main=NULL)
dev.off()

pdf("result1.pdf")
n <- c(1,3,4,5,6,7,8)
plot(data$gamma[n], data$vmax[n], xlim=c(0,2), ylim=c(0,250), xlab=expression(gamma), ylab=expression(v[max]), pch=seq(length(n)))
points(data$gamma_a[n], data$vmax[n], pch=seq(length(n)), col="red")
points(data$gamma_b[n], data$vmax[n], pch=seq(length(n)), col="blue")
dev.off()
