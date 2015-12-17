args <- commandArgs(trailingOnly = TRUE)

png(paste(args[2], "/relate.png", sep=""))

GEN=read.table(args[1], header=T, as.is=T)
hist(GEN$PI_HAT,50)

dev.off()