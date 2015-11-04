args <- commandArgs(trailingOnly = TRUE)

GEN=read.table(args[1], header=T, as.is=T)
hist(GEN$PI_HAT,50)
