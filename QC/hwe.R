args <- commandArgs(trailingOnly = TRUE)

png(paste(args[2], "/hwe.png", sep=""))

HWE=read.table(args[1], header=T, as.is=T)
P = HWE$P
n = length(P)
plot( qchisq((1:n)/(n+1),2), sort(-2*log(P)), main="Q-Q plot of log(control HWE P-values)", xlab="Expected quantile", ylab="Observed quantile" )
lines( c(0,50), c(0,50) ); grid()

dev.off()