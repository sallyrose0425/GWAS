args <- commandArgs(trailingOnly = TRUE)

png(paste(args[3], "/missing.png", sep=""))


IMISS=read.table(args[1], header=T, as.is=T)

LMISS=read.table(args[2], header=T, as.is=T)

oldpar=par(mfrow=c(1,2))

plot( (1:dim(IMISS)[1])/(dim(IMISS)[1]-1), sort(1-IMISS$F_MISS), main="Individual 

call rate cumulative distribution", xlab="Quantile", ylab="Call Rate" ); 

grid()

plot( (1:dim(LMISS)[1])/(dim(LMISS)[1]-1), sort(1-LMISS$F_MISS), main="SNP coverage 

cumulative distribution", xlab="Quantile", ylab="Coverage" ); grid()

par(oldpar)
dev.off()