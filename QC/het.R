args <- commandArgs(trailingOnly = TRUE)

HET=read.table(args[1], header=T, as.is=T)
H = (HET$N.NM.-HET$O.HOM.)/HET$N.NM.
oldpar=par(mfrow=c(1,2))
hist(H,50)
hist(HET$F,50)
par(oldpar)
HET[order(HET$F),]
