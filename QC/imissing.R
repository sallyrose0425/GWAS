PMISS=read.table("your2_phenomiss.missing", header=T, as.is=T)
GMISS=read.table("your2_genomiss.missing.hap", header=T, as.is=T)
LMISS=read.table("your2a_miss.lmiss", header=T, as.is=T, row.names=2)
FREQ=read.table("your2a_freq.frq", header=T, as.is=T, row.names=2)
oldpar=par(mfrow=c(2,2))
OK = (LMISS$F_MISS>0)&(FREQ$MAF>0); OK[is.na(OK)]=FALSE
plot( LMISS$F_MISS[OK], PMISS$P[OK], main=NULL, xlab="Missingness", ylab="P (miss vs pheno)" ); abline(lm(PMISS$P[OK]~LMISS$F_MISS[OK])); grid()
plot( FREQ$MAF[OK], PMISS$P[OK], main=NULL, xlab="MAF", ylab="P (miss vs pheno)" ); 
abline(lm(PMISS$P[OK]~LMISS$F_MISS[OK])); grid()
minP = tapply( GMISS$P, as.factor(GMISS$SNP), min, na.rm=T )
minP=minP[is.finite(minP)]
OK=names(minP)
plot( LMISS[OK,"F_MISS"], minP, main=NULL, xlab="Missingness", ylab="P (miss vs geno)" ); abline(lm(minP~LMISS[OK,"F_MISS"])); grid()
plot( FREQ[OK,"MAF"], minP, main=NULL, xlab="MAF", ylab="P (miss vs geno)" ); 
abline(lm(minP~FREQ[OK,"MAF"])); grid()
par(oldpar)
#Alternatively, try displaying coplots against different MAF frequency intervals as follows
OK = (LMISS$F_MISS>0)&(FREQ$MAF>0); OK[is.na(OK)]=FALSE
coplot(PMISS$P[OK]~LMISS$F_MISS[OK] | FREQ$MAF[OK] )
OK=names(minP)
coplot(minP~LMISS[OK,"F_MISS"] | FREQ[OK,"MAF"] )
