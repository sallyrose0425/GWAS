library('GWASTools')

args <- commandArgs(trailingOnly = TRUE)

data=read.table(args[1], header = TRUE)

chr=data[,1]

p=data[,9]

manhattanPlot(p, chr, thinThreshold=4, pointsPerBin=1000)

