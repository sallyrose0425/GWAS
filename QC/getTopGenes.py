#########################
#
#  Usage: getTopGenes.py PHENO_top P
#
#########################

import sys

topFile = open(sys.argv[1], 'r')

P=float(sys.argv[2])

uniqueGenes = []

first = 1
for line in topFile:
  if first == 1:
    first = 0
  else:
    linesplit=line.split()
    thisP = float(linesplit[8])
    genes = linesplit[9]
    if thisP > P:
      break
    else:
      if genes != '.':
        genelist = genes.split('|')
        for gene in genelist:
          if gene.split('(')[0] not in uniqueGenes:
            uniqueGenes.append(gene.split('(')[0])

topFile.close()
outfile=open(sys.argv[1].split('_')[0]+'.genelist', 'w')
for unique in uniqueGenes:
  if '=' not in unique:
    outfile.write(unique+'\n')
outfile.close()
