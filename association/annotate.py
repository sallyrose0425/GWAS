import sys

f=open("~/dave/glist-hg19_sorted", 'r')

genes={}

for line in f:
  chr, start, stop, gene = line.split()
  if chr not in genes:
    genes[chr] = []
  start = int(start)
  stop=int(stop)
  genes[chr]=append(start)
  genes[chr][start]=[stop, gene]
f.close()

f2=open(sys.argv[1], 'r')
f3=open(sys.argv[1]+'.new', 'w')

f2.close()
f3.close()
  
