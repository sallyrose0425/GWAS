### USAGE
# python Add_PCs.py .evec pheno.txt new_pheno.txt
####

import sys

#get PC file
pc_file=open(sys.argv[1], 'r')

pc={}
first = 1
for line in pc_file:
  if first == 1:
    first = 0
  else:
    id, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10, status = line.split()
    if id not in pc:
      pc[id]=''
    pc[id]=pc1+' '+pc2+' '+pc3+' '+pc4+' '+pc5+' '+pc6+' '+pc7+' '+pc8+' '+pc9+' '+pc10
pc_file.close

#get pheno file
pheno_file=open(sys.argv[2], 'r')

#open new pheno file
new_pheno_file=open(sys.argv[3], 'w')

overlap=0
first = 1
for line in pheno_file:
  if first == 1:
    new_pheno_file.write(line.strip()+' pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10\n')
    first = 0
  else:
    linesplit=line.split()
    if linesplit[1] in pc:
      new_pheno_file.write(line.strip()+' '+pc[linesplit[1]]+'\n')
      overlap=overlap+1
    else:
      new_pheno_file.write(line.strip()+' -9 -9 -9 -9 -9 -9 -9 -9 -9 -9\n')
pheno_file.close()
new_pheno_file.close()
print('\nOverlap of IDs: '+str(overlap)+'\n')
