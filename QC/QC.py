import sys, os, subprocess, statistics

pfile=open(sys.argv[1], 'r')

params = {}
for line in pfile:
  paramname, paramvalue = line.split()
  params[paramname[1:]]=paramvalue

bfile=params['WORK']+'/'+params['INPUT'].split('/')[-1]
logfile=open(params['WORK']+'/final_QC.log', 'w')
logfile.write('Input parameters: '+sys.argv[1]+'\n')
for p in params:
  logfile.write(p+'\t'+params[p]+'\n')
logfile.write('************************************\n\n')

#CHECK SEX
if params['SEX'] == 'yes':
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --check-sex --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+' --out '+bfile+'_sexcheck')
  #number of samples
  proc = subprocess.Popen(['grep', 'people', bfile+'_sexcheck.log'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()
  n = int(tmp.split()[0])
  logfile.write('n = '+str(n)+'\n')

  #number of variants
  proc = subprocess.Popen(['grep', 'variants', bfile+'_sexcheck.log'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()
  v = int(tmp.split()[0])
  logfile.write('variants = '+str(v)+'\n\n')

  #problems with sex check
  proc = subprocess.Popen(['grep', 'PROBLEM', bfile+'_sexcheck.sexcheck'], stdout=subprocess.PIPE)
  tmp, err = proc.communicate()
  tmp = tmp.decode('ascii')
  tmp_split=tmp.split()
  print(tmp)
  print(tmp_split)
  nproblems=int(len(tmp_split)/6)

  #remove problems
  logfile.write('Check Sex Problems: '+str(nproblems)+'\n')
  if nproblems > 0:
    rfile=open(bfile+'_remove.txt', 'w')
    for k in range(nproblems):
      rfile.write(tmp_split[k*6]+'\t'+tmp_split[1+k*6]+'\n')
    rfile.close()
    os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --remove '+bfile+'_remove.txt --make-bed --out '+bfile+'_sex')
    params['INPUT'] = bfile+'_sex'
    bfile=bfile+'_sex'
    proc = subprocess.Popen(['grep', 'remove:', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    n = int(tmp.split()[1])
    logfile.write('n =  '+str(n)+'\n\n')

#investigate missingness
os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --missing --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+'  --out '+bfile+'_miss')
os.system(params['RPATH']+' '+params['SCRIPTS']+'/missing.R '+bfile+'_miss.imiss '+bfile+'_miss.lmiss '+params['WORK']+ ' ')
#os.system('mv Rplots.pdf '+bfile+'_missing.pdf')
imiss = bfile+'_miss.imiss'

#if n and v not obtained during sex check, do here
if params['SEX'] == 'no':
  #number of samples
  proc = subprocess.Popen(['grep', 'people', bfile+'_miss.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  n = int(tmp.split()[0])
  logfile.write('n = '+str(n)+'\n')

  #number of variants
  proc = subprocess.Popen(['grep', 'variants', bfile+'_miss.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  v = int(tmp.split()[0])
  logfile.write('variants = '+str(v)+'\n\n')

#apply thresholds
if params['MAF'] == 'na':
  params['MAF'] = str(10.0/n)
  logfile.write('MAF = 10/n = '+params['MAF']+'\n')
if 'MIND1' not in params:
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --geno '+params['GENO']+' --maf '+params['MAF']+' --mind '+params['MIND']+' --make-bed --out '+bfile+'_thresh')
  #now bfile and input should always be same so just use bfile
  bfile=bfile+'_thresh'
  #report on QC
  proc = subprocess.Popen(['grep', 'mind)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind = '+str(mind)+'\n')
  if (float(params['GENO']) < 1):
    proc = subprocess.Popen(['grep', 'geno)', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    geno = int(tmp.split()[0])
  else:
    geno = 0
  logfile.write('Variants removed after --geno = '+str(geno)+'\n')
  proc = subprocess.Popen(['grep', 'maf)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  maf = int(tmp.split()[0])
  logfile.write('Variants removed after --maf = '+str(maf)+'\n')
  proc = subprocess.Popen(['grep', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldv = v
  oldn = n
  v = int(tmp.split()[0])
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\nvariants = '+str(v)+'\n')
  #check
  if oldv != (v+maf+geno):
    logfile.write('Check variant QC\n')
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')
#2-tier MIND theshold
else:
#first MIND threshold
  os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --mind '+params['MIND1']+' --make-bed --out '+bfile+'_thresh')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh'
  #report on QC
  proc = subprocess.Popen(['grep', 'mind)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind(1) = '+str(mind)+'\n')
  proc = subprocess.Popen(['grep', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldn = n
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\n')
  #check
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')
#other thresholds
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --geno '+params['GENO']+' --maf '+params['MAF']+' --make-bed --out '+bfile+'_thresh2')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh2'
  #report on QC
  if (float(params['GENO']) < 1):
    proc = subprocess.Popen(['grep', 'geno)', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    geno = int(tmp.split()[0])
  else:
    geno = 0
  logfile.write('Variants removed after --geno = '+str(geno)+'\n')
  proc = subprocess.Popen(['grep', 'minor allele threshold', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  maf = int(tmp.split()[0])
  logfile.write('Variants removed after --maf = '+str(maf)+'\n')
  proc = subprocess.Popen(['grep', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldv = v
  v = int(tmp.split()[0])
  logfile.write('\nvariants = '+str(v)+'\n')
  #check
  if oldv != (v+maf+geno):
    logfile.write('Check variant QC\n')
#2nd tier mind
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --mind '+params['MIND']+' --make-bed --out '+bfile+'_thresh3')
  #now bfile and input should always be same sp just use bfile
  bfile=bfile+'_thresh3'
  #report on QC
  proc = subprocess.Popen(['grep', 'mind)', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  mind = int(tmp.split()[0])
  logfile.write('Samples removed after --mind = '+str(mind)+'\n')
  proc = subprocess.Popen(['grep', 'pass', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  oldn = n
  n = int(tmp.split()[3])
  logfile.write('\nn = '+str(n)+'\n')
  #check
  if oldn != (n+mind):
    logfile.write('Check sample QC\n')


#heterozygosity
if 'HET' in params:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --het --out '+bfile+'_het')
  os.system(params['RPATH']+' '+params['SCRIPTS']+'/het.R '+bfile+'_het.het > '+bfile+'_f.txt '+params['WORK'])
  #os.system('mv Rplots.pdf '+bfile+'_het.pdf')

  logfile.write('\nSamples sorted by f written to '+bfile+'_f.txt\n')
  hfile=open(bfile+'_hremove.txt', 'w')
  #if QC by F
  if 'FMIN' in params and 'FMAX' in params:
    fmin=float(params['FMIN'])
    fmax=float(params['FMAX'])
    fout=open(bfile+'_f.txt', 'r')
    first = 0
    hunder = 0
    hover = 0
    for line in fout:
      if first < 1:
        first = first + 1
      else:
        linesplit=line.split()
        if float(linesplit[6]) < fmin:
          hunder = hunder + 1
          hfile.write(linesplit[1]+'\t'+linesplit[2]+'\n')
        elif float(linesplit[6]) > fmax:
          hover = hover + 1
          hfile.write(linesplit[1]+'\t'+linesplit[2]+'\n')
    fout.close()
    logfile.write(str(hunder)+' sample(s) under FMIN\n'+str(hover)+' sample(s) over FMAX\n')
  #if QC by H
  else:
    hout=open(bfile+'_het.het', 'r')
    first = 0
    hunder = 0
    hover = 0
    hvector = []
    for line in hout:
      if first == 0:
        first = first + 1
      else:
        fid, iid, o_hom, e_hom, n_nm, f = line.split()
        hvector.append((float(n_nm)-float(o_hom))/float(n_nm))
    hout.close()
    hmean = statistics.mean(hvector)
    hstd = statistics.stdev(hvector)
    hout=open(bfile+'_het.het', 'r')
    first = 0
    for line in hout:
      if first == 0:
        first = first + 1
      else:
        fid, iid, o_hom, e_hom, n_nm, f = line.split()
        if ((float(n_nm)-float(o_hom))/float(n_nm)) < (hmean - 3*hstd):
          hunder = hunder +1
          hfile.write(fid+'\t'+iid+'\n')
        if ((float(n_nm)-float(o_hom))/float(n_nm)) > (hmean + 3*hstd):
          hover = hover +1
          hfile.write(fid+'\t'+iid+'\n')
    hout.close()
    logfile.write('H mean: '+str(hmean)+'\nH stdev: '+str(hstd)+'\n'+str(hunder)+' sample(s) under 3 stdev of H\n'+str(hover)+' sample(s) over 3 stdev of H\n')
  hfile.close()
  if hunder+hover > 0:
    os.system(params['PLINKPATH']+' --bfile '+bfile+' --remove '+bfile+'_hremove.txt --make-bed --out '+bfile+'_h')
    bfile=bfile+'_h'
    proc = subprocess.Popen(['grep', 'remove:', bfile+'.log'], stdout=subprocess.PIPE)
    tmp = proc.stdout.read()
    n = int(tmp.split()[1])
    logfile.write('n =  '+str(n)+'\n\n')

#Hardy-Weinberg Equilibrium
os.system(params['PLINKPATH']+' --bfile '+bfile+'  --hardy --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+' --out '+bfile+'_hwe')
os.system('grep "O(HET)\|UNAFF" '+bfile+'_hwe.hwe > '+bfile+'_hwe_ctrl.hwe')
os.system(params['RPATH']+' '+params['SCRIPTS']+'/hwe.R '+bfile+'_hwe_ctrl.hwe '+params['WORK'])
#os.system('mv Rplots.pdf '+bfile+'_hwe.pdf')

hwe=float(params['HWE'])
hweout=open(bfile+'_hwe_ctrl.hwe', 'r')
hwefile=open(bfile+'_hweremove.txt', 'w')
first = 0
hweunder = 0
for line in hweout:
  if first < 1:
    first = first + 1
  else:
    linesplit=line.split()

    try:
      if float(linesplit[8]) < hwe:
        hweunder = hweunder + 1
        hwefile.write(linesplit[1]+'\n')
    except ValueError:
        logfile.write('Value Error in HWE: '+line)
hwefile.close()
hweout.close()
logfile.write('\nVariants with p-value lower than HWE: '+str(hweunder)+'\n')
if hweunder > 0:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --exclude '+bfile+'_hweremove.txt --make-bed --out '+bfile+'_hwe')
  bfile=bfile+'_hwe'
  proc = subprocess.Popen(['grep', 'exclude:', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  v = int(tmp.split()[1])
  logfile.write('variants =  '+str(v)+'\n\n')

#cryptic relatedness
#if list of assayed SNPs given, use this for pruning
if 'ASSAYED' in params:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --extract '+params['ASSAYED']+' --indep-pairwise 1500 150 0.2 --out '+bfile+'_thin')
else:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --indep-pairwise 1500 150 0.2 --out '+bfile+'_thin')
os.system(params['PLINKPATH']+' --bfile '+bfile+' --extract '+bfile+'_thin.prune.in --genome --min 0.05 --out '+bfile+'_genome')
os.system(params['RPATH']+' '+params['SCRIPTS']+'/relate.R '+bfile+'_genome.genome '+params['WORK'])
#os.system('mv Rplots.pdf '+bfile+'_rel.pdf')
relout=open(bfile+'_genome.genome', 'r')
relfile=open(bfile+'_relremove.txt', 'w')
first=0
nrel = []
for line in relout:
  if first < 1:
    first = first + 1
  else:
    linesplit=line.split()
    if float(linesplit[9]) < float(params['PI_HAT']):
      #nrel = nrel + 1
      one = linesplit[1]
      two = linesplit[3]
      proc = subprocess.Popen(['grep', one, imiss], stdout=subprocess.PIPE)
      tmp = proc.stdout.read()
      FMISS1 = float(tmp.split()[5])
      proc = subprocess.Popen(['grep', two, imiss], stdout=subprocess.PIPE)
      tmp = proc.stdout.read()
      FMISS2 = float(tmp.split()[5])
      print(one+' '+str(FMISS1)+'\n'+two+' '+str(FMISS2)+'\n')
      if FMISS1 > FMISS2:
        relfile.write(linesplit[0]+'\t'+one+'\n')
        print(one+' written to file\n')
        if one not in nrel:
          nrel.append(one)
      else:
        relfile.write(linesplit[2]+'\t'+two+'\n')
        print(two+' written to file\n')
        if two not in nrel:
          nrel.append(two)
relfile.close()
relout.close()
logfile.write('Related pairs: '+str(len(nrel))+'\n')
if len(nrel) > 0:
  os.system(params['PLINKPATH']+' --bfile '+bfile+' --remove '+bfile+'_relremove.txt --make-bed --out '+bfile+'_rel')
  bfile=bfile+'_rel'
  proc = subprocess.Popen(['grep', 'remove:', bfile+'.log'], stdout=subprocess.PIPE)
  tmp = proc.stdout.read()
  n = int(tmp.split()[1])
  logfile.write('n =  '+str(n)+'\n\n')

os.system('tar czvf '+params['WORK']+'/final_graphs.tgz '+params['WORK']+'/*.png')
logfile.write('Final bfile: '+bfile+'\nGraphs in final_graphs.tgz\n')

if 'CLEAN' in params:
  #os.system('mv '+bfile+'+.log final.log && mv '+bfile+'.bed final.bed.t && mv '+bfile+'.bim final.bim.t && mv '+bfile+'.fam final.fam.t && rm *.bed && rm *.bim && rm *.fam && mv final.bed.t final.bed && mv final.bim.t final.bim && mv final.fam.t final.fam ')
  os.system('mv '+bfile+'.log '+params['WORK']+'/final.log')
  os.system('mv '+bfile+'.bed '+params['WORK']+'/final.bed')
  os.system('mv '+bfile+'.bim '+params['WORK']+'/final.bim')
  os.system('mv '+bfile+'.fam '+params['WORK']+'/final.fam')
  os.system('rm '+params['INPUT']+'_*')
  print('\nFiles cleaned\n')
