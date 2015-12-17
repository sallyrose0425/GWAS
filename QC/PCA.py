import sys, os

pfile=open(sys.argv[1], 'r')

params = {}
for line in pfile:
  paramname, paramvalue = line.split()
  params[paramname[1:]]=paramvalue

os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --extract '+params['ASSAYED']+' --indep-pairwise 1500 150 0.2 --out '+params['WORK']+'/thin')
os.system(params['PLINKPATH']+' --bfile '+params['INPUT']+' --extract '+params['WORK']+'/thin.prune.in --out '+params['WORK']+'/smartpca --pheno '+params['PFILE']+' --pheno-name '+params['PHENO']+' --recode')

parf=open(params['WORK']+'/par.txt', 'w')
parf.write('genotypename:    '+params['WORK']+'/smartpca.ped\n')
parf.write('snpname:         '+params['WORK']+'/smartpca.map\n')
parf.write('indivname:       '+params['WORK']+'/smartpca.ped\n')
parf.write('evecoutname:     '+params['WORK']+'/smartpca.evec\n')
parf.write('evaloutname:     '+params['WORK']+'/smartpca.eval\n')
parf.write('altnormstyle:    NO\n')
parf.write('familynames:     NO\n')
parf.write('grmoutname:      '+params['WORK']+'/grmjunk\n')
parf.close()

os.system(params['SMARTPCA']+' -p '+params['WORK']+'/par.txt > '+params['WORK']+'/smartpca.out')
os.system("grep REMOVED "+params['WORK']+"/smartpca.out | awk '{print $3}' | sed 's/^/0 /g' > "+params['WORK']+"/remove.txt")

