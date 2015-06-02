import os, sys

plink='/home/sel228/dave/plink'

pheno=sys.argv[1]

#gwas='/home/sel228/dave/QC/1000G/adc_combined_notQCfirst/combine_thresh_thresh2_thresh3_h_hwe_rel'
#pfile='/home/sel228/dave/pheno_dave.txt'

gwas = sys.argv[2]
pfile=sys.argv[3]

#outliers=sys.argv[4]

os.system(plink+' --bfile '+gwas+' --out '+pheno+'_assoc --assoc --pheno '+pfile+' --pheno-name '+pheno)

os.system('/home/sel228/R-3.1.2/bin/Rscript /home/sel228/dave/association/scripts/man.R '+pheno+'_assoc.*assoc')

os.system('mv Rplots.pdf '+pheno+'.pdf')
