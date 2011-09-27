#!/bin/ksh

# This script submits jobs to the farm to calculate the various density
# features for a particular core database.

# If an output dir is not specified, output from bjobs will be written to PWD.

# Default values for command line switches:

host='ens-staging'  # -h
port='3306'         # -P
user='ensadmin'     # -u
pass=               # -p
dbname=             # -d
species=            # -s
outdir=             # -o


while getopts 'h:P:d:s:u:p:o:' opt; do
  case ${opt} in
    h)  host=${OPTARG} ;;
    P)  port=${OPTARG} ;;
    d)  dbname=${OPTARG} ;;
    s)  species=${OPTARG};;
    u)  user=${OPTARG} ;;
    p)  pass=${OPTARG} ;;
    o)  outdir=${OPTARG} ;;
  esac
done

if [[ -z ${host} || -z ${port} || -z ${dbname} || -z ${species} || -z ${user} || -z ${pass} ]]
then
  print -u2 "Usage:\n\t$0 -h host -P port -d database -s species -u user -p password"
  exit 1
fi

# Enter an output dir so results don't get written in the same place as your script
if [[ -z ${outdir} || ! -d ${outdir} ]]
then
  outdir=$PWD
fi

print "Output dir is '${outdir}'"

# Make sure this is a core database.
if [[ -n ${dbname##*_core_*} ]]; then
  print -u2 "The database '${dbname}' is not a core database"
  exit 1
fi

print "Submitting percent GC calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_gc.out'"
bsub -q normal -J gc_calc \
  -oo ${outdir}/${dbname}_gc.out \
  -eo ${outdir}/${dbname}_gc.err \
  perl ./percent_gc_calc.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -d ${dbname}

print "Submitting gene density calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_gene.out'"
bsub -q normal -J gene_density \
  -oo ${outdir}/${dbname}_gene.out \
  -eo ${outdir}/${dbname}_gene.err \
  perl ./gene_density_calc.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -d ${dbname}

print "Submitting repeat coverage calculation to queue 'long'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_repeat.out'"
bsub -q long -J repeat_cov \
  -oo ${outdir}/${dbname}_repeat.out \
  -eo ${outdir}/${dbname}_repeat.err \
  perl ./repeat_coverage_calc.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -d ${dbname}

print "Submitting variation density calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_var.out'"
bsub -q normal -J var_density \
  -oo ${outdir}/${dbname}_var.out \
  -eo ${outdir}/${dbname}_var.err \
  perl ./variation_density.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -s ${species}

print "Submitting seq region gene stats calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_seqreg.out'"
bsub -q normal -J seqreg_stats_gene \
  -oo ${outdir}/${dbname}_seqreg_gene.out \
  -eo ${outdir}/${dbname}_seqreg_gene.err \
  perl ./seq_region_stats.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -d ${dbname} \
  -s gene

print "Submitting seq region snp stats calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_seqreg.out'"
bsub -q normal -J seqreg_stats_snp \
  -oo ${outdir}/${dbname}_seqreg_snp.out \
  -eo ${outdir}/${dbname}_seqreg_snp.err \
  perl ./seq_region_stats.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -d ${dbname} \
  -s snp

print "Submitting gene_gc content calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_genegc.out'"
bsub -q normal -J genegc_stats \
  -oo ${outdir}/${dbname}_genegc.out \
  -eo ${outdir}/${dbname}_genegc.err \
  perl ../gene_gc.pl \
  -h ${host} \
  -port ${port} \
  -u ${user} \
  -p ${pass} \
  -pattern ${dbname}

print "All jobs submitted."

# $Id$
