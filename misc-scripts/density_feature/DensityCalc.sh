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
outdir=             # -o


while getopts 'h:P:d:u:p:o:' opt; do
  case ${opt} in
    h)  host=${OPTARG} ;;
    P)  port=${OPTARG} ;;
    d)  dbname=${OPTARG} ;;
    u)  user=${OPTARG} ;;
    p)  pass=${OPTARG} ;;
    o)  outdir=${OPTARG} ;;
  esac
done

if [[ -z ${host} || -z ${port} || -z ${dbname} || -z ${user} || -z ${pass} ]]
then
  print -u2 "Usage:\n\t$0 -h host -P port -d database -u user -p password"
  exit 1
fi

# Enter an output dir so results don't get written in the same place as your script
if [[ -z ${outdir} || ! -d ${outdir} ]]
then
  outdir=$PWD
fi

echo "\nOutput dir is '${outdir}'\n"

# Make sure this is a core database.
if [[ -n ${dbname##*_core_*} ]]; then
  print -u2 "The database '${dbname}' is not a core database"
  exit 1
fi

print "Submitting percent GC calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_gc.out'"
bsub -q normal -J gc_calc -oo ${outdir}/${dbname}_gc.out \
  perl ./percent_gc_calc.pl \
  -host ${host} \
  -port ${port} \
  -user ${user} \
  -pass ${pass} \
  -dbname ${dbname}

print "Submitting gene density calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_gene.out'"
bsub -q normal -J gene_density -oo ${outdir}/${dbname}_gene.out \
  perl ./gene_density_calc.pl \
  -host ${host} \
  -port ${port} \
  -user ${user} \
  -pass ${pass} \
  -dbname ${dbname}

print "Submitting repeat coverage calculation to queue 'long'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_repeat.out'"
bsub -q long -J repeat_cov -oo ${outdir}/${dbname}_repeat.out \
  perl ./repeat_coverage_calc.pl \
  -host ${host} \
  -port ${port} \
  -user ${user} \
  -pass ${pass} \
  -dbname ${dbname}

print "Submitting variation density calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_var.out'"
bsub -q normal -J var_density -oo ${outdir}/${dbname}_var.out \
  perl ./variation_density.pl \
  -host ${host} \
  -port ${port} \
  -user ${user} \
  -pass ${pass} \
  -dbname ${dbname}

print "Submitting seq region stats calculation to queue 'normal'"
print "\tThe output from this job goes to the file"
print "\t'${dbname}_seqreg.out'"
bsub -q normal -J seqreg_stats -oo ${outdir}/${dbname}_seqreg.out \
  perl ./seq_region_stats.pl \
  -host ${host} \
  -port ${port} \
  -user ${user} \
  -pass ${pass} \
  -dbname ${dbname}

print "All jobs submitted."

# $Id$
