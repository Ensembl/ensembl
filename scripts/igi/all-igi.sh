#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$
# 
# run the complete IGI show. May not be all that useful (since often I have
# to repair things by hand ... in fact, a Makefile would be more convenient),
# but at least can serve as some kind of documentation of what's supposed to
# happen. See also README in this directory, of course.

# prepare the sources into a form that can be merged; results go to
# subdir out/collated
all-collates.sh > all-collates.out 2>all-collates.log

# do the merge; results go to subdir out/merged
all-merges.sh > all-merges.out 2>all-merges.log

# run statistics on the merges, and also produce the summary and final gtf
# file(s). Results go to  out/{stats,mapping,summary,final}
all-stats.sh > all-stats.out 2>all-stats.log

# do the peptide business; first, collate them:
pep-collate.sh 2> pep-collate.log

# Now pull out the longest peptide of each original igi
# peptide file, using the 'valid' igi's from the summary files (which are
# proudced by all-stats.sh.
summary=out/summary/ens_affy_fgenesh.summary
peptidefiles="ensembl/ensembl.pep affymetrix/affymetrix.pep fgenesh/fgenesh.pep"
outdir=out/pep
outfile=$outdir/igi3.pep
logfile=$outdir/igi3.log
gtfsummary2pep.pl $summary $peptidefiles > $outfile 2> $logfile
gzip -c $outfile > $outfile.gz
# ftp stuff, and we're done.

