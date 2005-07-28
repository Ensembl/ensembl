# Parse regulatory feature data from a CisRED file into a file
# in a format suitable for reading into an ensembl database using
# load_regulatory.pl

# Cisred data is available under http://www.cisred.org/files/Databases
# e.g. http://www.cisred.org/files/Databases/cisREDdb-Hsap-1.1b/SQL/data.zip

# File format of cisred features.txt file (tab separated)
# id batch_id seqname source feature start end score strand frame ensembl_gene_id


use strict;

use DBI;
use Getopt::Long;

my ($infile, $outfile);

GetOptions( "infile=s",   \$infile,
	    "outfile=s",  \$outfile,
	    "help",       \&usage);

usage() if (!$infile || !$outfile);

print "Reading from $infile\n";
my $count = 0;

open (INFILE,  "<$infile") || die "Can't open $infile";
open (OUTFILE, ">$outfile") || die "Can't open $outfile";

while (<INFILE>) {

  my ($id,$batch_id,$seqname,$source,$feature,$start,$end,$score,$strand,$frame,$ensembl_gene_id) = split;

  print OUTFILE join ("\t", "Similarity", $feature, $source, "transcription_factor", $seqname, $start, $end, $strand, $frame,$score,"gene","id","\"$ensembl_gene_id\"" ) . "\n";

  # XXX transcription_factor?

  $count++;

}

close INFILE;
close OUTFILE;

print "Wrote $count regulatory features to $outfile\n";

sub usage {

  print <<EOF;
Usage: perl parse_cisred.pl
         -infile    : Cisred features.txt file
         -outfile   : File to be written for load_regulatory.pl

Data files available from e.g.

http://www.cisred.org/files/Databases/cisREDdb-Hsap-1.1b/SQL/data.zip

EOF

  exit(0);

}
