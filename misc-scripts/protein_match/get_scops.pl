use strict;

#Some doc will come
#The scop mapping which should be used here can be found at:
#http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.dom.scop.txt_1.53  (for the 1.53 release)

#perl ../../../src/ensembl-live/misc-scripts/protein_match/get_scops.pl -scop ../secondary/dir.dom.scop.txt_1.53 -dbmap mapdb.map -out scops.map

use Getopt::Long;

my ($scop,$dbmap,$out);

my %map;


&GetOptions(
	    'scop:s'=>\$scop,
	    'out:s'=>\$out,
	    'dbmap:s'=>\$dbmap
            );

open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n"; 
open (SCOP,"$scop") || die "Can't open file $scop\n";
open (OUT,">$out") || die "Can't open file $out\n";
    open (ERROR, ">scop.err") || die "Can't open file scop.err\n";

print STDERR "Reading dbmap\n";

while (<DBMAP>) {
    chomp;
    my ($mapac,$mapdb) = split(/\t/,$_);
    
    $map{$mapac} = $mapdb;
}

print STDERR "Reading Scop file\n";

while (<SCOP>) {
    chomp;
    my ($scopac, $pdb, $chain, $scopnb) = split(/\t/,$_);
    
	if (!defined $map{$scopac}) {
	    print ERROR "can't map $scopac\n";
	}

    print OUT "$map{$scopac}\t$scopac\tSCOP1\t$pdb\|\|$chain\n";
    print OUT "$map{$scopac}\t$scopac\tSCOP1\t$scopnb\n";
}
