use strict;

#Some doc will come

use Getopt::Long;

my ($refseq,$dbmap,$out);

my %map;


&GetOptions(
	    'refseq:s'=>\$refseq,
	    'out:s'=>\$out,
	    'dbmap:s'=>\$dbmap
            );

open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n"; 
open (REFSEQ,"$refseq") || die "Can't open file $refseq\n";
open (OUT,">$out") || die "Can't open file $out";

print STDERR "Reading dbmap\n";

while (<DBMAP>) {
    chomp;
    my ($mapac,$mapdb) = split(/\t/,$_);
    
    $map{$mapac} = $mapdb;
}


#Separate by entry (each entry goes into $_)
$/ = "\/\/\n";

print STDERR "Reading Refseq file\n";

while (<REFSEQ>) {
    my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
    my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;
    
       
    my ($mim) = $_ =~ /\/db_xref=\"MIM:(\d+)/;
    my ($locus) = $_ =~ /\/db_xref=\"LocusID:(\d*)/;

    if ($mim) {
	if (!defined $map{$dna_ac}) {
	    die "can't map $dna_ac\n";
	}
	print OUT "$map{$dna_ac}\t$dna_ac\tOMIM\t$mim\n";
    }

    if ($locus) {
	if (!defined $map{$dna_ac}) {
	    die "can't map $dna_ac\n";
	}
	print OUT "$map{$dna_ac}\t$dna_ac\tLOCUS\t$locus\n";
    }
}







