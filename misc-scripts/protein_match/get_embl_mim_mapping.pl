use strict;

#Some doc will come ...

#perl ../../../src/ensembl-live/misc-scripts/protein_match/get_embl_mim_mapping.pl -sp ../primary/hum_sp_sptrembl.pep -dbmap mapdb.map -output sp_embl_mim.map

use Getopt::Long;
use Bio::SeqIO;

my ($sp,$dbmap,$out);

my %map;


&GetOptions(
	    'sp:s'=>\$sp,
	    'output:s'=>\$out,
	    'dbmap:s'=>\$dbmap
            );

open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n"; 
open (OUT,">$out") || die "Can't open file\n";

print STDERR "Reading DBmap\n";

while (<DBMAP>) {
    chomp;
     my ($mapac,$mapdb) = split(/\t/,$_);
     $map{$mapac} = $mapdb;
}

print STDERR "Reading SP\n";

my $in1  = Bio::SeqIO->new(-file => $sp, '-format' =>'swiss');

while ( my $seq1 = $in1->next_seq() ) {
    my $ac = $seq1->accession;
    my @dblink = $seq1->annotation->each_DBLink;
    
    
    foreach my $link(@dblink) {
	if ($link->database eq "EMBL") {
	    if (!defined $map{$ac}) {
		die "Can't map $ac\n";
	    }
	    
	    print OUT "$map{$ac}\t$ac\tEMBL_AC\t".$link->primary_id,"\n";

	    my ($protac) = $link->optional_id =~ /^(\w+).\S+/;
	    print OUT "$map{$ac}\t$ac\tEMBL_PROT_AC\t$protac\n";
	}
	if  ($link->database eq "MIM") {
	    if (!defined $map{$ac}) {
		die "Can't map $ac\n";
	    }
	    
	    print OUT "$map{$ac}\t$ac\tOMIM\t".$link->primary_id,"\n";
	}
	
	
    }


}














