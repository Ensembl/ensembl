use strict;

=head1 get_all_external

=head2 Description

This script will set up a mapping of each primary database (SP, SPTREMBL, Refseq and PDB) to their respective name. This mapping will be then used to give a database name to each ACs.

=head2 Options

-sp
-sptrembl
-refseq
-pdb

=head2 Contact

mongin@ebi.ac.uk
birney@ebi.ac.uk

=cut

#perl ../../../src/ensembl-live/misc-scripts/protein_match/get_dbmap.pl -sp ../primary/01_04_hum_sprot.txt -sptrembl ../primary/07_04_hum_sptrembl.txt -refseq ../primary/hs.gnp -scop ../primary/scop.fas

use Getopt::Long;
use Bio::SeqIO;

my ($sp, $sptrembl, $refseq, $scop, $out);

&GetOptions(
	    'sp:s'=>\$sp,
	    'sptrembl:s'=>\$sptrembl,
	    'refseq:s'=>\$refseq,
	    'scop:s'=>\$scop,
	    'out:s'=>\$out
	    );

open (OUT,">mapdb.map");

print STDERR "Getting SP mapping\n";

my $in  = Bio::SeqIO->new(-file => $sp, '-format' =>'swiss');

while ( my $seq = $in->next_seq() ) {
    print OUT $seq->accession."\tSP\n";
    
    my @secs = $seq->get_secondary_accessions;
    if (@secs) {
	foreach my $sec (@secs) {
	    print OUT "$sec\tSP\n";
	}
    }
    
}

print STDERR "Getting SPTREMBL mapping\n";

my $in1  = Bio::SeqIO->new(-file => $sptrembl, '-format' =>'swiss');

while ( my $seq1 = $in1->next_seq() ) {
    print OUT $seq1->accession."\tSPTREMBL\n";
    
    my @secs = $seq1->get_secondary_accessions;
    if (@secs) {
	foreach my $sec (@secs) {
	    print OUT "$sec\tSPTREMBL\n";
	}
    }
    
}

print "Getting refseq mapping\n";

open (REFSEQ,"$refseq") || die "Can't open file\n";

$/ = "\/\/\n";

while (<REFSEQ>) {
    my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
    my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;
    
    if ($dna_ac) {
	print OUT "$dna_ac\tREFSEQ\n";
    }
}
close (REFSEQ);

$/ = "\n";

print "Getting Scop mapping\n";

open (SCOP,"$scop") || die "Can't open file\n";

while (<SCOP>) {
    my ($scop_ac) = $_ =~ /^>(\S+)/;
    if ($scop_ac) {
	print OUT "$scop_ac\tSCOP\n";
    }
}
close (SCOP);
    
close (OUT);









