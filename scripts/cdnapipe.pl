#!/usr/local/bin/perl

use strict;

=head1 NAME

cdnapipe.pl

=head1 DESCRIPTION

script that interacts with Tims directories to
run a sensible set of cdna-s against genomic DNA
in the system

=cut

use Bio::EnsEMBL::Analysis::CdnaResolver;
use Bio::SeqIO;

my $file = shift;

$file || system("perldoc $0");


my ($dir,$clone,@files);

if( $file=~ /(.*)\/([^\/]+)\.seq/) {
    $dir=$1;
    $clone=$2;
    local *TMPDIR;
    opendir(TMPDIR,$dir);
    @files=readdir(TMPDIR);
    closedir(TMPDIR);
} else {
    die  "cannot parse $file as timdb seq file";
}

my $resolver = Bio::EnsEMBL::Analysis::CdnaResolver->new();

my $seqin = Bio::SeqIO->new( -format => 'Fasta' , -file => $file );

while( my $seq = $seqin->next_seq() ) {
    my @embl_acc;
    my $contig_id = $seq->id();
    
    # we have already extracted all the files for this directory. Pickout the 
    # guys which hit peptides for this contig.
    
    my %embl_acc_hash;

    foreach my $tfile ( @files ) {
	if( $tfile =~ /$contig_id\.(\d+)\.tblastn_vert\.msptmp$/ ) {
	    open(MSP,"$dir/$tfile");
	    while( <MSP> ) {
		!/Homo/ && next; # skip non human entries.
		s/^\s+//g; # trim leading white space.
		my @fields = split;
		$embl_acc_hash{$fields[8]} = 1;
		#print STDERR "Adding $fields[8]\n";
	    }
	    close(MSP);
	}
    }

    @embl_acc = keys %embl_acc_hash;

    # bug out if no embl accessions

    if( scalar @embl_acc == 0 ) { next; }

    # call into cdna resolver.

    my $aseq = Bio::AnnSeq->new();
    $aseq->seq($seq);
    
    foreach my $e ( @embl_acc ) {
	print "seen $e\n";
    }

    my @sqfeatures = $resolver->resolve($aseq,shift @embl_acc);

    $aseq->add_SeqFeature(@sqfeatures);

    $aseq->write_GFF;
}


    
    
    
