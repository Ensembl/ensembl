
#
# Ensembl module for Bio::EnsEMBL::Utils::GTF_Merge
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Utils::GTF_Merge - Module having the GTF merge subroutines

=head1 SYNOPSIS

    use Bio::EnsEMBL::Utils::GTF_Merge('gtf_merge');

    &gtf_merge($input_file1,$input_file2,$output_file,$new_gene_prefix);
 
=head1 DESCRIPTION

Module containing the sub routine gtf_merge, which works on files to
provide merging between two gtf files. It has to work file based as
there needs to be a 2 pass approach.

The merge is done on overlapping exon lines

The end result is a new file in $output file which has the gene name
of input1 and input2 replaced by the "new_gene_prefix_unique number"
appropiate to the merge.

This only works if the transcript identifiers in the two files are
distinct, but does handle unsorted files fine.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Utils::GTF_Merge;
use vars qw(@ISA @EXPORT_OK);
use strict;
use Exporter;
use Carp;

@ISA= ('Exporter');

@EXPORT_OK = ( 'gtf_merge');


sub gtf_merge {
    my ($infile1,$infile2,$outfile,$new_prefix) = @_;


    if( !defined $new_prefix ) {
	&confess("Did not pass in new_prefix - must be infile1,infile2,outfile,new_prefix");
    }

    if( !-e $infile1 ) { &confess("$infile1 does not exist!") }
    if( !-e $infile2 ) { &confess("$infile2 does not exist!") }
    open(OUT,">$outfile") || &confess("Unable to open outfile $!");
    
    # this sorts all the exon lines first by ctg, then strand, then start position
    open(SORT,"cat $infile1 $infile2 | grep 'exon' | sort +0 +6 +4n |") || &confess("Unable to spawn pipe!");

    my ($prevctg,$prevstrand,$prevtrans,$prevend);
    my %thash;
    my $unique_number = 1;
    

    # main loop. Only process things with transcript_id lines
    # see whether previous line is contig or strand different, in which case,
    # simply move to the next line.
    
    # otherwise see if there is an overlap. %thash is a hash of transcript_ids
    # pointing to new gene ids. The gene_ids start being called as 'single'.
    # and get reassigned when overlaps occur

    while( <SORT> ) {
	/transcript_id\s+"(\S+)"/ || next;
	my $trans = $1;
	if( !defined $thash{$trans} ) {
	    $thash{$trans} = 'single';
	}

	my ($ctg,$source,$tag,$start,$end,$score,$strand) = split;
	if( !defined $prevctg || $ctg ne $prevctg || $strand ne $prevstrand) {
	    # new contig/strand/file
	    $prevctg   = $ctg;
	    $prevend = $end;
	    $prevtrans   = $trans;
	    $prevstrand = $strand;
	    next;
	} 
	# else potential merge
	  
	if( $start <= $prevend ) {
	    # merge candidate
	    if( $thash{$prevtrans} eq 'single' ) {
		# new cluster between prevtrans and this one!
		my $new_gene = $new_prefix . "_" . $unique_number;
		$unique_number++;
		$thash{$prevtrans} = $new_gene;
		$thash{$trans}     = $new_gene;
	    } else {
		# existing cluster
		$thash{$trans} = $thash{$prevtrans};
	    }
	}

	$prevctg    = $ctg;
	$prevend    = $end;
	$prevtrans  = $trans;
	$prevstrand = $strand;    

    }


    close(SORT) || &confess("Pipe did not open successfully, $!");

    # assign all 'single' transcripts their own separate gene

    foreach my $t ( keys %thash ) {
	if( $thash{$t} eq 'single' ) {
	    my $new_gene = $new_prefix . "_" . $unique_number;
	    $unique_number++;
	    $thash{$t} = $new_gene;
	}
    }

    # dump results for the moment

    foreach my $t ( keys %thash ) {
	print STDERR "$thash{$t}\t$t\n";
    }

    # process final output

    open(CAT,"cat $infile1 $infile2 |") || &confess("Could not spawn pipe");

    while( <CAT> ) {
	/^#/ && do { print OUT $_; next; };

	my ($trans,$gene);
	/gene_id\s+"(\S+)"/ && do { $gene = $1; };
	/transcript_id\s+"(\S+)"/ && do { $trans = $1; };

	if( !defined $trans || !defined $gene ) {
	    print STDERR "Line has no transcript or gene id. Cannot merge";
	    next;
	}

	my $new_gene = $thash{$trans};

	s/gene_id\s+"$gene"/gene_id "$new_gene"/;

	print OUT $_;
    }
	
    close(CAT) || &confess("cat pipe did not open successfully $!");


}

    






