
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

=head2 gtf_merge

 Title   : gtf_merge
 Usage   : &gtf_merge(\*INPUT,\*OUTPUT,'IGI_M3_');
 Function: Merges transcripts on the basis of overlapping
           exons and writes out the merged file to output
           The genes are given the new numbering scheme based on
           the prefix written above
 Example :
 Returns : nothing
 Args    : read stream, MUST BE SORTED by contig, strand, start
           (sort -k1,1 -k7,7 -k4,4n in UNIX commands)
           output stream
           prefix


=cut

sub gtf_merge {
    my ($sort,$out,$new_prefix) = @_;


    if( !defined $new_prefix ) {
	&confess("Did not pass in new_prefix - must be sort,out,new_prefix");
    }

    my ($prevctg,$prevstrand,$prevtrans,$prevend,$prevstart);
    my %thash;
    my %seenctg;
    my %ghash;
    my $unique_number = 1;
    

    # main loop. Only process things with transcript_id lines
    # see whether previous line is contig or strand different, in which case,
    # simply move to the next line.

    # save each line in @lines
    
    # each transcript creates an igi when first seen
    # when an overlap occurs, all the igis assigned to the later
    # transcript (may well be more than one, though one is likely to
    # to be the most often case) get transfered to previous overlapped igi
    
    my @lines;

    while( <$sort> ) {
	push(@lines,$_);
	my ($ctg,$source,$tag,$start,$end,$score,$strand) = split;
	if( $tag ne 'exon' ) { next; }

	/transcript_id\s+(\S+)/ || next;
	my $trans = $1;
	$trans =~ s/\"//g;
	$trans =~ s/;$//g;
	
	
	if( !defined $thash{$trans} ) {
	    my $new_igi = "temp_" . $unique_number;
	    $unique_number++;
	    $thash{$trans} = $new_igi;
	    $ghash{$new_igi} = [];
	    push(@{$ghash{$new_igi}},$trans);
	}

	if( !defined $prevctg || $ctg ne $prevctg || $strand ne $prevstrand) {
	    # new contig/strand/file
	    if( $strand eq $prevstrand && $seenctg{$ctg} ) {
		&confess("Not sorted GTF file - Seen $ctg already but in a between contig move!");
	    }
	    $seenctg{$ctg} = 1;
	    $prevctg   = $ctg;
	    $prevend = $end;
	    $prevtrans   = $trans;
	    $prevstrand = $strand;
	    $prevstart  = $start;
	    next;
	} 
	# else potential merge
	  
	# non sorted catcher
	if( $start < $prevstart ) {
	    &confess("Not sorted GTF file - start $ctg:$start is smaller than previous line of $ctg:$prevstart");
	}


	if( $start <= $prevend ) {
	    # merge $trans into $prevtrans

	    my $combined_igi  =  $thash{$prevtrans};
	    my $dead_igi      = $thash{$trans};

	    # of course, this could be the second time we
	    # see this merge ;)
	    if( $combined_igi ne $dead_igi ) {
		
		foreach my $t ( @{$ghash{$dead_igi}} ) {
	
		    # move it across to the prevtrans igi
		    push(@{$ghash{$combined_igi}},$t);
		    # reset this guys thash
		    $thash{$t} = $combined_igi;
		} 
		# delete igi from $trans
		delete $ghash{$dead_igi};
	    }
	}

	if( $prevend <= $end ) {
	    $prevend    = $end;
	    $prevtrans  = $trans;
	}
	$prevstart = $start;

    }


    close($sort) || &confess("Unable to close sorted GTF file stream $!");


    # dump results for the moment
    #
    #foreach my $t ( keys %thash ) {
    #	print STDERR "$thash{$t}\t$t\n";
    #}

    # reassign unique numbers
    $unique_number = 1;
    my %finaligi;
    
    foreach my $igi ( keys %ghash ) {
	my $new_igi = $new_prefix . "_" . $unique_number;
	$unique_number++;
	# point final to old array
	$finaligi{$new_igi} = $ghash{$igi};
	# loop over old array, switching thash to new
	foreach my $t ( @{$finaligi{$new_igi}} ) {
	    $thash{$t} = $new_igi;
	}
    }


    
    foreach ( @lines ) {
	/^#/ && do { print $out $_; next; };
	
	my ($trans,$gene);
	/transcript_id\s+(\S+)/ && do { $trans = $1; };


	if( !defined $trans ) {
	    print STDERR "Line starting [",substr($_,0,50),"] has no transcript id. Cannot provide IGI\n";
	    print $out $_;
	    next;
	}

	$trans =~ s/\"//g;
	$trans =~ s/;$//g;

	my $igi = $thash{$trans};
	s/transcript_id/igi_id "$igi"; transcript_id/;

	print $out $_;
    }



}

    






