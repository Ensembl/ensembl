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

    &gtf_merge($input_file1, $output_file, $new_gene_prefix);
 
=head1 DESCRIPTION

Module containing the sub routine gtf_merge, which works on one sorted stream
containing gtf features from several different sources.

The merge is done on overlapping exon lines

The end result is a new file in $output file which has the gene name
of input1 and input2 replaced by the "new_gene_prefix_unique number"
appropiate to the merge.

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
 Function: Merges genes on the basis of overlapping
           exons and writes out the merged file to output
           The genes are given the new numbering scheme based on
           the prefix written above
 Example :
 Returns : nothing
 Args    : read stream, MUST BE SORTED by contig, strand, start
             (sort -k1,1 -k7,7 -k4,4n in UNIX commands); if not sorted,
             module will gtf_merge will die with an error message.  
           output stream
           prefix (optional)
           log output stream for intermediate data

=cut

sub gtf_merge {
    my ($sort,$out,$new_prefix,$log) = @_;

    if( !defined $new_prefix ) {
	&confess("Did not pass in new_prefix - must be sort,out,new_prefix");
    }

    my ($prevctg,$prevstrand,$prevgene,$prevend,$prevstart);
    my %thash;
    my %seenctg;
    my %ghash;
    my $unique_number = 1;
    

    # main loop. Only process things with gene_id lines
    # see whether previous line is contig or strand different, in which case,
    # simply move to the next line.

    # save each line in @lines
    
    # each gene creates an igi when first seen
    # when an overlap occurs, all the igis assigned to the later
    # gene (may well be more than one, though one is likely to
    # to be the most often case) get transfered to previous overlapped igi
    
    my @lines;

    while( <$sort> ) {
	push(@lines,$_);
	my ($ctg,$source,$tag,$start,$end,$score,$strand) = split;
	if( $tag ne 'exon' ) { next; }
#	push(@lines,$_);

	if( $strand != '-' && $strand != '+' ) {
	    &confess("Strand is neither + or -. A GTF file error? At $_");
	}

	/gene_id\s+(\S+)/ || next;
	my $gene = $1;
	$gene =~ s/\"//g;
	$gene =~ s/;$//g;
	
	if( !defined $thash{$gene} ) {
	    my $new_igi = "temp_" . $unique_number;
	    $unique_number++;
	    $thash{$gene} = $new_igi; # one igi of gene
	    $ghash{$new_igi} = [];     
	    push(@{$ghash{$new_igi}},$gene); # all genes of igi
	}

	if( !defined $prevctg || $ctg ne $prevctg || $strand ne $prevstrand) {
	    # new contig or strand
	    if( $strand eq $prevstrand && $seenctg{$ctg} ) {
		&confess("Not sorted GTF file - Seen $ctg already but in a between contig move!");
	    }
	    $seenctg{$ctg} = 1;
	    $prevctg   = $ctg;
	    $prevend = $end;
	    $prevgene   = $gene;
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
	    # merge $gene into $prevgene
	    if( $log ) {
		print $log "$prevgene\t$gene\t$prevstart:$prevend\t$start:$end\n";
	    }
	    my $combined_igi  =  $thash{$prevgene};
	    my $dead_igi      = $thash{$gene};

	    # of course, this could be the second time we
	    # see this merge ;)
	    if( $combined_igi ne $dead_igi ) {
		
		foreach my $t ( @{$ghash{$dead_igi}} ) {
	
		    # move it across to the prevgene igi
		    push(@{$ghash{$combined_igi}},$t);
		    # reset this guys thash
		    $thash{$t} = $combined_igi;
		} 
		# delete igi from $gene
		delete $ghash{$dead_igi};
	    }
	}

	if( $prevend <= $end ) {
	    $prevend    = $end;
	    $prevgene  = $gene;
	}
	$prevstart = $start;

    }


    close($sort) || &confess("Unable to close sorted GTF file stream $!");


    # dump results for the moment
    #
    foreach my $t ( keys %thash ) {
    	print STDERR "$thash{$t}\t$t\n";
    }

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
	
	my ($gene);
	/gene_id\s+(\S+)/ && do { $gene = $1; };


	if( !defined $gene ) {
	    print STDERR "Line with no gene_id: $_\n";
	    print $out $_;
	    next;
	} 
	$gene =~ s/\"//g;
	$gene =~ s/;$//g;

	my $igi = $thash{$gene};
	s/gene_id/igi_id "$igi"; gene_id/;

 	print $out $_;
#	print $out $_ if 
    }



}

    






