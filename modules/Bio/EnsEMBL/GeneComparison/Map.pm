
#
# Ensembl module for Bio::EnsEMBL::GeneComparison::Counts
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GeneComparison::Map - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::GeneComparison::Map;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

    my $self = {};
    bless $self,$class;

    return $self;
}


=head2 process_gene_list

 Title   : process_gene_list
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub process_gene_list{
   my ($self,$standard,$predicted) = @_;

   # loop over standard and predicted, calculating start/end points
   # to store in hash. This is going to make subsequent processing much
   # faster. At the same time, sum up all totals

   # parallel hashes 
   my (%start,%end,%strand);

   foreach my $gene ( @$standard ) {
       my $start = 1000000000000000;
       my $end   = 1;
       if( !defined $gene->id ) {
	   $self->throw("Cannot process genes without ids in the genes - and they have to be distinct!");
       }

       $self->total_standard_genes($self->total_standard_gene+1);
       foreach my $trans ( $gene->each_Transcript ) {
	   foreach my $exon ( $trans->each_Exon ) {
	       $self->total_standard_exons($self->total_standard_exon+1);
	       $self->total_standard_basepairs($self->total_standard_basepairs + $exon->length);

	       if( $exon->start < $start ) {
		   $start = $exon->start;
	       }
	       if( $exon->end > $end ) {
		   $end = $exon->end;
	       }

	       if( !defined $strand{$gene->id} ) {
		   $strand{$gene->id} = $exon->strand;
	       }
	   }
       }
       $start{$gene->id} = $start;
       $end{$gene->id} = $end;
   }

   # do the same for predicted

   foreach my $gene ( @$predicted ) {
       my $start = 1000000000000000;
       my $end   = 1;
       if( !defined $gene->id ) {
	   $self->throw("Cannot process genes without ids in the genes - and they have to be distinct!");
       }

       $self->total_predicted_genes($self->total_standard_gene+1);
       foreach my $trans ( $gene->each_Transcript ) {
	   foreach my $exon ( $trans->each_Exon ) {
	       $self->total_predicted_exons($self->total_standard_exon+1);
	       $self->total_predicted_basepairs($self->total_standard_basepairs + $exon->length);

	       if( $exon->start < $start ) {
		   $start = $exon->start;
	       }
	       if( $exon->end > $end ) {
		   $end = $exon->end;
	       }

	       if( !defined $strand{$gene->id} ) {
		   $strand{$gene->id} = $exon->strand;
	       }
	   }
       }
       $start{$gene->id} = $start;
       $end{$gene->id} = $end;
   }


   # set up exon exact/non overlap


   my (%std_exact,%pre_exact,%std_over,%pre_over,%std_gene_exact,%std_gene_over,%pre_gene_exact,%pre_gene_over);

   foreach my $gene ( @$standard) {
       $std_gene_exact{$gene->id} = 0;
       $std_gene_over{$gene->id} = 0;

       foreach my $exon ( $gene->each_unique_Exon ) {
	   $std_exact{$exon->id} = 0;
	   $std_over{$exon->id} = 0;
       }
   }

   foreach my $gene ( @$predicted) {
       $pre_gene_exact{$gene->id} = 0;
       $pre_gene_over{$gene->id} = 0;

       foreach my $exon ( $gene->each_unique_Exon ) {
	   $pre_exact{$exon->id} = 0;
	   $pre_over{$exon->id} = 0;
       }
   }



   # ok. Ready to rock and roll
   # foreach predicted gene, see if there is a standard gene
   # of choice. Note the overlapping/missing gene flags.

 GENE:

   foreach my $gene ( @$predicted ) {
       foreach my $std ( @$standard ) {
	   if( $start{$gene->id} > $end{$std->id} || 
	       $end{$gene->id} < $start{$std->id} ) { 
	       next;
	   }
	   
	   TRANSCRIPT :
	   foreach my $trans ( $gene->each_Transcript ) {
	       my $not_exact = 0;
	       my $seen_base_pairs = 0;
	       foreach my $std_trans ( $std->each_Transcript ) {
		   # see if we have an exact exon match
		   my @p = $trans->each_Exon;
		   my @s = $std_trans->each_Exon;

		   if( scalar(@p) != scalar(@s) ) {
		       $not_exact = 1;
		   } else {
		       my $i;
		       for($i=0 ; $i <= $#p;$i++) {
			   if( $p[$i]->start  != $s[$i]->start ||
			       $p[$i]->end    != $s[$i]->end ||
			       $p[$i]->strand != $s[$i]->strand ) { 
			       
			       $not_exact = 1;
			       last;
			   } else {
			       $seen_base_pairs += $p[$i]->length;
			   }
		       }
		   }
		   if( $not_exact == 0 ) {
		       # we have an exact match!
		       $std_gene_exact{$std->id} = 1;
		       $pre_gene_exact{$gene->id} = 1;
		   }


		   # ok. Process exons now

		   foreach my $exon ( @p ) {
		       foreach my $stde ( @s ) {
			   if( $exon->overlaps($stde) ) {
			       my ($ex,$overlap,$basepair) = $self->exon_overlap($stde,$exon);
			       
			       $self->correct_base_pairs($self->correct_base_pairs+$basepair);
			       
			       if( $ex == 1 ) {
				   $std_exact{$stde->id} = 1;
				   $pre_exact{$exon->id} = 1;
			       }
			       
			       if( $overlap == 1 ) {
				   $std_over{$stde->id} = 1;
				   $pre_over{$exon->id} = 1;
				   $std_gene_over{$std->id} = 1;
				   $pre_gene_over{$gene->id} = 1;
			       }
			       
			   }
		       }
		       
		   }
	       } # end of each standard Transcript
	   } # end of each predicted Transcript
	   
       }
   }# end of foreach predicted gene
       
   #foreach $gene 
       

   # process the hashes for exact/overlap genes

   # process the hashes for exact/overlap exons
   
   # process overlapping gene hash to make merged gene object	       
   
}


#
# add/each missing_standard_gene
#

#
# add merged gene, each merged gene
#


#
# get sets on the counts
#



=head2 correct_base_pairs

 Title   : correct_base_pairs
 Usage   : $obj->correct_base_pairs($newval)
 Function: 
 Returns : value of correct_base_pairs
 Args    : newvalue (optional)


=cut

sub correct_base_pairs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'correct_base_pairs'} = $value;
    }
    return $obj->{'correct_base_pairs'};

}

=head2 missing_base_pairs

 Title   : missing_base_pairs
 Usage   : $obj->missing_base_pairs($newval)
 Function: 
 Returns : value of missing_base_pairs
 Args    : newvalue (optional)


=cut

sub missing_base_pairs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'missing_base_pairs'} = $value;
    }
    return $obj->{'missing_base_pairs'};

}

=head2 overpredicted_base_pairs

 Title   : overpredicted_base_pairs
 Usage   : $obj->overpredicted_base_pairs($newval)
 Function: 
 Returns : value of overpredicted_base_pairs
 Args    : newvalue (optional)


=cut

sub overpredicted_base_pairs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overpredicted_base_pairs'} = $value;
    }
    return $obj->{'overpredicted_base_pairs'};

}

=head2 exactly_correct_exons

 Title   : exactly_correct_exons
 Usage   : $obj->exactly_correct_exons($newval)
 Function: 
 Returns : value of exactly_correct_exons
 Args    : newvalue (optional)


=cut

sub exactly_correct_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_correct_exons'} = $value;
    }
    return $obj->{'exactly_correct_exons'};

}

=head2 exactly_missed_exons

 Title   : exactly_missed_exons
 Usage   : $obj->exactly_missed_exons($newval)
 Function: 
 Returns : value of exactly_missed_exons
 Args    : newvalue (optional)


=cut

sub exactly_missed_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_missed_exons'} = $value;
    }
    return $obj->{'exactly_missed_exons'};

}

=head2 exactly_overpredicted_exons

 Title   : exactly_overpredicted_exons
 Usage   : $obj->exactly_overpredicted_exons($newval)
 Function: 
 Returns : value of exactly_overpredicted_exons
 Args    : newvalue (optional)


=cut

sub exactly_overpredicted_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_overpredicted_exons'} = $value;
    }
    return $obj->{'exactly_overpredicted_exons'};

}

=head2 overlapping_correct_exons

 Title   : overlapping_correct_exons
 Usage   : $obj->overlapping_correct_exons($newval)
 Function: 
 Returns : value of overlapping_correct_exons
 Args    : newvalue (optional)


=cut

sub overlapping_correct_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_correct_exons'} = $value;
    }
    return $obj->{'overlapping_correct_exons'};

}

=head2 overlapping_missing_exons

 Title   : overlapping_missing_exons
 Usage   : $obj->overlapping_missing_exons($newval)
 Function: 
 Returns : value of overlapping_missing_exons
 Args    : newvalue (optional)


=cut

sub overlapping_missing_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_missing_exons'} = $value;
    }
    return $obj->{'overlapping_missing_exons'};

}

=head2 overlapping_overpredicted_exons

 Title   : overlapping_overpredicted_exons
 Usage   : $obj->overlapping_overpredicted_exons($newval)
 Function: 
 Returns : value of overlapping_overpredicted_exons
 Args    : newvalue (optional)


=cut

sub overlapping_overpredicted_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_overpredicted_exons'} = $value;
    }
    return $obj->{'overlapping_overpredicted_exons'};

}

=head2 exactly_correct_genes

 Title   : exactly_correct_genes
 Usage   : $obj->exactly_correct_genes($newval)
 Function: 
 Returns : value of exactly_correct_genes
 Args    : newvalue (optional)


=cut

sub exactly_correct_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_correct_genes'} = $value;
    }
    return $obj->{'exactly_correct_genes'};

}

=head2 exactly_missing_genes

 Title   : exactly_missing_genes
 Usage   : $obj->exactly_missing_genes($newval)
 Function: 
 Returns : value of exactly_missing_genes
 Args    : newvalue (optional)


=cut

sub exactly_missing_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_missing_genes'} = $value;
    }
    return $obj->{'exactly_missing_genes'};

}


=head2 exactly_overpredicted_genes

 Title   : exactly_overpredicted_genes
 Usage   : $obj->exactly_overpredicted_genes($newval)
 Function: 
 Returns : value of exactly_overpredicted_genes
 Args    : newvalue (optional)


=cut

sub exactly_overpredicted_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'exactly_overpredicted_genes'} = $value;
    }
    return $obj->{'exactly_overpredicted_genes'};

}

=head2 overlapping_correct_genes

 Title   : overlapping_correct_genes
 Usage   : $obj->overlapping_correct_genes($newval)
 Function: 
 Returns : value of overlapping_correct_genes
 Args    : newvalue (optional)


=cut

sub overlapping_correct_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_correct_genes'} = $value;
    }
    return $obj->{'overlapping_correct_genes'};

}

=head2 overlapping_missing_genes

 Title   : overlapping_missing_genes
 Usage   : $obj->overlapping_missing_genes($newval)
 Function: 
 Returns : value of overlapping_missing_genes
 Args    : newvalue (optional)


=cut

sub overlapping_missing_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_missing_genes'} = $value;
    }
    return $obj->{'overlapping_missing_genes'};

}

=head2 overlapping_overpredicted_genes

 Title   : overlapping_overpredicted_genes
 Usage   : $obj->overlapping_overpredicted_genes($newval)
 Function: 
 Returns : value of overlapping_overpredicted_genes
 Args    : newvalue (optional)


=cut

sub overlapping_overpredicted_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'overlapping_overpredicted_genes'} = $value;
    }
    return $obj->{'overlapping_overpredicted_genes'};

}

=head2 total_standard_basepairs

 Title   : total_standard_basepairs
 Usage   : $obj->total_standard_basepairs($newval)
 Function: 
 Returns : value of total_standard_basepairs
 Args    : newvalue (optional)


=cut

sub total_standard_basepairs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_standard_basepairs'} = $value;
    }
    return $obj->{'total_standard_basepairs'};

}

=head2 total_predicted_basepairs

 Title   : total_predicted_basepairs
 Usage   : $obj->total_predicted_basepairs($newval)
 Function: 
 Returns : value of total_predicted_basepairs
 Args    : newvalue (optional)


=cut

sub total_predicted_basepairs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_predicted_basepairs'} = $value;
    }
    return $obj->{'total_predicted_basepairs'};

}

=head2 total_standard_exons

 Title   : total_standard_exons
 Usage   : $obj->total_standard_exons($newval)
 Function: 
 Returns : value of total_standard_exons
 Args    : newvalue (optional)


=cut

sub total_standard_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_standard_exons'} = $value;
    }
    return $obj->{'total_standard_exons'};

}

=head2 total_predicted_exons

 Title   : total_predicted_exons
 Usage   : $obj->total_predicted_exons($newval)
 Function: 
 Returns : value of total_predicted_exons
 Args    : newvalue (optional)


=cut

sub total_predicted_exons{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_predicted_exons'} = $value;
    }
    return $obj->{'total_predicted_exons'};

}

=head2 total_standard_genes

 Title   : total_standard_genes
 Usage   : $obj->total_standard_genes($newval)
 Function: 
 Returns : value of total_standard_genes
 Args    : newvalue (optional)


=cut

sub total_standard_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_standard_genes'} = $value;
    }
    return $obj->{'total_standard_genes'};

}

=head2 total_predicted_genes

 Title   : total_predicted_genes
 Usage   : $obj->total_predicted_genes($newval)
 Function: 
 Returns : value of total_predicted_genes
 Args    : newvalue (optional)


=cut

sub total_predicted_genes{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_predicted_genes'} = $value;
    }
    return $obj->{'total_predicted_genes'};

}

