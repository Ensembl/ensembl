
#
# Ensembl module for Bio::EnsEMBL::Assembly::Slice
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Slice - Arbitary Slice of a genome

=head1 SYNOPSIS


   foreach $gene ( $slice->get_all_Genes ) {
      # do something with a gene
   }
       

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Slice;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

# new() is written here 

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;

  my ($chr,$start,$end,$type,$adaptor) = $self->_rearrange([qw(CHR_NAME CHR_START CHR_END ASSEMBLY_TYPE ADAPTOR)],@args);

  if( !defined $chr || !defined $start || !defined $end || !defined $type ) {
    $self->throw("Do not have all the parameters for slice");
  }

  $self->chr_name($chr);
  $self->chr_start($start);
  $self->chr_end($end);
  $self->assembly_type($type);
  $self->adaptor($adaptor);

# set stuff in self from @args
  return $self;
}


=head2 First pass implementation

The first pass implementation tries to mimic precisely the important
parts of the old VirtualContig system for the pipeline. These
are the methods to implement

=cut

=head2 get_all_SimilarityFeatures_above_score

 Title   : get_all_SimilarityFeatures_above_score
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_score{
   my ($self,$score) = @_;

   my @out;

   if( !defined $score ) {
     $self->throw("No defined score.");
   }

   push(@out,$self->adaptor->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice_and_score($self,$score));

   return @out;
}



=head2 get_all_SimilarityFeatures_above_pid

 Title   : get_all_SimilarityFeatures_above_pid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_pid{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");
}


=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures{
   my ($self,@args) = @_;


   my @repeats = $self->repeat_adaptor->fetch_by_Slice($self);

   foreach my $repeat ( @repeats ) {
       $repeat->transform_location($self->start);
   }

   return @repeats;
}


=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PredictionFeatures{
   my ($self,@args) = @_;


   $self->throw("Ewan has not implemented this function! Complain!!!!");
}


=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ExternalFeatures{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}


=head2 get_Genes_by_type

 Title   : get_Genes_by_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_type{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}

=head2 invert

 Title   : invert
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub invert{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}


=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_seq{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}


=head2 convert_Gene_to_raw_contig

 Title   : convert_Gene_to_raw_contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub convert_Gene_to_raw_contig{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");
}


=head2 Backward Compatibility functions

=cut

=head2 get_all_Genes_exononly

 Title   : get_all_Genes_exononly
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes_exononly{
   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l get_all_Genes_exononly has been deprecated. get_all_Genes called");

   return $self->get_all_Genes();
}



=head2 Internal functions

=cut

=head2 chr_name

 Title   : chr_name
 Usage   : $obj->chr_name($newval)
 Function: 
 Example : 
 Returns : value of chr_name
 Args    : newvalue (optional)


=cut

sub chr_name{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'chr_name'} = $value;
    }
    return $self->{'chr_name'};

}

=head2 chr_start

 Title   : chr_start
 Usage   : $obj->chr_start($newval)
 Function: 
 Example : 
 Returns : value of chr_start
 Args    : newvalue (optional)


=cut

sub chr_start{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'chr_start'} = $value;
    }
    return $self->{'chr_start'};

}

=head2 chr_end

 Title   : chr_end
 Usage   : $obj->chr_end($newval)
 Function: 
 Example : 
 Returns : value of chr_end
 Args    : newvalue (optional)


=cut

sub chr_end{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'chr_end'} = $value;
    }
    return $self->{'chr_end'};

}

=head2 assembly_type

 Title   : assembly_type
 Usage   : $obj->assembly_type($newval)
 Function: 
 Example : 
 Returns : value of assembly_type
 Args    : newvalue (optional)


=cut

sub assembly_type{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'assembly_type'} = $value;
    }
    return $self->{'assembly_type'};

}


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: 
 Example : 
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}



1;

