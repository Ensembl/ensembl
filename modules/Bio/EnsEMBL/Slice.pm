
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
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeqI;

@ISA = qw(Bio::EnsEMBL::Root Bio::PrimarySeqI);

# new() is written here 

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;

  my ($chr,$start,$end,$strand,$type,$adaptor, $dbID) = $self->_rearrange([qw(CHR_NAME CHR_START CHR_END STRAND ASSEMBLY_TYPE ADAPTOR DBID)],@args);

  if( !defined $chr || !defined $start || !defined $end || !defined $type ) {
    print "Chr: " . $chr . "\t" . "Start: " . $start . "\t" . "End: " . $end . "\t" . "Type: " . $type . "\n";
    $self->throw("Do not have all the parameters for slice");
  }

  $self->chr_name($chr);
  $self->chr_start($start);
  $self->chr_end($end);
  $self->id("$chr.$start-$end");

  #set strand to a default of 1 if it is not set
  if ( undef $strand) {
    $self->strand($strand);
  }
  else {
    $self->strand('1');
  }

  $self->assembly_type($type);
  $self->adaptor($adaptor);
  $self->dbID( $dbID );
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



=head2 seq

  Args      : none
  Function  : returns the entire sequence string for this Slice
              needs the adaptor to be set.
  Returntype: txt
  Exceptions: none
  Caller    : general

=cut

sub seq {
  my $self = shift;
  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  my $seq = $seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1, -1, 1 );

  return $seq;
}


=head2 subseq

  Arg  1    : int $startBasePair
              relative to start of slice, which is 1.
  Arg  2    : int $endBasePair
              relative to start of slice.
  Arg  3    : int $strand
  Function  : returns string of dna sequence
  Returntype: txt
  Exceptions: end should be at least as big as start
              strand must be set
  Caller    : general

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end < $start ) {
    $self->throw("End coord is less then start coord to call on Slice subseq.");
  }

  if ( !defined $strand || ( $strand != -1 && $strand != 1 )) {
#    $self->throw("Incorrect strand information set to call on Slice subseq.");
    $strand = 1;
  }

  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  my $seq = $seqAdaptor->fetch_by_Slice_start_end_strand( $self, $start, $end, $strand );

  return $seq;
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

=head2 repeat_adaptor

 Title   : repeat_adaptor
 Usage   : $obj->repeat_adaptor
 Function: For getting hold of the repeat feature adaptor within a slice obj
 Example :
 Returns : RepeatFeatureAdaptor
 Args    : none


=cut


sub repeat_adaptor{
    my ($self) = @_;

    my $db = $self->adaptor->{'db'};

    my $rep_adaptor = $db->get_RepeatFeatureAdaptor;

    return $rep_adaptor;
}


=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   : $obj->get_all_PredictionFeatures
 Function: Use to derive a list of prediction features specific to the analysis type specified by the logic name.
 Example : my @pred_rm_feat = $obj->get_all_PredictionFeatures('RepeatMasker');
 Returns : a list of Bio::EnsEMBL::PredictionTranscript objects
 Args    : a logic name - the name of the analysis that created or returned the prediction feature.


=cut

sub get_all_PredictionFeatures{
   my ($self,@args) = @_;
   
   my @pred_feat = $self->adaptor->fetch_all_prediction_transcripts($self);
   
   return @pred_feat;
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

   my $gene_adaptor = $self->adaptor->db->get_GeneAdaptor();
   my @genes = $gene_adaptor->fetch_by_Slice($self);

   return @genes;
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

    if( $self->{'_virtual_primary_seq'} ){
	return $self->{'_virtual_primary_seq'};
    }

   my $seq = $self->seq();
    my $slice_seq = Bio::PrimarySeq->new( 
					 -id    =>$self->id,
					 -'seq' =>$seq
					);

    $self->{'_virtual_primary_seq'} = $slice_seq;
    return $slice_seq;
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

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Example : 
 Returns : value of strand
 Args    : newvalue (optional)


=cut

sub strand{
   my ($self,$value) = @_;

   if( defined $value) {
      $self->{'strand'} = $value;
    }
    return $self->{'strand'};

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


=head2 dbID

  Arg [1]   : int databaseInternalId
              A slice might exist in the database and will than have this
              internal id.
  Function  : attribute function
  Returntype: int
  Exceptions: none
  Caller    : DBSQL::SliceAdaptor

=cut

sub dbID {
   my ( $self, $value ) = @_;
   if( defined $value ) {
     $self->{'dbID'} = $value;
   }
   return $self->{'dbID'};
}

sub id {
   my ( $self, $value ) = @_;
   if( defined $value ) {
     $self->{'id'} = $value;
   }
   return $self->{'id'};
}



1;

