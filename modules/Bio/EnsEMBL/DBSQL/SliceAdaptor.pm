
#
# Ensembl module for Bio::EnsEMBL::Assembly::SliceFactory
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SliceAdaptor - Adaptors for slices

=head1 SYNOPSIS
  



=head1 DESCRIPTION

Factory for getting out slices of assemblies. WebSlice is the highly
accelerated version for the web site.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::SliceAdaptor;
use vars qw(@ISA);
use strict;


# Object preamble - inherits from Bio::EnsEMBL::Root
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;


@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


# new is inherieted from BaseAdaptor

=head2 new_slice

 Title   : new_slice
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_slice{
    my ($self,$chr,$start,$end,$type) = @_;


    my $slice = Bio::EnsEMBL::Slice->new( -chr_name  => $chr,
					  -chr_start => $start,
					  -chr_end   => $end,
					  -assembly_type      => $type,
					-adaptor => $self);

    return $slice;
}


=head2 new_web_slice

 Title   : new_web_slice
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_web_slice{
    my ($self,$chr,$start,$end,$type) = @_;
    
    die "Not implemented new slice yet";
    
}


sub fetch_all_repeat_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all repeat features if con't have a slice to fetch them for\n");
  }

  my @repeats = $self->db->get_RepeatFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  return @repeats;

}

sub fetch_all_simple_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }

  my @simple = $self->db->get_SimpleFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  return @simple;

}

sub fetch_all_prediction_transcripts{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }

  my @prediction = $self->db->get_PredictionTranscriptAdaptor->fetch_by_Slice($slice, $logic_name);

  return @prediction;

}





sub fetch_all_similarity_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice($slice, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}

sub fetch_all_similarity_features_above_score{
  my($self, $slice, $score, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  if(!$score){
    $self->throw("need score even if it 0\n");
  }
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice_and_score($slice, $score, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice_and_score($slice, $score, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}


sub fetch_all_similarity_features_above_pid{
  my($self, $slice, $pid, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  if(!$pid){
    $self->throw("need percent_id even if it 0\n");
  }
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice_and_pid($slice, $pid, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice_and_pid($slice, $pid, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}
