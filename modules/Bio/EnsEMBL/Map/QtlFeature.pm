# EnsEMBL module for QtlFeature
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::QtlFeature

=head1 SYNOPSIS


=head1 AUTHOR

Arne Stabenau stabenau@ebi.ac.uk

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 DESCRIPTION

Represents a QtlFeature in the EnsEMBL database. QtlFeatures are generally very
long and its not clear wether a representation in Contig coordinates
actually makes sense. In the database they will have chromosomal coordinates.


 
=cut

package Bio::EnsEMBL::Map::QtlFeature;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::SeqFeature);



=head2 new

  Arg [1]    : Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor $adaptor
  Example    : none
  Description: Create a QtlFeature
  Returntype : Bio::EnsEMBL::Map::QtlFeature
  Exceptions : none
  Caller     : general, DBSQL::QtlFeatureAdaptor

=cut

sub new {
  my ( $class, $adaptor, $contig, $start, $end, $qtl, $analysis ) = @_;

  $class = ref( $class ) ||$class;
  my $self = bless( {
		     'adaptor' => $adaptor,
		     '_gsf_seq' => $contig,
		     '_gsf_start' => $start,
		     '_gsf_end' => $end,
		     'qtl' => $qtl,
		     'analysis' => $analysis
		    }, $class );
  $self->is_splittable( 1 );
  return $self;
}



=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor $adaptor
  Example    : none
  Description: Getter/Setter attribute adaptor
  Returntype : Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor
  Exceptions : none
  Caller     : DBSQL::QtlFeatureAdaptor

=cut

sub adaptor {
  my $self = shift;

  if(@_) {
    $self->{'adaptor'} = shift;
  }

  return $self->{'adaptor'};
}


=head2 qtl

  Arg [1]    : Bio::EnsEMBL::Map::Qtl $qtl
               the qtl object for this feature
  Example    : none
  Description: return the Qtl object associated with this location
  Returntype : Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general

=cut

sub qtl {
  my $self = shift;

  if(@_) {
    $self->{'qtl'} = shift;
  }

  return $self->{'qtl'};
}

=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysis $analysis
  Example    : none
  Description: getter/setter for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut

sub analysis {
  my $self = shift;

  if(@_) {
    $self->{'analysis'} = shift;
  }

  return $self->{'analysis'};
}


=head2 strand

  Arg [1]    : none
	Example    : $strand = $qtl_feat->strand();
  Description: Overrides the SeqFeature strand method to always return a 
               value of 0 for qtl features (they are unstranded features)
  Returntype : int (always 0)
  Exceptions : none
  Caller     : general

=cut

sub strand {
	my $self = shift;
  
  return 0;
}


1;


