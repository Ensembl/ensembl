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

use Bio::EnsEMBL::Feature;

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [1]    : Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor $adaptor
  Example    : none
  Description: Create a QtlFeature
  Returntype : Bio::EnsEMBL::Map::QtlFeature
  Exceptions : none
  Caller     : general, DBSQL::QtlFeatureAdaptor
  Status     : Stable

=cut

sub new {
  my ( $class, $adaptor, $slice, $start, $end, $qtl, $analysis ) = @_;

  $class = ref( $class ) ||$class;
  my $self = bless( {
		     'adaptor'  => $adaptor,
		     'slice'    => $slice,
		     'start'    => $start,
		     'end'      => $end,
		     'qtl'      => $qtl,
		     'analysis' => $analysis,
         'strand'   => 0
		    }, $class );

  return $self;
}


=head2 qtl

  Arg [1]    : Bio::EnsEMBL::Map::Qtl $qtl
               the qtl object for this feature
  Example    : none
  Description: return the Qtl object associated with this location
  Returntype : Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub qtl {
  my $self = shift;

  if(@_) {
    $self->{'qtl'} = shift;
  }

  return $self->{'qtl'};
}



=head2 strand

  Arg [1]    : none
	Example    : $strand = $qtl_feat->strand();
  Description: Overrides the Feature strand method to always return a
               value of 0 for qtl features (they are unstranded features)
  Returntype : int (always 0)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub strand {
	my $self = shift;
  return 0;
}



=head2 move

  Arg [1]    : $start - The new end of this qtl feature
  Arg [2]    : $end - The new start of this qtl feature
  Arg [3]    : $strand - ignored always set to 0
  Example    : $qtl_feat->move(1, 10_000);
  Description: Overrides superclass move() method to ensure strand is always 0.
               See Bio::EnsEMBL::Feature::move
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub move {
  my ($self, $start, $end, $strand) = @_;

  #maintain a strandedness of 0
  return $self->SUPER::move($start,$end,0);
}



1;


