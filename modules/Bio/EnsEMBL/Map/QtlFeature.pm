=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Map::QtlFeature

=head1 SYNOPSIS

=head1 DESCRIPTION

Represents a QtlFeature in the EnsEMBL database. QtlFeatures are
generally very long and its not clear wether a representation in Contig
coordinates actually makes sense. In the database they will have
chromosomal coordinates.

=head1 METHODS

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
		     'slice'    => $slice,
		     'start'    => $start,
		     'end'      => $end,
		     'qtl'      => $qtl,
		     'analysis' => $analysis,
         'strand'   => 0
		    }, $class );

  $self->adaptor($adaptor);
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


