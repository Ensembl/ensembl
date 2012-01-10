
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

Bio::EnsEMBL::SlicingEvent - Object representing an alternative splicing event

=head1 SYNOPSIS

  my $ase =
    Bio::EnsEMBL::SplicingEvent->new( -START   => 123,
                                      -END     => 1045,
                                      -STRAND  => 1,
                                      -GENE_ID => $gene->dbID,
                                      -SLICE   => $slice );

  # set some additional attributes
  $ase->name('ENSG00000000003-CNE-3');
  $ase->type('CNE');

=head1 DESCRIPTION

A representation of an Alternative Splicing Event within the Ensembl system.

=head1 METHODS

=cut

package Bio::EnsEMBL::SplicingEvent;

use strict;

use POSIX;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);

## Add to gene get_all_splicing_events

sub gene_id {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'gene_id'} = $value;
  }

  return $self->{'gene_id'};
}

sub name {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'name'} = $value;
  }

  return $self->{'name'};
}

sub type {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'type'} = $value;
  }

  return $self->{'type'};
}

sub add_Feature {
  my ( $self, $feature ) = @_;

  if (    !ref($feature)
       || !$feature->isa("Bio::EnsEMBL::SplicingEventFeature") )
  {
    throw("$feature is not a Bio::EnsEMBL::SplicingEventFeature!");
  }

  $self->{'_feature_array'} ||= [];

  push( @{ $self->{'_feature_array'} }, $feature );
}

sub get_all_Features {
  my ($self) = @_;

  if ( !exists( $self->{'_feature_array'} ) ) {
    if ( defined( $self->adaptor() ) ) {
      my $fta =
        $self->adaptor()->db()->get_SplicingEventFeatureAdaptor();
      my $features = $fta->fetch_all_by_SplicingEvent($self);
      $self->{'_feature_array'} = $features;
    }
  }

  return $self->{'_feature_array'};
}

sub add_Pair {
  my ( $self, $feature ) = @_;

  if (    !ref($feature)
       || !$feature->isa("Bio::EnsEMBL::SplicingEventPair") )
  {
    throw("$feature is not a Bio::EnsEMBL::SplicingEventPair!");
  }

  $self->{'_pair_array'} ||= [];

  push( @{ $self->{'_pair_array'} }, $feature );
}

sub get_all_Pairs {
  my ($self) = @_;

  if ( !exists( $self->{'_pair_array'} ) ) {
    if ( defined( $self->adaptor() ) ) {
      my $pa =
        $self->adaptor()->db()->get_SplicingTranscriptPairAdaptor();
      my $pairs = $pa->fetch_all_by_SplicingEvent($self);
      $self->{'_pair_array'} = $pairs;
    }
  }

  return $self->{'_pair_array'};
}

1;
