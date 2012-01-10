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

Bio::EnsEMBL::Map::Qtl

=head1 SYNOPSIS

=head1 DESCRIPTION

Represents a Qtl in the EnsEMBL database. A quantitative trait locus is
defined by three markers, two flanking and one peak (optional) marker.
Its a region (or more often a group of regions) which is likely to
affect the phenotype (trait) described in this Qtl.

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::Qtl;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

@ISA = qw(Bio::EnsEMBL::Storable);



=head2 new

  Arg [1]    : int $dbID
  Arg [2]    : Bio::EnsEMBL::Map::DBSQL::QtlAdaptor $adaptor
  Arg [3]    : Bio::EnsEMBL::Map::Marker $flank_marker_1
  Arg [4]    : Bio::EnsEMBL::Map::Marker $peak_marker
  Arg [5]    : Bio::EnsEMBL::Map::Marker $flank_marker_2
  Arg [6]    : string $trait
  Arg [7]    : float $lod_score
  Arg [8]    : hashref $synonyms
               A hashref with source keys and identifier values
  Example    : none
  Description: Creates a new Qtl object. Usually done by Adaptor
  Returntype : Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general, DBSQL::QtlAdaptor, DBSQL::QtlFeatureAdaptor
  Status     : stable

=cut

sub new {
  my ( $class, $dbID, $adaptor, $flank_marker_1, $peak_marker,
       $flank_marker_2, $trait, $lod_score,
       $synonyms ) = @_;

  $class = ref( $class ) ||$class;
  my $self = bless( {
		     'dbID'           => $dbID,
		     'flank_marker_1' => $flank_marker_1,
		     'flank_marker_2' => $flank_marker_2,
		     'peak_marker'    => $peak_marker,
		     'trait'          => $trait,
		     'lod_score'      => $lod_score,
		     'synonyms'       => $synonyms
		    }, $class );
  $self->adaptor($adaptor);
  return $self;
}


=head2 add_synonym

  Arg [1]    : string $source
               The source of the synonym
  Arg [2]    : string $identifier
               The identifier from this source
  Example    : $qtl->add_synonym('rat genome database', '65516');
  Description: Adds a synonym to this qtl
  Returntype : none
  Exceptions : thrown if arguments are not provided
  Caller     : general
  Status     : stable

=cut

sub add_synonym {
  my $self = shift;
  my $source = shift;
  my $identifier = shift;

  unless($source && $identifier) {
    throw('Source and identifier arguments are required');
  }

  $self->{'synonyms'}->{$source} = $identifier;
}


=head2 get_synonyms

  Arg [1]    : none
  Example    :
     foreach my $source ($keys %{$qtl->get_synonyms}) {
       print $source . ':'. $qtl->get_synonyms->{$source};
     }
  Description: Returns a hashref of synonyms keyed on their source name 
  Returntype : hashref of synonyms keyed on their source name
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub get_synonyms {
  my $self = shift;

  return $self->{'synonyms'} || {};
}



=head2 trait

  Arg [1]    : string $trait
               Phenotype of this Qtl
  Example    : none
  Description: Getter/Setter for the trait attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub trait {
  my $self = shift;

  if(@_) {
    $self->{'trait'} = shift;
  }

  return $self->{'trait'};
}


=head2 lod_score

  Arg [1]    : float $lod_score
               A score for the Qtl
  Example    : none
  Description: Getter/Setter for attribute lod_score
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub lod_score {
  my $self = shift;

  if(@_) {
    $self->{'lod_score'} = shift;
  }

  return $self->{'lod_score'};
}


=head2 peak_marker

  Arg [1]    : Bio::EnsEMBL::Map::Marker $peak_marker
               an optional Marker which has the peak probablitity
               for this traits occurence
  Example    : none
  Description: Getter/Setter for attribute peak_marker
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub peak_marker {
  my $self = shift;

  if(@_) {
    $self->{'peak_marker'} = shift;
  }

  return $self->{'peak_marker'};
}


=head2 flank_marker_1

  Arg [1]    : Bio::EnsEMBL::Map::Marker $flank_marker_1
               One flanking marker of the interest region, the two flanking
               markers define the region
  Example    : none
  Description: Getter/Setter attribute flanking_marker_1
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub flank_marker_1 {
  my $self = shift;

  if(@_) {
    $self->{'flank_marker_1'} = shift;
  }

  return $self->{'flank_marker_1'};
}



=head2 flank_marker_2

  Arg [1]    : Bio::EnsEMBL::Map::Marker $flank_marker_2
               One flanking marker of the interest region, the two flanking
               markers define the region
  Example    : none
  Description: Getter/Setter attribute flanking_marker_2
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : none
  Caller     : general
  Status     : stable

=cut


sub flank_marker_2 {
  my $self = shift;

  if(@_) {
    $self->{'flank_marker_2'} = shift;
  }

  return $self->{'flank_marker_2'};
}



=head2 get_QtlFeatures

  Args       : none
  Example    : none
  Description: return the qtl feature which is associated with this
               Qtl. It comes in chromosomal slice coordinates. There can 
               only be one.
  Returntype : Bio::EnsEMBL::Map::QtlFeature
  Exceptions : only works with adaptored Qtls
  Caller     : general
  Status     : stable

=cut

sub get_QtlFeature {
  my $self = shift;

  my $adaptor = $self->adaptor();
  return undef unless $adaptor;
  my $result = $adaptor->db()->get_QtlFeatureAdaptor()->
    fetch_all_by_Qtl( $self );

  if( @$result ) {
    return $result->[0];
  } else {
    return;
  }
}





=head2 source_database

This method is deprecated.  Use get_synonyms or add_synonym instead.

=cut

sub source_database {
  my $self = shift;

  deprecate('Use get_synonyms or add_synonym instead');

  my $syns = $self->get_synonyms;
  my ($source) = keys %$syns;

  return $source || '';
}


=head2 source_primary_id

This method is deprecated. Use get_synonyms or add_synonym instead.

=cut

sub source_primary_id {
  my $self = shift;

  deprecate('Use get_synonyms or add_synonym instead');

  my $syns = $self->get_synonyms;
  my ($source) = keys %$syns;

  if($source) {
    return $syns->{$source};
  }

  return '';
}


1;
