# EnsEMBL module for Qtl
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::Qtl

=head1 SYNOPSIS


=head1 AUTHOR

Arne Stabenau stabenau@ebi.ac.uk

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=head1 DESCRIPTION

Represents a Qtl in the EnsEMBL database. A quantitative trait locus is 
defined by three markers, two flanking and one peak (optional) marker.
Its a region (or more often a group of regions) which is likely to
affect the phenotype (trait) described in this Qtl.
 
=cut

package Bio::EnsEMBL::Map::Qtl;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Arg [1]    : int $dbID
  Arg [2]    : Bio::EnsEMBL::Map::DBSQL::QtlAdaptor $adaptor
  Arg [3]    : Bio::EnsEMBL::Map::Marker $flank_marker_1
  Arg [4]    : Bio::EnsEMBL::Map::Marker $peak_marker
  Arg [5]    : Bio::EnsEMBL::Map::Marker $flank_marker_2
  Arg [6]    : string $trait
  Arg [7]    : float $lod_score
  Arg [8]    : string $source_database
  Arg [9]    : string $source_primary_id
  Example    : none
  Description: Creates a new Qtl object. Usually done by Adaptor
  Returntype : Bio::EnsEMBL::Map::Qtl
  Exceptions : none
  Caller     : general, DBSQL::QtlAdaptor, DBSQL::QtlFeatureAdaptor

=cut

sub new {
  my ( $class, $dbID, $adaptor, $flank_marker_1, $peak_marker, 
       $flank_marker_2, $trait, $lod_score,
       $source_database, $source_primary_id ) = @_;

  $class = ref( $class ) ||$class;
  my $self = bless( {
		     'dbID' => $dbID,
		     'adaptor' => $adaptor,
		     'flank_marker_1' => $flank_marker_1,
		     'flank_marker_2' => $flank_marker_2,
		     'peak_marker' => $peak_marker,
		     'trait' => $trait,
		     'lod_score' => $lod_score,
		     'source_database' => $source_database,
		     'source_primary_id' => $source_primary_id
		    }, $class );
  return $self;
}



=head2 dbID

  Arg  [1]   : int $dbID
  Example    : none
  Description: get/set/clear attribute dbID
  Returntype : int
  Exceptions : none
  Caller     : DBSQL::QtlAdaptor

=cut

sub dbID {
  my $self = shift;
  
  if(@_) {
    $self->{'dbID'} = shift;
  }

  return $self->{'dbID'};
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::Map::DBSQL::QtlAdaptor $adaptor
  Example    : none
  Description: Getter/Setter attribute adaptor
  Returntype : Bio::EnsEMBL::Map::DBSQL::QtlAdaptor
  Exceptions : none
  Caller     : DBSQL::QtlAdaptor

=cut

sub adaptor {
  my $self = shift;

  if(@_) {
    $self->{'adaptor'} = shift;
  }

  return $self->{'adaptor'};
}


=head2 source_database

  Arg [1]    : string $source_database
               Name for the database that provided the QTL. It should
               give possibility to URL link to the QTL.
  Example    : none
  Description: Getter/Setter for the source_database attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub source_database {
  my $self = shift;

  if(@_) {
    $self->{'source_database'} = shift;
  }

  return $self->{'source_database'};
}


=head2 source_primary_id

  Arg [1]    : string $source_primary_id
               primary_id for this qtl in the source database. Should enable
               URL construction to get additional information about Qtl
  Example    : none
  Description: Getter/Setter for attribute source_primary_id
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub source_primary_id {
  my $self = shift;

  if(@_) {
    $self->{'source_primary_id'} = shift;
  }

  return $self->{'source_primary_id'};

}


=head2 trait

  Arg [1]    : string $trait
               Phenotype of this Qtl
  Example    : none
  Description: Getter/Setter for the trait attribute
  Returntype : string
  Exceptions : none
  Caller     : general

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

