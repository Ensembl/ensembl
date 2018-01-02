=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Map::Marker

=head1 SYNOPSIS

=head1 DESCRIPTION

Represents a marker in the EnsEMBL database.  The marker object
is unpositioned on the genome.  Markers which are positioned are
represented by the MarkerFeature object.

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::Marker;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw);

@ISA = qw(Bio::EnsEMBL::Storable);



=head2 new

  Arg [1]    : (optional) int $dbID
  Arg [2]    : (optional) Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor $adaptor
  Arg [3]    : (optional) string $left_primer
  Arg [4]    : (optional) string $right_primer
  Arg [5]    : (optional) int $primer_distance
  Arg [6]    : (optional) int $priority
  Arg [7]    : (optional) string $type
  Arg [8]    : (optional) Bio::EnsEMBL::Map::MarkerSynonym $display_synonym
  Arg [9]    : (optional) listref of Bio::EnsEMBL::Map::MarkerSynonyms $syns
  Arg [10]   : (optional) listref of Bio::EnsEMBL::Map::MapLocations $locs
  Example    : $marker = Bio::EnsEMBL::Map::MarkerSynonym->new
                            (123, $adaptor,
			     $left_primer, $right_primer, 400,
			     80, $ms1, [$ms1, $ms2], [$mloc1, $mloc2]);
  Description: Creates a new Marker
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub new {
  my ($caller, $dbID, $adaptor, $left_primer, $right_primer,
      $min_primer_dist, $max_primer_dist, $priority, $type, $display_synonym,
      $syns, $mlocs) = @_;

  my $class = ref($caller) || $caller;

  my $self = bless( {'dbID'            => $dbID,
                     'left_primer'     => $left_primer,
                     'right_primer'    => $right_primer,
                     'min_primer_dist' => $min_primer_dist,
                     'max_primer_dist' => $max_primer_dist,
                     'priority'        => $priority,
                     'type'            => $type,
                     'display_marker_synonym' => $display_synonym
                    }, $class);

  $self->adaptor($adaptor);

  #only load the marker synononyms if they were supplied, otherwise they 
  # will be lazy-loaded
  if($syns && @$syns) {
    $self->{'marker_synonyms'} = $syns;
  }

  #only load the map_locations if they were supplied, otherwise they will
  # be lazy-loaded
  if($mlocs) {
    foreach my $ml (@$mlocs) {
      $self->add_MapLocation($ml);
    }
  }

  return $self;
}


=head2 left_primer

  Arg [1]    : (optional) string $left_primer
  Example    : $left_primer = $marker->left_primer;
  Description: Getter/Setter for the left primer sequence of this marker
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub left_primer {
  my $self = shift;
  
  if(@_) {
    $self->{'left_primer'} = shift;
  }

  return $self->{'left_primer'};
}



=head2 right_primer

  Arg [1]    : (optional) string $right_primer
  Example    : $right_primer = $marker->right_primer;
  Description: Getter/Setter for the right primer sequence of this marker
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub right_primer {
  my $self = shift;
  
  if(@_) {
    $self->{'right_primer'} = shift;
  }

  return $self->{'right_primer'};
}



=head2 min_primer_dist

  Arg [1]    : (optional) string $min
  Example    : $dist = $marker->min_primer_dist;
  Description: Getter/Setter for the minimum seperation distance between the 
               left and right primers of this marker
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub min_primer_dist {
  my $self = shift;

  if(@_) {
    $self->{'min_primer_dist'} = shift;
  }

  return $self->{'min_primer_dist'};
}


=head2 max_primer_dist

  Arg [1]    : (optional) string $max
  Example    : $dist = $marker->max_primer_dist;
  Description: Getter/Setter for the maximum seperation distance between the 
               left and right primers of this marker
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub max_primer_dist {
  my $self = shift;

  if(@_) {
    $self->{'max_primer_dist'} = shift;
  }

  return $self->{'max_primer_dist'};
}



=head2 priority

  Arg [1]    : (optional) int $priority
  Example    : $priority = $marker->priority;
  Description: Getter/Setter for priority of this marker which can be used to 
               determine which markers are displayed.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub priority {
  my $self = shift;

  if(@_) {
    $self->{'priority'} = shift;
  }
  
  return $self->{'priority'};
}




=head2 type

  Arg [1]    : (optional) string $type
  Example    : $type = $marker->type;
  Description: Getter/Setter for type of this marker. Rat markers are typed
               as 'est' or 'microsatellite'.  Other markers may not have 
               defined types.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub type {
  my $self = shift;

  if(@_) {
    $self->{'type'} = shift;
  }
  
  return $self->{'type'};
}




=head2 get_all_MarkerSynonyms

  Arg [1]    : none
  Example    : @synonyms = @{$marker->get_all_MarkerSynonyms};
  Description: Retrieves a list of marker synonyms associated with this
               marker.  If this marker is connected to the datbase (i.e. it
               has an adaptor and 
  Returntype : listref of Bio::EnsEMBL::Map::MarkerSynonyms
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub get_all_MarkerSynonyms {
  my $self = shift;
  
  #lazy-load the marker synonyms if they haven't been retrieved
  if(!exists $self->{'marker_synonyms'} && 
     $self->adaptor && $self->{'dbID'}) {
    $self->adaptor->fetch_attributes($self);
  }
  
  return $self->{'marker_synonyms'} || [];
}



=head2 add_MarkerSynonyms

  Arg [1]    : Bio::EnsEMBL::MarkerSynonym $ms
  Example    : $marker->add_MarkerSynonym($ms);
  Description: Associates a new synonym with this marker.  Adding marker 
               synonyms to a marker which has not yet retrieved its
               synonyms from the database will prevent the loading of these
               from the database at request time (unless flush_MarkerSynonyms
               is called first).
  Returntype : none
  Exceptions : thrown if incorrect argument is passed
  Caller     : general
  Status     : stable

=cut

sub add_MarkerSynonyms {
  my ($self, @ms) = @_;

  #create the array if it does not exist it
  $self->{'marker_synonyms'} ||= [];

  push(@{$self->{'marker_synonyms'}}, @ms);
}



=head2 flush_MarkerSynonyms

  Arg [1]    : none
  Example    : $marker->flush_MarkerSynonyms;
  Description: clears all of the marker sysnonyms which have been added to 
               this marker.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub flush_MarkerSynonyms {
  my $self = shift;

  delete $self->{'marker_synonyms'};
}



=head2 display_MarkerSynonym

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::MarkerSynonym $ms
  Example    : none
  Description: Getter/Setter for the 'display' synonym of this marker
  Returntype : Bio::EnsEMBL::Map::MarkerSynonym
  Exceptions : thrown if the argument is invalid
  Caller     : general
  Status     : stable

=cut

sub display_MarkerSynonym {
  my $self = shift;
  
  if(@_) {
    my $ms = shift;
    if($ms && !(ref $ms && $ms->isa('Bio::EnsEMBL::Map::MarkerSynonym'))) {
     throw("ms arg must be Bio::EnsEMBL::Map::MarkerSynonym not [$ms]");
    }    
    $self->{'display_marker_synonym'} = $ms;
  } 


  return $self->{'display_marker_synonym'};
}



=head2 get_all_MarkerFeatures

  Arg [1]    : none
  Example    : @marker_features = @{$marker->get_all_MarkerFeatures};
  Description: Retrieves the marker features which are associated with this
               marker.  I.e. locations where this marker has been mapped to
               the genome via e-PCR
  Returntype : listref of Bio::EnsEMBL::Map::MarkerFeatures
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub get_all_MarkerFeatures {
  my $self = shift;

  my $marker_feature_adaptor = $self->adaptor->db->get_MarkerFeatureAdaptor;

  #these results are not cached to avoid a circular reference loop
  return $marker_feature_adaptor->fetch_all_by_Marker($self);
}



=head2 get_all_MapLocations

  Arg [1]    : none
  Example    : @map_locations = @{$marker->get_all_MapLocations};
  Description: Retrieves all map locations which are associated with this 
               marker. 
  Returntype : listref of Bio::EnsEMBL::Map::MapLocations
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub get_all_MapLocations {
  my $self = shift;

  #lazy-load the map locations if they have not been fetched yet
  if(!exists $self->{'map_locations'} && 
     $self->adaptor && $self->{'dbID'}) {
    $self->adaptor->fetch_attributes($self);
  }

  my @out = values %{$self->{'map_locations'}};

  return \@out;
}



=head2 get_MapLocation

  Arg [1]    : string $map_name 
  Example    : $map_location = $marker->get_MapLocation('genethon');
  Description: Retrieves the location of this marker in a specified map.
               If this marker is not defined in the specified map then 
               undef is returned.
  Returntype : Bio::EnsEMBL::Map::MapLocation
  Exceptions : thrown if the map_name arg is not provided
  Caller     : general
  Status     : stable

=cut

sub get_MapLocation {
  my $self = shift;
  my $map_name = shift;

  #lazy-load the map locations if they have not been fetched yet
  if(!exists $self->{'map_locations'} && 
     $self->adaptor && $self->{'dbID'}) {
    $self->adaptor->fetch_attributes($self);
  }

  unless($map_name) {
    throw('map_name argument is required, or use get_all_MapLocations');
  }

  return $self->{'map_locations'}->{$map_name}; 
}



=head2 add_MapLocations

  Arg [1..n] : @mlocs list of Bio::EnsEMBL::MapLocations
  Example    : $marker->add_MapLocations(@mlocs);
  Description: Associates 1 or more map locations with this marker
               using this function to manually load map locations will prevent
               lazy-loading of locations from the database. 
  Returntype : listref of Bio::EnsEMBL::MapLocations
  Exceptions : throws if map location has no name
  Caller     : general
  Status     : stable

=cut

sub add_MapLocations {
  my ($self, @mlocs) = @_;

  foreach my $ml (@mlocs) {
    unless($ml && ref $ml && $ml->isa('Bio::EnsEMBL::Map::MapLocation')) {
      throw("args must be Bio::EnsEMBL::Map::MapLocations not [$ml]");
    }

    my $mname = $ml->map_name;
    unless($mname) {
      throw("map location arg [$ml] does not define a map name");
    }

    $self->{'map_locations'}->{$mname} = $ml;  
  }
}




=head2 flush_MapLocations

  Arg [1]    : none
  Example    : $marker->get_all_MapLocations;
  Description: Removes map locations associated with this marker.  Markers may
               be lazy-loaded from the database (again) after this.
  Returntype : none
  Exceptions : 
  Caller     : 
  Status     : stable

=cut

sub flush_MapLocations{
  my $self = shift;

  delete $self->{'map_locations'};
}


1;




