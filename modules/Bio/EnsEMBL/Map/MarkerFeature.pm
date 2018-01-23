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

Bio::EnsEMBL::Map::MarkerFeature

=head1 SYNOPSIS

=head1 DESCRIPTION

Represents a marker feature in the EnsEMBL database.  A marker feature
is a marker which has been mapped to the genome by ePCR.  Each marker
has one marker feature per mapped location.

=head1 METHODS

=cut

package Bio::EnsEMBL::Map::MarkerFeature;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Feature;

@ISA = qw(Bio::EnsEMBL::Feature);



=head2 new

  Arg [1]    : (optional) int $dbID
  Arg [2]    : (optional) Bio::EnsEMBL::Adaptor $adaptor
  Arg [3]    : (optional) int $start
  Arg [4]    : (optional) int $end
  Arg [5]    : (optional) Bio::EnsEMBL::Slice $slice
  Arg [6]    : (optional) Bio::EnsEMBL::Analysis
  Arg [7]    : (optional) int $marker_id
  Arg [8]    : (optional) int $map_weight
  Arg [9]    : (optional) Bio::EnsEMBL::Map::Marker $marker 
  Example    : $marker = Bio::EnsEMBL::Map::MarkerFeature->new(123, $adaptor,
							       100, 200, 
							       $ctg, 123);
  Description: Creates a new MarkerFeature
  Returntype : Bio::EnsEMBL::Map::MarkerFeature
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub new {
  my ($caller, $dbID, $adaptor, $start, $end, $slice, $analysis,
      $marker_id, $map_weight, $marker) = @_;

  my $class = ref($caller) || $caller;

  my $self = bless( {
		 'dbID'        => $dbID,
		 'start'       => $start,
		 'end'         => $end,
		 'strand'      => 0,
		 'slice'       => $slice,
		 'analysis'    => $analysis,
		 'marker_id'   => $marker_id,
		 'marker'      => $marker,
     'map_weight'  => $map_weight }, $class);

  $self->adaptor($adaptor);
  return $self;
}



=head2 _marker_id

  Arg [1]    : (optional) int $marker_id
  Example    : none
  Description: PRIVATE Getter/Setter for the internal identifier of the marker
               associated with this marker feature
  Returntype : int
  Exceptions : none
  Caller     : internal
  Status     : stable

=cut

sub _marker_id {
  my $self = shift;

  if(@_) {
    $self->{'marker_id'} = shift;
  }

  return $self->{'marker_id'};
}



=head2 marker

  Arg [1]    : (optional) Bio::EnsEMBL::Map::Marker $marker
  Example    : $marker = $marker_feature->marker;
  Description: Getter/Setter for the marker associated with this marker feature
               If the marker has not been set and the database is available
               the marker will automatically be retrieved (lazy-loaded).
  Returntype : Bio::EnsEMBL::Map::Marker
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub marker {
  my $self = shift;

  if(@_) {
    $self->{'marker'} = shift;
  } elsif(!$self->{'marker'} && $self->adaptor && $self->{'marker_id'}) {
    #lazy load the marker if it is not already loaded
    my $ma = $self->adaptor->db->get_MarkerAdaptor;
    $self->{'marker'} = $ma->fetch_by_dbID($self->{'marker_id'});
  }

  return $self->{'marker'};
}



=head2 map_weight

  Arg [1]    : (optional) int $map_weight
  Example    : $map_weight = $marker_feature->map_weight;
  Description: Getter/Setter for the map weight of this marker.  This is the
               number of times that this marker has been mapped to the genome.
               E.g.  a marker iwth map weight 3 has been mapped to 3 locations
               in the genome.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub map_weight {
  my $self = shift;

  if(@_) {
    $self->{'map_weight'} = shift;
  }

  return $self->{'map_weight'};
}



=head2 display_id

  Arg [1]    : none
  Example    : print $mf->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For marker features this is the
               name of the display synonym or '' if it is not defined.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : stable

=cut

sub display_id {
  my $self = shift;
  my $marker = $self->{'marker'};

  return '' if(!$marker);
  my $ms = $marker->display_MarkerSynonym();
  return '' if(!$ms);
  return $ms->name() || '';
}


1;
