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

Bio::EnsEMBL::Mapper::Coordinate

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a mapped region in a sequence; returned from Mapper.pm
when the target region maps on to valid sequence.

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::Coordinate;

use strict;

=head2 new

  Arg [1]     char|int id of object in mapped region
  Arg [2]     int start of region   
  Arg [3]     int end of region
  Arg [4]     int strand if region
  Arg [5]     Bio::EnsEMBL::CoordSystem  coordsytem of the region
  Function    creates a new Coordinate object
  Returntype  Bio::EnsEMBL::Mapper::Coordinate
  Exceptions  none
  Status      Stable

=cut

sub new {
  my ( $proto, $id, $start, $end, $strand, $coord_system, $rank ) = @_;

  my $class = ref($proto) || $proto;

  return
    bless( { 'id'           => $id,
             'start'        => $start,
             'end'          => $end,
             'strand'       => $strand,
             'coord_system' => $coord_system,
             'rank'         => $rank || 0
           },
           $class );
}


=head2 start

  Arg  1      int $start
              start coordinate of object in mapped region
  Function    getter/setter method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start'} = $value;
  }

  return $self->{'start'};
}


=head2 end

  Arg  1      int $end
              end coordinate of object in mapped region
  Function    getter/setter method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}

=head2 strand

  Arg  1      int $strand
              strand of object in mapped region
  Function    getter/setter method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub strand {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'strand'} = $value;
  }

  return $self->{'strand'};
}

=head2 id

  Arg  1      char|int $id
              id of object in mapped region
	      e.g. seq_region_id
  Function    getter/setter method
  Returntype  char|int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub id {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'id'} = $value;
  }

  return $self->{'id'};
}

=head2 name

  Arg  1      char name
              name of object in mapped region
	      e.g. seq_region_name
  Function    getter/setter method
  Returntype  char
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub name {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'name'} = $value;
  }

  return $self->{'name'};
}

=head2 coord_system

  Arg  1      Bio::EnsEMBL::CoordSystem
  Function    getter/setter method
  Returntype  Bio::EnsEMBL::CoordSystem
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Coordinate
  Status      Stable

=cut

sub coord_system {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'coord_system'} = $value;
  }

  return $self->{'coord_system'};
}

=head2 length

  Function    getter method
  Returntype  int
  Exceptions  none
  Caller      ?
  Status      Stable

=cut

sub length {
  my ($self) = @_;

  return $self->{'end'} - $self->{'start'} + 1;
}

sub rank {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'rank'} = $value;
  }

  return $self->{'rank'};
}

1;
