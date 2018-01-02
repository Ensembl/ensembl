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

Bio::EnsEMBL::Mapper::Unit - One side of a map pair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are each
represented by a Bio::EnsEMBL::Mapper::Unit and joined together as a
Bio::EnsEMBL::Mapper::Pair.

=head1 METHODS

=cut


package Bio::EnsEMBL::Mapper::Unit;

use strict;

sub new {
  my ( $proto, $id, $start, $end ) = @_;

  my $class = ref($proto) || $proto;

  return
    bless( { 'id' => $id, 'start' => $start, 'end' => $end }, $class );
}

=head2 id

  Arg  1      int|char $id
	      the id of the object (e.g. seq_region_name) which is mapped
  Function    accessor method
  Returntype  int|char
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
  Status      Stable

=cut

sub id {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'id'} = $value;
  }

  return $self->{'id'};
}

=head2 start

  Arg  1      int $start
	      the start coordinate of the mapped
	      region which this object represents
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
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
	      the end coordinate of the mapped
	      region which this object represents
  Function    accessor method
  Returntype  int
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Unit
  Status      Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}

1;
