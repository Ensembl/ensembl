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

Bio::EnsEMBL::Mapper::Gap

=head1 SYNOPSIS

=head1 DESCRIPTION

Representation of a gap in a sequence; returned from Mapper.pm when the
target region is in a gap.

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::Gap;

use strict;

=head2 new

  Arg [1]    : int $start
  Arg [2]    : int $end
  Example    : $gap = Bio::EnsEMBL::Mapper::Gap($start, $end);
  Description: Creates a new Gap object.
  Returntype : Bio::EnsEMBL::Mapper::Gap
  Exceptions : none
  Caller     : Bio::EnsEMBL::Mapper
  Status     : Stable

=cut

sub new {
  my ( $proto, $start, $end, $rank ) = @_;

  my $class = ref($proto) || $proto;

  return bless( { 'start' => $start, 'end' => $end, 'rank' => $rank  || 0 }, $class );
}

=head2 start

  Arg [1]    : (optional) int $start
               start coordinate of gap region
  Example    : $start = $gap->start();
  Description: Getter/Setter for the start attribute
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'start'} = $value;
  }

  return $self->{'start'};
}

=head2 end

  Arg [1]    : (optional) int $newval
               The new value to set the end coordinate to
  Example    : $end = $gap->end()
  Description: Getter/Setter for the end coordinate of the gap region
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'end'} = $value;
  }

  return $self->{'end'};
}

=head2 length

  Arg [1]    : none
  Example    : $len = $gap->length();
  Description: Getter for the length of this gap region
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

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
