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

Bio::EnsEMBL::Mapper::Pair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are each
represented by a Bio::EnsEMBL::Mapper::Unit and joined together as a
Bio::EnsEMBL::Mapper::Pair.

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::Pair;

use strict;

sub new {
  my ( $proto, $from, $to, $ori ) = @_;

  my $class = ref($proto) || $proto;

  return
    bless( { 'from' => $from, 'to' => $to, 'ori' => $ori }, $class );
}

=head2 to

  Arg  1      Bio::EnsEMBL::Mapper::Unit $seqobj
	      from and to represent the two regions
	      which are mapped to each other
  Function    accessor method
  Returntype  Bio::EnsEMBL::Mapper::Unit
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Pair
  Status     : Stable

=cut

sub to {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'to'} = $value;
  }

  return $self->{'to'};
}

=head2 from

  Arg  1      Bio::EnsEMBL::Mapper::Unit $seqobj
	      from and to represent the two regions
	      which are mapped to each other
  Function    accessor method
  Returntype  Bio::EnsEMBL::Mapper::Unit
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Pair
  Status     : Stable

=cut
sub from {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'from'} = $value;
  }

  return $self->{'from'};
}

=head2 ori

  Arg  1      Bio::EnsEMBL::Mapper::Unit $ori
  Function    accessor method
	      relative orientation of the the
	      two mapped regions
  Returntype  Bio::EnsEMBL::Mapper::Unit
  Exceptions  none
  Caller      Bio::EnsEMBL::Mapper::Pair
  Status     : Stable

=cut

sub ori {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    $self->{'ori'} = $value;
  }

  return $self->{'ori'};
}

1;
