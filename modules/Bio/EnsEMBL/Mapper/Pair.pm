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
