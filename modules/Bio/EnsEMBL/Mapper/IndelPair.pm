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

Bio::EnsEMBL::Mapper::IndelPair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are each
represented by a Bio::EnsEMBL::Mapper::Unit and joined together as a
Bio::EnsEMBL::Mapper::Pair, when one of the regions is an indel.

=head1 METHODS

=cut

package Bio::EnsEMBL::Mapper::IndelPair;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::Mapper::Pair);

sub new {
  my ($proto, @args) = @_;

  my $class = ref($proto) || $proto;

  my $self = $class->SUPER::new(@args);    # create the Pair object
  $self->{'indel'} = 1;                    # and add the Indel flag

  return $self;
}

1;
