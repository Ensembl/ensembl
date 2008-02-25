#
# Ensembl module for Bio::EnsEMBL::Mapper::IndelPair
#
# Written by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::IndelPair

=head1 SYNOPSIS

=head1 DESCRIPTION

Two regions mapped between different coordinate systems are
each represented by a Bio::EnsEMBL::Mapper::Unit and joined
together as a Bio::EnsEMBL::Mapper::Pair, when one of the 
regions is an indel.

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Post general queries to B<ensembl-dev@ebi.ac.uk>

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
