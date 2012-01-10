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

Bio::EnsEMBL::Utils::BitString - bitstring object implementation

=head1 DESCRIPTION

This is an implementation of a bitstring object, taken from Damian
Convey's book "Object Oriented Perl".

=head1 METHODS

=cut


package Bio::EnsEMBL::Utils::BitString;

use strict;
use warnings;
no warnings 'uninitialized';


sub new {
  my $class = ref($_[0])||$_[0];
  my $initbits = join '', map {$_?1:0} @_[1..$#_];
  my $bs = pack 'b*', $initbits;
  bless \$bs, $class;
}


sub get {
  my ($self, $bitnum) = @_;
  return vec($$self,$bitnum,1);
}


sub set {
  my ($self, $bitnum, $newval) = @_;
  vec($$self,$bitnum,1) = $newval?1:0;
}


sub bitcount {
  8 * length ${$_[0]};
}


sub complement {
  my ($self) = @_;
  my $complement = ~$$self;
  bless \$complement, ref($self);
}


sub print_me {
  my ($self) = @_;
  for (my $i=0; $i < $self->bitcount(); $i++)
  {
    print $self->get($i);
    print ' ' unless ($i+1)%8;
    print "\n" unless ($i+1)%64;
  }
  print "\n";
}


1;

