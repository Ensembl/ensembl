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

