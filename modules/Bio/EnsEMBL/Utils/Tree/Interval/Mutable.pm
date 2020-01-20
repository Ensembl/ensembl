=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Utils::Tree::Interval::Mutable

=head1 SYNOPSIS

  # start with an empty tree
  my $tree = Bio::EnsEMBL::Utils::Tree::Interval::Mutable->new();

  # add a few intervals (i.e. Bio::EnsEMBL::Utils::Interval)
  $tree->insert($i1);
  $tree->insert($i2);
  $tree->insert($i3);

  # query the tree
  my $result = $tree->search(85, 100);
  if (scalar @{$result}) {
    print "Found overlapping interval: [", $result->[0]->start, ', ', $result->[0]->end, "\n";
  }

=head1 DESCRIPTION

Class representing a dynamic, i.e. mutable, interval tree implemented as an augmented AVL balanced binary tree.

This module is a wrapper around two possible implementations: one using the Perl extension (XS) mechanisms, and
a pure Perl (PP) one. 

The module is capable of detecting whether the XS module is available and it loads it in that
case; it falls back to the PP implementation otherwise.
 
=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::Mutable;

use strict;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);

# the modules providing the underlying implementation,
# either XS or pure perl fallback
my $XS = 'Bio::EnsEMBL::XS::Utils::Tree::Interval::Mutable';
my $PP = 'Bio::EnsEMBL::Utils::Tree::Interval::Mutable::PP';

# if XS is used, version at least 1.3.1 is required (provides the interval tree library)
my $VERSION_XS = '1.3.1';

my @public_methods = qw/ insert search remove size /;

# import either XS or PP methods into namespace
unless ($Bio::EnsEMBL::Utils::Tree::Interval::Mutable::IMPL) {
  # first check if XS is available and try to load it,
  # otherwise fall back to PP implementation
  _load_xs() or _load_pp() or throw "Couldn't load implementation: $@";
}


=head2 new

  Arg []      : none
  Example     : my $tree = Bio::EnsEMBL::Utils::Tree::Mutable();
  Description : Constructor. Creates a new mutable tree instance
  Returntype  : Bio::EnsEMBL::Utils::Tree::Interval::Mutable
  Exceptions  : none
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  # for ($XS|$PP)::new(0);
  return eval qq| $Bio::EnsEMBL::Utils::Tree::Interval::Mutable::IMPL\::new( \$caller ) | unless $caller; ## no critic

  if (my $self = $Bio::EnsEMBL::Utils::Tree::Interval::Mutable::IMPL->new(@_)) {
    $self->{_IMPL} = $Bio::EnsEMBL::Utils::Tree::Interval::Mutable::IMPL;
    bless($self, $class);
    return $self
  }

  return;
}

sub _load_xs {
  _load($XS, $VERSION_XS);
}

sub _load_pp {
  _load($PP);
}

sub _load {
  my ($module, $version) = @_;
  $version ||= '';

  eval qq| use $module $version |; ## no critic
  info(sprintf("Cannot load %s interval tree implementation", $module eq $XS?'XS':'PP'), 2000)
    and return if $@;

  push @Bio::EnsEMBL::Utils::Tree::Interval::Mutable::ISA, $module;
  $Bio::EnsEMBL::Utils::Tree::Interval::Mutable::IMPL = $module;

  local $^W;
  no strict qw(refs); ## no critic

  for my $method (@public_methods) {
    *{"Bio::EnsEMBL::Utils::Tree::Interval::Mutable::$method"} = \&{"$module\::$method"};
  }
  
  return 1;
}

1;

