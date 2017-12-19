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

Bio::EnsEMBL::XS::Utils::Tree::Interval::PP

=head1 SYNOPSIS


=head1 DESCRIPTION

Pure Perl fall back implementation of a mutable interval tree, uses
augmented AVL binary balanced trees.
 
=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Tree::Interval::PP;

use strict;

=head2 new

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  # mmhh, should probably return just the hash ref
  return bless({}, $class);
}

sub insert {
  my $self = shift;
  return "Bio::EnsEMBL::Utils::Tree::Interval::PP::insert";
}

sub find {
  my $self = shift;
  return "Bio::EnsEMBL::Utils::Tree::Interval::PP::find";
}

1;

