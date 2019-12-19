=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

package Xref::Schema::Result::AltAllele;


=head1 NAME

Xref::Schema::Result::AltAllele

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<alt_allele>

=cut

__PACKAGE__->table("alt_allele");

=head1 ACCESSORS

=head2 alt_allele_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 is_reference

  data_type: 'integer'
  default_value: 0
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "alt_allele_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "gene_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "is_reference",
  {
    data_type => "integer",
    default_value => 0,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<allele_idx>

=over 4

=item * L</alt_allele_id>

=item * L</gene_id>

=back

=cut

__PACKAGE__->add_unique_constraint("allele_idx", ["alt_allele_id", "gene_id"]);

=head2 C<gene_idx>

=over 4

=item * L</gene_id>

=back

=cut

# __PACKAGE__->add_unique_constraint("gene_idx", ["gene_id"]);

__PACKAGE__->has_one('gene', 'Xref::Schema::Result::GeneStableId', { 'foreign.internal_id' => 'self.gene_id'});

1;
