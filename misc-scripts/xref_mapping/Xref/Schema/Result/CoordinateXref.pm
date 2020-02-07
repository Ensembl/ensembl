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

package Xref::Schema::Result::CoordinateXref;

=head1 NAME

Xref::Schema::Result::CoordinateXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<coordinate_xref>

=cut

__PACKAGE__->table("coordinate_xref");

=head1 ACCESSORS

=head2 coord_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 species_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 accession

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 chromosome

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 strand

  data_type: 'tinyint'
  is_nullable: 0

=head2 txstart

  data_type: 'integer'
  is_nullable: 0

=head2 txend

  data_type: 'integer'
  is_nullable: 0

=head2 cdsstart

  data_type: 'integer'
  is_nullable: 1

=head2 cdsend

  data_type: 'integer'
  is_nullable: 1

=head2 exonstarts

  data_type: 'text'
  is_nullable: 0

=head2 exonends

  data_type: 'text'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "coord_xref_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "species_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "accession",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "chromosome",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "strand",
  { data_type => "tinyint", is_nullable => 0 },
  "txstart",
  { data_type => "integer", is_nullable => 0 },
  "txend",
  { data_type => "integer", is_nullable => 0 },
  "cdsstart",
  { data_type => "integer", is_nullable => 1 },
  "cdsend",
  { data_type => "integer", is_nullable => 1 },
  "exonstarts",
  { data_type => "text", is_nullable => 0 },
  "exonends",
  { data_type => "text", is_nullable => 0 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<coord_xref_idx>

=over 4

=item * L</coord_xref_id>

=back

=cut

__PACKAGE__->set_primary_key('coord_xref_id');

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.coord_xref_id' });
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
1;
