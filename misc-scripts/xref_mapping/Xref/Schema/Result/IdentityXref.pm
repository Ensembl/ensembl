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

package Xref::Schema::Result::IdentityXref;

=head1 NAME

Xref::Schema::Result::IdentityXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<identity_xref>

=cut

__PACKAGE__->table("identity_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 query_identity

  data_type: 'integer'
  is_nullable: 1

=head2 target_identity

  data_type: 'integer'
  is_nullable: 1

=head2 hit_start

  data_type: 'integer'
  is_nullable: 1

=head2 hit_end

  data_type: 'integer'
  is_nullable: 1

=head2 translation_start

  data_type: 'integer'
  is_nullable: 1

=head2 translation_end

  data_type: 'integer'
  is_nullable: 1

=head2 cigar_line

  data_type: 'text'
  is_nullable: 1

=head2 score

  data_type: 'double precision'
  is_nullable: 1

=head2 evalue

  data_type: 'double precision'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "object_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "query_identity",
  { data_type => "integer", is_nullable => 1 },
  "target_identity",
  { data_type => "integer", is_nullable => 1 },
  "hit_start",
  { data_type => "integer", is_nullable => 1 },
  "hit_end",
  { data_type => "integer", is_nullable => 1 },
  "translation_start",
  { data_type => "integer", is_nullable => 1 },
  "translation_end",
  { data_type => "integer", is_nullable => 1 },
  "cigar_line",
  { data_type => "text", is_nullable => 1 },
  "score",
  { data_type => "double precision", is_nullable => 1 },
  "evalue",
  { data_type => "double precision", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</object_xref_id>

=back

=cut

__PACKAGE__->set_primary_key("object_xref_id");

__PACKAGE__->has_one('object_xref', 'Xref::Schema::Result::ObjectXref', 'object_xref_id' );
1;
