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

package Xref::Schema::Result::DisplayXrefPriority;

=head1 NAME

Xref::Schema::Result::DisplayXrefPriority

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<display_xref_priority>

=cut

__PACKAGE__->table("display_xref_priority");

=head1 ACCESSORS

=head2 ensembl_object_type

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 0
  size: 100

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 priority

  data_type: 'smallint'
  extra: {unsigned => 1}
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  'display_xref_source_id',
  { data_type => 'integer', extra => { unsigned => 1}, is_nullable => 0, is_auto_increment => 1},
  "ensembl_object_type",
  { data_type => "varchar", default_value => "", is_nullable => 0, size => 100 },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "priority",
  { data_type => "smallint", extra => { unsigned => 1 }, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</ensembl_object_type>

=item * L</source_id>

=back

=cut

__PACKAGE__->set_primary_key('display_xref_source_id');
__PACKAGE__->add_unique_constraint('unique_idx', ["ensembl_object_type", "source_id"]);

__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
# Can we link ensembl_object_type?

1;
