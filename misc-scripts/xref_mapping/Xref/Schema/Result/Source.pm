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

package Xref::Schema::Result::Source;

=head1 NAME

Xref::Schema::Result::Source

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source>

=cut

__PACKAGE__->table("source");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 source_release

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 priority

  data_type: 'integer'
  default_value: 1
  extra: {unsigned => 1}
  is_nullable: 1

=head2 priority_description

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 1
  size: 40

=cut

__PACKAGE__->add_columns(
  "source_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "source_release",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "priority",
  {
    data_type => "integer",
    default_value => 1,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
  "priority_description",
  { data_type => "varchar", is_nullable => 1, size => 40 },
);

=head1 PRIMARY KEY

=over 4

=item * L</source_id>

=back

=cut

__PACKAGE__->set_primary_key('source_id');

__PACKAGE__->has_many('xrefs', 'Xref::Schema::Result::Xref', 'source_id');
__PACKAGE__->might_have(
  'dependent_source',
  'Xref::Schema::Result::DependentSource',
  { 'foreign.master_source_id' => 'self.source_id' }
);
__PACKAGE__->might_have(
  'source_url',
  'Xref::Schema::Result::SourceUrl',
  { 'foreign.source_id' => 'self.source_id' }
);

1;
