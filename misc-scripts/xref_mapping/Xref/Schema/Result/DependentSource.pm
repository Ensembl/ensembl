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

package Xref::Schema::Result::DependentSource;

=head1 NAME

Xref::Schema::Result::DependentSource

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<dependent_source>

=cut

__PACKAGE__->table("dependent_source");

=head1 ACCESSORS

=head2 master_source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 dependent_name

  data_type: 'varchar'
  is_nullable: 0
  description: Name of the source that must be imported before this source can be processed
  size: 255

=cut

__PACKAGE__->add_columns(
  "dependent_relation_id",
  { data_type => "integer", extra => { unsigned => 1, auto_increment => 1}, is_nullable => 0},
  "master_source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "dependent_name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
);

=head1 PRIMARY KEY

=over 4

=item * L</master_source_id>

=item * L</dependent_name>

=back

=cut

__PACKAGE__->set_primary_key('dependent_relation_id');
# __PACKAGE__->add_unique_constraint("master_source_id", ["dependent_name"]);

__PACKAGE__->belongs_to('master_source', 'Xref::Schema::Result::Source', { 'foreign.source_id' => 'self.master_source_id' } );

1;
