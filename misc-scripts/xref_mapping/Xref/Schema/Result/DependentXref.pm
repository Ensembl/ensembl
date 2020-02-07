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

package Xref::Schema::Result::DependentXref;


=head1 NAME

Xref::Schema::Result::DependentXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<dependent_xref>

=cut

__PACKAGE__->table("dependent_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 master_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 dependent_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 linkage_annotation

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 linkage_source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  'dependency_id',
  { data_type => 'integer', extra => { unsigned => 1 }, is_nullable => 0, is_auto_increment => 1 },
  "object_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "master_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "dependent_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "linkage_annotation",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "linkage_source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
);

__PACKAGE__->set_primary_key('dependency_id');


__PACKAGE__->has_one('object_xref', 'Xref::Schema::Result::ObjectXref', 'object_xref_id');
__PACKAGE__->has_one('dependent_xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.dependent_xref_id'} );
__PACKAGE__->has_one('master_xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.master_xref_id'} );
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', { 'foreign.source_id' => 'self.linkage_source_id' } );

1;
