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

package Xref::Schema::Result::Pair;

=head1 NAME

Xref::Schema::Result::Pair

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<pairs>

=cut

__PACKAGE__->table("pairs");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 accession1

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 accession2

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=cut

__PACKAGE__->add_columns(
  "pair_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 , is_auto_increment => 1},
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "accession1",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "accession2",
  { data_type => "varchar", is_nullable => 0, size => 255 },
);


__PACKAGE__->set_primary_key('pair_id');
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );

1;
