=head1 LICENSE

See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
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

package Xref::Schema::Result::GeneStableId;

=head1 NAME

Xref::Schema::Result::GeneStableId

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<gene_stable_id>

=cut

__PACKAGE__->table("gene_stable_id");

=head1 ACCESSORS

=head2 internal_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 0
  size: 128

=head2 display_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 desc_set

  data_type: 'integer'
  default_value: 0
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "internal_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "stable_id",
  { data_type => "varchar", is_nullable => 0, size => 128 },
  "display_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "desc_set",
  {
    data_type => "integer",
    default_value => 0,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</stable_id>

=back

=cut

__PACKAGE__->set_primary_key("internal_id");

# __PACKAGE__->might_have('display_xref', 'Xref::Schema::Result::Xref', {'foreign.xref_id' => 'self.display_xref_id'} );

1;
