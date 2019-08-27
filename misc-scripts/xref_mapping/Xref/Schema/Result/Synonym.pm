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

package Xref::Schema::Result::Synonym;

=head1 NAME

Xref::Schema::Result::Synonym

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<synonym>

=cut

__PACKAGE__->table("synonym");

=head1 ACCESSORS

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 synonym

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "synonym_relation_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0, is_auto_increment => 1 },
  "xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "synonym",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<xref_id>

=over 4

=item * L</xref_id>

=item * L</synonym>

=back

=cut

__PACKAGE__->add_unique_constraint("xref_id", ["xref_id", "synonym"]);

__PACKAGE__->set_primary_key('synonym_relation_id');

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', 'xref_id' );
1;
