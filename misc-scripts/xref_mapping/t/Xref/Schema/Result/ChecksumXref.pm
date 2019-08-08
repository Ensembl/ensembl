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

package Xref::Schema::Result::ChecksumXref;

=head1 NAME

Xref::Schema::Result::ChecksumXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<checksum_xref>

=cut

__PACKAGE__->table("checksum_xref");

=head1 ACCESSORS

=head2 checksum_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 accession

  data_type: 'char'
  is_nullable: 0
  size: 14

=head2 checksum

  data_type: 'char'
  is_nullable: 0
  size: 32

=cut

__PACKAGE__->add_columns(
  "checksum_xref_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "accession",
  { data_type => "char", is_nullable => 0, size => 14 },
  "checksum",
  { data_type => "char", is_nullable => 0, size => 32 },
);

=head1 PRIMARY KEY

=over 4

=item * L</checksum_xref_id>

=back

=cut

__PACKAGE__->set_primary_key("checksum_xref_id");

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.checksum_xref_id' });
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );

1;
