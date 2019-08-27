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

package Xref::Schema::Result::PrimaryXref;

=head1 NAME

Xref::Schema::Result::PrimaryXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<primary_xref>

=cut

__PACKAGE__->table("primary_xref");

=head1 ACCESSORS

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 sequence

  accessor: undef
  data_type: 'mediumtext'
  is_nullable: 1

=head2 sequence_type

  data_type: 'enum'
  extra: {list => ["dna","peptide"]}
  is_nullable: 1

=head2 status

  data_type: 'enum'
  extra: {list => ["experimental","predicted"]}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "sequence",
  { data_type => "mediumtext", is_nullable => 1 },
  "sequence_type",
  {
    data_type => "enum",
    extra => { list => ["dna", "peptide"] },
    is_nullable => 1,
  },
  "status",
  {
    data_type => "enum",
    extra => { list => ["experimental", "predicted"] },
    is_nullable => 1,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</xref_id>

=back

=cut

__PACKAGE__->set_primary_key("xref_id");

__PACKAGE__->add_column('+sequence' => {accessor => 'raw_seq'});
__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', 'xref_id' );
1;
