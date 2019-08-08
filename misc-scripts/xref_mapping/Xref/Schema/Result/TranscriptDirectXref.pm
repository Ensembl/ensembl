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

package Xref::Schema::Result::TranscriptDirectXref;

=head1 NAME

Xref::Schema::Result::TranscriptDirectXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<transcript_direct_xref>

=cut

__PACKAGE__->table("transcript_direct_xref");

=head1 ACCESSORS

=head2 general_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 ensembl_stable_id

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 linkage_xref

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "general_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "ensembl_stable_id",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "linkage_xref",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.general_xref_id' });
1;
