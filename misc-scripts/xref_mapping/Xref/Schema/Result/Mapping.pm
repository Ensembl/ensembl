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

package Xref::Schema::Result::Mapping;

=head1 NAME

Xref::Schema::Result::Mapping

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<mapping>

=cut

__PACKAGE__->table("mapping");

=head1 ACCESSORS

=head2 job_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 type

  data_type: 'enum'
  extra: {list => ["dna","peptide","UCSC"]}
  is_nullable: 1

=head2 command_line

  data_type: 'text'
  is_nullable: 1

=head2 percent_query_cutoff

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 percent_target_cutoff

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 method

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 array_size

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "job_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "type",
  {
    data_type => "enum",
    extra => { list => ["dna", "peptide", "UCSC"] },
    is_nullable => 1,
  },
  "command_line",
  { data_type => "text", is_nullable => 1 },
  "percent_query_cutoff",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "percent_target_cutoff",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "method",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "array_size",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);

1;
