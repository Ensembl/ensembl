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

package Xref::Schema::Result::SourceMappingMethod;

=head1 NAME

Xref::Schema::Result::SourceMappingMethod

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source_mapping_method>

=cut

__PACKAGE__->table("source_mapping_method");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 method

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "method",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<source_id>

=over 4

=item * L</source_id>

=item * L</method>

=back

=cut

__PACKAGE__->add_unique_constraint("source_id", ["source_id", "method"]);

1;
