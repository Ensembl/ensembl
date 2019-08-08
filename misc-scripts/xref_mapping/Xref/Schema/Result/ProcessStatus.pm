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

package Xref::Schema::Result::ProcessStatus;

=head1 NAME

Xref::Schema::Result::ProcessStatus

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<process_status>

=cut

__PACKAGE__->table("process_status");

=head1 ACCESSORS

=head2 id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 status

  data_type: 'enum'
  extra: {list => ["xref_created","parsing_started","parsing_finished","alt_alleles_added","xref_fasta_dumped","core_fasta_dumped","core_data_loaded","mapping_submitted","mapping_finished","mapping_processed","direct_xrefs_parsed","prioritys_flagged","processed_pairs","biomart_test_finished","source_level_move_finished","alt_alleles_processed","official_naming_done","checksum_xrefs_started","checksum_xrefs_finished","coordinate_xrefs_started","coordinate_xref_finished","tests_started","tests_failed","tests_finished","core_loaded","display_xref_done","gene_description_done"]}
  is_nullable: 1

=head2 date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "status",
  {
    data_type => "enum",
    extra => {
      list => [
        "xref_created",
        "parsing_started",
        "parsing_finished",
        "alt_alleles_added",
        "xref_fasta_dumped",
        "core_fasta_dumped",
        "core_data_loaded",
        "mapping_submitted",
        "mapping_finished",
        "mapping_processed",
        "direct_xrefs_parsed",
        "prioritys_flagged",
        "processed_pairs",
        "biomart_test_finished",
        "source_level_move_finished",
        "alt_alleles_processed",
        "official_naming_done",
        "checksum_xrefs_started",
        "checksum_xrefs_finished",
        "coordinate_xrefs_started",
        "coordinate_xref_finished",
        "tests_started",
        "tests_failed",
        "tests_finished",
        "core_loaded",
        "display_xref_done",
        "gene_description_done",
      ],
    },
    is_nullable => 1,
  },
  "date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 0,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</id>

=back

=cut

__PACKAGE__->set_primary_key("id");

1;
