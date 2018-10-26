use utf8;
package Xref::Schema::Result::ProcessStatus;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::ProcessStatus

=cut

use strict;
use warnings;

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


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:LaqSMJg2ID1qohxQHiNRYQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
