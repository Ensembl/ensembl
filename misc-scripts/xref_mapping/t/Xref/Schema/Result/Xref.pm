use utf8;
package Xref::Schema::Result::Xref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Xref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<xref>

=cut

__PACKAGE__->table("xref");

=head1 ACCESSORS

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 accession

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 version

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 label

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 description

  data_type: 'text'
  is_nullable: 1

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 species_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 info_type

  data_type: 'enum'
  default_value: 'NONE'
  extra: {list => ["NONE","PROJECTION","MISC","DEPENDENT","DIRECT","SEQUENCE_MATCH","INFERRED_PAIR","PROBE","UNMAPPED","COORDINATE_OVERLAP","CHECKSUM"]}
  is_nullable: 0

=head2 info_text

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 0
  size: 255

=head2 dumped

  data_type: 'enum'
  extra: {list => ["MAPPED","NO_DUMP_ANOTHER_PRIORITY","UNMAPPED_NO_MAPPING","UNMAPPED_NO_MASTER","UNMAPPED_MASTER_FAILED","UNMAPPED_NO_STABLE_ID","UNMAPPED_INTERPRO"]}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "xref_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "accession",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "version",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "label",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "description",
  { data_type => "text", is_nullable => 1 },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "species_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "info_type",
  {
    data_type => "enum",
    default_value => "NONE",
    extra => {
      list => [
        "NONE",
        "PROJECTION",
        "MISC",
        "DEPENDENT",
        "DIRECT",
        "SEQUENCE_MATCH",
        "INFERRED_PAIR",
        "PROBE",
        "UNMAPPED",
        "COORDINATE_OVERLAP",
        "CHECKSUM",
      ],
    },
    is_nullable => 0,
  },
  "info_text",
  { data_type => "varchar", default_value => "", is_nullable => 0, size => 255 },
  "dumped",
  {
    data_type => "enum",
    extra => {
      list => [
        "MAPPED",
        "NO_DUMP_ANOTHER_PRIORITY",
        "UNMAPPED_NO_MAPPING",
        "UNMAPPED_NO_MASTER",
        "UNMAPPED_MASTER_FAILED",
        "UNMAPPED_NO_STABLE_ID",
        "UNMAPPED_INTERPRO",
      ],
    },
    is_nullable => 1,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</xref_id>

=back

=cut

__PACKAGE__->set_primary_key("xref_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<acession_idx>

=over 4

=item * L</accession>

=item * L</label>

=item * L</source_id>

=item * L</species_id>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "acession_idx",
  ["accession", "label", "source_id", "species_id"],
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:p9mn5/VyM92UXNPGE2DjXQ

__PACKAGE__->might_have('object_xref', 'Xref::Schema::Result::ObjectXref','xref_id'); # Not true until after the mapping is completed
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id');
__PACKAGE__->has_many('synonym', 'Xref::Schema::Result::Synonym', 'xref_id' );

1;
