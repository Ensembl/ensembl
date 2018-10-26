use utf8;
package Xref::Schema::Result::SourceUrl;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::SourceUrl

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source_url>

=cut

__PACKAGE__->table("source_url");

=head1 ACCESSORS

=head2 source_url_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 species_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 url

  data_type: 'mediumtext'
  is_nullable: 1

=head2 release_url

  data_type: 'mediumtext'
  is_nullable: 1

=head2 checksum

  data_type: 'varchar'
  is_nullable: 1
  size: 1025

=head2 file_modified_date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 upload_date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 1

=head2 parser

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "source_url_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "species_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "url",
  { data_type => "mediumtext", is_nullable => 1 },
  "release_url",
  { data_type => "mediumtext", is_nullable => 1 },
  "checksum",
  { data_type => "varchar", is_nullable => 1, size => 1025 },
  "file_modified_date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "upload_date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 1,
  },
  "parser",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 PRIMARY KEY

=over 4

=item * L</source_url_id>

=back

=cut

__PACKAGE__->set_primary_key("source_url_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:trrjT7dUapyKlNbMafZUXQ


__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
1;
