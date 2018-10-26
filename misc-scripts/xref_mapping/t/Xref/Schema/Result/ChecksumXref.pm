use utf8;
package Xref::Schema::Result::ChecksumXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::ChecksumXref

=cut

use strict;
use warnings;

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


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:eaHDk7+LBk96mKw29RH3/A

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.checksum_xref_id' });
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );

1;
