use utf8;
package Xref::Schema::Result::TranscriptStableId;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::TranscriptStableId

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<transcript_stable_id>

=cut

__PACKAGE__->table("transcript_stable_id");

=head1 ACCESSORS

=head2 internal_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 0
  size: 128

=head2 display_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 biotype

  data_type: 'varchar'
  is_nullable: 0
  size: 40

=cut

__PACKAGE__->add_columns(
  "internal_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "stable_id",
  { data_type => "varchar", is_nullable => 0, size => 128 },
  "display_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "biotype",
  { data_type => "varchar", is_nullable => 0, size => 40 },
);

=head1 PRIMARY KEY

=over 4

=item * L</stable_id>

=back

=cut

__PACKAGE__->set_primary_key("stable_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:dly7ooI2KAJxP+mDXIo8dw


# __PACKAGE__->might_have('display_xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.display_xref_id' });
1;
