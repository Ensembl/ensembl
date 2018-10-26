use utf8;
package Xref::Schema::Result::CoordinateXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::CoordinateXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<coordinate_xref>

=cut

__PACKAGE__->table("coordinate_xref");

=head1 ACCESSORS

=head2 coord_xref_id

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

=head2 accession

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 chromosome

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 strand

  data_type: 'tinyint'
  is_nullable: 0

=head2 txstart

  data_type: 'integer'
  is_nullable: 0

=head2 txend

  data_type: 'integer'
  is_nullable: 0

=head2 cdsstart

  data_type: 'integer'
  is_nullable: 1

=head2 cdsend

  data_type: 'integer'
  is_nullable: 1

=head2 exonstarts

  data_type: 'text'
  is_nullable: 0

=head2 exonends

  data_type: 'text'
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "coord_xref_id",
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
  "accession",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "chromosome",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "strand",
  { data_type => "tinyint", is_nullable => 0 },
  "txstart",
  { data_type => "integer", is_nullable => 0 },
  "txend",
  { data_type => "integer", is_nullable => 0 },
  "cdsstart",
  { data_type => "integer", is_nullable => 1 },
  "cdsend",
  { data_type => "integer", is_nullable => 1 },
  "exonstarts",
  { data_type => "text", is_nullable => 0 },
  "exonends",
  { data_type => "text", is_nullable => 0 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<coord_xref_idx>

=over 4

=item * L</coord_xref_id>

=back

=cut

# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:ti2S4oUA/BNLH96uYmMe6A
__PACKAGE__->set_primary_key('coord_xref_id');
__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.coord_xref_id' });
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
1;
