use utf8;
package Xref::Schema::Result::Species;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Species

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<species>

=cut

__PACKAGE__->table("species");

=head1 ACCESSORS

=head2 species_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 taxonomy_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 aliases

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "species_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "taxonomy_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "aliases",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<species_taxonomy_idx>

=over 4

=item * L</species_id>

=item * L</taxonomy_id>

=back

=cut

__PACKAGE__->add_unique_constraint("species_taxonomy_idx", ["species_id", "taxonomy_id"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:NMUKqMinB64ZXCMrhdWHJg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
