use utf8;
package Xref::Schema::Result::Meta;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Meta

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<meta>

=cut

__PACKAGE__->table("meta");

=head1 ACCESSORS

=head2 meta_id

  data_type: 'integer'
  is_auto_increment: 1
  is_nullable: 0

=head2 species_id

  data_type: 'integer'
  default_value: 1
  extra: {unsigned => 1}
  is_nullable: 1

=head2 meta_key

  data_type: 'varchar'
  is_nullable: 0
  size: 40

=head2 meta_value

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 date

  data_type: 'datetime'
  datetime_undef_if_invalid: 1
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  "meta_id",
  { data_type => "integer", is_auto_increment => 1, is_nullable => 0 },
  "species_id",
  {
    data_type => "integer",
    default_value => 1,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
  "meta_key",
  { data_type => "varchar", is_nullable => 0, size => 40 },
  "meta_value",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "date",
  {
    data_type => "datetime",
    datetime_undef_if_invalid => 1,
    is_nullable => 0,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</meta_id>

=back

=cut

__PACKAGE__->set_primary_key("meta_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<species_key_value_idx>

=over 4

=item * L</meta_id>

=item * L</species_id>

=item * L</meta_key>

=item * L</meta_value>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "species_key_value_idx",
  ["meta_id", "species_id", "meta_key", "meta_value"],
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:zf12byEOvyRayGNCaO1Nvw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
