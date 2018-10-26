use utf8;
package Xref::Schema::Result::TranslationStableId;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::TranslationStableId

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<translation_stable_id>

=cut

__PACKAGE__->table("translation_stable_id");

=head1 ACCESSORS

=head2 internal_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 stable_id

  data_type: 'varchar'
  is_nullable: 0
  size: 128

=cut

__PACKAGE__->add_columns(
  "internal_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "stable_id",
  { data_type => "varchar", is_nullable => 0, size => 128 },
);

=head1 PRIMARY KEY

=over 4

=item * L</internal_id>

=back

=cut

__PACKAGE__->set_primary_key("internal_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:K9FnDaWMbxu2P4XMNs61xQ


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
