use utf8;
package Xref::Schema::Result::TranslationDirectXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::TranslationDirectXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<translation_direct_xref>

=cut

__PACKAGE__->table("translation_direct_xref");

=head1 ACCESSORS

=head2 general_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 ensembl_stable_id

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 linkage_xref

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "general_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "ensembl_stable_id",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "linkage_xref",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:4foz7mMkOoK2uFPtiTfrWw


__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.general_xref_id' });
1;
