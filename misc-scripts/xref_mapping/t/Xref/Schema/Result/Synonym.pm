use utf8;
package Xref::Schema::Result::Synonym;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Synonym

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<synonym>

=cut

__PACKAGE__->table("synonym");

=head1 ACCESSORS

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 synonym

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "synonym_relation_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0, is_auto_increment => 1 },
  "xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "synonym",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<xref_id>

=over 4

=item * L</xref_id>

=item * L</synonym>

=back

=cut

__PACKAGE__->add_unique_constraint("xref_id", ["xref_id", "synonym"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:i3+AZOCmv6yvIe1X5ZBrlg
__PACKAGE__->set_primary_key('synonym_relation_id');

__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', 'xref_id' );
1;
