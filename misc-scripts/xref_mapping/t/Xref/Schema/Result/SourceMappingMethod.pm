use utf8;
package Xref::Schema::Result::SourceMappingMethod;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::SourceMappingMethod

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source_mapping_method>

=cut

__PACKAGE__->table("source_mapping_method");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 method

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "method",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<source_id>

=over 4

=item * L</source_id>

=item * L</method>

=back

=cut

__PACKAGE__->add_unique_constraint("source_id", ["source_id", "method"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:L6QvTS7L+O4rueKfkRiNyg


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
