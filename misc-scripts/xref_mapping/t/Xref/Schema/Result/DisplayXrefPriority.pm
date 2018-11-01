use utf8;
package Xref::Schema::Result::DisplayXrefPriority;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::DisplayXrefPriority

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<display_xref_priority>

=cut

__PACKAGE__->table("display_xref_priority");

=head1 ACCESSORS

=head2 ensembl_object_type

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 0
  size: 100

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 priority

  data_type: 'smallint'
  extra: {unsigned => 1}
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  'display_xref_source_id',
  { data_type => 'integer', extra => { unsigned => 1}, is_nullable => 0, is_auto_increment => 1},
  "ensembl_object_type",
  { data_type => "varchar", default_value => "", is_nullable => 0, size => 100 },
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "priority",
  { data_type => "smallint", extra => { unsigned => 1 }, is_nullable => 0 },
);

=head1 PRIMARY KEY

=over 4

=item * L</ensembl_object_type>

=item * L</source_id>

=back

=cut

__PACKAGE__->set_primary_key('display_xref_source_id');
__PACKAGE__->add_unique_constraint('unique_idx', ["ensembl_object_type", "source_id"]);

# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:z9CuvQ6YrJOYr/ictX9n1Q


__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
# Can we link ensembl_object_type?

1;
