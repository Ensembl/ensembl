use utf8;
package Xref::Schema::Result::GoXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::GoXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<go_xref>

=cut

__PACKAGE__->table("go_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  default_value: 0
  extra: {unsigned => 1}
  is_nullable: 0

=head2 linkage_type

  data_type: 'char'
  is_nullable: 0
  size: 3

=head2 source_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "go_xref_link_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0, is_auto_increment => 1 },
  "object_xref_id",
  {
    data_type => "integer",
    default_value => 0,
    extra => { unsigned => 1 },
    is_nullable => 0,
  },
  "linkage_type",
  { data_type => "char", is_nullable => 0, size => 3 },
  "source_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<object_xref_id_2>

=over 4

=item * L</object_xref_id>

=item * L</source_xref_id>

=item * L</linkage_type>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "object_xref_id_2",
  ["object_xref_id", "source_xref_id", "linkage_type"],
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:7zjw9MqtfPIsKbvVlRiBUw

__PACKAGE__->set_primary_key('go_xref_link_id');

__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
__PACKAGE__->has_one('object_xref', 'Xref::Schema::Result::ObjectXref', 'object_xref_id' );

1;
