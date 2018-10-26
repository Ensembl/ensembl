use utf8;
package Xref::Schema::Result::Source;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Source

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source>

=cut

__PACKAGE__->table("source");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 status

  data_type: 'enum'
  default_value: 'NOIDEA'
  extra: {list => ["KNOWN","XREF","PRED","ORTH","PSEUDO","LOWEVIDENCE","NOIDEA"]}
  is_nullable: 0

=head2 source_release

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 download

  data_type: 'enum'
  default_value: 'Y'
  extra: {list => ["Y","N"]}
  is_nullable: 1

=head2 ordered

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 priority

  data_type: 'integer'
  default_value: 1
  extra: {unsigned => 1}
  is_nullable: 1

=head2 priority_description

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 1
  size: 40

=cut

__PACKAGE__->add_columns(
  "source_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "status",
  {
    data_type => "enum",
    default_value => "NOIDEA",
    extra => {
      list => ["KNOWN", "XREF", "PRED", "ORTH", "PSEUDO", "LOWEVIDENCE", "NOIDEA"],
    },
    is_nullable => 0,
  },
  "source_release",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "download",
  {
    data_type => "enum",
    default_value => "Y",
    extra => { list => ["Y", "N"] },
    is_nullable => 1,
  },
  "ordered",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "priority",
  {
    data_type => "integer",
    default_value => 1,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
  "priority_description",
  { data_type => "varchar", default_value => "", is_nullable => 1, size => 40 },
);

=head1 PRIMARY KEY

=over 4

=item * L</source_id>

=back

=cut

__PACKAGE__->set_primary_key("source_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:B7yvA5bB9elVK5YKqevV5A
__PACKAGE__->has_many('xrefs', 'Xref::Schema::Result::Xref', 'source_id');
1;
