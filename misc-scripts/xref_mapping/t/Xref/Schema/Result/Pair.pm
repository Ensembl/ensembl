use utf8;
package Xref::Schema::Result::Pair;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Pair

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<pairs>

=cut

__PACKAGE__->table("pairs");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 accession1

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 accession2

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=cut

__PACKAGE__->add_columns(
  "pair_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 , is_auto_increment => 1},
  "source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "accession1",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "accession2",
  { data_type => "varchar", is_nullable => 0, size => 255 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:Eotmov9SiCBJmcX5nGOEEA

__PACKAGE__->set_primary_key('pair_id');
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );

1;
