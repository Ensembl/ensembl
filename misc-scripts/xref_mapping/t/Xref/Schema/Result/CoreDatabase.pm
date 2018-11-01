use utf8;
package Xref::Schema::Result::CoreDatabase;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::CoreDatabase

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<core_database>

=cut

__PACKAGE__->table("core_database");

=head1 ACCESSORS

=head2 port

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 user

  data_type: 'varchar'
  is_nullable: 1
  size: 16

=head2 pass

  data_type: 'varchar'
  is_nullable: 1
  size: 16

=head2 dbname

  data_type: 'varchar'
  is_nullable: 1
  size: 16

=head2 xref_dir

  data_type: 'text'
  is_nullable: 1

=head2 core_dir

  data_type: 'text'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "port",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "user",
  { data_type => "varchar", is_nullable => 1, size => 16 },
  "pass",
  { data_type => "varchar", is_nullable => 1, size => 16 },
  "dbname",
  { data_type => "varchar", is_nullable => 1, size => 16 },
  "xref_dir",
  { data_type => "text", is_nullable => 1 },
  "core_dir",
  { data_type => "text", is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:9JycspWuFjAceCPVuusEAA


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
