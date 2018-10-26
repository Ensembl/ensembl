use utf8;
package Xref::Schema::Result::DependentSource;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::DependentSource

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<dependent_source>

=cut

__PACKAGE__->table("dependent_source");

=head1 ACCESSORS

=head2 master_source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 dependent_name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=cut

__PACKAGE__->add_columns(
  "dependent_relation_id",
  { data_type => "integer", extra => { unsigned => 1, auto_increment => 1}, is_nullable => 0},
  "master_source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "dependent_name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
);

=head1 PRIMARY KEY

=over 4

=item * L</master_source_id>

=item * L</dependent_name>

=back

=cut

__PACKAGE__->set_primary_key('dependent_relation_id');
__PACKAGE__->add_unique_constraint("master_source_id", ["dependent_name"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:GpVgVNpaNwS4lYuJVOB7FA

__PACKAGE__->has_one('master_source', 'Xref::Schema::Result::Source', { 'foreign.source_id' => 'self.master_source_id' } );

1;
