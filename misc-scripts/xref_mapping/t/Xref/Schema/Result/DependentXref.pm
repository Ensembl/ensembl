use utf8;
package Xref::Schema::Result::DependentXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::DependentXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<dependent_xref>

=cut

__PACKAGE__->table("dependent_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 master_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 dependent_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 linkage_annotation

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 linkage_source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=cut

__PACKAGE__->add_columns(
  'dependency_id',
  { data_type => 'integer', extra => { unsigned => 1 }, is_nullable => 0, is_auto_increment => 1 },
  "object_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "master_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "dependent_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "linkage_annotation",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "linkage_source_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:1BZDfZzGW5D2eb8uEkfBtA

__PACKAGE__->set_primary_key('dependency_id');


__PACKAGE__->has_one('object_xref', 'Xref::Schema::Result::ObjectXref', 'object_xref_id');
__PACKAGE__->has_one('dependent_xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.dependent_xref_id'} );
__PACKAGE__->has_one('master_xref', 'Xref::Schema::Result::Xref', { 'foreign.xref_id' => 'self.master_xref_id'} );
__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', { 'foreign.source_id' => 'self.linkage_source_id' } );

1;
