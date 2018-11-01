use utf8;
package Xref::Schema::Result::MappingJob;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::MappingJob

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<mapping_jobs>

=cut

__PACKAGE__->table("mapping_jobs");

=head1 ACCESSORS

=head2 root_dir

  data_type: 'text'
  is_nullable: 1

=head2 map_file

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 status

  data_type: 'enum'
  extra: {list => ["SUBMITTED","FAILED","SUCCESS"]}
  is_nullable: 1

=head2 out_file

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 err_file

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 array_number

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 job_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 failed_reason

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 object_xref_start

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 object_xref_end

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "root_dir",
  { data_type => "text", is_nullable => 1 },
  "map_file",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "status",
  {
    data_type => "enum",
    extra => { list => ["SUBMITTED", "FAILED", "SUCCESS"] },
    is_nullable => 1,
  },
  "out_file",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "err_file",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "array_number",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "job_id",
  { data_type => "integer", extra => { unsigned => 1 } },
  "failed_reason",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "object_xref_start",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "object_xref_end",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:HbB5THeDV9Bf/Cu0ia3xZA
__PACKAGE__->set_primary_key('job_id');

__PACKAGE__->has_one('job', 'Xref::Schema::Result::Mapping', 'job_id' );
1;
