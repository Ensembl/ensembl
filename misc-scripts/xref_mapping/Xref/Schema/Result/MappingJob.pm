=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Xref::Schema::Result::MappingJob;

=head1 NAME

Xref::Schema::Result::MappingJob

=cut

use strict;
use warnings;
use utf8;

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

__PACKAGE__->set_primary_key('job_id');

__PACKAGE__->has_one('job', 'Xref::Schema::Result::Mapping', 'job_id' );
1;
