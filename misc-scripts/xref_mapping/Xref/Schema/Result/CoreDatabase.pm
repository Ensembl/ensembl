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

package Xref::Schema::Result::CoreDatabase;

=head1 NAME

Xref::Schema::Result::CoreDatabase

=cut

use strict;
use warnings;
use utf8;

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


1;
