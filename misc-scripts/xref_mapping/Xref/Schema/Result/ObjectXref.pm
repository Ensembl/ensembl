=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Xref::Schema::Result::ObjectXref;

=head1 NAME

Xref::Schema::Result::ObjectXref

=cut

use strict;
use warnings;
use utf8;

use base 'DBIx::Class::Core';

=head1 TABLE: C<object_xref>

=cut

__PACKAGE__->table("object_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 ensembl_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 ensembl_object_type

  data_type: 'enum'
  extra: {list => ["RawContig","Transcript","Gene","Translation"]}
  is_nullable: 0

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 linkage_annotation

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 linkage_type

  data_type: 'enum'
  extra: {list => ["PROJECTION","MISC","DEPENDENT","DIRECT","SEQUENCE_MATCH","INFERRED_PAIR","PROBE","UNMAPPED","COORDINATE_OVERLAP","CHECKSUM"]}
  is_nullable: 1

=head2 ox_status

  data_type: 'enum'
  default_value: 'DUMP_OUT'
  extra: {list => ["DUMP_OUT","FAILED_PRIORITY","FAILED_CUTOFF","NO_DISPLAY","MULTI_DELETE"]}
  is_nullable: 0

=head2 unused_priority

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 master_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "object_xref_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "ensembl_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "ensembl_object_type",
  {
    data_type => "enum",
    extra => { list => ["RawContig", "Transcript", "Gene", "Translation"] },
    is_nullable => 0,
  },
  "xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "linkage_annotation",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "linkage_type",
  {
    data_type => "enum",
    extra => {
      list => [
        "PROJECTION",
        "MISC",
        "DEPENDENT",
        "DIRECT",
        "SEQUENCE_MATCH",
        "INFERRED_PAIR",
        "PROBE",
        "UNMAPPED",
        "COORDINATE_OVERLAP",
        "CHECKSUM",
      ],
    },
    is_nullable => 1,
  },
  "ox_status",
  {
    data_type => "enum",
    default_value => "DUMP_OUT",
    extra => {
      list => [
        "DUMP_OUT",
        "FAILED_PRIORITY",
        "FAILED_CUTOFF",
        "NO_DISPLAY",
        "MULTI_DELETE",
      ],
    },
    is_nullable => 0,
  },
  "unused_priority",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "master_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</object_xref_id>

=back

=cut

__PACKAGE__->set_primary_key("object_xref_id");

=head1 UNIQUE CONSTRAINTS

=head2 C<ensembl_object_type>

=over 4

=item * L</ensembl_object_type>

=item * L</ensembl_id>

=item * L</xref_id>

=item * L</ox_status>

=item * L</master_xref_id>

=back

=cut

__PACKAGE__->add_unique_constraint(
  "ensembl_object_type",
  [
    "ensembl_object_type",
    "ensembl_id",
    "xref_id",
    "ox_status",
    "master_xref_id",
  ],
);

__PACKAGE__->has_one('source', 'Xref::Schema::Result::Source', 'source_id' );
# Cannot express the gene/transcript/translation relationship here due to non-normalised format
__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', 'xref_id' );

1;
