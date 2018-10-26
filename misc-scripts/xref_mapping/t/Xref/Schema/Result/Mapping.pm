use utf8;
package Xref::Schema::Result::Mapping;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::Mapping

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<mapping>

=cut

__PACKAGE__->table("mapping");

=head1 ACCESSORS

=head2 job_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 type

  data_type: 'enum'
  extra: {list => ["dna","peptide","UCSC"]}
  is_nullable: 1

=head2 command_line

  data_type: 'text'
  is_nullable: 1

=head2 percent_query_cutoff

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 percent_target_cutoff

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=head2 method

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 array_size

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "job_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "type",
  {
    data_type => "enum",
    extra => { list => ["dna", "peptide", "UCSC"] },
    is_nullable => 1,
  },
  "command_line",
  { data_type => "text", is_nullable => 1 },
  "percent_query_cutoff",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "percent_target_cutoff",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
  "method",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "array_size",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:xOYxmTgAepF22rgzODMIgw


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
