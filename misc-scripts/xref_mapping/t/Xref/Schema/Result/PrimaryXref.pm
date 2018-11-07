use utf8;
package Xref::Schema::Result::PrimaryXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::PrimaryXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<primary_xref>

=cut

__PACKAGE__->table("primary_xref");

=head1 ACCESSORS

=head2 xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 sequence

  accessor: undef
  data_type: 'mediumtext'
  is_nullable: 1

=head2 sequence_type

  data_type: 'enum'
  extra: {list => ["dna","peptide"]}
  is_nullable: 1

=head2 status

  data_type: 'enum'
  extra: {list => ["experimental","predicted"]}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "sequence",
  { accessor => undef, data_type => "mediumtext", is_nullable => 1 },
  "sequence_type",
  {
    data_type => "enum",
    extra => { list => ["dna", "peptide"] },
    is_nullable => 1,
  },
  "status",
  {
    data_type => "enum",
    extra => { list => ["experimental", "predicted"] },
    is_nullable => 1,
  },
);

=head1 PRIMARY KEY

=over 4

=item * L</xref_id>

=back

=cut

__PACKAGE__->set_primary_key("xref_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:j8WfubBbnx0pSFcQcAak6w

__PACKAGE__->add_column('+sequence' => {accessor => 'raw_seq'});
__PACKAGE__->has_one('xref', 'Xref::Schema::Result::Xref', 'xref_id' );
1;
