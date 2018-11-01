use utf8;
package Xref::Schema::Result::AltAllele;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::AltAllele

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<alt_allele>

=cut

__PACKAGE__->table("alt_allele");

=head1 ACCESSORS

=head2 alt_allele_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 gene_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 is_reference

  data_type: 'integer'
  default_value: 0
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "alt_allele_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "gene_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "is_reference",
  {
    data_type => "integer",
    default_value => 0,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<allele_idx>

=over 4

=item * L</alt_allele_id>

=item * L</gene_id>

=back

=cut

__PACKAGE__->add_unique_constraint("allele_idx", ["alt_allele_id", "gene_id"]);

=head2 C<gene_idx>

=over 4

=item * L</gene_id>

=back

=cut

# __PACKAGE__->add_unique_constraint("gene_idx", ["gene_id"]);


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:FwgD3o32xb489GhixaUXEw

__PACKAGE__->has_one('gene', 'Xref::Schema::Result::GeneStableId', { 'foreign.internal_id' => 'self.gene_id'});

1;
