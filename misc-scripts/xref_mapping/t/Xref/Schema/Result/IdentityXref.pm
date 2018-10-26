use utf8;
package Xref::Schema::Result::IdentityXref;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::IdentityXref

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<identity_xref>

=cut

__PACKAGE__->table("identity_xref");

=head1 ACCESSORS

=head2 object_xref_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 query_identity

  data_type: 'integer'
  is_nullable: 1

=head2 target_identity

  data_type: 'integer'
  is_nullable: 1

=head2 hit_start

  data_type: 'integer'
  is_nullable: 1

=head2 hit_end

  data_type: 'integer'
  is_nullable: 1

=head2 translation_start

  data_type: 'integer'
  is_nullable: 1

=head2 translation_end

  data_type: 'integer'
  is_nullable: 1

=head2 cigar_line

  data_type: 'text'
  is_nullable: 1

=head2 score

  data_type: 'double precision'
  is_nullable: 1

=head2 evalue

  data_type: 'double precision'
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "object_xref_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "query_identity",
  { data_type => "integer", is_nullable => 1 },
  "target_identity",
  { data_type => "integer", is_nullable => 1 },
  "hit_start",
  { data_type => "integer", is_nullable => 1 },
  "hit_end",
  { data_type => "integer", is_nullable => 1 },
  "translation_start",
  { data_type => "integer", is_nullable => 1 },
  "translation_end",
  { data_type => "integer", is_nullable => 1 },
  "cigar_line",
  { data_type => "text", is_nullable => 1 },
  "score",
  { data_type => "double precision", is_nullable => 1 },
  "evalue",
  { data_type => "double precision", is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</object_xref_id>

=back

=cut

__PACKAGE__->set_primary_key("object_xref_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:eJCbX5UcdedzsjKcXEhn+Q


__PACKAGE__->has_one('object_xref', 'Xref::Schema::Result::ObjectXref', 'object_xref_id' );
1;
