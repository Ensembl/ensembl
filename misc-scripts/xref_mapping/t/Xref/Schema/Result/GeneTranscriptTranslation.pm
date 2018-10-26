use utf8;
package Xref::Schema::Result::GeneTranscriptTranslation;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Xref::Schema::Result::GeneTranscriptTranslation

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<gene_transcript_translation>

=cut

__PACKAGE__->table("gene_transcript_translation");

=head1 ACCESSORS

=head2 gene_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 transcript_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 translation_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 1

=cut

__PACKAGE__->add_columns(
  "gene_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "transcript_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "translation_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 1 },
);

=head1 PRIMARY KEY

=over 4

=item * L</transcript_id>

=back

=cut

__PACKAGE__->set_primary_key("transcript_id");


# Created by DBIx::Class::Schema::Loader v0.07042 @ 2018-10-23 11:58:10
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:6EVVf9WeLQW+LnZ1wks/1g


# You can replace this text with custom code or comments, and it will be preserved on regeneration
1;
