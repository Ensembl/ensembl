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

package Xref::Schema::Result::GeneTranscriptTranslation;

=head1 NAME

Xref::Schema::Result::GeneTranscriptTranslation

=cut

use strict;
use warnings;
use utf8;

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

1;
