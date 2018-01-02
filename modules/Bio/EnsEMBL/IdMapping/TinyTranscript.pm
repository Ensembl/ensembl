=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::IdMapping::TinyTranscript - lightweight transcript object

=head1 SYNOPSIS

  # fetch a transcript from the db and create a lightweight
  # transcript object from it
  my $tr = $transcript_adaptor->fetch_by_stable_id('ENST000345437');
  my $lightweight_tr =
    Bio::EnsEMBL::IdMapping::TinyTranscript->new_fast( [
      $tr->dbID,          $tr->stable_id,
      $tr->version,       $tr->created_date,
      $tr->modified_date, $tr->start,
      $tr->end,           $tr->strand,
      $tr->length,        md5_hex( $tr->spliced_seq ),
    ] );

=head1 DESCRIPTION

This is a lightweight transcript object for the stable Id mapping. See
the documentation in TinyFeature for general considerations about its
design.

=head1 METHODS

  start
  end
  strand
  length
  seq_md5_sum
  add_Translation
  translation
  add_Exon
  get_all_Exons

=cut

package Bio::EnsEMBL::IdMapping::TinyTranscript;

# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  start
#  6  end
#  7  strand
#  8  length
#  9  seq_md5_sum
# 10  translation
# 11  [exons]
# 12  biotype
# 13  slice


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 start

  Arg[1]      : (optional) Int - the transcript's start coordinate
  Description : Getter/setter for the transcript's start coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub start {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


=head2 end

  Arg[1]      : (optional) Int - the transcript's end coordinate
  Description : Getter/setter for the transcript's end coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub end {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


=head2 strand

  Arg[1]      : (optional) Int - the transcript's strand
  Description : Getter/setter for the transcript's strand.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub strand {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


=head2 length

  Arg[1]      : (optional) Int - the transcript's length
  Description : Getter/setter for the transcript's length. Note that this is
                *not* the distance between start and end, but rather the sum of
                the lengths of all exons.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub length {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


=head2 seq_md5_sum

  Arg[1]      : (optional) String - the md5 digest of the transcript's sequence
  Description : Getter/setter for the md5 digest of the transcript's sequence.
                Note that when used as a setter, you are expected to pass a
                digest, not the raw sequence (i.e. the digest is not created for
                you).
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq_md5_sum {
  my $self = shift;
  $self->[9] = shift if (@_);
  return $self->[9];
}



=head2 add_Translation

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyTranslation $tl - the translation
                to add
  Example     : $tiny_transcript->add_Translation($tiny_translation);
  Description : Adds a translation to this transcript.
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_Translation {
  my $self = shift;
  my $tl = shift;

  unless ($tl && $tl->isa('Bio::EnsEMBL::IdMapping::TinyTranslation')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranslation.');
  }

  $self->[10] = $tl;
}


=head2  translation

  Description : Getter for the transcript's translation.
  Return type : Bio::EnsEMBL::IdMapping::TinyTranslation
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub translation {
  return $_[0]->[10];
}


=head2 add_Exon

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyExon $exon - the exon to add
  Example     : $tiny_transcript->add_Exon($tiny_exon);
  Description : Adds an exon to this transcript.
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_Exon {
  my $self = shift;
  my $exon = shift;

  unless ($exon && $exon->isa('Bio::EnsEMBL::IdMapping::TinyExon')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyExon.');
  }

  push @{ $self->[11] }, $exon;
}


=head2 get_all_Exons

  Example     : foreach my $exon (@{ $tiny_transcript->get_all_Exons }) {
                  # do something with exon
                }
  Description : Returns all exons attached to that transcript.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::TinyExon objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_Exons {
  return $_[0]->[11] || [];
}

=head2 biotype

  Arg[1]      : (optional) String - the gene's biotype
  Description : Getter/setter for the gene's biotype.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub biotype {
  my $self = shift;
  $self->[12] = shift if (@_);
  return $self->[12];
}

=head2 seq_region_name

  Arg[1]      : (optional) String - seq_region name
  Description : Getter/setter for the seq_region name of the slice the gene is
                on.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq_region_name {
  my $self = shift;
  $self->[13] = shift if (@_);
  return $self->[13];
}


1;

