=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
      ( $tr->is_known ? 1 : 0 ),
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
  is_known
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
# 10  is_known
# 11  translation
# 12  [exons]
# 13  biotype


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


=head2 is_known

  Arg[1]      : (optional) Boolean - the transcript's "known" status
  Description : Getter/setter for the transcript's "known" status.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_known {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
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

  $self->[11] = $tl;
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
  return $_[0]->[11];
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

  push @{ $self->[12] }, $exon;
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
  return $_[0]->[12] || [];
}

sub biotype {
  my $self = shift;
  $self->[13] = shift if (@_);
  return $self->[13];
}


1;

