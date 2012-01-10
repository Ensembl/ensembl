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

Bio::EnsEMBL::IdMapping::TinyGene - lightweight gene object

=head1 SYNOPSIS

 # fetch a gene from the db and create a lightweight gene object from it
  my $gene = $gene_adaptor->fetch_by_stable_id('ENSG000345437');
  my $lightweight_gene = Bio::EnsEMBL::IdMapping::TinyGene->new_fast( [
      $gene->dbID,                   $gene->stable_id,
      $gene->version,                $gene->created_date,
      $gene->modified_date,          $gene->start,
      $gene->end,                    $gene->strand,
      $gene->slice->seq_region_name, $gene->biotype,
      $gene->status,                 $gene->analysis->logic_name,
      ( $gene->is_known ? 1 : 0 ),
  ] );

=head1 DESCRIPTION

This is a lightweight gene object for the stable Id mapping. See the
documentation in TinyFeature for general considerations about its
design.

=head1 METHODS

  start
  end
  strand
  seq_region_name
  biotype
  status
  logic_name
  is_known
  add_Transcript
  get_all_Transcripts
  length

=cut

package Bio::EnsEMBL::IdMapping::TinyGene;

# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  start
#  6  end
#  7  strand
#  8  seq_region_name
#  9  biotype
# 10  status
# 11  logic_name
# 12  is_known
# 13  [transcripts]


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 start

  Arg[1]      : (optional) Int - the gene's start coordinate
  Description : Getter/setter for the gene's start coordinate.
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

  Arg[1]      : (optional) Int - the gene's end coordinate
  Description : Getter/setter for the gene's end coordinate.
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

  Arg[1]      : (optional) Int - the gene's strand
  Description : Getter/setter for the gene's strand.
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
  $self->[8] = shift if (@_);
  return $self->[8];
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
  $self->[9] = shift if (@_);
  return $self->[9];
}


=head2 strand

  Arg[1]      : (optional) String - the gene's status
  Description : Getter/setter for the gene's status.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub status {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
}


=head2 logic_name

  Arg[1]      : (optional) String - the gene's analysis' logic_name
  Description : Getter/setter for the gene's analysis' logic_name.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub logic_name {
  my $self = shift;
  $self->[11] = shift if (@_);
  return $self->[11];
}


=head2 is_known

  Arg[1]      : (optional) Boolean - the gene's "known" status
  Description : Getter/setter for the gene's "known" status.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_known {
  my $self = shift;
  $self->[12] = shift if (@_);
  return $self->[12];
}


=head2 add_Transcript

  Arg[1]      : Bio::EnsEMBL::IdMapping::TinyTranscript $tr - the transcript to
                add
  Example     : $tiny_gene->add_Transcript($tiny_transcript);
  Description : Adds a transcript to a gene.
  Return type : none
  Exceptions  : thrown on wrong or missing argument
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub add_Transcript {
  my $self = shift;
  my $tr = shift;

  unless ($tr && $tr->isa('Bio::EnsEMBL::IdMapping::TinyTranscript')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranscript.');
  }

  push @{ $self->[13] }, $tr;
}


=head2 get_all_Transcripts

  Example     : foreach my $tr (@{ $tiny_gene->get_all_Transcripts }) {
                  # do something with transcript
                }
  Description : Returns all transcripts attached to that gene.
  Return type : Arrayref of Bio::EnsEMBL::IdMapping::TinyTranscript objects
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub get_all_Transcripts {
  return $_[0]->[13] || [];
}


=head2 length

  Description : Returns the gene length (distance between start and end).
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub length {
  my $self = shift;
  return ($self->end - $self->start + 1);
}


1;

