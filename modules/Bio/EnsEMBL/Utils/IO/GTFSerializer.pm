=pod

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

Bio::EnsEMBL::Utils::IO::GTFSerializer - Transcript to GTF converter

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::IO::GTFSerializer;

my $serializer = Bio::EnsEMBL::Utils::IO::GTFSerializer->new($output_fh);

=head1 DESCRIPTION

Subclass of Serializer that can turn a transcript into a series of lines 
for the GTF format.

=cut

package Bio::EnsEMBL::Utils::IO::GTFSerializer;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

my %strand_conversion = ( '1' => '+', '0' => '?', '-1' => '-');

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Transcript
    Example    : $serializer->print_feature($transcript)
    Description: 
    Returntype : none

=cut

sub print_feature {
  my $self = shift;
  my $feature = shift;

  throw(sprintf "Feature is of type %s. Cannot write non transcripts to GTF", ref($feature))
    unless check_ref($feature, "Bio::EnsEMBL::Transcript");

  my $text_buffer = "";

  #filehandle is inherited
  my $fh = $self->{'filehandle'};
  print $fh $text_buffer;

}

1;
