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

Bio::EnsEMBL::SplicingTranscriptPair - Object representing an alternative splicing transcript pair

=head1 SYNOPSIS

  my $ase = Bio::EnsEMBL::SplicingTranscriptPair->new(
    -START  => 123,
    -END    => 1045,
    -TRANSCRIPT_ID_1 => $tran1->dbID,
    -TRANSCRIPT_ID_2 => %tran2->dbID
  );

=head1 DESCRIPTION

A representation of an Alternative Splicing Transcrript Pair within the Ensembl system.

=head1 METHODS

=cut

package Bio::EnsEMBL::SplicingTranscriptPair;

use strict;

use POSIX;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);




sub transcript_id_1{
  my $self = shift;
  $self->{'transcript_id_1'} = shift if (@_);

  if (defined $self->{'transcript_id_1'}) {
    return $self->{'transcript_id_1'};
  }

  return undef;
}

sub transcript_id_2{
  my $self = shift;
  $self->{'transcript_id_2'} = shift if (@_);

  if (defined $self->{'transcript_id_2'}) {
    return $self->{'transcript_id_2'};
  }

  return undef;
}


sub start{
  my $self = shift;
  $self->{'start'} = shift if (@_);

  if (defined $self->{'start'}) {
    return $self->{'start'};
  }

  return undef;
}

sub end{
  my $self = shift;
  $self->{'end'} = shift if (@_);

  if (defined $self->{'end'}) {
    return $self->{'end'};
  }

  return undef;
}





1;
