package Bio::EnsEMBL::IdMapping::TinyGene;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub start {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


sub end {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


sub strand {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


sub seq_region_name {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


sub biotype {
  my $self = shift;
  $self->[9] = shift if (@_);
  return $self->[9];
}


sub status {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
}


sub logic_name {
  my $self = shift;
  $self->[11] = shift if (@_);
  return $self->[11];
}


sub is_known {
  my $self = shift;
  $self->[12] = shift if (@_);
  return $self->[12];
}


sub add_Transcript {
  my $self = shift;
  my $tr = shift;

  unless ($tr && $tr->isa('Bio::EnsEMBL::IdMapping::TinyTranscript')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranscript.');
  }

  push @{ $self->[13] }, $tr;
}


sub get_all_Transcripts {
  return $_[0]->[13] || [];
}


sub length {
  my $self = shift;
  return ($self->end - $self->start + 1);
}


1;

