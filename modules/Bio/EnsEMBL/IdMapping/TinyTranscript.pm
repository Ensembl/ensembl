package Bio::EnsEMBL::IdMapping::TinyTranscript;

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


sub length {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


sub seq_md5_sum {
  my $self = shift;
  $self->[9] = shift if (@_);
  return $self->[9];
}


sub is_known {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
}


sub add_Translation {
  my $self = shift;
  my $tl = shift;

  unless ($tl && $tl->isa('Bio::EnsEMBL::IdMapping::TinyTranslation')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranslation.');
  }

  $self->[11] = $tl;
}


sub translation {
  return $_[0]->[11];
}


sub add_Exon {
  my $self = shift;
  my $exon = shift;

  unless ($exon && $exon->isa('Bio::EnsEMBL::IdMapping::TinyExon')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyExon.');
  }

  push @{ $self->[12] }, $exon;
}


sub get_all_Exons {
  return $_[0]->[12] || [];
}


1;

