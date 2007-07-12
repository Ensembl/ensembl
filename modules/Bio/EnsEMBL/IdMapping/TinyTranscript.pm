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
  $self->[2] = shift if (@_);
  return $self->[2];
}


sub end {
  my $self = shift;
  $self->[3] = shift if (@_);
  return $self->[3];
}


sub strand {
  my $self = shift;
  $self->[4] = shift if (@_);
  return $self->[4];
}


sub length {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


sub add_Translation {
  my $self = shift;
  my $tl = shift;

  unless ($tl && $tl->isa('Bio::EnsEMBL::IdMapping::TinyTranslation')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranslation.');
  }

  $self->[6] = $tl;
}


sub add_Exon {
  my $self = shift;
  my $exon = shift;

  unless ($exon && $exon->isa('Bio::EnsEMBL::IdMapping::TinyExon')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyExon.');
  }

  push @{ $self->[7] }, $exon;
}


sub get_all_Exons {
  return $_[0]->[7] || [];
}


1;

