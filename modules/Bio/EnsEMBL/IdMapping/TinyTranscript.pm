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

use Bio::EnsEMBL::IDMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IDMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub add_Translation {
  my $self = shift;
  my $tl = shift;

  unless ($self->[0] eq 'tr' and $tl->[0] eq 'tl') {
    throw('You can only add a translation to a transcript.');
  }

  $self->[10] = $tl;
}


sub add_Exon {
  my $self = shift;
  my $exon = shift;

  unless ($exon && $exon->isa('Bio::EnsEMBL::IdMapping::TinyExon')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyExon.');
  }

  push @{ $self->[11] }, $exon;
}


sub get_all_Exons {
  return $_->[11] || [];
}


1;

