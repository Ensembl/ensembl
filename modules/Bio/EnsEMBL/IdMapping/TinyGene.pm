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

use Bio::EnsEMBL::IDMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IDMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub add_Transcript {
  my $self = shift;
  my $tr = shift;

  unless ($tr && $tr->isa('Bio::EnsEMBL::IdMapping::TinyTranscript')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TinyTranscript.');
  }

  push @{ $self->[9] }, $tr;
}


sub get_all_Transcripts {
  return $_->[9] || [];
}


1;

