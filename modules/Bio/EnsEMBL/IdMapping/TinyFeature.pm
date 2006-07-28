package Bio::EnsEMBL::IdMapping::TinyFeature;

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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub new_fast {
  my $class = shift;
  my $array_ref = shift;
  return bless $array_ref, $class;
}


sub id {
  my $self = shift;
  $self->[1] = shift if (@_);
  return $self->[1];
}


sub stable_id {
  my $self = shift;
  $self->[2] = shift if (@_);
  return $self->[2];
}


sub add_Transcript {
  my $self = shift;
  my $tr = shift;

  unless ($self->[0] eq 'g' and $tr->[0] eq 'tr') {
    throw('You can only add transcripts to a gene.');
  }

  push @{ $self->[9] }, $tr;
}


sub get_all_Transcripts {
  my $self = shift;

  throw('Only genes have transcripts.') unless ($self->[0] eq 'g'); 
  
  return $self->[9] || [];
}


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

  unless ($self->[0] eq 'tr' and $exon->[0] eq 'e') {
    throw('You can only add exons to a transcript.');
  }

  push @{ $self->[11] }, $exon;
}


sub get_all_Exons {
  my $self = shift;

  throw('Only transcripts have exons.') unless ($self->[0] eq 'tr'); 
  
  return $self->[11] || [];
}


1;

