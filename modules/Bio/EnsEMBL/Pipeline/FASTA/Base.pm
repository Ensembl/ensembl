package Bio::EnsEMBL::Pipeline::FASTA::Base;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::Base/;

use File::Spec;

sub fasta_path {
  my ( $self, @extras ) = @_;
  return $self->get_dir('fasta', $self->param('species'), @extras);
}

sub old_path {
  my ($self, $species) = @_;
  my $base = $self->param('ftp_dir');
  my $prod = $self->production_name($species);
  my $release = $self->param('previous_release');
  my $dir = File::Spec->catdir($base, "release-$release", 'fasta', $prod, 'dna');
}

1;
