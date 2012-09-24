=pod

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

=head1 NAME

Bio::EnsEMBL::Pipeline::FASTA::CreatePrimaryAssembly

=head1 DESCRIPTION

Scans the given file and attempts to create a primary assembly file which
is a file containing only those regions we believe to be part of the 
core assembly e.g. in Human this is 1-22, X, Y & MT

Allowed parameters are:

=over 8

=item species - Required to indicate which species we are working with

=item file - The file to filter

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::CreatePrimaryAssembly;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::FASTA::Base/;
use Bio::EnsEMBL::Utils::IO qw/gz_work_with_file/;

sub fetch_input {
  my ($self) = @_;
  foreach my $key (qw/species file/) {
    $self->throw("Cannot find the required parameter $key") unless $self->param($key);
  }
  return;
}

# Creates the file only if required i.e. has non-reference toplevel sequences
sub run {
  my ($self) = @_;
  if($self->has_non_refs()) {
    my $source = $self->param('file');
    my $target = $self->target_file();
    $self->info('Decompressing to %s', $target);
    gz_work_with_file($target, 'w', sub {
      my ($trg_fh) = @_;
      $self->filter_fasta_for_nonref($source, $trg_fh);
      return;      
    });
  }
  $self->cleanup_DBAdaptor();
  return;
}

sub target_file {
  my ($self) = @_;
  # File name format looks like:
  # <species>.<assembly>.<release>.<sequence type>.<id type>.<id>.fa.gz
  # e.g. Homo_sapiens.GRCh37.64.dna_rm.toplevel.fa.gz -> Homo_sapiens.GRCh37.64.dna_rm.primary_assembly.fa.gz
  my $file = $self->param('file');
  $file =~ s/\.toplevel\./.primary_assembly./;
  return $file;
}

1;
