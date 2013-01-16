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

Bio::EnsEMBL::Pipeline::Flatfile::CheckFlatfile

=head1 DESCRIPTION

Takes in a file and passes it through BioPerl's SeqIO parser code. This
is just a smoke test to ensure the files are well formatted.

Allowed parameters are:

=over 8

=item file - The file to parse

=item type - Passed into SeqIO; the format to parse

=back

=cut

package Bio::EnsEMBL::Pipeline::Flatfile::CheckFlatfile;

use strict;
use warnings;

use Bio::SeqIO;

use base qw/Bio::EnsEMBL::Pipeline::Flatfile::Base/;

sub fetch_input {
  my ($self) = @_;
  $self->throw("No 'file' parameter specified") unless $self->param('file');
  $self->throw("No 'type' parameter specified") unless $self->param('type');
  return;
}

sub run {
  my ($self) = @_;
  my $fh = $self->get_fh();
  my $type = $self->param('type');
  my $stream = Bio::SeqIO->new(-FH => $fh, -FORMAT => $type);
  my $count = 0;
  while ( (my $seq = $stream->next_seq()) ) {
    $self->fine("Found the record %s", $seq->accession());
    $count++;
  }
  $self->info("Processed %d record(s)", $count);
  close $fh;
  return;
}

sub get_fh {
  my ($self) = @_;
  my $file = $self->param('file');
  $self->throw("Cannot find file $file") unless -f $file;
  my $fh;
  if($file =~ /\.gz$/) {
    open $fh, '-|', 'gzip -c -d '.$file or die "Cannot open $file for gunzip: $!";
  }
  else {
    open $fh, '<', $file or die "Cannot open file $file: $!";
  }
  return $fh;
}

1;