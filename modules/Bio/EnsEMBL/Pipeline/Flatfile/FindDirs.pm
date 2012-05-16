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

Bio::EnsEMBL::Pipeline::Flatfile::FindDirs

=head1 DESCRIPTION

Finds all directories under the given species directory. This is used to
flow any further processing only dependent on the directory

Allowed parameters are:

=over 8

=item species - The species to work with

=back

=cut

package Bio::EnsEMBL::Pipeline::FASTA::FindDirs;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::FindDirs Bio::EnsEMBL::Pipeline::Flatfile::Base/;

use File::Spec;

sub fetch_input {
  my ($self) = @_;
  $self->throw("No 'species' parameter specified") unless $self->param('species');
  $self->param('path', $self->data_path());
  $self->SUPER::fetch_input();
  return;
}

1;
