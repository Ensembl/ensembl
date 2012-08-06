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

Bio::EnsEMBL::Pipeline::Flatfile::DumpTypeFactory

=head1 DESCRIPTION

Small extension of the job factory to do default type submission as otherwise
we get every type of file being produced.

Allowed parameters are:

=over 8

=item types - The types to use; defaults to embl and genbank

=back

=cut

package Bio::EnsEMBL::Pipeline::Flatfile::DumpTypeFactory;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Hive::RunnableDB::JobFactory/;

sub param_defaults {
  my ($self) = @_;
  return {
    column_names => ['type'],
    default_types => [qw/embl genbank/],
  };
}

sub fetch_input {
  my ($self) = @_;
  my $user_types = $self->param('types');
  if($user_types && @{$user_types}) {
    $self->param('inputlist', $user_types);
  }
  else {
    $self->param('inputlist', $self->param('default_types'));
  }
  return;
}

1;