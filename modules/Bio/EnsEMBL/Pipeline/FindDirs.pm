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

Bio::EnsEMBL::Pipeline::FindDirs

=head1 DESCRIPTION

Finds all directories under the given path.

Allowed parameters are:

=over 8

=item path - The path to search

=back

=cut

package Bio::EnsEMBL::Pipeline::FindDirs;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Hive::RunnableDB::JobFactory/;

use File::Spec;

sub fetch_input {
  my ($self) = @_;
  $self->throw("No 'path' parameter specified") unless $self->param('path');
  my $dirs = $self->dirs();
  $self->param('inputlist', $dirs);
  return;
}

sub dirs {
  my ($self) = @_;
  
  my @dirs;
  
  my $dir = $self->param('path');
  $self->info('Searching directory %s', $dir);

  opendir(my $dh, $dir) or die "Cannot open directory $dir";
  my @files = sort { $a cmp $b } readdir($dh);
  closedir($dh) or die "Cannot close directory $dir";

  foreach my $file (@files) {
    next if $file =~ /^\./;         #hidden file or up/current dir
    my $path = File::Spec->catdir($dir, $file);
    if(-d $path) {
      $self->fine('Adding %s to the list of found dirs', $path);
      push(@dirs, $path);
    }
  }
  
  return \@dirs;
}

1;
