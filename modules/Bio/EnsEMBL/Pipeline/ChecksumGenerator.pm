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

Bio::EnsEMBL::Pipeline::ChecksumGenerator

=head1 DESCRIPTION

Creates a CHECKSUMS file in the given directory which is produced from running
the sum command over every file in the directory. This excludes the CHECKSUMS
file, parent directory or any hidden files.

Allowed parameters are:

=over 8

=item dir - The directory to generate checksums for

=item gzip - If the resulting file should be gzipped. Defaults to false

=back

=cut

package Bio::EnsEMBL::Pipeline::ChecksumGenerator;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Base/;

use File::Spec;
use Bio::EnsEMBL::Utils::IO qw/work_with_file gz_work_with_file/;

sub param_defaults {
  my ($self) = @_;
  return {
    gzip => 0
  };
}

sub fetch_input {
  my ($self) = @_;
  my $dir = $self->param('dir');
  $self->throw("No 'dir' parameter specified") unless $dir;
  $self->throw("Dir $dir does not exist") unless -d $dir;
  return;
}

sub run {
  my ($self) = @_;
  my @checksums;
  
  my $dir = $self->param('dir');
  $self->info('Checksumming directory %s', $dir);

  opendir(my $dh, $dir) or die "Cannot open directory $dir";
  my @files = sort { $a cmp $b } readdir($dh);
  closedir($dh) or die "Cannot close directory $dir";

  foreach my $file (@files) {
    next if $file =~ /^\./;         #hidden file or up/current dir
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($dir, $file);
    my $checksum = $self->checksum($path);
    push(@checksums, [$checksum, $file])
  }
  
  $self->param('checksums', \@checksums);
  return;
}

sub write_output {
  my ($self) = @_;
  my $dir = $self->param('dir');
  my $checksum = File::Spec->catfile($dir, 'CHECKSUMS');
  $checksum .= '.gz' if $self->param('gzip');
  if(-f $checksum) {
    $self->info('Checksum file already exists. Removing');
    unlink $checksum;
  }
  
  my @checksums = @{$self->param('checksums')};
  
  return unless @checksums;
  
  my $writer = sub {
    my ($fh) = @_;
    foreach my $entry (@checksums) {
      my $line = join(qq{\t}, @{$entry});
      print $fh $line;
      print $fh "\n"; 
    }
    return;
  };
  my @params = ($checksum, 'w', $writer);
  
  
  if($self->param('gzip')) {
    gz_work_with_file(@params);
  } 
  else {
    work_with_file(@params);
  } 
  
  $self->permissions($checksum);
  return;
}

sub checksum {
  my ($self, $path) = @_;
  my $checksum = `sum $path`;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

sub permissions {
  my ($self, $file) = @_;
  my $mode = 0666;
  chmod($mode, $file) or $self->throw("Cannot perform the chmod to mode $mode for file $file");
  return;
}

1;
