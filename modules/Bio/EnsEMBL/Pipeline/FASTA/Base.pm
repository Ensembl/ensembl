package Bio::EnsEMBL::Pipeline::FASTA::Base;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Hive::Process/;

use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/check_ref/;
use File::Find;
use File::Spec;
use File::Path qw/mkpath/;
use POSIX qw/strftime/;

# Takes in a key, checks if the current $self->param() was an empty array
# and replaces it with the value from $self->param_defaults()
sub reset_empty_array_param {
  my ($self, $key) = @_;
  my $param_defaults = $self->param_defaults();
  my $current = $self->param($key); 
  my $replacement = $self->param_defaults()->{$key};
  if(check_ref($current, 'ARRAY') && check_ref($replacement, 'ARRAY')) {
    if(! @{$current}) {
      $self->fine('Restting param %s because the given array was empty', $key);
      $self->param($key, $replacement);
    }
  }
  return;
}

sub get_Slices {
  my ($self, $type) = @_;
  my $sa = $self->get_DBAdaptor($type)->get_SliceAdaptor();
  return [ sort { $a->length() <=> $b->length() }  @{$sa->fetch_all('toplevel', undef, 1, undef, undef)} ];
}

# Registry is loaded by Hive (see beekeeper_extra_cmdline_options() in conf)
sub get_DBAdaptor {
  my ($self, $type) = @_;
  my $species = $self->param('species');
  $type ||= 'core';
  return Bio::EnsEMBL::Registry->get_DBAdaptor($species, $type);
}

sub cleanup_DBAdaptor {
  my ($self, $type) = @_;
  my $dba = $self->get_DBAdaptor($type);
  $dba->clear_caches;
  $dba->dbc->disconnect_if_idle;
  return;
}

sub get_dir {
  my ($self, @extras) = @_;
  my $base_dir = $self->param('base_path');
  my $dir = File::Spec->catdir($base_dir, @extras);
  mkpath($dir);
  return $dir;
}

sub fasta_path {
  my ( $self, @extras ) = @_;
  return $self->get_dir('fasta', $self->param('species'), @extras);
}

sub web_name {
  my ($self) = @_;
#  my $mc = $self->get_DBAdaptor()->get_MetaContainer();
#  my $name = $mc->single_value_by_key('species.url'); # change back
  my $name = ucfirst($self->production_name());
  return $name;
}

sub assembly {
  my ($self) = @_;
  my $dba = $self->get_DBAdaptor();
  return $dba->get_CoordSystemAdaptor()->fetch_all()->[0]->version();
}

sub production_name {
  my ($self, $name) = @_;
  my $dba;
  if($name) {
    $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($name, 'core');
  }
  else {
    $dba = $self->get_DBAdaptor();
  }
  my $mc = $dba->get_MetaContainer();
  my $prod = $mc->get_production_name();
  $dba->dbc()->disconnect_if_idle();
  return $prod;
}

sub old_path {
  my ($self, $species) = @_;
  my $base = $self->param('ftp_dir');
  my $prod = $self->production_name($species);
  my $version = $self->param('previous_version');
  my $dir = File::Spec->catdir($base, "release-$version", 'fasta', $prod, 'dna');
}

# Closes file handle, and deletes the file stub if it contains no data
# Returns success type

sub tidy_file_handle {
  my ($self, $fh, $path) = @_;
  if ($fh->tell() == 0) {
    $fh->close;
    unlink($path);
    return 1;
  }
  close($fh);
  return 0;
}

sub info {
  my ($self, $msg, @params) = @_;
  if ($self->debug() > 1) {
    my $formatted_msg;
    if(scalar(@params)) {
      $formatted_msg = sprintf($msg, @params);
    } 
    else {
      $formatted_msg = $msg;
    }
    printf STDERR "INFO: %s %s\n", strftime('%c',localtime()), $formatted_msg;
  }
  return
}

sub fine {
  my ($self, $msg, @params) = @_;
  if ($self->debug() > 2) {
    my $formatted_msg;
    if(scalar(@params)) {
      $formatted_msg = sprintf($msg, @params);
    } 
    else {
      $formatted_msg = $msg;
    }
    printf STDERR "FINE: %s %s\n", strftime('%c',localtime()), $formatted_msg;
  }
  return
}

sub find_files {
  my ($self, $dir, $boolean_callback) = @_;
  $self->throw("Cannot find path $dir") unless -d $dir;
  my @files;
  find(sub {
    my $path = $File::Find::name;
    if($boolean_callback->($_)) {
      push(@files, $path);
    }
  }, $dir);
  return \@files;
}

1;
