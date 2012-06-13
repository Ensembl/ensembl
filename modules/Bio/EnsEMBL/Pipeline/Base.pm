package Bio::EnsEMBL::Pipeline::Base;

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

=head2 get_Slices

	Arg[1]      : String type of DB to use (defaults to core)
	Arg[2]      : Boolean should we filter the slices if it is human
  Example     : my $slices = $self->get_Slices('core', 1);
  Description : Basic get_Slices() method to return all distinct slices
                for a species but also optionally filters for the 
                first portion of Human Y which is a non-informative region
                (composed solely of N's). The code will only filter for 
                GRCh37 forcing the developer to update the test for other 
                regions. 
  Returntype  : ArrayRef[Bio::EnsEMBL::Slice] 
  Exceptions  : Thrown if you are filtering Human but also are not on GRCh37

=cut

sub get_Slices {
  my ($self, $type, $filter_human) = @_;
  my $dba = $self->get_DBAdaptor($type);
  throw "Cannot get a DB adaptor" unless $dba;
  
  my $sa = $dba->get_SliceAdaptor();
  my @slices = @{$sa->fetch_all('toplevel', undef, 1, undef, undef)};
  
  if(!$filter_human) {
    my $production_name = $self->production_name();
    if($production_name eq 'homo_sapiens') {
      my ($cs) = @{$dba->get_CoordSystem()->fetch_all()};
      my $expected = 'GRCh37';
      if($cs->version() ne $expected) {
        throw sprintf(q{Cannot continue as %s's coordinate system %s is not the expected %s }, $production_name, $cs->version(), $expected);
      }
      @slices = grep { $_->seq_region_name() eq 'Y' && $_->end() < 2649521 } @slices;
    }
  }
  
  return [ sort { $a->length() <=> $b->length() }  @slices ];
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

sub web_name {
  my ($self) = @_;
#  my $mc = $self->get_DBAdaptor()->get_MetaContainer();
#  my $name = $mc->single_value_by_key('species.url'); # change back
  my $name = ucfirst($self->production_name());
  return $name;
}

sub scientific_name {
  my ($self) = @_;
  my $dba = $self->get_DBAdaptor();
  my $mc = $dba->get_MetaContainer();
  my $name = $mc->get_scientific_name();
  $dba->dbc()->disconnect_if_idle();
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

# Closes file handle, and deletes the file stub if it contains no data
# Returns success type

sub tidy_file_handle {
  my ($self, $fh, $path) = @_;
  if ($fh->tell() == 0) {
    $fh->close() if $fh->opened();
    unlink($path) if -f $path;
    return 1;
  }
  $fh->close() if $fh->opened();
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

sub unlink_all_files {
  my ($self, $dir) = @_;
  $self->info('Removing files from the directory %s', $dir);
  #Delete anything which is a file & not the current or higher directory
  my $boolean_callback = sub {
    return ( $_[0] =~ /^\.\.?$/) ? 0 : 1;
  };
  my $files = $self->find_files($dir, $boolean_callback);
  foreach my $file (@{$files}) {
    $self->fine('Unlinking %s', $file);
    unlink $file;
  }
  $self->info('Removed %d file(s)', scalar(@{$files}));
  return;
}

1;
