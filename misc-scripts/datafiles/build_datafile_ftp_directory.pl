#!/usr/bin/env perl

package Script;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use File::Find;
use File::Spec;
use File::Path qw/mkpath/;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;

my $rcsid = '$Revision$';
our ($VERSION) = $rcsid =~ /(\d+\.\d+)/;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->check();
  $self->setup();
  $self->process();
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {
    port => 3306
  };
  GetOptions(
    $opts, qw/
      release|version=i
      dry
      host|hostname|h=s
      port|P=i
      username|user|u=s
      password|pass|p=s
      datafile_dir=s
      ftp_dir=s
      verbose 
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub opts {
  my ($self) = @_;
  return $self->{'opts'};
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params = qw/host username datafile_dir ftp_dir release/;

  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message => "-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  return;
}

sub setup {
  my ($self) = @_;
  my $o = $self->opts();
  
  $self->v('Using the database server %s@%s:%d', map { $o->{$_} } qw/username host port/);
  
  ##SETTING UP REGISTRY
  my %args = (
    -HOST => $o->{host}, -PORT => $o->{port}, 
    -USER => $o->{username}, -DB_VERSION => $o->{release},
  ); 
  $args{-PASS} = $o->{password} if $o->{password};
#  $args{-VERBOSE} = 1 if $o->{verbose};
  my $loaded = Bio::EnsEMBL::Registry->load_registry_from_db(%args);
  $self->v('Loaded %d DBAdaptor(s)', $loaded);
  
  return;
}

sub process {
  my ($self) = @_;
  my $dbas = $self->_get_core_like_dbs();
  foreach my $dba (@$dbas) {
    $self->_process_dba($dba);
  }
  return;
}

sub _process_dba {
  my ($self, $dba) = @_;
  $self->v('Working with species %s', $dba->species());
  my $datafiles = $dba->get_DataFileAdaptor()->fetch_all();
  foreach my $df (@{$datafiles}) {
    $self->_process_datafile($df);
  }
  return;
}

sub _process_datafile {
  my ($self, $datafile) = @_;
  return if $datafile->absolute();
  my $target_dir = $self->_target_root($datafile);
  if(! -d $target_dir) {
    if($self->opts->{dry}) {
      $self->v("Would have created directory '%s'", $target_dir);
    }
    else {
      $self->v("Creating directory '%s'", $target_dir);
      mkpath($target_dir) or die "Cannot create the directory $target_dir: $!";
    }
  }
  my $files = $self->_files($datafile);
  foreach my $filepath (@{$files}) {
    my ($file_volume, $file_dir, $name) = File::Spec->splitpath($filepath);
    my $target = File::Spec->catfile($target_dir, $name);
    if($self->opts()->{dry}) {
      $self->v("Would have linked '%s' -> '%s'", $filepath, $target);
    }
    else {
      if(-e $target) {
        if(-f $target) {
          my $id = $datafile->dbID();
          die "Cannot unlink $target as it is a file and not a symbolic link. Datafile ID was $id";
        }
        elsif(-h $target) {
          unlink $target;
        }
      }
      $self->v('Linking %s -> %s', $filepath, $target);
      symlink($filepath, $target) or die "Cannot symbolically link $filepath to $target: $!";
    }
  }
  return;
}

# Expected path: /root/RELEASE/DATATYPE/TYPE/SPECIES/files
# e.g. pub/66/bam/genebuild/pan_trogladytes/chimp_1.bam
sub _target_root {
  my ($self, $datafile) = @_;
  my $base = $self->opts()->{ftp_dir};
  my $file_type = $self->_datafile_to_type($datafile); 
  my $ftp_type = $self->_dba_to_ftp_type($datafile->adaptor()->db());
  my $species = $datafile->adaptor()->db()->get_MetaContainer()->get_production_name();
  return File::Spec->catdir($base, $file_type, $ftp_type, $species);
}

sub _datafile_to_type {
  my ($self, $datafile) = @_;
  return lc($datafile->file_type());
}

sub _dba_to_ftp_type {
  my ($self, $dba) = @_;
  my $group = $dba->group();
  my $type = {
    core => 'genebuild',
    rnaseq => 'genebuild',
    otherfeatures => 'genebuild',
    
    variation => 'variation',
    
    funcgen => 'funcgen',
  }->{$group};
  die "No way to convert from $group to a type" unless $type;
  return $type;
}

sub _get_core_like_dbs {
  my ($self) = @_;
  my $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors();
  my @final_dbas;
  while(my $dba = shift @{$dbas}) {
    next if $dba->species() eq 'multi';
    next if lc($dba->species()) eq 'ancestral sequences';
    
    my $type = $dba->get_MetaContainer()->single_value_by_key('schema_type');
    $dba->dbc()->disconnect_if_idle();
    next unless $type;
    push(@final_dbas, $dba) if $type eq 'core';
  }
  $self->v('Found %d core like database(s)', scalar(@final_dbas));
  return \@final_dbas;
}

sub _files {
  my ($self, $datafile) = @_;
  my $source_file = $datafile->path($self->opts()->{datafile_dir});
  my ($volume, $dir, $name) = File::Spec->splitpath($source_file);
  my $regex = qr/^$name.*/;
  my @files;
  find(sub {
    push(@files, $File::Find::name) if $_ =~ $regex;
  }, $dir);
  return \@files;
}

sub v {
  my ($self, $msg, @params) = @_;
  return unless $self->opts()->{verbose};
  printf(STDERR $msg."\n", @params);
  return;
}

Script->run();