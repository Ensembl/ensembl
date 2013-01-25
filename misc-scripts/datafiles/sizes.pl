#!/usr/bin/env perl

package Script;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Fcntl ':mode';
use File::Basename qw/dirname/;
use File::Spec::Functions qw/:ALL/;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use Pod::Usage;
use Test::More;

my $Test = Test::Builder->new();

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
      host|hostname|h=s
      port|P=i
      username|user|u=s
      password|pass|p=s
      datafile_dir|dir=s
      species=s
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

  my @required_params = qw/host username datafile_dir/;

  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message => "-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  foreach my $key (qw/datafile_dir/) {
    my $dir = $o->{$key};
    if(! -d $dir) {
      pod2usage(
        -message => "-${key} given location '${dir}' does not exist",
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
    -USER => $o->{username}
  );
  $args{-DB_VERSION} = $o->{release};
  $args{-PASS} = $o->{password} if $o->{password};
  my $loaded = Bio::EnsEMBL::Registry->load_registry_from_db(%args);
  $self->v('Loaded %d DBAdaptor(s)', $loaded);
  
  return;
}

sub process {
  my ($self) = @_;
  my $dbas = $self->_get_core_like_dbs();
  my $total_size = 0;
  while (my $dba = shift @{$dbas}) {
    my $size = $self->_process_dba($dba);
    $total_size += 0;
    my $size_in_gb = $size / 1024 / 1024 /1024;
    $self->v('Species size is %dGB', $size_in_gb);
  }
  my $total_size_in_gb = $total_size / 1024 / 1024 /1024;
  $self->v('Total size is %dGB', $total_size_in_gb);
  return;
}

sub _process_dba {
  my ($self, $dba) = @_;
  $self->v('Working with species %s', $dba->species());
  my $size = 0;
  my $datafiles = $dba->get_DataFileAdaptor()->fetch_all();
  if(! @{$datafiles}) {
    $self->v("No datafiles found");
  }
  else {
    foreach my $data_file (@{$datafiles}) {
      my $paths = $data_file->get_all_paths($self->opts->{datafile_dir});
      foreach my $path (@{$paths}) {
        if(-f $path) {
          $size += -s $path;
        }
      }
    }
  }
  $dba->dbc()->disconnect_if_idle();
  return $size;
}

sub _get_core_like_dbs {
  my ($self) = @_;
  my $dbas;
  if($self->opts->{species}) {
    $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-SPECIES => $self->opts->{species});
  }
  else {
     $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors();
  }
  my @final_dbas;
  while(my $dba = shift @{$dbas}) {
    next if $dba->species() eq 'multi';
    next if lc($dba->species()) eq 'ancestral sequences';
    next if $dba->dbc()->dbname() =~ /^.+_userdata$/xms;
    
    my $type = $dba->get_MetaContainer()->single_value_by_key('schema_type');
    $dba->dbc()->disconnect_if_idle();
    next unless $type;
    push(@final_dbas, $dba) if $type eq 'core';
  }
  $self->v('Found %d core like database(s)', scalar(@final_dbas));
  return \@final_dbas;
}

sub v {
  my ($self, $msg, @params) = @_;
  note sprintf($msg, @params);
  return;
}

sub _get_gid {
  my ($self, $group) = @_;
  my $group_uid;
  if ($group =~ /^\d+/) {
    $group_uid = $group;
    $group = ( getgrgid $group )[0];
  }
  else {
    $group_uid = (getgrnam($group))[2];
  }
  return $group_uid;
}

Script->run();
done_testing();

1;
__END__

=pod

=head1 NAME

sizes.pl

=head1 SYNOPSIS

  #BASIC
  ./sizes.pl -release VER -user USER -pass PASS -host HOST [-port PORT] -datafile_dir DIR [-species SPECIES] [-help | -man]
  
  #EXAMPLE
  ./sizes.pl -release 69 -host ensembdb.ensembl.org -port 5306 -user anonymous -datafile_dir /my/datafile
  
=head1 DESCRIPTION

A script which says how much space is used in total

=head1 OPTIONS

=over 8

=item B<--username | --user | -u>

REQUIRED. Username of the connecting account

=item B<--password | -pass | -p>

REQUIRED. Password of the connecting user.

=item B<--release | --version>

REQUIRED. Indicates the release of Ensembl to process

=item B<--host | --host | -h>

REQUIRED. Host name of the database to connect to

=item B<--port | -P>

Optional integer of the database port. Defaults to 3306.

=item B<--datafile_dir | --dir>

  -datafile_dir /datafile/dir

REQUIRED. Source directory which is the intended root of the datafiles.

=item B<--species>

Specify the tests to run over a single species' set of core databases

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item Bio::EnsEMBL

=item Post 66 databases

=back

=end

