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

#Optionally bring in Bio::DB::Sam for BAM querying
my $NO_SAM_PERL;
eval 'require Bio::DB::Sam';
$NO_SAM_PERL = $@ if $@;

#Optionally look for samtools
my $NO_SAMTOOLS;
if($NO_SAM_PERL) {
  my $output = `which samtools`;
  $output = `locate samtools` unless $output;
  if(!$output) {
    $NO_SAMTOOLS = 'Unavailable after searching using `which` and `locate`. Add to your PATH';
  }
}

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
      unix_group=s
      species=s
      group=s
      no_bamchecks!
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
  
  $o->{unix_group} = 'ensembl' unless $o->{unix_group};
  
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

sub test_path {
  my ($self, $path, $data_file) = @_;
  
  #Just record the dir and file for later use when looking for extra files
  $path = rel2abs($path);
  my ($vol, $dir, $file) = splitpath($path);
  push(@{$self->{dirs}->{$dir}}, $file);
  
  my $name = $data_file->name();
  my $prefix = "Data file $name | File [$path]";
  
  #First pass. Check file & if not there then bail
  my $file_ok = ok(-f $path, "$prefix exists");
  return unless $file_ok;
  
  #File attributes now we know it's here
  my @stat = stat($path);
  my $mode = $stat[2];
  my $user_r = ($mode & S_IRUSR) >> 6;
  my $user_w = ($mode & S_IWUSR) >> 6;
  my $group_r = ($mode & S_IRGRP) >> 3;
  my $other_r = ($mode & S_IROTH) >> 0;
  
  my $user_rwx = ($mode & S_IRWXU) >> 6;
  my $group_rwx = ($mode & S_IRWXG) >> 3;
  my $other_rwx = ($mode & S_IRWXO);
  
  my $file_gid  = $stat[5];
  
  #Now do the tests
  ok(-s $path, "$prefix has data");
  is($user_rwx, 6, "$prefix is ReadWrite (mode 6) by user");
  is($group_rwx, 4, "$prefix is Read (mode 4) by group");
  is($other_rwx, 4, "$prefix is Read (mode 4) by owner");
  
  if($self->opts->{unix_group}) {
    my $unix_group = $self->opts->{unix_group};
    my $group_gid = $self->_get_gid($unix_group);
    if($group_gid) {
      my $group_ok = ok($group_gid == $file_gid, "$prefix is owned by group $unix_group");
      if(!$group_ok) {
        my $real_group = ( getgrgid $file_gid )[0];
        diag("$prefix belongs to $real_group ($file_gid) not $unix_group ($group_gid)");
      }
    }
    else {
      fail("The UNIX group $unix_group is not known on this system");
    }
  }
  
  return;
}

sub test_dirs {
  my ($self) = @_;
  foreach my $dir (keys %{$self->{dirs}}) {
    my $files = $self->{dirs}->{$dir};
    my %lookup = map { $_ => 1 } @{$files};
    my @all_file_paths = grep { $_ =~ /\/\.$/ && $_ =~ /\/\.{2}$/ } glob catfile($dir, "*.*");
    foreach my $path (@all_file_paths) {
      my ($vol, $dir, $file) = splitpath($path);
      ok($lookup{$file}, "$dir | $file was an expected file");
    }
  }
  return;
}

sub process {
  my ($self) = @_;
  my $dbas = $self->_get_core_like_dbs();
  while (my $dba = shift @{$dbas}) {
    $self->_process_dba($dba);
  }
  $self->test_dirs();
  return;
}

sub _process_dba {
  my ($self, $dba) = @_;
  $self->v('Working with species %s and group %s', $dba->species(), $dba->group());
  my $datafiles = $dba->get_DataFileAdaptor()->fetch_all();
  if(! @{$datafiles}) {
    $self->v("No datafiles found");
  }
  else {
    foreach my $data_file (@{$datafiles}) {
      my $paths = $data_file->get_all_paths($self->opts->{datafile_dir});
      foreach my $path (@{$paths}) {
        $self->test_path($path, $data_file);
      }
      if(!$self->opts->{no_bamchecks} && $data_file->file_type() eq 'BAM') {
        $self->_process_bam($dba, $data_file);
      }
    }
  }
  $dba->dbc()->disconnect_if_idle();
  return;
}

sub _process_bam {
  my ($self, $dba, $data_file) = @_;
  my $bam_names = $self->_get_bam_region_names($data_file);
  return unless $bam_names;
  my $toplevel_names = $self->_get_toplevel_slice_names($dba);
  my @missing_names;
  foreach my $bam_seq_name (@{$bam_names}) {
    if(! exists $toplevel_names->{$bam_seq_name}) {
      push(@missing_names, $bam_seq_name);
    }
  }
  
  my $missing_count = scalar(@missing_names);
  if($missing_count > 0) {
    fail("We have names in the BAM file which are not part of this assembly. Please see note output for more details");
    note(sprintf('Missing names: [%s]', join(q{,}, @missing_names)));
  }
  else {
    pass("All names in the BAM file are toplevel sequence region names");
  }
  
  return;
}

sub _get_core_like_dbs {
  my ($self) = @_;
  my $dbas;
  my %args;
  $args{-SPECIES} = $self->opts->{species};
  $args{-GROUP} = $self->opts->{group};
  if(%args) {
    $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(%args);
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
  return [sort { $a->species() cmp $b->species() } @final_dbas];
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

sub _get_bam_region_names {
  my ($self, $data_file) = @_;
  my $path = $data_file->path($self->opts->{datafile_dir});
  my $names;
  
  if(!$NO_SAM_PERL) {
    $names = $self->_get_bam_region_names_from_perl($path);
  }
  else {
    diag "Cannot use Bio::DB::Bam as it is not installed. Falling back to samtools compiled binary";
    if($NO_SAMTOOLS) {
      diag $NO_SAMTOOLS;
    }
    else {
      $names = $self->_get_bam_region_names_from_samtools($path);
    }
  }
  
  return $names;
}

sub _get_bam_region_names_from_perl {
  my ($self, $path) = @_;
  my @names;
  my $bam = Bio::DB::Bam->open($path);
  my $header = $bam->header();
  my $target_names = $header->target_name;
  my $length_arrayref = $header->target_len;
  for(my $i = 0; $i < $length_arrayref; $i++) {
    my $seq_id = $target_names->[$i];
    push(@names, $seq_id);
  }
  return \@names;
}

sub _get_bam_region_names_from_samtools {
  my ($self, $path) = @_;
  return if $NO_SAMTOOLS;
  my @names;  
  my $output = `samtools idxstats $path`;
  my @lines = split(/\n/, $output);
  foreach my $line (@lines) {
    my ($name) = split(/\s/, $line);
    next if $name eq '*';
    push(@names, $name);
  }
  return \@names;
}

sub _get_toplevel_slice_names {
  my ($self, $dba) = @_;
  my $species = $dba->species();
  if(! exists $self->{toplevel_names}->{$species}) {
    delete $self->{toplevel_names};
    my $core = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
    my $slices = $core->get_SliceAdaptor()->fetch_all('toplevel');
    my %lookup = map { $_->seq_region_name() => 1 } @{$slices};
    $self->{toplevel_names}->{$species} = \%lookup;
  }
  return $self->{toplevel_names}->{$species};
}

Script->run();
done_testing();

1;
__END__

=pod

=head1 NAME

check_datafiles.pl

=head1 SYNOPSIS

  #BASIC
  ./check_datafiles.pl -release VER -user USER -pass PASS -host HOST [-port PORT] -datafile_dir DIR [-unix_group UNIX_GROUP] \
                      [-species SPECIES] [-group GROUP] [-no_bamchecks] \
                      [-help | -man]
  
  #EXAMPLE
  ./check_datafiles.pl -release 69 -host ensembdb.ensembl.org -port 5306 -user anonymous -datafile_dir /my/datafile
  
  #Skipping BAM checks
  ./check_datafiles.pl -release 71 -host ensembdb.ensembl.org -port 5306 -user anonymous -datafile_dir /my/datafile -no_bamchecks
  
=head1 DESCRIPTION

A script which ensures a data files directory works for a single release.

The code outputs the results in Perl's TAP output making the format compatible with any TAP compatible binary such as prove. You should look
for all tests to pass. Scan for 'not ok' if there is an issue.

This code will also check the integrity of a BAM file e.g. that all BAM regions are toplevel sequences in the given Ensembl species. If
you want to disable this check then do so using the C<-no_bamchecks> flag.

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

=item B<--unix_group>

Specify the UNIX group these files should be readable by. Defaults to ensembl

=item B<--datafile_dir | --dir>

  -datafile_dir /datafile/dir

REQUIRED. Source directory which is the intended root of the datafiles.

=item B<--species>

Specify the tests to run over a single species' set of core databases

=item B<--group>

Only run tests on a single type of database registry group

=item B<--no_bamchecks>

Specify to stop any BAM file checks from occuring. This is normally sanity checks such
as are all BAM regions toplevel core regions

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

=item Bio::DB::Bam or samtools binary on ENV $PATH

=back

=end

