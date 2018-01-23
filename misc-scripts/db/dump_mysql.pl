#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


package Script;

use strict;
use warnings;

use Carp;
use DBI;
use File::Spec;
use File::Path qw/mkpath/;
use Getopt::Long qw/:config no_ignore_case auto_version bundling_override/;
use IO::Compress::Gzip qw/gzip $GzipError/;
use Pod::Usage;
use Sys::Hostname;

my $PIGZ_BINARY = 'pigz';
my $PIGZ_PROCESSORS = 2; #some machines only have 4 cores so do not go mad
my $MAX_FILE_SIZE = 1 * 1024 * 1024; #anything greater than 1MB farm out

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->logging();
  $self->check();
  $self->defaults();
  $self->dry() if $self->opts()->{dry};
  $self->process() if ! $self->opts()->{dry};

  if ($self->{oldfh}) {
    select($self->{oldfh});
  }
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {};
  GetOptions(
    $opts, qw/
      defaults=s
      version|release=i
      dry
      host|hostname|h=s
      port|P=i
      username|user|u=s
      password|pass|p=s
      directory|dir=s
      databases|database|db=s@
      groups=s@
      species=s@
      tables|table=s@
      pattern=s
      sql
      perlgzip
      testcompatible
      verbose|v
      log=s
      help
      man
      /
  ) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub logging {
  my ($self) = @_;
  my $o = $self->opts();
  if ($o->{log}) {
    $o->{verbose} = 1;
    my $file = $o->{log};
    open my $fh, '>', $file or die "Cannot open log file '${file}' for writing: $!";
    my $oldfh = select($fh);
    $self->{oldfh} = $oldfh;
  }
  return;
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params;

  if (! $o->{defaults}) {
    @required_params = qw/host username/;
    pod2usage(
              -message => '-pattern is not supported with -databases mode',
              -verbose => 1,
              -exitval => 1
    ) if $o->{pattern} && $o->{databases};
  }

  foreach my $r (@required_params) {
    if (!$o->{$r}) {
      pod2usage(
        -message =>
"-${r} has not been given at the command line but is a required parameter",
        -verbose => 1,
        -exitval => 1
      );
    }
  }
  
  #Check if gzip command is available
  if($o->{perlgzip}) {
    $self->{pigz_binary} = 0;
    $self->v('Forcing Perl based GZip compression');
  }
  else {
    `$PIGZ_BINARY --version >/dev/null 2>/dev/null`;
    $self->{pigz_binary} = ($? == 0) ? 1 : 0;
    my $feedback = ($self->{pigz_binary}) ? 'available' : 'not available';
    $self->v(q{pigz binary '%s' is %s}, $PIGZ_BINARY, $feedback);
  }

  return;
}

sub defaults {
  my ($self) = @_;
  my $o = $self->opts();

  #Processing -opt 1 -opt 2,3 into opt => [1,2,3]
  $self->_cmd_line_to_array('databases') if $o->{databases};
  $self->_cmd_line_to_array('groups')    if $o->{groups};
  $self->_cmd_line_to_array('species')   if $o->{species};
  
  my $original_databases_args = $o->{databases}; 

  #Tables
  if ($o->{tables} && !$o->{sql}) {
    $self->_cmd_line_to_array('tables');
    $self->v(q{Will work with the tables [%s]}, join(q{,}, @{ $o->{tables} }));
  }
  
  if ($o->{defaults}) {
    $self->_set_opts_from_hostname();
  } else {
    $o->{port} = 3306 if !$o->{port};
    if ($o->{pattern}) {
      my $p = $o->{pattern};
      $p = qr/$p/;
      $o->{databases} = $self->_all_dbs($p);
    }
    $o->{directory} = File::Spec->rel2abs($o->{directory});
  }
  
    if(! $o->{username}) {
    pod2usage(
      -msg     => 'No -username given on the command line or in the configuration file',
      -exitval => 1,
      -verbose => 0
    );
  }

  $self->v(q{Using the database server %s@%s:%d},
           map { $o->{$_} } qw/username host port/);

  #Filter for those on the specified server; sometimes redundant
  my %dbs = map { $_ => 1 } @{ $self->_all_dbs() };
  my @final_dbs;
  foreach my $db (@{$o->{databases}}) {
    if($dbs{$db}) {
      push(@final_dbs, $db);
    }
    else {
      $self->v('DB %s is not available from the specified server', $db);
    }
  }
  $o->{databases} = \@final_dbs;

  #Filtering DBs based on groups & species
  if ($o->{groups}) {
    my %dbs;
    foreach my $group (@{ $o->{groups} }) {
      $self->v('Filtering for group %s', $group);
      %dbs = map { $_ => 1 } grep { / _ $group _ /xms } @{ $o->{databases} };
    }
    $o->{databases} = [ keys %dbs ];
  }
  if ($o->{species}) {
    my %dbs;
    foreach my $species (@{ $o->{species} }) {
      $self->v('Filtering for species %s', $species);
      %dbs = map { $_ => 1 } grep { / $species _ /xms } @{ $o->{databases} };
    }
    $o->{databases} = [ keys %dbs ];
  }

  #Do we have any DBs left to process?
  my $db_count = scalar(@{ $o->{databases} });
  if ($db_count == 0) {
    my $msg = 'No databases found on the server ' . $o->{host};
    if($original_databases_args && $o->{defaults}) {
      my $version  = $o->{version} || '-NONE';
      $msg .= qq{. You specified the -database arg and -defaults. Are you on the correct server or did  -version '$version' excluded this DB?};
    }
    pod2usage(-msg => $msg, -exitval => 1, -verbose => 0);
  }
  $self->v(q{Working %d database(s)}, $db_count);

  $o->{databases} = [ sort { $a cmp $b } @{ $o->{databases} } ];
  
  $o->{verbose} = 1 if $o->{dry};

  return;
}

sub dry {
  my ($self) = @_;
  my $databases = $self->opts()->{databases};
  my $list = join(q{,}, @{$databases});
  $self->v(q{The following databases would have been dumped [%s]}, $list);
  return;
}

sub process {
  my ($self) = @_;
  
  
  my $test_case = $self->opts()->{testcompatible}; 
  
  $self->v('Producing test case compatible dumps') if $test_case;
  
  my $databases = $self->opts()->{databases};
  foreach my $db (@{$databases}) {
    $self->v('Working with database %s', $db);

    #Setup connection
    $self->dbh($db);

    $self->_setup_dir($db);

    #Get all tables
    my @tables = keys %{ $self->tables() };

    #Do data dumps if we didn't ask for just SQL
    if (!$self->opts()->{sql}) {
      my @tables_to_process;
      if ($self->opts()->{tables}) {
        my %lookup;
        %lookup = map { $_ => 1 } @{ $self->opts()->{tables} };
        @tables_to_process = grep { $lookup{$_} } @tables;
      } else {
        @tables_to_process = @tables;
      }
      foreach my $table (sort { $a cmp $b } @tables_to_process) {
        next if $self->is_view($table);
        $self->data($table);
      }
    } else {
      $self->v('-sql mode is on so no data dumping will occur');
    }

    #Do SQL
    my $sql_file;
    my $fh;
    if($test_case) {
      $sql_file = $self->file('table.sql');
      unlink $sql_file if -f $sql_file;
      open $fh, '>', $sql_file or croak "Cannot open filehandle to $sql_file: $!"; 
    }
    else {
      $sql_file = $self->file($db . '.sql.gz');
      unlink $sql_file if -f $sql_file;
      $fh = IO::Compress::Gzip->new($sql_file) or croak "Cannot create gzip stream to $sql_file: $GzipError";
    }
    
    my $writer = sub {
      my (@tabs) = @_;
      foreach my $table (sort { $a cmp $b } @tabs) {
        my $sql = $self->sql($table);
        print $fh $sql, ';', "\n" x 2;
      }
    };
    $writer->(grep { !$self->is_view($_) } @tables);
    $writer->(grep { $self->is_view($_) } @tables);
    $fh->close();
    $self->permissions($sql_file);
    
    #Checksum the DB's files
    $self->checksum() if ! $test_case;

    #Reset everything
    $self->clear_dbh();
    $self->clear_current_dir();
    $self->clear_tables();

    $self->v('Finished with database %s', $db);
  }
  return;
}

sub sql {
  my ($self, $table) = @_;
  my $q_table = $self->dbh()->quote_identifier($table);
  my $array =
    $self->dbh()
    ->selectcol_arrayref(qq{SHOW CREATE TABLE $q_table}, { Columns => [2] });
  my $sql = $array->[0];
  return $self->modify_sql($sql, $table);
}

sub modify_sql {
  my ($self, $sql, $table) = @_;
  if ($self->is_view($table)) {
    $sql =~ s/DEFINER=.+ \s+ SQL/DEFINER=CURRENT_USER() SQL/xms;
    $sql =~ s/SQL \s+ SECURITY \s+ DEFINER/SQL SECURITY INVOKER/xms;
  }
  if($self->opts()->{testcompatible}) {
    $sql =~ s/DEFAULT\s+CHARSET=latin1//xms;
    $sql =~ s/COLLATE=latin1_bin//xms;
    $sql =~ s/AUTO_INCREMENT=\d+//xms;
    $sql =~ s/CHARACTER SET latin1//g;
    $sql =~ s/COLLATE latin1_bin//g;
  }
  return $sql;
}

sub data {
  my ($self, $table) = @_;
  return if $self->is_view($table);
  $self->v('Dumping table %s', $table);
  my $q_table      = $self->dbh()->quote_identifier($table);
  my $file         = $self->file($table . '.txt');
  my $force_escape = q{FIELDS ESCAPED BY '\\\\'};
  my $sql          = sprintf(q{SELECT * FROM %s INTO OUTFILE '%s' %s},
                    $q_table, $file, $force_escape);
  unlink $file if -f $file;
  $self->dbh()->do($sql);
  $self->compress($file) if ! $self->opts()->{testcompatible};
  return;
}

sub checksum {
  my ($self) = @_;
  my $dir = $self->current_dir();

  $self->v('Checksumming directory %s', $dir);

  opendir(my $dh, $dir) or die "Cannot open directory $dir";
  my @files = sort { $a cmp $b } readdir($dh);
  closedir($dh) or die "Cannot close directory $dir";

  my $checksum = $self->file('CHECKSUMS');
  unlink $checksum if -f $checksum;

  open my $fh, '>', $checksum or croak "Cannot open filehandle to $checksum: $!";
  foreach my $file (@files) {
    next if $file =~ /^\./;         #hidden file or up/current dir
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($dir, $file);
    my $sum = `sum $path`;
    $sum =~ s/\s* $path//xms;
    chomp($sum);
    print $fh "${sum}\t${file}\n";
  }
  $fh->close();
  $self->permissions($checksum);

  return;
}

sub dbh {
  my ($self, $database) = @_;
  if (!exists $self->{'dbh'}) {
    my $o = $self->opts();
    my %args = (host => $o->{host}, port => $o->{port});
    $args{database} = $database if defined $database;

    my $dsn =
      'DBI:mysql:' . join(q{;}, map { $_ . '=' . $args{$_} } keys %args);
    $self->v('DBI connection URI %s', $dsn);
    my $dbh =
      DBI->connect($dsn, $o->{username}, $o->{password}, { RaiseError => 1 });

    $self->{dbh} = $dbh;
  }
  return $self->{'dbh'};
}

sub clear_dbh {
  my ($self) = @_;
  if (exists $self->{dbh}) {
    $self->{dbh}->disconnect();
    delete $self->{dbh};
  }
  return;
}

sub is_view {
  my ($self, $table) = @_;
  return ($self->tables()->{$table} eq 'VIEW') ? 1 : 0;
}

sub tables {
  my ($self) = @_;
  if (!exists $self->{tables}) {
    my $array =
      $self->dbh()->selectcol_arrayref(
'select TABLE_NAME, TABLE_TYPE from information_schema.TABLES where TABLE_SCHEMA = DATABASE()',
      { Columns => [ 1, 2 ] }
      );
    my %hits = @{$array};
    $self->{tables} = \%hits;
  }
  return $self->{tables};
}

sub clear_tables {
  my ($self) = @_;
  delete $self->{tables};
  return;
}

sub file {
  my ($self, $filename) = @_;
  return File::Spec->catfile($self->current_dir(), $filename);
}

sub current_dir {
  my ($self, $current_dir) = @_;
  $self->{'current_dir'} = $current_dir if defined $current_dir;
  return $self->{'current_dir'};
}

sub clear_current_dir {
  my ($self) = @_;
  delete $self->{current_dir};
  return;
}

sub opts {
  my ($self) = @_;
  return $self->{'opts'};
}

sub compress {
  my ($self, $file) = @_;
  my $target_file = $file . '.gz';

  $self->v(q{Compressing '%s' to '%s'}, $file, $target_file);

  if (-f $target_file) {
    unlink $target_file
      or die "Cannot remove the existing gzip file $target_file: $!";
  }
  
  my @stats = stat($file);
  my $size = $stats[7];
  if($self->{pigz_binary} && $size >= $MAX_FILE_SIZE) {
    system ("$PIGZ_BINARY --processes $PIGZ_PROCESSORS -4 $file") and confess "Could not pigz $file using $PIGZ_BINARY";
  }
  else {
    gzip $file => $target_file
      or die "gzip failed from $file to $target_file : $GzipError\n";
  }
  if (-f $target_file && -f $file) {
    unlink $file or die "Cannot remove the file $file: $!";
  }
  $self->permissions($target_file);
  return $target_file;
}

sub permissions {
  my ($self, $file) = @_;
  my $mode = 0666;
  chmod($mode, $file) or die "Cannot perform the chmod to mode $mode for file $file";
  return;
}

sub v {
  my ($self, $msg, @args) = @_;
  return unless $self->opts()->{verbose};
  my $s_msg = sprintf($msg, @args);
  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
    localtime(time());
  print sprintf("[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
                $mday, $mon, $year + 1900,
                $hour, $min, $sec, $s_msg);
  return;
}

sub _setup_dir {
  my ($self, $db) = @_;
  my @path = ($self->opts()->{directory});
  if($self->opts()->{testcompatible}) {
    if( $db =~ /^(?:\w+_test_db_)?([a-zA-Z0-9_]+)_([a-z]+)_\d+/) {
      push(@path, $1, $2);
    }
    else {
      $self->v("Cannot decipher name and group from $db. Using the database name");
      push(@path, $db);
    }
  }
  else {
    push(@path, $db);
  }
  my $dir = File::Spec->catdir(@path);
  $self->current_dir($dir);
  if (!-d $dir) {
    mkpath($dir) or die "Cannot create directory $dir: $!";
    chmod(0777, $dir)
      or die "Cannot change permissions on dir for everyone to write: $!";
  }
  return $dir;
}

sub _set_opts_from_hostname {
  my ($self) = @_;
  my $o = $self->opts();
  my $defaults = $o->{defaults};
  return unless $o->{defaults};
  confess "The given location '$defaults' does not exist" if ! -f $defaults; 

  my $host     = $self->_host();
  my $settings = $self->_hostname_opts_from_config()->{$host};
  confess "Specified -defaults but $host is not known. Check your $defaults ini file"
    if !$settings;

  #Setup default connection params
  $o->{host}      = $settings->{host} || $host; # use a configured host otherwise use hostname
  
  #only use if specified
  $o->{port}      = $settings->{port} if $settings->{port};
  $o->{username}  = $settings->{username} if $settings->{username};
  $o->{password}  = $settings->{password} if $settings->{password};
  $o->{sql}       = $settings->{sql} if $settings->{sql};

  if (!$o->{databases}) {
    my $opts_pattern = $o->{pattern};
    $opts_pattern = qr/$opts_pattern/ if $opts_pattern; 
    my $settings_pattern = $settings->{pattern};
    my $pattern = (defined $opts_pattern) ? $opts_pattern : $settings_pattern;
    $o->{databases} = $self->_all_dbs($pattern);
  }

  #Set default dir
  $o->{directory} = $settings->{dir};

  return;
}

#Assume normal ini-file format
sub _hostname_opts_from_config {
  my ($self) = @_;
  my $hostname_opts = {};
  my $target_dir;
  if($self->opts()->{version}) {
    $target_dir = 'release-' . $self->opts()->{version};
  }
  else {
    $target_dir = 'dumps';
  }
  
  my $defaults = $self->opts()->{defaults};
  open my $fh, '<', $defaults or confess "Cannot open defaults file '$defaults' for reading: $!";
  
  my $current_section;
  my $counter = 0;
  while(my $line = <$fh>) {
    $counter++;
    next if $line =~ /^\s*(?:\#|\;|$)/; #next for comments & empty lines
    $line =~ s/\s\;\s.+$//xmsg; #remove inline comments
    #Section [sec]
    if(my ($section) = $line =~ /^\s*\[\s*(.+?)\s*\]\s*$/xms) {
      $current_section = $section;
      $hostname_opts->{$current_section} = {};
      next;
    }
    # key = value
    if( my ($key, $value) = $line =~ /^\s*([^=]+?)\s*=\s*(.*?)\s*$/) {
      
      #Compile into a regex
      if($key eq 'pattern') {
        $value = qr/$value/;
      }
      #Change into the correct location
      elsif($key eq 'dir') {
        $value = File::Spec->catdir($value, $target_dir);
      }
      $hostname_opts->{$current_section}->{$key} = $value;
      next;
    }
    confess "Error in ini file '$defaults' at line $counter: '$line'";
  }
  close $fh;
  return $hostname_opts;
}

sub _host {
  my ($self) = @_;
  my $host = hostname();
  return $host;
}

#Always filter by version if it was given
sub _all_dbs {
  my ($self, $pattern) = @_;
  my $like = '%';
  if($self->opts()->{version}) {
    $like = '%\\_'.$self->opts()->{version}.'%';
    $self->v(q{Looking for databases with the pattern '%s' }, $like);
  }
  my $dbh       = $self->dbh();
  my $databases = $dbh->selectcol_arrayref('show databases like ?',
                                           { Columns => [1] }, $like);
  $self->clear_dbh();
  return [ grep { $_ =~ $pattern } @{$databases}] if $pattern;
  return $databases;
}

sub _cmd_line_to_array {
  my ($self, $key) = @_;
  my $array = $self->opts()->{$key};
  $array = (ref($array) && ref($array) eq 'ARRAY') ? $array : [$array];
  my $string = join(q{,}, @{$array});
  my @new_array = split(/,/, $string);
  $self->opts()->{$key} = \@new_array;
  return;
}

Script->run();

1;
__END__

=pod

=head1 NAME

dump_mysql.pl

=head1 SYNOPSIS

  #Basic
  ./dump_mysql.pl (-version VER | -release VER) [-defaults] | [ -username USER -password PASS -host HOST [-port PORT] [-pattern 'REGEX' | -databases DB] [-tables TABLE] -directory DIR] [-verbose] [-help | -man]
  
  #Test Case compatbile dumps
  ./dump_mysql.pl -username root -password pass -host 127.0.0.1 -testcompatible -verbose -directory /tmp/test-genome-DBs -database homo_sapiens_core_testdb
  
  #Using defaults ini file
  ./dump_mysql.pl --defaults my.ini --username root --password p --version 64
  
  ./dump_mysql.pl --defaults my.ini --username root --password p --release 64 -dry
  
  ./dump_mysql.pl --defaults my.ini --username root --password p --version 64 --tables dna
  
  ./dump_mysql.pl --defaults my.ini --username root --password p --version 64 --tables meta,meta_coord --tables analysis --groups core,otherfeatures --groups vega
  
  ./dump_mysql.pl --defaults my.ini --username root --password p --version 64 --tables meta,meta_coord --tables analysis --groups core,otherfeatures --groups vega --sql
  
  #Using host
  ./dump_mysql.pl --host srv --username root --password p --pattern '.+_64.+' --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db --databases other_db --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db,toto_db --databases other_db --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db --tables dna,dnac --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --db my_db --tables dna --tables dnac --directory $PWD/dumps

=head1 DESCRIPTION

A script which is used to generate MySQL dumps which take into account issues
surrounding BLOB handling, VIEWS and other oddities of the Ensembl MySQL dump
process.

You B<MUST> be on the database server the dumps are going to be generated
from. If not then all file manipulations will fail. 

As a pose to normal scripts this version is aware of webteam database setup
and therefore can automatically configure itself to a set of useful 
parameters rather than having to manually configure the setup.

=head1 OPTIONS

=over 8

=item B<--username | --user | -u>

REQUIRED. Username of the connecting account. Must be able to perform 
C<SELECT INTO OUTFILE> calls.

=item B<--password | -pass | -p>

REQUIRED. Password of the connecting user.

=item B<--defaults>

Uses the default mechanism which involves looking at the host the script
is executing on and setting a number of options for databases to look for
as well as port settings. C<-defaults> can be used in conjunction with
C<-groups>, C<-species> and C<-tables> but not with parameters like <--host>.

The options set are specified by your custom ini-file.

=item B<--version | --release>

If you are using C<--defaults> then you must also specify the version
of the databases you are dumping. Once specified the program will only
consider databases with the version number in there (specifically the 
occurance of C<%\_VERSION%>). C<--release> can also be used.

=item B<--host | --hostname | -h>

Host name of the database to connect to. Cannot be used with <--defaults>.

=item B<--port | -P>

Optional integer of the database port. Defaults to 3306. Cannot be used 
with <--defaults>.

=item B<--pattern>

Allows the specification of a regular expression to select databases with. 
Cannot be used in conjunction with the C<--databases> argument.

=item B<--databases | --database | --db>

Allows database name specification and can be used more than once. Cannot
be used in conjunction with C<--pattern>. Comma separated values are 
supported.

=item B<--tables | --table>

Allows you to specify a table to perform the dumps for. This will be applied
to all databases matching the given pattern or the list of databases. Be
warned that this will cause a full SQL re-dump and checksum re-calculation.

=item B<--directory | --dir>

Target directory to place all dumps. A sub-directory will be created here;
one per database dump. Cannot be used with <--defaults>.

=item B<--sql>

Force a dump of the SQL for the database and nothing else.

=item B<--testcompatible>

If specified will create MySQL dumps compatible with the Ensembl test
framework. This creates 2 levels of directory based on the database
species and the group it belongs to e.g. homo_sapiens and otherfeatures. We
also avoid gzipping any file and produce a single table.sql file.

=item B<--verbose>

Makes the program give more information about what is going on. Otherwise
the program is silent.

=item B<--log>

If given the script will write all logs to output. Switches on C<--verbose>

=item B<--perlgzip>

Force the use of Perl's GZip libraries rather than using external zipping 
like pigz.

=item B<--dry>

If specified the script will list all databases which have been found and
will be dumped but will not run any dumping process.

=item B<--help>

Help message

=item B<--man>

Man page

=back

=head1 DEFAULTS FILE FORMAT

The defaults file format is an ini file which respects all basic rules about
comments (preceeded by a ;), section headers and key/value pairs. The basic
form of an entry is

  ; basic format
  [server-name]
  port = 3306               ; port of the DB
  pattern = ^homo_sap\w+$   ; regular expression to filter DBs by
  dir = /path/to/dump/dir   ;
  
  ;more complex
  [other-server-name]
  port = 3306               ; port of the DB
  pattern = ^web\w+$        ; regular expression to filter DBs by
  dir = /path/to/dump/dir   ;
  sql = 1                   ; dump just the SQL for these databases
  
  ;if your host isn't the same as the server you are running the script on
  [mydumpserver]
  host = my-real-server     ; host you want the script to connect to
  port = 3306               ; port of the DB
  pattern = ^web\w+$        ; regular expression to filter DBs by
  dir = /path/to/dump/dir   ;
  sql = 1                   ; dump just the SQL for these databases
  
  ; and if you wanted everything in the config file so our cmd line becomes
  ; ./dump_mysql.pl -defaults mydbcfg.ini
  [myserver]
  port = 5306
  username = uberuser
  password = uberpassword
  pattern = ^homo_sap.+var.+$
  dir = /path/to/dump/dir

As an example of one which grabs all core dbs from a-m and puts it in /dumps

  [genome.mysql.server]
  port = 3306
  pattern = ^[a-m]\w*_core_\d+_\d+[a-z]?$
  dir = /dumps

The server name should be the same as what is emitted from

  perl -MSys::Hostname -e 'print hostname(), "\n"'

Also using line bounded regular expressions i.e. C<^> and C<$> will improve
the accuracy of the databases you are looking for.

=head1 REQUIREMENTS

=over 8

=item Perl 5.8+

=item DBI

=item DBD::mysql

=item IO::Compress::Gzip

=back

=head1 MAKING IT FASTER

=over 8

=item pigz L<http://www.zlib.net/pigz/> a parallel GZip compressor

=back

=end
