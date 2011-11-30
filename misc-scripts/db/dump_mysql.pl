#!/usr/bin/env perl

package Script;

use strict;
use warnings;

use Carp;
use DBI;
use File::Spec;
use File::Path qw/mkpath/;
use Getopt::Long;
use IO::Compress::Gzip qw/gzip $GzipError/;
use Pod::Usage;
use Sys::Hostname;
use IO::File;

my $PIGZ_BINARY = 'pigz';
my $PIGZ_PROCESSORS = 2; #ensdb-1-* only have 4 cores
my $MAX_FILE_SIZE = 1_048_576; #anything greater than 1MB farm out

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->logging();
  $self->check();
  $self->defaults();
  $self->process();

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
      defaults
      version=i
      host=s
      port=i
      username=s
      password=s
      directory=s
      databases=s@
      groups=s@
      species=s@
      tables=s@
      pattern=s
      sql
      perlgzip
      verbose
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
    my $fh = IO::File->new($file, 'w');
    $fh->autoflush(1);
    my $oldfh = select($fh);
    $self->{oldfh} = $oldfh;
  }
  return;
}

sub check {
  my ($self) = @_;
  my $o = $self->opts();

  my @required_params;

  if ($o->{defaults}) {
    @required_params = qw/version/;
    pod2usage(
              -message => '-pattern is not supported with -defaults mode',
              -verbose => 1,
              -exitval => 1
    ) if $o->{pattern};
  } else {
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
  
  if(! $o->{username}) {
    pod2usage(
      -msg     => 'No -username given on the command line',
      -exitval => 1,
      -verbose => 0
    );
  }
  
  if(! $o->{password}) {
    pod2usage(
      -msg     => 'No -password given on the command line',
      -exitval => 1,
      -verbose => 0
    );
  }
  
  if ($o->{defaults}) {
    $self->_set_opts_from_hostname();
  } else {
    $o->{port} = 3306 if !$o->{port};
    if ($o->{pattern}) {
      $o->{databases} = $self->_all_dbs($o->{pattern});
    }
    $o->{directory} = File::Spec->rel2abs($o->{directory});
  }

  $self->v(q{Using the database server %s@%s:%d},
           map { $o->{$_} } qw/username host port/);

  #Filter for those on the specified server; sometimes redundant
  my %dbs = map { $_ => 1 } @{ $self->_all_dbs() };
  my @final_dbs = grep { $dbs{$_} } @{ $o->{databases} };
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
      $msg .= '. You specified the -database arg and -defaults. Are you on the correct server?';
    }
    pod2usage(-msg => $msg, -exitval => 1, -verbose => 0);
  }
  $self->v(q{Working %d database(s)}, $db_count);

  $o->{databases} = [ sort { $a cmp $b } @{ $o->{databases} } ];

  return;
}

sub process {
  my ($self) = @_;
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
    my $sql_file = $self->file($db . '.sql.gz');
    unlink $sql_file if -f $sql_file;
    my $fh = IO::Compress::Gzip->new($sql_file) or croak "Cannot create gzip stream to $sql_file: $GzipError";
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
    $self->compress($sql_file);

    #Checksum the DB's files
    $self->checksum();

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
  $self->compress($file);
  return;
}

sub checksum {
  my ($self) = @_;
  my $dir = $self->current_dir();

  $self->v('Checksumming directory %s', $dir);

  opendir(my $dh, $dir) or die "Cannot open directory $dir";
  my @files = sort { $a cmp $b } readdir($dh);
  closedir($dh) or die "Cannot close directory $dir";

  my $checksum = $self->file('CHECKSUMS.gz');
  unlink $checksum if -f $checksum;

  my $fh = IO::Compress::Gzip->new($checksum) or croak "Cannot create gzip stream to $checksum: $GzipError";
  foreach my $file (@files) {
    next if $file =~ /^\./;         #hidden file or up/current dir
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($dir, $file);
    my $sum = `sum $path`;
    $sum =~ s/\s* $path//xms;
    print $fh $file, "\t", $sum;
  }
  $fh->close();

  $self->compress($checksum);

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
  return $target_file;
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
  my $dir = File::Spec->catdir($self->opts()->{directory}, $db);
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
  return unless $o->{defaults};

  my $host     = $self->_host();
  my $settings = $self->_hostname_opts()->{$host};
  confess
"Specified -defaults but $host is not known to this script. Please edit the subroutine _hostname_opts() if you think it should be"
    if !$settings;

  #Setup default connection params
  $o->{host}     = $host;
  $o->{port}     = $settings->{port};

  if (!$o->{databases}) {
    $o->{databases} = $self->_all_dbs_regex($settings->{pattern});
  }

  #Set default dir
  $o->{directory} = $settings->{dir};

  return;
}

sub _hostname_opts {
  my ($self) = @_;

  my $version = $self->opts()->{version};

  my $target_dir = 'release-' . $version;
  my $default_dir =
    File::Spec->catdir(File::Spec->rootdir(), qw/mysql dumps/, $target_dir);

  return {
    'ensdb-1-01' => {
                      port     => 5306,
                      pattern  => qr/[a-m]\w*_ $version _\d+[a-z]?/xms,
                      dir      => $default_dir
    },
    'ensdb-1-02' => {
                      port     => 5306,
                      pattern  => qr/[n-z]\w*_ $version _\d+[a-z]?/xms,
                      dir      => $default_dir
    },
    'ensdb-1-03' => {
                      port     => 5303,
                      pattern  => qr/ensembl_compara_ $version/xms,
                      dir      => $default_dir
    },
    'ensdb-1-04' => {
                      port     => 5303,
                      pattern  => qr/ensembl_(ancestral|ontology)_ $version/xms,
                      dir      => $default_dir
    },
    'ensdb-1-05' => {
                      port     => 5316,
                      pattern  => qr/[efvg]\w+_mart_\w* $version/xms,
                      dir      => $default_dir
    },
    'ensdb-1-06' => {
                      port     => 5316,
                      pattern  => qr/[os]\w+_*mart_\w* $version/xms,
                      dir      => $default_dir
    },
    'ensdb-1-13' => {
                      port    => 5307,
                      pattern => qr/ensembl_website_ $version|ensembl_production_ $version/xms,
                      dir     => $default_dir
    },
  };
}

sub _host {
  my ($self) = @_;
  my $host = hostname();
  return $host;
}

sub _all_dbs_regex {
  my ($self, $pattern) = @_;
  confess 'No pattern given' if !$pattern;
  my $databases = $self->_all_dbs();
  return [ grep { $_ =~ $pattern } @{$databases} ];
}

sub _all_dbs {
  my ($self, $pattern) = @_;
  $pattern ||= '%';
  my $dbh       = $self->dbh();
  my $databases = $dbh->selectcol_arrayref('show databases like ?',
                                           { Columns => [1] }, $pattern);
  $self->clear_dbh();
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
  ./dump_mysql.pl --username USER --password PASS [--defaults] | [--host HOST [--port PORT] [-pattern '%' | -databases DB] [-tables TABLE] -directory DIR] [-help | -man]
  
  #Using defaults
  ./dump_mysql.pl --defaults --username root --password p --version 64
  
  ./dump_mysql.pl --defaults --username root --password p --version 64 --tables dna
  
  ./dump_mysql.pl --defaults --username root --password p --version 64 --tables meta,meta_coord --tables analysis --groups core,otherfeatures --groups vega
  
  ./dump_mysql.pl --defaults --username root --password p --version 64 --tables meta,meta_coord --tables analysis --groups core,otherfeatures --groups vega --sql 
  
  #Using host
  ./dump_mysql.pl --host srv --username root --password p --pattern '%_64%' --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db --databases other_db --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db,toto_db --databases other_db --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db --tables dna,dnac --directory $PWD/dumps
  
  ./dump_mysql.pl --host srv --username root --password p --databases my_db --tables dna --tables dnac --directory $PWD/dumps

=head1 DESCRIPTION

A script which is used to generate MySQL dumps which take into account issues
surrounding BLOB handling, VIEWS and other oddities of the Ensembl MySQL dump
process.

As a pose to normal scripts this version is aware of webteam database setup
and therefore can automatically configure itself to a set of useful 
parameters rather than having to manually configure the setup.

=head1 OPTIONS

=over 8

=item B<--username>

REQUIRED. Username of the connecting account. Must be able to perform 
C<SELECT INTO OUTFILE> calls.

=item B<--password>

REQUIRED. Password of the connecting user.

=item B<--defaults>

Uses the default mechanism which involves looking at the host the script
is executing on and setting a number of options for databases to look for
as well as port settings. C<-defaults> can be used in conjunction with
C<-groups>, C<-species> and C<-tables> but not with parameters like <--host>
and C<--pattern>.

=item B<--version>

If you are using C<--defaults> then you must also specify the version
of the databases you are dumping.

=item B<--host>

Host name of the database to connect to. Cannot be used with <--defaults>.

=item B<--port>

Optional integer of the database port. Defaults to 3306. Cannot be used 
with <--defaults>.

=item B<--pattern>

Allows the specification of a LIKE pattern to select databases with. Cannot
be used in conjunction with the C<--databases> argument. Cannot be used 
with <--defaults>.

=item B<--databases>

Allows database name specification and can be used more than once. Cannot
be used in conjunction with C<--pattern>. Comma separated values are 
supported.

=item B<--tables>

Allows you to specify a table to perform the dumps for. This will be applied
to all databases matching the given pattern or the list of databases. Be
warned that this will cause a full SQL re-dump and checksum re-calculation.

=item B<--directory>

Target directory to place all dumps. A sub-directory will be created here;
one per database dump. Cannot be used with <--defaults>.

=item B<--sql>

Force a dump of the SQL for the database and nothing else.

=item B<--verbose>

Makes the program give more information about what is going on. Otherwise
the program is silent.

=item B<--log>

If given the script will write all logs to output. Switches on C<--verbose>

=item B<--perlgzip>

Force the use of Perl's GZip libraries rather than using external zipping 
like pigz.

=item B<--help>

Help message

=item B<--man>

Man page

=back

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
