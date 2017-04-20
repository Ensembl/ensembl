#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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


use strict;
use warnings;

use DBI;
use Socket;
use English qw( -no_match_vars );
use File::Copy;
use File::Spec::Functions
  qw( rel2abs curdir canonpath updir catdir catfile );
use File::Temp;
use Getopt::Long;
use IO::File;
use Sys::Hostname;
use Tie::File;

$OUTPUT_AUTOFLUSH = 1;

my $start_time = time();

sub short_usage {
  print <<SHORT_USAGE_END;
Usage:
  $0 --pass=XXX \\
  \t[--user=XXX] \\
  \t[--noflush] [--nocheck] [--notargetflush]\\
  \t[--noopt] [--noinnodb] [--skip_views] [--force] \\
  \t[ --only_tables=XXX,YYY | --skip_tables=XXX,YYY ] \\
  \t[ --source_dir=/data/dir ] [ --target_dir=/data/dir ] \\
  \t[ --update] \\
  \t[ input_file |
  \t  --source=db\@host[:port] \\
  \t  --target=db\@host[:port] ]

  $0 --help

  Use --help to get a much longer help text.

SHORT_USAGE_END
}

sub long_usage {
  short_usage();

  print <<LONG_USAGE_END;
Description:

  Safetly copy MySQL databases between different servers and run
  myisamchk on the indexes when done.

  The script A) transfers the database files to a local staging
  directory, B) checks the files for errors using myisamchk, and C)
  moves them into place in the database server directory.


Command line switches:

  --pass=XXX        (Required)
                    The password for the MySQL user to
                    connect to the source database server.

  --user=XXX        (Optional)
                    MySQL user to connect to the source database server.
                    Default will be 'ensadmin'.

  --target_pass=XXX (Required)
                    The password for the MySQL user to
                    connect to the target database server.
                    If not specified, will be the same as --pass

  --target_user=XXX (Optional)
                    MySQL user to connect to the target database server.
                    If not specified, will be the same as --user

  --noflush         (Optional)
                    Skips table flushing completely.  Use very
                    sparsingly as copying un-flushed databases
                    potentially means missing data not yet flushed to
                    disk.  Use only after due consideration.

  --nocheck         (Optional)
                    Skip running myisamchk on the copied table files.
                    Use this only if you are *absolutely* sure that the
                    table files are ok and that you really do not have
                    time to wait for the check to run.  Use only after
                    due consideration.

  --noopt           (Optional)
                    Skip the optimization step.  The database tables are
                    optimized (as with "OPTIMIZE TABLE") to flush the
                    index data onto disk.  This may be a time consuming
                    operation for very large tables.  The --noopt flag
                    disables the optimization.

                    DO NOT RUN OPTIMIZE WHEN YOU ARE ON THE STAGING
                    MACHINES. THESE ARE OPTIMIZED BY THE PRODUCTION
                    TEAM DURING THE RELEASE CYCLE. ATTEMPTING TO DO
                    THIS WILL CAUSE THE SCRIPT TO DIE.

  --notargetflush   (Optional)
                    Skips table flushing on the target machine.

  --noinnodb        (Optional)
                    Skip the copy of any InnoDB table encountered.
                    Default is to include and fail when InnoDB seen.

  --skip_views      (Optional)
                    Exclude 'view' tables in any copy.  Default is to
                    process them.

  --force           (Optional)
                    Ordinarily, the script refuses to overwrite an
                    already existing staging directory.  This switch
                    will force it to re-use the directory if it exists.
                    This is useful for continuing an aborted copy.

                    NOTE: This switch will not ever force overwriting of
                    a server database directory.

                    NOTE: Using --force will bypass the table
                    optimization step.

  --update          (Optional)
                    Only sync tables that have changed on source database. This option
                    uses rsync with update and checksum to find the changes.
                    This will not replace any newer tables on target database.

  --drop            (Optional)
                     Drop the database on the target server before the copy.

  --only_tables=XXX,YYY
                    (Optional)
                    Only copy the tables specified in the
                    comma-delimited list.

  --skip_tables=XXX,YYY
                    (Optional)
                    Copy all tables except the ones specified in the
                    comma-delimited list.

  --tmpdir=TMP      (Optional)
                    Allows for a non-standard staging location to be
                    used during the copy process.  Normally it will
                    create a directory called 'tmp' in the directory
                    above the target data directory.

  --routines        (Optional)
                    Also copies functions and procedures

  --source_dir=/data/dir
                    (Optional)
                    MySQL server database source directory if different
                    from /mysql/data_3306/databases/

  --target_dir=/data/dir
                    (Optional)
                    MySQL server database target directory if different
                    from /mysql/data_3306/databases/. This option will 
                    also create a symlink for the database
                    from target dir to MySQL directory.

  --help            (Optional)
                    Displays this text.

The script takes either an input file of a certain format (see below),
or two additional switches:

  --source          (Optional, but required if there is no input file)
                    The source database to copy on the form
                    "db\@host:port" where the ":port" part is optional.

  --target          (Optional, but required if --source is used)
                    The target database to copy to on the same form as
                    with the --source switch.  The target database must
                    not already exist.

Input file format:

  The input file should be a tab/space-separated file with six fields
  per line:

    1. Source server
    2. Source server port
    3. Source database name

    4. Target server
    5. Target server port
    6. Target database name
    7. Target location (optional)

  For example:

  genebuild1 3306 at6_gor2_wga ens-staging1 3306 gorilla_gorilla_core_57_1

  (with each field being separated by a tab or space).

  This corresponds to the following command line options:

      --source=at6_gor2_wga\@genebuild1 \\
      --target=gorilla_gorilla_core_57_1\@ens-staging1

  (note that the port 3306 is the default port, and thus does not need
  to be mentioned)

  Blank lines, lines containing only whitespaces, and lines starting
  with '#', are silently ignored.

  Column 7 is used only when you need to copy the database to a location
  which is not the MySQL server's data directory.  The same rules
  applies; the mysqlens user at Sanger or other user at EBI must have write access to this directory
  and to the one above to create any temporary staging directory
  structures.  There is currently no way to specify this on the
  command-line.


Script restrictions:

  1. You must run the script on the destination server.

  2. The script must be run as the 'mysqlens' Unix user at Sanger. Talk to a
     recent release coordinator for access.

  3. The script will only copy MYISAM tables.  Databases with InnoDB
     tables will have to be copied manually using mysqldump.  InnoDB
     tables will make the script throw an error in the table checking
     stage.


LONG_USAGE_END
} ## end sub long_usage

my ( $opt_password, $opt_user, $opt_password_tgt, $opt_user_tgt, $opt_only_tables, $opt_skip_tables, $opt_help );

my $opt_flush    = 1;    # Flush by default.
my $opt_check    = 1;    # Check tables by default.
my $opt_optimize = 1;    # Optimize the tables by default.
my $opt_force = 0; # Do not reuse existing staging directory by default.
my $opt_update = 0; # Only copy tables that have changed
my $opt_skip_views = 0;    # Process views by default
my $opt_innodb     = 1;    # Don't skip InnoDB by default
my $opt_flushtarget = 1;
my $opt_drop =0 ;          # We don't drop database on target server by default
my $opt_tmpdir;
my $opt_routines = 0;
my ( $opt_source, $opt_target, $opt_source_dir, $opt_target_dir);

if ( !GetOptions( 
                  'pass=s'        => \$opt_password,
                  'user=s'        => \$opt_user,
                  'target_pass=s' => \$opt_password_tgt,
                  'target_user=s' => \$opt_user_tgt,
                  'flush!'        => \$opt_flush,
                  'flushtarget!'  => \$opt_flushtarget,
                  'check!'        => \$opt_check,
                  'opt!'          => \$opt_optimize,
                  'force!'        => \$opt_force,
                  'update!'       => \$opt_update,
                  'drop!'         => \$opt_drop,
                  'only_tables=s' => \$opt_only_tables,
                  'skip_tables=s' => \$opt_skip_tables,
                  'innodb!'       => \$opt_innodb,
                  'skip_views!'   => \$opt_skip_views,
                  'tmpdir=s'      => \$opt_tmpdir,
                  'help!'         => \$opt_help,
                  'source=s'      => \$opt_source,
                  'target=s'      => \$opt_target,
                  'routines'      => \$opt_routines,
                  'source_dir=s'  => \$opt_source_dir,
                  'target_dir=s'  => \$opt_target_dir,
     ) ||
     ( !defined($opt_password) && !defined($opt_help) ) )
{
  short_usage();
  exit 1;
}

if ( defined($opt_help) ) {
  long_usage();
  exit 0;
}

if ( ( defined($opt_source) && !defined($opt_target) ) ||
     ( !defined($opt_source) && defined($opt_target) ) )
{
  die("You must use --source and --target together.\n");
}

if ( defined($opt_skip_tables) && defined($opt_only_tables) ) {
  die("Can't use both --only_tables and --skip_tables\n");
}

if ( $opt_force && $opt_optimize ) {
  print("Note: Using --force will bypass optimization.\n");
  $opt_optimize = 0;
}

if ( $opt_force && $opt_update ) {
  die("Can't use both --force and --opt_update\n");
}

my %only_tables;
if ( defined($opt_only_tables) ) {
  %only_tables = map( { $_ => 1 } split( /,/, $opt_only_tables ) );
}
my %skip_tables;
if ( defined($opt_skip_tables) ) {
  %skip_tables = map( { $_ => 1 } split( /,/, $opt_skip_tables ) );
}

if ( scalar( getpwuid($<) ) ne 'ensmysql' and scalar( getpwuid($<) ) ne 'mysqlens' ) {
  warn("You are not running this script as the 'ensmysql' user at the EBI or as the 'mysqlens' at Sanger.\n");
}

if (!defined($opt_user))
{
  $opt_user='ensadmin';
}

$opt_user_tgt ||= $opt_user;
$opt_password_tgt ||= $opt_password;

my $input_file;

if ( !defined($opt_source) ) {
  $input_file = shift(@ARGV);

  if ( !defined($input_file) ) {
    short_usage();
    exit 1;
  }
}

my @executables =('myisamchk','rsync');

# Make sure we can find all executables.
foreach my $key (@executables) {
  my $output = `which $key`;
  my $rc     = $? >> 8;
  
  if ( !$opt_check && $key eq 'myisamchk' ) {
      print( "Can not find 'myisamchk' " .
      "but --nocheck was specified so skipping\n" ),;
  }
  elsif($rc != 0) {
    chomp $output;
    die(
      sprintf(
      "Can not find '%s' and 'which %s' "
      . "yields anything useful. Check your \$PATH",
      $key, $key
      ) );
  }
} ## end foreach my $key ( @executables...)

my $run_hostaddr = hostname_to_ip(hostname());
my $working_dir = rel2abs( curdir() );

##====================================================================##
##  Read the configuration file line by line and try to validate all  ##
##  parts of each line.  Store the validated information in the @todo ##
##  list (a list of hashes) for later processing.                     ##
##====================================================================##

my $do_unlink_tmp_file = 0;

if ( !defined($input_file) ) {
  # The user uses --source and --target rather than an input file.
  # Create our own temporary input file with the info from these options
  # and then parse that.  #FIXME: this is a hack.
  my $out = File::Temp->new( UNLINK => 0 ) or
    die("Can not create temporary file.\n");

  my ( $source_db, $source_server, $source_port ) =
    ( $opt_source =~ /^([^@]+)@([^:]+):?(\d+)?$/ );
  my ( $target_db, $target_server, $target_port ) =
    ( $opt_target =~ /^([^@]+)@([^:]+):?(\d+)?$/ );

  $out->printf( "%s %d %s %s %d %s\n",
                $source_server,
                defined($source_port) ? $source_port : 3306,
                $source_db,
                $target_server,
                defined($target_port) ? $target_port : 3306,
                $target_db );

  $input_file = $out->filename();

  # This temporary file will be unlinked further down.
  $do_unlink_tmp_file = 1;
} ## end if ( !defined($input_file...))

my $in = IO::File->new( '<' . $input_file ) or
  die( sprintf( "Can not open '%s' for reading", $input_file ) );

my @todo;    # List of verified databases etc. to copy.

my $lineno = 0;
while ( my $line = $in->getline() ) {
  ++$lineno;

  $line =~ s/^\s+//;    # Strip leading whitespace.
  $line =~ s/\s+$//;    # Strip trailing whitespace.

  if ( $line =~ /^\#/ )   { next }    # Comment line.
  if ( $line =~ /^\s*$/ ) { next }    # Empty line.

  my $failed = 0;                     # Haven't failed so far...

  my ( $source_server, $source_port, $source_db,
       $target_server, $target_port, $target_db,
       $target_location
  ) = split( /\s+/, $line );

  my $source_hostaddr = hostname_to_ip($source_server);
  my $target_hostaddr = hostname_to_ip($target_server);

  # Verify source server and port.
  if ( !defined($source_hostaddr) || $source_hostaddr eq '' ) {
    warn( sprintf( "line %d: Source server '%s' is not valid.\n",
                   $lineno, $source_server
          ) );
    $failed = 1;
  }

  if ( !defined($source_port) || $source_port =~ /\D/ ) {
    warn( sprintf( "line %d: Source port '%s' is not a number.\n",
                   $lineno, $source_port || ''
          ) );
    $failed = 1;
  }

  # Verify target server and port.
  if ( !defined($target_hostaddr) || $target_hostaddr eq '' ) {
    warn( sprintf( "line %d: Target server '%s' is not valid.\n",
                   $lineno, $target_server
          ) );
    $failed = 1;
  }

  if ( !defined($target_port) || $target_port =~ /\D/ ) {
    warn( sprintf( "line %d: Target port '%s' is not a number.\n",
                   $lineno, $target_port || ''
          ) );
    $failed = 1;
  }

  # Make sure we running on the target server.
  if ( !$failed && $run_hostaddr ne $target_hostaddr ) {
    warn( sprintf(
            "line %d: " .
              "This script needs to be run on the destination server " .
              "'%s' ('%s').\n",
            $lineno, $target_server, $target_hostaddr
          ) );
    $failed = 1;
  }

  if(! $failed && $source_server =~ /^ens-staging\d?$/ && $opt_optimize) {
    my $tmpl = q{line %d: Source server '%s' is an ens-staging machine. Do not optimize on this server. Rerun with -noopt};
    warn(sprintf($tmpl, $lineno, $source_server));
    $failed = 1;
  }

  if ( !$failed ) {
    push( @todo, {
            'source_server'   => $source_server,
            'source_hostaddr' => $source_hostaddr,
            'source_port'     => $source_port,
            'source_db'       => $source_db,
            'target_server'   => $target_server,
            'target_hostaddr' => $target_hostaddr,
            'target_port'     => $target_port,
            'target_db'       => $target_db,
            'target_location' => $target_location,
         } );
  }
} ## end while ( my $line = $in->getline...)

$in->close();

if ($do_unlink_tmp_file) {
  # Unlink the temporary input file.
  unlink($input_file);
}

##====================================================================##
##  Take the copy specifications from the @todo list and for each     ##
##  specification copy the database to a staging area using rsync,    ##
##  check it with myisamchk, and move it in place in the database     ##
##  directory.                                                        ##
##====================================================================##

TODO:
foreach my $spec (@todo) {
  my $source_server   = $spec->{'source_server'};
  my $source_hostaddr = $spec->{'source_hostaddr'};
  my $source_port     = $spec->{'source_port'};
  my $source_db       = $spec->{'source_db'};
  my $target_server   = $spec->{'target_server'};
  my $target_hostaddr = $spec->{'target_hostaddr'};
  my $target_port     = $spec->{'target_port'};
  my $target_db       = $spec->{'target_db'};
  my $target_location = $spec->{'target_location'};

  my $label = sprintf( "{ %s -> %s }==", $source_db, $target_db );
  print( '=' x ( 80 - length($label) ), $label, "\n" );

  print( "CONNECTING TO SOURCE AND TARGET DATABASES " .
         "TO GET 'datadir'\n" );

  my $source_dsn = sprintf( "DBI:mysql:database=%s;host=%s;port=%d",
                            $source_db, $source_hostaddr,
                            $source_port );

  my $source_dbh = DBI->connect( $source_dsn,
                                 $opt_user,
                                 $opt_password, {
                                   'PrintError' => 1,
                                   'AutoCommit' => 0
                                 } );

  if ( !defined($source_dbh) ) {
    warn( sprintf(
               "Failed to connect to the source database '%s@%s:%d'.\n",
               $source_db, $source_server, $source_port
          ) );

    $spec->{'status'} =
      sprintf( "FAILED: can not connect to source database '%s@%s:%d'.",
               $source_db, $source_server, $source_port );
    next TODO;
  }

  my $target_dsn;
  # If using the update option, connect to the target database  
  if ($opt_update){
     $target_dsn = sprintf( "DBI:mysql:database=%s;host=%s;port=%d",
                            $target_db, $target_hostaddr, $target_port );
  }
  else{
    $target_dsn = sprintf( "DBI:mysql:host=%s;port=%d",
                             $target_hostaddr, $target_port );
  }

  my $target_dbh = DBI->connect( $target_dsn,
                                 $opt_user_tgt,
                                 $opt_password_tgt, {
                                   'PrintError' => 1,
                                   'AutoCommit' => 0
                                 } );

  if ( !defined($target_dbh) ) {
    warn( sprintf( "Failed to connect to the target server '%s:%d'.\n",
                   $target_server, $target_port
          ) );

    $spec->{'status'} =
      sprintf( "FAILED: can not connect to target server '%s:%d'.",
               $target_server, $target_port );

    $source_dbh->disconnect();
    next TODO;
  }

  # Assigning Source and target directory.
  my $source_dir;
  my $target_dir;

  # Get source and target server data directories.
  if ( defined($opt_source_dir) ) {
    $source_dir = $opt_source_dir;
  }
  else {
    $source_dir =
      $source_dbh->selectall_arrayref("SHOW VARIABLES LIKE 'datadir'")
      ->[0][1];
  }
  if ( defined($opt_target_dir) ) {
    $target_dir = $opt_target_dir;
  }
  else {
    $target_dir =
      $target_dbh->selectall_arrayref("SHOW VARIABLES LIKE 'datadir'")
      ->[0][1];
  }
  # Getting staging machine MySQL database directory.
  my $mysql_database_dir =
    $target_dbh->selectall_arrayref("SHOW VARIABLES LIKE 'datadir'")
    ->[0][1];
  
  if ( !defined($source_dir) ) {
    warn(
      sprintf(
        "Failed to find data directory for source server at '%s:%d'.\n",
        $source_server, $source_port
      ) );

    $spec->{'status'} = sprintf(
        "FAILED: can not find data directory on source server '%s:%d'.",
        $source_server, $source_port );

    $source_dbh->disconnect();
    next TODO;
  }

  if ( defined($target_location) ) { $target_dir = $target_location }

  if ( !defined($target_dir) ) {
    warn(
      sprintf(
        "Failed to find data directory for target server at '%s:%d'.\n",
        $target_server, $target_port
      ) );

    $spec->{'status'} = sprintf(
        "FAILED: can not find data directory on target server '%s:%d'.",
        $target_server, $target_port );

    $source_dbh->disconnect();
    next TODO;
  }

  printf( "SOURCE 'datadir' = '%s'\n", $source_dir );
  printf( "TARGET 'datadir' = '%s'\n", $target_dir );

  my $destination_dir = catdir( $target_dir, $target_db );
  my $staging_dir;
  my $tmp_dir;

  # If we update the database, we don't need a tmp dir
  # We will use the dest dir instead of staging dir.
  if (!$opt_update) {
    if ( defined($opt_tmpdir) ) { $tmp_dir = $opt_tmpdir }
    else {
      if ((scalar( getpwuid($<) ) eq 'mysqlens')) {
        $tmp_dir = canonpath( catdir( $target_dir, updir(), 'tmp' ));
      }
      elsif ((scalar( getpwuid($<) ) eq 'ensmysql')){
       $tmp_dir = canonpath( catdir( $target_dir, updir(), 'temp' ));
      }
    }

    if ( !-d $tmp_dir ) {
      die(sprintf( "Can not find the temporary directory '%s'", $tmp_dir )
      );
    }

    printf( "TMPDIR = %s\n",             $tmp_dir );

    $staging_dir = catdir( $tmp_dir, sprintf( "tmp.%s", $target_db ) );

    # Try to make sure the temporary directory and the final destination
    # directory actually exists, and that the staging directory within the
    # temporary directory does *not* exist.  Allow the staging directory
    # to be reused when the --force switch is used.



    if ( !$opt_drop && -d $destination_dir ) {
      warn( sprintf( "Destination directory '%s' already exists.\n",
                   $destination_dir ) );

      $spec->{'status'} =
        sprintf( "FAILED: database destination directory '%s' exists.".
                 " You can use the --drop option to drop the database".
                 " on the target server" ,
               $destination_dir );

      $source_dbh->disconnect();
      next TODO;
    }

    if ($opt_drop){
      printf ("Dropping database %s on %s \n", $target_db, $target_server);
      $target_dbh->do("DROP DATABASE $target_db;");
      $target_dbh->disconnect();
    }
    else{
      $target_dbh->disconnect();
    }

    if ( !$opt_force && -d $staging_dir ) {
      warn( sprintf( "Staging directory '%s' already exists.\n",
                   $staging_dir ) );

      $spec->{'status'} =
        sprintf( "FAILED: staging directory '%s' exists.".
          " You can use the --force option to overwrite database in staging dir", $staging_dir);

      $source_dbh->disconnect();
      next TODO;
    }

    if ( !mkdir($staging_dir) ) {
      if ( !$opt_force || !-d $staging_dir || !$opt_update) {
        warn( sprintf( "Failed to create staging directory '%s'.\n",
                     $staging_dir ) );

        $spec->{'status'} =
          sprintf( "FAILED: can not create staging directory '%s'.",
                 $staging_dir );

        $source_dbh->disconnect();
        next TODO;
      }
    }
  }
  else{
    $staging_dir=$destination_dir;
  }

  $spec->{'status'} = 'SUCCESS';    # Assume success until failure.

  my @tables;
  my @views;

  my $table_sth = $source_dbh->prepare('SHOW TABLE STATUS');

  $table_sth->execute();

  my %row;
  # Fancy magic from DBI manual.
  $table_sth->bind_columns( \( @row{ @{ $table_sth->{'NAME_lc'} } } ) );

TABLE:
  while ( $table_sth->fetch() ) {
    my $table  = $row{'name'};
    my $engine = $row{'engine'};

    if ( defined($opt_only_tables) && !exists( $only_tables{$table} ) )
    {
      next TABLE;
    }
    elsif ( defined($opt_skip_tables) &&
            exists( $skip_tables{$table} ) )
    {
      next TABLE;
    }

    if ( defined($engine) ) {
      if ( $engine eq 'InnoDB' ) {
        if ( !$opt_innodb ) {
          warn( sprintf( "SKIPPING InnoDB table '%s.%s'\n",
                         $source_db, $table
                ) );
          next TABLE;
        }
        else {
          $spec->{'status'} =
            sprintf( "FAILED: can not copy InnoDB tables.  " .
                       "Please convert '%s.%s' to MyISAM " .
                       "or run with --noinnodb.",
                     $source_db, $table );
          next TODO;
        }
      }
    }
    else {
      if ($opt_skip_views) {
        warn( sprintf( "SKIPPING view '%s'\n", $table ) );
      }
      else {
        push( @views, $table );
      }
      next TABLE;
    }

    push( @tables, $table );
  } ## end while ( $table_sth->fetch...)
  ## For MySQL version 5.6 and above,  FLUSH TABLES is not permitted when there is an active READ LOCK.
  if ($source_dbh->{mysql_serverversion} >= "50600"){
    if ($opt_flush) {
      # Flush and Lock tables with a read lock.
      print("FLUSHING AND LOCKING TABLES SOURCE...\n");
      $source_dbh->do(
                  sprintf( "FLUSH TABLES %s WITH READ LOCK", join( ', ', @tables ) ) );
    }
    else {
      warn( sprintf("You are running the script with --noflush and MySQL source server version is $source_dbh->{mysql_serverversion}. " .
                       "The database will not be locked during the copy. " .
                       "This is not recomended!!!"));
    }
  }
  else {
    # Lock tables with a read lock.
    print("LOCKING TABLES...\n");
    $source_dbh->do(
         sprintf( "LOCK TABLES %s READ", join( ' READ, ', @tables ) ) );

    if ($opt_flush) {
      # Flush tables.

      print("FLUSHING TABLES...\n");
      $source_dbh->do(
                  sprintf( "FLUSH TABLES SOURCE %s", join( ', ', @tables ) ) );
    }
  }

  # If we use the update option, also flush and lock target server
    if ($opt_update and ($target_dbh->{mysql_serverversion} >= "50600")){
      if ($opt_flush) {
        print("FLUSHING AND LOCKING TABLES TARGET...\n");
        $target_dbh->do(
                  sprintf( "FLUSH TABLES %s WITH READ LOCK", join( ', ', @tables ) ) );
      }
      else {
      warn( sprintf("You are running the script with --noflush and MySQL target server version is $target_dbh->{mysql_serverversion}. " .
                       "The database will not be locked during the copy. " .
                       "This is not recomended!!!"));
      }
    }
  # If update, also flush target server
    elsif ($opt_update and ($target_dbh->{mysql_serverversion} < "50600")){
          # Lock tables with a read lock.
      print("LOCKING TABLES...\n");
      $target_dbh->do(
         sprintf( "LOCK TABLES %s READ", join( ' READ, ', @tables ) ) );
      if ($opt_flush) {
        print("FLUSHING TABLES...\n");
        $target_dbh->do(
                  sprintf( "FLUSH TABLES TARGET %s", join( ', ', @tables ) ) );
      }
    }
  ##------------------------------------------------------------------##
  ## OPTIMIZE                                                         ##
  ##------------------------------------------------------------------##

  if ($opt_optimize) {
    print( '-' x 35, ' OPTIMIZE ', '-' x 35, "\n" );
    print 'OPTIMIZING TABLES...', "\n";
    foreach my $table (@tables) {
      $source_dbh->do( sprintf( "OPTIMIZE TABLE %s", $table ) );
    }
  }

  ##------------------------------------------------------------------##
  ## COPY                                                             ##
  ##------------------------------------------------------------------##

  print( '-' x 37, ' COPY ', '-' x 37, "\n" );

  # Set up database copying.  We're using rsync for this because it's
  # using SSH for network transfers, because it may be used for local
  # copy too, and because it has good inclusion/exclusion filter
  # options.

  my @copy_cmd;
  
  @copy_cmd = ('rsync');

  push(@copy_cmd, '--whole-file', '--archive', '--progress' );

  if ($opt_force) {
    push( @copy_cmd, '--delete', '--delete-excluded' );
  }
  
  # If we use the update option, we want to make sure that
  # We don't overwrite newer tables on the target db
  # We want to use checksum to find the changes
  # We want to exclude the index files MYI as we will do an 
  # optimize on the target database after copy.
  if ($opt_update) {
    push( @copy_cmd, '--update', '--checksum', '--exclude=*.MYI' );
  }

  # Set files permission to 755 (rwxr-xr-x)
  push (@copy_cmd, '--chmod=Du=rwx,go=rx,Fu=rwx,go=rx');  

  # Add TCP with arcfour encryption, TCP does go pretty fast (~110 MB/s) and is a better choice in LAN situation.
  push(@copy_cmd, '-e', q{ssh -c arcfour} );

  if ( defined($opt_only_tables) ) {
    push( @copy_cmd, '--include=db.opt' );

    push( @copy_cmd,
          map { sprintf( '--include=%s.*', $_ ) } keys(%only_tables) );

    # Partitioned tables:
    push( @copy_cmd,
          map { sprintf( '--include=%s#P#*.*', $_ ) }
            keys(%only_tables) );

    push( @copy_cmd, "--exclude=*" );
  }
  elsif ( defined($opt_skip_tables) ) {
    push( @copy_cmd, '--include=db.opt' );

    push( @copy_cmd,
          map { sprintf( '--exclude=%s.*', $_ ) } keys(%skip_tables) );

    # Partitioned tables:
    push( @copy_cmd,
          map { sprintf( '--exclude=%s#P#*.*', $_ ) }
            keys(%skip_tables) );

    push( @copy_cmd, "--include=*" );
  }

  if ( $source_hostaddr eq $target_hostaddr ) {
    # Local copy.
    push( @copy_cmd,
          sprintf( "%s/", catdir( $source_dir, $source_db ) ) );
  }
  else {
    # Copy from remote server.
    push( @copy_cmd,
          sprintf( "%s:%s/",
                   $source_hostaddr, catdir( $source_dir, $source_db ) )
    );
  }

  push( @copy_cmd, sprintf( "%s/", $staging_dir ) );

  # Perform the copy and make sure it succeeds.

  printf( "COPYING '%s:%d/%s' TO STAGING DIRECTORY '%s'\n",
          $source_server, $source_port, $source_db, $staging_dir );

  # For debugging:
  # print( join( ' ', @copy_cmd ), "\n" );

  my $copy_failed = 0;
  if ( system(@copy_cmd) != 0 ) {
    warn( sprintf( "Failed to copy database.\n" .
                     "Please clean up '%s' (if needed).",
                   $staging_dir
          ) );
    $copy_failed = 1;
  }

  # Unlock tables.
  print("UNLOCKING TABLES SOURCE...\n");
  $source_dbh->do('UNLOCK TABLES');

  $source_dbh->disconnect();

  if ($opt_update){
    print("UNLOCKING TABLES TARGET...\n");
    $target_dbh->do('UNLOCK TABLES');
    $target_dbh->disconnect();
  }

  if ($copy_failed) {
    $spec->{'status'} =
      sprintf( "FAILED: copy failed (cleanup of '%s' may be needed).",
               $staging_dir );
    next TODO;
  }

  ##------------------------------------------------------------------##
  ## VIEW REPAIR                                                      ##
  ##------------------------------------------------------------------##
  print( '-' x 36, ' VIEW REPAIR ', '-' x 37, "\n" );

  if ($opt_skip_views) {
    print 'SKIPPING VIEWS...', "\n";
  }
  else {
    print 'PROCESSING VIEWS...', "\n";

    if ( $source_db eq $target_db ) {
      print( "Source and target names (${source_db}) are the same; " .
             "views do not need repairing\n" );
    }
    else {
      my $ok = 1;

    VIEW:
      foreach my $current_view (@views) {
        print "Processing $current_view\n";

        my $view_frm_loc =
          catfile( $staging_dir, "${current_view}.frm" );

        if ( tie my @view_frm, 'Tie::File', $view_frm_loc ) {
          for (@view_frm) {
            s/`$source_db`/`$target_db`/g;
          }
          untie @view_frm;
        }
        else {
          warn( sprintf( q{Cannot tie file '%s' for VIEW repair. Error},
                         $view_frm_loc ) );
          $ok = 0;
          next VIEW;
        }
      }

      if ( !$ok ) {
        $spec->{status} = sprintf(
                   "FAILED: view cleanup failed " .
                     "(cleanup of view frm files in '%s' may be needed",
                   $staging_dir );
      }
    } ## end else [ if ( $source_db eq $target_db)]
  } ## end else [ if ($opt_skip_views) ]

  # if we use the update option, optimize the target database
  if ($opt_update) {
    ##------------------------------------------------------------------##
    ## OPTIMIZE TARGET                                                  ##
    ##------------------------------------------------------------------##

    $target_dbh = DBI->connect( $target_dsn,
                                $opt_user_tgt,
                                $opt_password_tgt, {
                                'PrintError' => 1,
                                'AutoCommit' => 0
                             } );
    if ($opt_optimize) {
      print( '-' x 35, ' OPTIMIZE TARGET ', '-' x 35, "\n" );
      print 'OPTIMIZING TABLES...', "\n";
      foreach my $table (@tables) {
        $target_dbh->do( sprintf( "OPTIMIZE TABLE %s", $table ) );
      }
    }
    $target_dbh->disconnect;
  }

  ##------------------------------------------------------------------##
  ## CHECK                                                            ##
  ##------------------------------------------------------------------##

  print( '-' x 36, ' CHECK ', '-' x 37, "\n" );

  # Check the copied table files with myisamchk.  Let myisamchk
  # automatically repair any broken or un-closed tables.

  if ( !$opt_check ) {
    print("NOT CHECKING...\n");
  }
  else {
    print("CHECKING TABLES...\n");

    my $check_failed = 0;

    foreach my $table (@tables) {
      foreach my $index (
         glob( catfile( $staging_dir, sprintf( '%s*.MYI', $table ) ) ) )
      {
        my @check_cmd = (
          'myisamchk', '--check', '--check-only-changed',
          '--update-state', '--silent', '--silent',    # Yes, twice.
          $index );

        if ( system(@check_cmd) != 0 ) {
          $check_failed = 1;

          warn( sprintf( "Failed to check some tables. " .
                           "Please clean up '%s'.\n",
                         $staging_dir
                ) );

          $spec->{'status'} =
            sprintf( "FAILED: MYISAM table check failed " .
                       "(cleanup of '%s' may be needed).",
                     $staging_dir );

          last;
        }
      }
    } ## end foreach my $table (@tables)
  } ## end else [ if ( !$opt_check ) ]

  ##------------------------------------------------------------------##
  ## MOVE                                                             ##
  ##------------------------------------------------------------------##
  
  # Only move the database from the temp directory to live directory if
  # we are not using the update option
  if (!$opt_update) {
    print( '-' x 37, ' MOVE ', '-' x 37, "\n" );

    # Move table files into place in $destination_dir using
    # File::Copy::move(), and remove the staging directory.  We already
    # know that the destination directory does not exist.

    printf( "MOVING '%s' TO '%s'...\n", $staging_dir, $destination_dir );

    if ( !mkdir($destination_dir) ) {
      warn( sprintf( "Failed to create destination directory '%s'.\n" .
                     "Please clean up '%s'.\n",
                   $destination_dir, $staging_dir
          ) );

      $spec->{'status'} =
        sprintf( "FAILED: can not create destination directory '%s' " .
                   "(cleanup of '%s' may be needed)",
               $destination_dir, $staging_dir );
      next TODO;
    }

    move( catfile( $staging_dir, 'db.opt' ), $destination_dir );

    foreach my $table (@tables, @views) {
      my @files =
        glob( catfile( $staging_dir, sprintf( "%s*", $table ) ) );

      printf( "Moving %s...\n", $table );

    FILE:
      foreach my $file (@files) {
        if ( !move( $file, $destination_dir ) ) {
          warn( sprintf( "Failed to move database.\n" .
                         "Please clean up '%s' and '%s'.\n",
                       $staging_dir, $destination_dir
              ) );

          $spec->{'status'} =
            sprintf( "FAILED: can not move '%s' from staging directory " .
                     "(cleanup of '%s' and '%s' may be needed)",
                   $file, $staging_dir, $destination_dir );
          next FILE;
        }
      }
    }

    # Remove the now empty staging directory.
    if ( !rmdir($staging_dir) ) {
      warn( sprintf( "Failed to unlink the staging directory '%s'.\n" .
                     "Clean this up manually.\n",
                   $staging_dir
          ) );

      $spec->{'status'} =
        sprintf( "SUCCESS: cleanup of '%s' may be needed", $staging_dir );
    }
  }

  ##------------------------------------------------------------------##
  ## Symlinks                                                         ##
  ##------------------------------------------------------------------##

  print( '-' x 37, ' SYMLINK ', '-' x 37, "\n" );

  # Create symlink if target directory is not /mysql/data_3306/databases/
  # If given target_dir is not the staging machine MySQL data directory, create symlink for them.
  if ($target_dir ne $mysql_database_dir)
  {
    my $create_symlinks="ln -s ".$target_dir."/".$target_db." ".$mysql_database_dir;
    if ( system($create_symlinks) != 0 ) {

          warn( sprintf( "Failed to create symlink from '%s' to '%s'/'%s'.\n ",
                         $target_dir,$mysql_database_dir,$mysql_database_dir
                ) );

          $spec->{'status'} =
            sprintf( "FAILED: cannot create symlink from '%s' to '%s'/'%s'.\n ",
                     $target_dir,$mysql_database_dir,$mysql_database_dir );
    }
    else {
      printf("Symlink successfully created from '%s' to '%s'/'%s' for the copied databases.\n",$target_dir,$mysql_database_dir,$mysql_database_dir);
    }
  }

  ##------------------------------------------------------------------##
  ## Flush                                                            ##
  ##------------------------------------------------------------------##
 
  print( '-' x 37, ' FLUSH ', '-' x 37, "\n" );

  # Flush tables on target.
  if ($opt_flushtarget) {
    print("FLUSHING TABLES ON TARGET...\n");
    my $tdbh = DBI->connect( $target_dsn,
                               $opt_user_tgt,
                               $opt_password_tgt, {
                                 'PrintError' => 1,
                                 'AutoCommit' => 0
                               } );
    $tdbh->do("use $target_db");
    my $ddl = sprintf('FLUSH TABLES %s', join(q{, }, @tables));
    $tdbh->do($ddl);
    $tdbh->disconnect();
  }
  #------------------------------------------------------------------##
  ## COPYING FUNCTIONS AND PROCEDURES                                ##
  ##-----------------------------------------------------------------##


  if($opt_routines){
    print( '-' x 31, ' Procedures ', '-' x 31, "\n" );
    $target_dbh = DBI->connect( $target_dsn,
                                $opt_user_tgt,
                                $opt_password_tgt, {
                                'PrintError' => 1,
                                'AutoCommit' => 0
                             } );

    $source_dbh = DBI->connect( $source_dsn,
                                $opt_user,
                                $opt_password, {
                                'PrintError' => 1,
                                'AutoCommit' => 0
                            } );

    my $sql_select = "select name, type, returns, body from mysql.proc where db = '$source_db'";
    my $proc_funcs = $source_dbh->selectall_hashref($sql_select, 'name') or die $source_dbh->errstr ;

    foreach my $name (sort keys %{$proc_funcs}){

      my $type = $proc_funcs->{$name}->{type};
      if($type !~ /FUNCTION|PROCEDURE/){
        warn "Copying '$type' not implemted. Skipping....";
        next;
      }

      # Functions must return something, Procedures must not return anything
      my $returns = '';
      if($type eq 'FUNCTION') {
        $returns = "RETURNS $proc_funcs->{$name}->{returns}";
      }

      my $sql = "CREATE $type $target_db.$name()\n $returns\n $proc_funcs->{$name}->{body}";
      print "COPYING $proc_funcs->{$name}->{type} $name\n";
      $target_dbh->do($sql) or die $source_dbh->errstr;
    }
    print "\n\t\tFinished copying functions and procedures\n";

    $target_dbh->disconnect;
    $source_dbh->disconnect;

  }

} ## end foreach my $spec (@todo)

# Display summary.

my $label = '{ SUMMARY }~~';
print( "\n", '~' x ( 80 - length($label) ), $label, "\n" );

foreach my $spec (@todo) {
  printf( "%s -> %s\n  %s\n\n",
          $spec->{'source_db'}, $spec->{'target_db'},
          $spec->{'status'} );
}

print("DONE!\n");

# Provides a way to convert from server names into IP addresses with some error checking in
# there. Usage is
#   my $ip = hostname_to_ip('ens-staging1')
sub hostname_to_ip {
  my ($name) = @_;
  if(! defined $name) {
    die "Cannot convert an undefined name into an IP address. Please check your input";
  }
  my $packed_ip = gethostbyname($name);
  if(! defined $packed_ip) {
    die "Cannot convert '${name}' into an IP address. Please check for typos";
  }
  my $resolved_ip = inet_ntoa($packed_ip);
  return $resolved_ip;
}

END {
  my $seconds = time() - $start_time;

  my $hours = int( $seconds/( 60*60 ) );
  $seconds -= 60*60*$hours;

  my $minutes = int( $seconds/60 );
  $seconds -= 60*$minutes;

  printf( "Time taken: %s%dm%ds\n",
          ( $hours > 0 ? sprintf( "%dh", $hours ) : '' ),
          $minutes, $seconds );
}
