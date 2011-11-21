#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use File::Spec::Functions;
use Getopt::Long qw( :config no_ignore_case auto_version );
use IO::Dir;

use Data::Dumper;

my $rcsid = '$Revision$';
our ($VERSION) = $rcsid =~ /(\d+\.\d+)/;

sub usage {
  my $indent = ' ' x length($0);

  print <<USAGE_END;
Usage:

  $0 --host=dbhost [ --port=dbport ] \\
  $indent --user=dbuser [ --pass=dbpass ] \\
  $indent --type=schema-type | --database=dbname \\
  $indent --release=new-release [ --from=old-release ] \\
  $indent [ --species=dbspecies ] \\
  $indent [ --cvsdir=/some/path ] \\
  $indent [ --dryrun ] \\
  $indent [ --verbose ] [ --quiet ] \\
  $indent [ --fix ]

  $0 --help | --about

  $0 --version

  --host / -h\tdatabase host name (required)
  --port / -P\tdatabase port (optional, default=3306)
  --user / -u\tdatabase user (required)
  --pass / -p\tdatabase user password (optional, no default)

  --type / -t   restrict to database schema type
                (i.e. core, funcgen, or variation)
                (required if --database is not specified)

  --database / -d   full name of database, or database name pattern
                    (required if --type is not specified)

  --release / -r    release number (required)

  --from / -f       only consider databases from this release
                    (optional, no default)

  --species / -s    restrict to species (optional, no default)

  --cvsdir          the directory where the relevant Ensembl CVS modules
                    have been checked out (optional, default='../../')

  --dryrun / -n     do not actually modify databases
                    (optional, default=not set)

  --verbose / -v    display extra information

  --quiet / -q      do not display warnings

  --fix             also go through all old patches to find any missing
                    patch (patching starts at release equal to the
                    oldest patch in the database)

  --help        display this text
  --about       display further information
  --version     display version and quit

USAGE_END
} ## end sub usage

sub about {
  print <<ABOUT_END;

    This script patches one or several Ensembl databases from older
    releases to the release specified by the user on the command line
    using the --release=NN command line switch.  To only patch databases
    from a particular Ensembl release, the user may use the --from=NN
    command line switch.  In this case, the script will use the value
    from the 'schema_version' meta key or, failing that, from the
    database name, to determine what databases should be or shouldn't be
    patched.

    The script is able to patch databases that have Ensembl Core
    schemas, Ensembl Regulation schemas, and Ensembl Variation schemas,
    provided that the apropriate CVS modules have been checked out
    and are available to this script.  The CVS root directory where
    all Ensembl CVS modules are located may be specified using the
    --cvsdir=/some/path command line switch if the script is unable to
    determine it by itself.

    The user has to specify either a particular database to patch (using
    --database=XXX) or by giving a schema type (using --type=XXX).  In
    the case where a single database is provided, the script will try to
    figure out what schema type the database has from the value of the
    'schema_type' meta key or, failing that, from the database name.

    If the user gives only a schema type, the script will look for
    databases with that schema type, again using the value of the
    'schema_type' meta key or, failing that, using the database names.

    The --database=XXX command line switch may also be used to specify a
    pattern to match database names with.

    To further restrict the set of databases that will be patched,
    the --species=XXX command line switch may be used.  The value
    will be used to match against the various values of the meta key
    'species.alias' or, failing that, against the database name.

    When the --fix command line switch is used, the script will also
    make it possible to apply older patches that might have been
    skipped.  For example: A database at schema version 65, with patches
    since release 40, will not only be patched to the appropriate
    release (specified with --release=NN), but all patches since release
    40 (inclusively) will be tested.  If the patch identifier for an old
    patch is missing from the meta table of the database, that patch
    will be applied.

    Examples

      The release coordinator patches all Ensembl Core-like databases
      from release 65 to release 66:

        $0 -h host -u user -p password \\
          -t core -f 65 -r 66

      A genebuilder patches one of her databases to release 66, and
      wants to look at what the script proposes to do before actually
      running it for real:

        $0 -h host -u user -p password \\
          -d my_database_66 -r 66 --dryrun

      The release coordinator patches all mouse Core-like databases to
      release 66.  She has checked out the 'ensembl' CVS modules in her
      ~/cvs directory:

        $0 -h host -u user -p password \\
          -t core -s mouse -r 66 --cvsdir=~/cvs

      A genebuilder (username 'my') patches all her human databases to
      release 66:

        $0 -h host -u user -p password \\
          -s homo_sapiens -r 66 -d 'my_%'

      A genebuilder makes sure that all patches up to and including
      those for release 66 are included in her database:

        $0 -h host -u user -p password \\
          -r 66 -d my_database --fix

ABOUT_END
} ## end sub about

my ( $opt_host, $opt_port ) = ( undef, '3306' );
my ( $opt_user, $opt_pass ) = ( undef, undef );
my ( $opt_species, $opt_type, $opt_release ) = ( undef, undef, undef );
my $opt_database;

my $opt_cvsdir = '../../';

my $opt_dryrun;
my $opt_from;
my $opt_fix;

my ( $opt_verbose, $opt_quiet );

if ( !GetOptions( 'host|h=s'     => \$opt_host,
                  'port|P=i'     => \$opt_port,
                  'user|u=s'     => \$opt_user,
                  'pass|p=s'     => \$opt_pass,
                  'species|s=s'  => \$opt_species,
                  'type|t=s'     => \$opt_type,
                  'from|f=i'     => \$opt_from,
                  'release|r=i'  => \$opt_release,
                  'database|d=s' => \$opt_database,
                  'cvsdir=s'     => \$opt_cvsdir,
                  'dryrun|n!'    => \$opt_dryrun,
                  'fix!'         => \$opt_fix,
                  'verbose|v!'   => \$opt_verbose,
                  'quiet|q!'     => \$opt_quiet,
                  'help!'        => sub { usage(); exit(0); },
                  'about!'       => sub { about(); exit(0); } ) ||
     !defined($opt_host) ||
     !defined($opt_user) ||
     !defined($opt_release) ||
     ( !defined($opt_database) && !defined($opt_type) ) )
{
  usage();
  exit(1);
}

if ( defined($opt_type) &&
     $opt_type ne 'core' &&
     $opt_type ne 'funcgen' &&
     $opt_type ne 'variation' )
{
  die( sprintf( "Unknown schema type: %s\n", $opt_type ) );
}

my %patches;

# Get available patches.

foreach my $thing ( [ 'ensembl',               'core' ],
                    [ 'ensembl-functgenomics', 'funcgen' ],
                    [ 'ensembl-variation',     'variation' ] )
{
  my $cvs_module  = $thing->[0];
  my $schema_type = $thing->[1];

  if ( defined($opt_type) && $schema_type ne $opt_type ) { next }

  my $sql_dir = canonpath( catdir( $opt_cvsdir, $cvs_module, 'sql' ) );
  my $dh = IO::Dir->new($sql_dir);

  if ( !defined($dh) ) {
    if ( !$opt_quiet ) {
      warn(sprintf( "Unable to find SQL directory '%s'\n", $sql_dir ) );
    }
    next;
  }

  while ( my $file_name = $dh->read() ) {
    if ( $file_name =~ /^patch_\d+_(\d+)_?[a-z]?\.sql$/ ) {
      my $patch_release = $1;

      if ($opt_verbose) {
        printf( "Found %s patch file '%s' for release %d\n",
                $schema_type, $file_name, $patch_release );
      }

      my $full_file_name = catfile( $sql_dir, $file_name );

      push( @{ $patches{$schema_type}{$patch_release} },
            { 'patch' => $file_name, 'path' => $full_file_name } );
    }
  }

} ## end foreach my $thing ( [ 'ensembl'...])

my $dsn = sprintf( "DBI:mysql:host=%s;port=%d", $opt_host, $opt_port );

my $dbh = DBI->connect( $dsn, $opt_user, $opt_pass,
                        { 'RaiseError' => 0, 'PrintError' => 0 } );

# Loop through the databases on the server, patch the ones we want to
# patch and filter out the ones that we don't want to patch.

my $sth;

if ( defined($opt_database) ) {
  $sth = $dbh->prepare("SHOW DATABASES LIKE ?");
  $sth->bind_param( 1, $opt_database, SQL_VARCHAR );
}
else { $sth = $dbh->prepare("SHOW DATABASES") }

$sth->execute();

my $database;
$sth->bind_col( 1, \$database );

DATABASE:
while ( $sth->fetch() ) {
  if ( $database =~ /^(?:information_schema|mysql)$/ ) { next }

  # Figure out schema version, schema type, and species name from the
  # database by querying its meta table.

  my $sth2 = $dbh->prepare(
            sprintf(
              "SELECT meta_key, meta_value FROM %s WHERE meta_key IN " .
                "('schema_version', 'schema_type', " .
                "'species.alias', 'species.common_name', 'patch')",
              $dbh->quote_identifier( undef, $database, 'meta' ) ) );

  $sth2->execute();

  my ( $key, $value );
  $sth2->bind_columns( \( $key, $value ) );

  my ( $schema_version_ok, $schema_type_ok, $species_ok );
  my ( $schema_version,    $schema_type,    $species );
  my %dbpatches;

  while ( $sth2->fetch() ) {
    if ( $key eq 'schema_version' ) {
      $schema_version = $value;
      if ( defined($opt_from) ) {
        if   ( $schema_version eq $opt_from ) { $schema_version_ok = 1 }
        else                                  { $schema_version_ok = 0 }
      }
      else { $schema_version_ok = 1 }
    }
    elsif ( $key eq 'schema_type' ) {
      $schema_type = $value;
      if ( defined($opt_type) ) {
        if   ( $schema_type eq $opt_type ) { $schema_type_ok = 1 }
        else                               { $schema_type_ok = 0 }
      }
      else { $schema_type_ok = 1 }
    }
    elsif ( $key eq 'species.alias' ) {
      if ( defined($opt_species) ) {
        if ( $value eq $opt_species ) { $species_ok = 1 }
      }
      else { $species_ok = 1 }
    }
    elsif ( $key eq 'species.common_name' ) {
      $species = $value;
    }
    elsif ( $key eq 'patch' ) {
      $value =~ /^(patch_\d+_(\d+)_?[a-z]?\.sql)\|(.*)$/;
      my $patch_ident   = $1;
      my $patch_release = $2;
      my $patch_info    = $3;
      $dbpatches{$patch_release}{$patch_ident} = $patch_info;
    }
  } ## end while ( $sth2->fetch() )

  # If we haven't yet found out the schema version, schema type, or
  # species, look to the database name to provide clues.

  if ( !defined($schema_version) ) {
    if ( $database =~ /_(\d+)_\w+$/ ) {
      $schema_version = $1;
      if ( defined($opt_from) ) {
        if   ( $schema_version == $opt_from ) { $schema_version_ok = 1 }
        else                                  { $schema_version_ok = 0 }
      }
      else { $schema_version_ok = 1 }
    }
    elsif ( !$opt_quiet ) {
      warn( sprintf( "Can not determine schema version from '%s'\n",
                     $database ) );
    }
  }
  if ( !defined($schema_type) ) {
    if ( $database =~ /_(core|funcgen|variation)_/ ) {
      $schema_type = $1;
      if ( defined($opt_type) ) {
        if   ( $schema_type eq $opt_type ) { $schema_type_ok = 1 }
        else                               { $schema_type_ok = 0 }
      }
      else { $schema_type_ok = 1 }
    }
    elsif ( !$opt_quiet ) {
      warn( sprintf( "Can not determine schema type from '%s'\n",
                     $database ) );
    }
  }
  if ( !defined($species) ) {
    if ($database =~ /([a-z][a-z_]+[a-z])_(?:core|funcgen|variation)_/ )
    {
      $species = $1;
      if ( defined($opt_species) ) {
        if ( $species eq $opt_species ) { $species_ok = 1 }
      }
      else { $species_ok = 1 }
    }
    elsif ( ( $opt_species && !$opt_quiet ) || $opt_verbose ) {
      warn(
        sprintf( "Can not determine species from '%s'\n", $database ) );
    }
  }

  if ( $schema_version_ok &&
       $schema_type_ok &&
       $species_ok &&
       ( ( !$opt_fix && $schema_version < $opt_release ) ||
         ( $opt_fix && $schema_version <= $opt_release ) ) )
  {
    print( '-' x ( $ENV{COLUMNS} || 80 ), "\n" );
    printf( "Considering '%s' [%s,%s,%d]\n",
            $database, defined($species) ? $species : 'unknown',
            $schema_type, $schema_version );
  }
  else { next }

  # Now figure out what patches we need to apply to this database.

  my $start_version;

  if ($opt_fix) {
    $start_version = ( sort { $a <=> $b } keys %dbpatches )[0];
    printf( "Earliest patch in database '%s' is from release %d\n",
            $database, $start_version );
  }
  else { $start_version = $schema_version + 1 }

  my @apply_these;

  for ( my $r = $start_version; $r <= $opt_release; ++$r ) {
    foreach my $entry ( sort { $a->{'patch'} cmp $b->{'patch'} }
                        @{ $patches{$schema_type}{$r} } )
    {
      my $patch = $entry->{'patch'};
      my $path  = $entry->{'path'};

      if ( exists( $dbpatches{$r}{$patch} ) ) {
        if ($opt_verbose) {
          printf( "Patch '%s' (%s) already applied\n",
                  $patch, $schema_type );
        }
      }
      else {
        if ( !$opt_dryrun ) {
          printf( "Will apply patch '%s' (%s)\n", $patch,
                  $schema_type );
          push( @apply_these, $entry );
        }
        else {
          printf( "Would apply patch '%s' (%s)\n",
                  $patch, $schema_type );
        }
      }

    } ## end foreach my $entry ( sort { ...})
  } ## end for ( my $r = $start_version...)

  if ($opt_dryrun)     { next }
  if ( !@apply_these ) { next }

  local $| = 1;
  print("Proceed with applying these patches? (y/N): ");

  my $yesno = <STDIN>;
  chomp($yesno);

  if ( lc($yesno) =~ /^y(?:es)?$/ ) {
  PATCH:
    foreach my $entry (@apply_these) {
      my $patch = $entry->{'patch'};
      my $path  = $entry->{'path'};

      my @cmd_list = ( 'mysql',
                       "--host=$opt_host",
                       "--user=$opt_user",
                       "--password=$opt_pass",
                       "--database=$database",
                       "--verbose",
                       "--execute=source $path" );

      printf( "Executing the following command:\n%s\n",
              join( ' ', @cmd_list ) );

      if ( system(@cmd_list) ) {
        warn( sprintf( "Failed to apply patch '%s' to database '%s'!\n",
                       $patch, $database ) );

        print("Next patch, next database, or abort? (p/d/A): ");

        my $response = <STDIN>;
        chomp($response);

        if    ( lc($response) =~ /^p(?:atch)$/ )    { next PATCH }
        elsif ( lc($response) =~ /^d(?:atabase)$/ ) { next DATABASE }
        else                                        { exit(1) }
      }
    } ## end foreach my $entry (@apply_these)
  } ## end if ( lc($yesno) =~ /^y(?:es)?$/)

} ## end while ( $sth->fetch() )

$dbh->disconnect();
