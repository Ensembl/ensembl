#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use File::Spec::Functions qw/:ALL/;
use Getopt::Long qw( :config no_ignore_case auto_version );
use IO::Dir;

my $rcsid = '$Revision$';
our ($VERSION) = $rcsid =~ /(\d+\.\d+)/;

sub usage {
  my $indent = ' ' x length($0);

  print <<USAGE_END;
Usage:

  $0 --host=dbhost [ --port=dbport ] \\
  $indent --user=dbuser [ --pass=dbpass ] \\
  $indent --type=schema-type | --database=dbname \\
  $indent [ --release=new-release ] [ --from=old-release ] \\
  $indent [ --species=dbspecies ] \\
  $indent [ --cvsdir=/some/path ] \\
  $indent [ --dryrun ] \\
  $indent [ --interactive 0|1 ]\\
  $indent [ --verbose ] [ --quiet ] \\
  $indent [ --mysql=optional_path ]  \\
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

  --release / -r    release number (optional, default is the latest
                    release that we can find patches for)

  --from / -f       only consider databases from this release
                    (optional, no default)

  --species / -s    restrict to species (optional, no default)

  --cvsdir          the directory where the relevant Ensembl CVS modules
                    have been checked out (optional, default=misc-scripts/../..)

  --dryrun / -n     do not actually modify databases
                    (optional, default=not set)

  --verbose / -v    display extra information

  --quiet / -q      do not display warnings

  --fix             also go through all old patches to find any missing
                    patch (patching starts at release equal to the
                    oldest patch in the database) >>USE WITH CAUTION<<
  
  --mysql           specify the location of the mysql binary if it is not on
                    \$PATH. Otherwise we default this to mysql
  
  --nointeractive   specify if you want an non-interactive patching environment
                    (default false). >>USE WITH CAUTION<<

  --help        display this text
  --about       display further information
  --version     display version and quit

USAGE_END
} ## end sub usage

sub about {
  print <<ABOUT_END;

    This script patches one or several Ensembl databases from older
    releases to the release specified by the user on the command line
    using the --release=NN command line switch, or to the latest release
    for which the script is able to find a patch if the --release=NN
    switch is not used.  To only patch databases from a particular
    Ensembl release, the user may use the --from=NN command line
    switch.  In this case, the script will use the value from the
    'schema_version' meta key or, failing that, from the database name,
    to determine what databases should be or shouldn't be patched.

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
          
      A genebuilder wishes to patch the same set as specified above but
      without being prompted to apply patches
      
        $0 -h host -u user -p password \\
          -t core -f 65 -r 66 --nointeractive

      A genebuilder patches one of her databases to release 66, and
      wants to look at what the script proposes to do before actually
      running it for real:

        $0 -h host -u user -p password \\
          -d my_database_66 -r 66 --dryrun

      The release coordinator patches all mouse Core-like databases to
      the latest release.  She has checked out the 'ensembl' CVS modules
      in her ~/cvs directory:

        $0 -h host -u user -p password \\
          -t core -s mouse --cvsdir=~/cvs

      A genebuilder (username 'my') patches all her human databases to
      the latest release.

        $0 -h host -u user -p password \\
          -s homo_sapiens -d 'my_%'

      A genebuilder makes sure that all patches up to and including
      those for release 66 are included in her database, without
      actually applying any patch (any missing patches needs to be
      manually checked!):

        $0 -h host -u user -p password \\
          -r 66 -d my_database --fix --dryrun

ABOUT_END
} ## end sub about

my ( $opt_host, $opt_port ) = ( undef, '3306' );
my ( $opt_user, $opt_pass ) = ( undef, undef );
my ( $opt_species, $opt_type, $opt_release ) = ( undef, undef, undef );
my $opt_database;

my $opt_cvsdir;

my $opt_dryrun;
my $opt_from;
my $opt_fix;
my $opt_mysql = 'mysql';
my $opt_interactive = 1;

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
                  'mysql=s'      => \$opt_mysql,
                  'interactive|i!' => \$opt_interactive,
                  'verbose|v!'   => \$opt_verbose,
                  'quiet|q!'     => \$opt_quiet,
                  'help!'        => sub { usage(); exit(0); },
                  'about!'       => sub { about(); exit(0); } ) ||
     !defined($opt_host) ||
     !defined($opt_user) ||
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

my $latest_release;
my %patches;

# Get available patches.

foreach my $thing ( [ 'ensembl',               'core', 'table.sql' ],
                    [ 'ensembl-functgenomics', 'funcgen', 'efg.sql' ],
                    [ 'ensembl-variation',     'variation', 'table.sql' ] )
{
  my ($cvs_module, $schema_type, $schema_file) = @{$thing};

  if ( defined($opt_type) && $schema_type ne $opt_type ) { next }

  my $sql_dir = _sql_dir($cvs_module, $schema_file);
  if(! defined $sql_dir) {
    if ( !$opt_quiet ) {
      warn(sprintf("No SQL directory found for CVS module %s\n", $cvs_module));
    }
    next;
  }
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

      if ( !defined($latest_release) ||
           $latest_release < $patch_release )
      {
        $latest_release = $patch_release;
      }

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

if ( defined($opt_release) && $opt_release > $latest_release ) {
  die( sprintf( "Release %d is too new, " .
                  "last release with patches is release %d\n",
                $opt_release, $latest_release ) );
}

if ( !defined($opt_release) ) {
  if ($opt_verbose) {
    printf( "Latest release with patches is release %d\n",
            $latest_release );
  }
  $opt_release = $latest_release;
}

my $dsn = sprintf( "DBI:mysql:host=%s;port=%d", $opt_host, $opt_port );

my $dbh = DBI->connect( $dsn, $opt_user, $opt_pass,
                        { 'RaiseError' => 0, 'PrintError' => 0 } );

if(! $dbh) {
  my $pass = ($opt_pass) ? 'with a' : 'with no';
  warn(sprintf(q{Cannot connect to DSN '%s' with user %s %s password. Check your settings}, $dsn, $opt_user, $pass));
  exit 1;  
}

# Loop through the databases on the server, patch the ones we want to
# patch and filter out the ones that we don't want to patch.

my $sth;
my $found_databases = 0;

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
      if($value =~ /^(patch_\d+_(\d+)_?[a-z]?\.sql)\|(.*)$/) {
        my $patch_ident   = $1;
        my $patch_release = $2;
        my $patch_info    = $3;
        $dbpatches{$patch_release}{$patch_ident} = $patch_info;
      }
      else {
        warn "The patch value $value from database $database does not conform to the pattern of 'patch_from_to_tag|description'. Please fix";
      }
    }
  } ## end while ( $sth2->fetch() )

  # If we haven't yet found out the schema version, schema type, or
  # species, look to the database name to provide clues.


  if ( ! $schema_version ) {
	#remove defined as version maybe empty string

    if ( $database =~ /_(\d+)_\w+$/ ) {

      $schema_version = $1;

      if ( defined($opt_from) ) {
        if   ( $schema_version == $opt_from ) { $schema_version_ok = 1 }
        else                                  { $schema_version_ok = 0 }
      }
      else {
		$schema_version_ok = 1 }
    }
    elsif ( ! $opt_quiet ) {
	  $schema_version_ok = 0;
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
    elsif ( $opt_species && !$opt_quiet ) {
      warn(
        sprintf( "Can not determine species from '%s'\n", $database ) );
    }
  }
  
  if ( $schema_version_ok &&
       $schema_type_ok &&
       ( !defined($opt_species) ||
         ( defined($opt_species) && $species_ok ) ) &&
       ( ( !$opt_fix && $schema_version < $opt_release ) ||
         ( $opt_fix && $schema_version <= $opt_release ) ) )
  {
    $found_databases = 1;
    print( '-' x ( $ENV{COLUMNS} || 80 ), "\n" );
    printf( "Considering '%s' [%s,%s,%d]\n",
            $database, defined($species) ? $species : 'unknown',
            $schema_type, $schema_version );
  }
  else { 
    if($opt_verbose) {
      printf("Skipping database %s (type: %s | version: %d)\n", $database, ($schema_type||'-'), ($schema_version || 0));
      if( $schema_type_ok && $schema_type eq $opt_type && $schema_version_ok && $schema_version == $opt_release) {
        
        my $release_patches = join(q{, }, sort map { $_->{patch} } @{$patches{$schema_type}{$schema_version}});
        my $db_patches = join(q{, }, sort keys %{$dbpatches{$schema_version}});
        
        if($release_patches ne $db_patches) {
          printf("\t%s patches [%s] are not the same as release %i patches [%s]; rerun with --fix and --dryrun\n", 
            $database, $db_patches, $opt_release, $release_patches);
        }
      }
      if($schema_type_ok && ! exists $patches{$schema_type}) {
        printf("\t%s patches could not be found. Check your --cvsdir option and try again\n", $schema_type);
      }
    }
    next; 
  }

  # Now figure out what patches we need to apply to this database.

  my $start_version;

  if ($opt_fix) {
    $start_version = ( sort { $a <=> $b } keys %dbpatches )[0];
    if ( !defined($start_version) ) {
      warn( sprintf( "No patches in database, " .
                       "beginning fix from release %d\n",
                     $schema_version ) );
      $start_version = $schema_version;
    }
    else {
      printf( "Earliest patch in database '%s' is from release %d\n",
              $database, $start_version );
    }
  }
  else { $start_version = $schema_version + 1 }

  my @apply_these;
  my $schema_version_warning = 0;
  
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

          if ( $r < $opt_release && $patch =~ /a\.sql$/ ) {
            # Warn about possible setting schema_version with an 'a'
            # patch.
            $schema_version_warning = 1;
          }
        }
        else {
          printf( "Would apply patch '%s' (%s)\n",
                  $patch, $schema_type );
        }
      }
    } ## end foreach my $entry ( sort { ...})

  } ## end for ( my $r = $start_version...)

  if ( $opt_dryrun || !@apply_these ) { print("\n"); next }

  my $apply_patches;
  local $| = 1;
  if($opt_interactive) {
    print("Proceed with applying these patches? (y/N): ");
    my $yesno = <STDIN>;
    chomp($yesno);
    $apply_patches = (lc($yesno) =~ /^y(?:es)?$/) ? 1 : 0;
  }
  else {
    $apply_patches = 1;
    print "Enterning non-interative mode. Will apply patches\n";
  }

  if ( $apply_patches ) {
  PATCH:
    foreach my $entry (@apply_these) {
      my $patch = $entry->{'patch'};
      my $path  = $entry->{'path'};

      my @cmd_list = (  $opt_mysql,
                        "--host=$opt_host",
                        "--user=$opt_user");
      push(@cmd_list,   "--password=$opt_pass") if $opt_pass;
      push(@cmd_list,   "--port=$opt_port",
                        "--database=$database",
                        "--verbose",
                        "--execute=source $path" );

      printf( "Executing the following command:\n%s\n",
              join( ' ', @cmd_list ) );

      if ( system(@cmd_list) ) {
        warn( sprintf( "Failed to apply patch '%s' to database '%s'!\n",
                       $patch, $database ) );

        if(!$opt_interactive) {
          warn('In non-interative mode; aborting current run');
          exit(1);
        }
        print("Next patch, next database, or abort? (p/d/A): ");

        my $response = <STDIN>;
        chomp($response);

        if    ( lc($response) =~ /^p(?:atch)?$/ )    { next PATCH }
        elsif ( lc($response) =~ /^d(?:atabase)?$/ ) { next DATABASE }
        else                                         { exit(1) }
      }
    } ## end foreach my $entry (@apply_these)

    if ( !$opt_quiet && $schema_version_warning ) {
      warn( "Applied one or several 'a' patches, " .
            "schema_version might have been updated\n" );
    }

  } ## end if ( lc($yesno) =~ /^y(?:es)?$/)

  print("\n");

} ## end while ( $sth->fetch() )

if(!$found_databases) {
  printf(('-'x80)."\n");
  printf("No databases considered. Check your --database/--type/--release flags\n");
  printf(('-'x80)."\n");
}

$dbh->disconnect();

sub _sql_dir {
  my ($cvs_module, $schema_file) = @_;
  my $cvs_dir;
  if($opt_cvsdir) {
    $cvs_dir = $opt_cvsdir;
  }
  else {
    my ($volume, $directories, $file) = splitpath(__FILE__);
    $directories = curdir() unless $directories;
    $cvs_dir = catdir($directories, updir(), updir());
  }
  my $sql_dir = rel2abs(canonpath( catdir( $cvs_dir, $cvs_module, 'sql' ) ));
  my $schema_location = catfile($sql_dir, $schema_file);
  if(! -f $schema_location) {
    if($opt_verbose) {
      printf("Could not find the schema file '%s' for E! module %s", $schema_location, $cvs_module);
      printf("\tTry using --cvsdir if your checkouts are in a non-standard location\n") if $opt_cvsdir;
    }
    return;
  }
  printf("Using '%s' as our SQL directory\n", $sql_dir);
  return $sql_dir;
}
