#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::EnsEMBL::Registry;

use Getopt::Long qw( :config no_ignore_case );

#-----------------------------------------------------------------------

sub usage {
  print("Info:\n");
  print <<EOT;
    Script to populate the meta table with species aliases.

    The script reads the already existing aliases from the
    meta table (meta_key 'species.alias') and adds to this
    aliases computed from the species name.  It also uses the
    information stored for the meta_keys species.taxonomy_id,
    species.common_name, species.ensembl_common_name, and
    species.ensembl_alias_name as aliases.

    If the -n or --dryrun options are *not* specified, the existing
    list of aliases is deleted from the meta table and the new list
    is inserted.  In any case, the list of aliases will be displayed
    on the console.

    If the -d or --dbname options are *not* used, the script will
    iterate over all Core databases.  If the -d or --dbname option
    *is* used, only that Core database will be examined.

    This script assumes that the database is a single-species
    database.

    This script does not check for alias duplications between
    species.


EOT

  print("Usage:\n");
  printf( "\t%s\t[-n] -h dbhost [-P dbport] \\\n"
      . "\t%s\t-u dbuser [-p dbpass] \\\n"
      . "\t%2\$s\t[-d dbname]\n",
    $0, ' ' x length($0) );
  print("\n");
  printf( "\t%s\t-?\n", $0 );
  print("\n");
  print("Arguments:\n");
  print("\t-n/--dryrun\t\tDry run, don't write to database\n");
  print("\t-h/--host dbhost\tDatabase server host name\n");
  print("\t-P/--port dbport\tDatabase server port (optional)\n");
  print("\t-u/--user dbuser\tDatabase user name\n");
  print("\t-p/--pass dbpass\tUser password (optional)\n");
  print("\t-d/--name dbname\tDatabase name (optional)\n");
  print("\t-?/--help\t\tDisplays this information\n");
}

#-----------------------------------------------------------------------

my $dryrun;
my ( $dbhost, $dbport );
my ( $dbuser, $dbpass );
my $dbname;

if (
  !GetOptions(
    'dryrun|n'        => \$dryrun,
    'dbhost|host|h=s' => \$dbhost,
    'dbport|port|P=i' => \$dbport,
    'dbuser|user|u=s' => \$dbuser,
    'dbpass|pass|p=s' => \$dbpass,
    'dbname|name|d=s' => \$dbname,
    'help|?'          => sub { usage(); exit } )
  || !defined($dbhost)
  || !defined($dbuser) )
{
  usage();
  exit;
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  '-host' => $dbhost,
  '-port' => $dbport,
  '-user' => $dbuser,
  '-pass' => $dbpass,
);

my $select_stmt = qq(
SELECT DISTINCT LCASE(meta_value)
FROM    meta
WHERE meta_key IN (
  'species.alias', 'species.taxonomy_id', 'species.common_name',
  'species.ensembl_common_name', 'species.ensembl_alias_name'
)
  AND species_id = 1
);

my @dbas = @{ $registry->get_all_DBAdaptors( '-group' => 'Core' ) };

foreach my $dba (@dbas) {
  my $dbh = $dba->dbc()->db_handle();
  if ( defined($dbname) && $dbname ne $dba->dbc()->dbname() ) { next }

  my $species = $dba->species();
  if ( $species =~ /^Ancestral/ ) { next }

  my %aliases;

  my $alias = $species;
  $aliases{$alias} = 1;

  $alias =~ tr [_] [ ];
  $aliases{$alias} = 1;

  $species =~ /^(.)[^_]*_(.*)$/;
  $alias = $1 . $2;
  $aliases{$alias} = 1;

  $species =~ /^(.)[^_]*_(...).*$/;
  $alias = $1 . $2;
  $aliases{$alias} = 1;

  $species =~ /^(...)[^_]*_(...).*$/;
  $alias = $1 . $2;
  $aliases{$alias} = 1;

  my $select_sth = $dbh->prepare($select_stmt);

  $select_sth->execute();

  my $meta_value;

  $select_sth->bind_columns( \$meta_value );

  while ( $select_sth->fetch() ) {
    $aliases{$meta_value} = 1;
  }

  my @aliases =
    sort { length($a) <=> length($b) || $a cmp $b } keys(%aliases);

  my $insert_stmt = sprintf(
    "INSERT IGNORE INTO meta (species_id, meta_key, meta_value) "
      . "VALUES %s",
    join(
      ', ',
      map {
        sprintf( "( 1, 'species.alias', %s )", $dbh->quote( lc($_) ) )
        } @aliases
    ) );

  printf( "Database = %s\n", $dba->dbc()->dbname() );
  printf( "Aliases  = \n\t%s\n", join( "\n\t", @aliases ) );

  if ( !$dryrun ) {
    # Delete old aliases.
    $dbh->do( "DELETE FROM meta WHERE species_id = 1 "
        . "AND meta_key = 'species.alias'" );

    # Insert new aliases.
    $dbh->do($insert_stmt);
  } else {
    print("(not writing to database)\n");
  }

} ## end foreach my $dba (@dbas)
