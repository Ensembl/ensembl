#!/usr/bin/perl -w

# A simple OBO file reader/loader

use strict;
use warnings;

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );
use IO::File;

#-----------------------------------------------------------------------

sub usage {
  print("Usage:\n");
  printf( "\t%s\t-h dbhost [-P dbport] \\\n"
      . "\t%s\t-u dbuser [-p dbpass] \\\n"
      . "\t%2\$s\t-d dbname [-t] \\\n"
      . "\t%2\$s\t-f file\n",
    $0, ' ' x length($0) );
  print("\n");
  printf( "\t%s\t-?\n", $0 );
  print("\n");
  print("Arguments:\n");
  print("\t-h/--host dbhost\tDatabase server host name\n");
  print("\t-P/--port dbport\tDatabase server port (optional)\n");
  print("\t-u/--user dbuser\tDatabase user name\n");
  print("\t-p/--pass dbpass\tUser password (optional)\n");
  print("\t-d/--name dbname\tDatabase name\n");
  print(
    "\t-t/--truncate\t\tTruncate (empty) each table\n");
  print("\t\t\t\tbefore writing (optional)\n");
  print("\t-f/--file file\t\tThe OBO file to parse\n");
  print("\t-?/--help\t\tDisplays this information\n");
}

#-----------------------------------------------------------------------

sub write_ontology {
  my ( $dbh, $truncate, $namespaces ) = @_;

  print("Writing to 'ontology' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE ontology") }

  $dbh->do("LOCK TABLE ontology WRITE");

  my $statement = "INSERT INTO ontology (name, namespace) VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $namespace ( keys( %{$namespaces} ) ) {
    my $ontology = $namespaces->{$namespace};

    $sth->bind_param( 1, $ontology,  SQL_VARCHAR );
    $sth->bind_param( 2, $namespace, SQL_VARCHAR );

    $sth->execute();

    $namespaces->{$namespace} =
      { 'id' => ++$count, 'name' => $ontology };
  }

  $dbh->do("OPTIMIZE TABLE ontology");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_ontology

#-----------------------------------------------------------------------

sub write_term {
  my ( $dbh, $truncate, $terms, $namespaces ) = @_;

  print("Writing to 'term' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE term") }

  $dbh->do("LOCK TABLE term WRITE");

  my $statement =
      "INSERT INTO term "
    . "(ontology_id, accession, name, definition) "
    . "VALUES (?,?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $accession ( keys( %{$terms} ) ) {
    my $term = $terms->{$accession};

    $sth->bind_param( 1, $namespaces->{ $term->{'namespace'} }{'id'},
      SQL_INTEGER );
    $sth->bind_param( 2, $accession,            SQL_VARCHAR );
    $sth->bind_param( 3, $term->{'name'},       SQL_VARCHAR );
    $sth->bind_param( 4, $term->{'definition'}, SQL_VARCHAR );

    $sth->execute();

    $term->{'id'} = ++$count;
  }

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_term

#-----------------------------------------------------------------------

sub write_relation_type {
  my ( $dbh, $truncate, $relation_types ) = @_;

  print("Writing to 'relation_type' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE relation_type") }

  $dbh->do("LOCK TABLE relation_type WRITE");

  my $statement = "INSERT INTO relation_type (name) VALUES (?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $relation_type ( keys( %{$relation_types} ) ) {
    $sth->bind_param( 1, $relation_type, SQL_VARCHAR );

    $sth->execute();

    $relation_types->{$relation_type} = { 'id' => ++$count };
  }

  $dbh->do("OPTIMIZE TABLE relation_type");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation_type

#-----------------------------------------------------------------------

sub write_relation {
  my ( $dbh, $truncate, $terms, $relation_types ) = @_;

  print("Writing to 'relation' table...\n");

  if ($truncate) { $dbh->do("TRUNCATE TABLE relation") }

  $dbh->do("LOCK TABLE relation WRITE");

  my $statement =
      "INSERT INTO relation "
    . "(child_term_id, parent_term_id, relation_type_id) "
    . "VALUES (?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;
  foreach my $child_term ( values( %{$terms} ) ) {
    foreach my $relation_type ( keys( %{ $child_term->{'parents'} } ) )
    {
      foreach
        my $parent_acc ( @{ $child_term->{'parents'}{$relation_type} } )
      {
        $sth->bind_param( 1, $child_term->{'id'},         SQL_INTEGER );
        $sth->bind_param( 2, $terms->{$parent_acc}{'id'}, SQL_INTEGER );
        $sth->bind_param( 3, $relation_types->{$relation_type}{'id'},
          SQL_INTEGER );

        $sth->execute();

        ++$count;
      }
    }
  }

  $dbh->do("OPTIMIZE TABLE relation");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation

#-----------------------------------------------------------------------

my ( $dbhost, $dbport );
my ( $dbuser, $dbpass );
my ( $dbname, $truncate, $obo_file_name );

$dbport   = '3306';
$truncate = 0;

if (
  !GetOptions(
    'dbhost|host|h=s' => \$dbhost,
    'dbport|port|P=i' => \$dbport,
    'dbuser|user|u=s' => \$dbuser,
    'dbpass|pass|p=s' => \$dbpass,
    'dbname|name|d=s' => \$dbname,
    'truncate|t'      => \$truncate,
    'file|f=s'        => \$obo_file_name,
    'help|?'          => sub { usage(); exit } )
  || !defined($dbhost)
  || !defined($dbuser)
  || !defined($dbname)
  || !defined($obo_file_name) )
{
  usage();
  exit;
}

my $dsn = sprintf( 'dbi:mysql:database=%s;host=%s;port=%s',
  $dbname, $dbhost, $dbport );

my $dbh =
  DBI->connect( $dsn, $dbuser, $dbpass,
  { 'RaiseError' => 1, 'PrintError' => 2 } );

my $statement =
  "SELECT meta_value FROM meta WHERE meta_key = 'OBO_file_date'";

my $stored_obo_file_date = $dbh->selectall_arrayref($statement)->[0][0];
my $obo_file_date;

my $obo = IO::File->new( $obo_file_name, 'r' ) or die;

my ( $state, $accession, $name, $namespace, $definition, %parents );
my ( %terms, %namespaces, %relation_types );

printf( "Reading OBO file '%s'...\n", $obo_file_name );

my $date_is_checked = 0;

while ( defined( my $line = $obo->getline() ) ) {
  chomp($line);

  if ( !defined($state) ) {
    if ( $line =~ /^\[(\w+)\]$/ ) { $state = $1; next }

    if ( !$date_is_checked && $line =~ /^date: (.+)$/ ) {
      $obo_file_date = sprintf( "%s/%s", $obo_file_name, $1 );

      if ( defined($stored_obo_file_date)
        && $stored_obo_file_date eq $obo_file_date )
      {
        print("This OBO file has already been processed.\n");
        $obo->close();
        $dbh->disconnect();
        exit;
      }

      $date_is_checked = 1;

    }

    next;
  }

  if ( $state eq 'Term' ) {
    if ( $line eq '' ) {
      $namespaces{$namespace} = 'GO';
      $terms{$accession}      = {
        'namespace'  => $namespace,
        'name'       => $name,
        'definition' => $definition
      };

      foreach my $relation_type ( keys(%parents) ) {
        if ( !exists( $relation_types{$relation_type} ) ) {
          $relation_types{$relation_type} = 1;
        }
        $terms{$accession}{'parents'} = {%parents};
      }

      $state = 'clear';
    } elsif ( $line =~ /^(\w+): (.+)$/ ) {
      my $type = $1;
      my $data = $2;

      if    ( $type eq 'id' )        { $accession = $data }
      elsif ( $type eq 'name' )      { $name      = $data }
      elsif ( $type eq 'namespace' ) { $namespace = $data }
      elsif ( $type eq 'def' ) {
        ($definition) = $data =~ /"([^"]+)"/;
      } elsif ( $type eq 'is_a' ) {
        my ($parent_acc) = $data =~ /(\S+)/;
        push( @{ $parents{'is_a'} }, $parent_acc );
      } elsif ( $type eq 'relationship' ) {
        my ( $relation_type, $parent_acc ) = $data =~ /^(\w+) (\S+)/;
        push( @{ $parents{$relation_type} }, $parent_acc );
      } elsif ( $type eq 'is_obsolete' ) {
        if ( $data eq 'true' ) { $state = 'clear' }
      }

    }
  } ## end if ( $state eq 'Term' )

  if ( $state eq 'clear' ) {
    undef($accession);
    undef($name);
    undef($namespace);
    undef($definition);
    %parents = ();

    undef($state);
  }

} ## end while ( defined( my $line...

$obo->close();

print("Finished reading OBO file, now writing to database...\n");

write_ontology( $dbh, $truncate, \%namespaces );
write_term( $dbh, $truncate, \%terms, \%namespaces );
write_relation_type( $dbh, $truncate, \%relation_types );
write_relation( $dbh, $truncate, \%terms, \%relation_types );

print("Updating meta table...\n");

my $sth =
  $dbh->prepare( "DELETE FROM meta "
      . "WHERE meta_key = 'OBO_file_date' "
      . "AND meta_value LIKE ?" );
  $sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ),
    SQL_VARCHAR );
  $sth->execute();
  $sth->finish();

$sth = $dbh->prepare( "INSERT INTO meta (meta_key, meta_value)"
    . "VALUES ('OBO_file_date', ?)" );
$sth->bind_param( 1, $obo_file_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

my $obo_load_date =
  sprintf( "%s/%s", $obo_file_name, scalar( localtime() ) );

$sth =
  $dbh->prepare( "DELETE FROM meta "
    . "WHERE meta_key = 'OBO_load_date' "
    . "AND meta_value LIKE ?" );
$sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ), SQL_VARCHAR );
$sth->execute();
$sth->finish();

$sth = $dbh->prepare( "INSERT INTO meta (meta_key, meta_value)"
    . "VALUES ('OBO_load_date', ?)" );
$sth->bind_param( 1, $obo_load_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

$dbh->disconnect();

print("Done.\n");

# $Id$
