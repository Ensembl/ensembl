#!/usr/local/ensembl/bin/perl -w

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../../../ONTO-PERL-1.31/lib";

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );
use IO::File;
use OBO::Core::Ontology;
use OBO::Core::Term;
use OBO::Core::Relationship;
use OBO::Core::RelationshipType;
use OBO::Core::SynonymTypeDef;
use OBO::Parser::OBOParser;
use OBO::Util::TermSet;

#-----------------------------------------------------------------------

sub usage {
  print("Usage:\n");
  printf( "\t%s\t-h dbhost [-P dbport] \\\n"
            . "\t%s\t-u dbuser [-p dbpass] \\\n"
            . "\t%2\$s\t-d dbname [-t] \\\n"
            . "\t%2\$s\t-f file\\\n"
            . "\t%2\$s\t-o ontology\n",
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
  print("\t-f/--file file\t\tThe OBO file to parse\n");
  print("\t-o/--ontology \t\tOntology name\n");
  print("\t-?/--help\t\tDisplays this information\n");
}

#-----------------------------------------------------------------------

sub write_ontology {
  my ( $dbh, $namespaces ) = @_;

  print("Writing to 'ontology' table...\n");

  my $statement = "INSERT INTO ontology (name, namespace) VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $id;
  my $count = 0;

  local $SIG{ALRM} = sub {
    printf( "\t%d entries, %d to go...\n",
            $count, scalar( keys( %{$namespaces} ) ) - $count );
    alarm(10);
  };
  alarm(10);

  foreach my $namespace ( sort( keys( %{$namespaces} ) ) ) {
    my $ontology = $namespaces->{$namespace};

    $sth->bind_param( 1, $ontology,  SQL_VARCHAR );
    $sth->bind_param( 2, $namespace, SQL_VARCHAR );

    $sth->execute();

    if ( !defined($id) ) {
      $id =
        $dbh->last_insert_id( undef, undef, 'ontology', 'ontology_id' );
    } else {
      ++$id;
    }

    $namespaces->{$namespace} = { 'id' => $id, 'name' => $ontology };

    ++$count;
  }
  alarm(0);

  $dbh->do("OPTIMIZE TABLE ontology");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_ontology

#-----------------------------------------------------------------------

sub write_subset {
  my ( $dbh, $subsets ) = @_;

  print("Writing to 'subset' table...\n");

  $dbh->do("LOCK TABLE subset WRITE");

  my $statement =
    "INSERT INTO subset " . "(name, definition) " . "VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $id;
  my $count = 0;

  local $SIG{ALRM} = sub {
    printf( "\t%d entries, %d to go...\n",
            $count, scalar( keys( %{$subsets} ) ) - $count );
    alarm(10);
  };
  alarm(10);

  foreach my $subset_name ( sort( keys( %{$subsets} ) ) ) {
    my $subset = $subsets->{$subset_name};

    if (!( defined($subset->{'name'}) && defined($subset->{'definition'}) ))
    {
	print "Null value encountered: subset name " . $subset->{'name'}." subset definition ". $subset->{'definition'} . "\n";
	exit;
    }

    $sth->bind_param( 1, $subset->{'name'},       SQL_VARCHAR );
    $sth->bind_param( 2, $subset->{'definition'}, SQL_VARCHAR );

    $sth->execute();

    if ( !defined($id) ) {
      $id = $dbh->last_insert_id( undef, undef, 'subset', 'subset_id' );
    } else {
      ++$id;
    }
    $subset->{'id'} = $id;

    ++$count;
  }
  alarm(0);

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_subset

#-----------------------------------------------------------------------

sub write_term {
  my ( $dbh, $terms, $subsets, $namespaces ) = @_;

  print("Writing to 'term' and 'synonym' tables...\n");

  $dbh->do("LOCK TABLES term WRITE, synonym WRITE");

  my $statement =
      "INSERT INTO term "
    . "(ontology_id, subsets, accession, name, definition) "
    . "VALUES (?,?,?,?,?)";

  my $syn_stmt = "INSERT INTO synonym (term_id, name) VALUES (?,?)";

  my $sth     = $dbh->prepare($statement);
  my $syn_sth = $dbh->prepare($syn_stmt);

  my $id;
  my $count     = 0;
  my $syn_count = 0;

  local $SIG{ALRM} = sub {
    printf( "\t%d entries, %d to go...\n",
            $count, scalar( keys( %{$terms} ) ) - $count );
    alarm(10);
  };
  alarm(10);

  foreach my $accession ( sort( keys( %{$terms} ) ) ) {
    my $term = $terms->{$accession};

    my $term_subsets;

    if ( exists( $term->{'subsets'} ) ) {
      $term_subsets = join( ',',
                            map { $subsets->{$_}{'name'} }
                              @{ $term->{'subsets'} } );
    }

    if ( !( defined($term->{'accession'}) && defined($term->{'name'}) ) )
    {
	print "Null value encountered: term accession " . $term->{'accession'}." term name ". $term->{'name'} . "\n";
	exit;
    }

    $sth->bind_param( 1, $namespaces->{ $term->{'namespace'} }{'id'},
                      SQL_INTEGER );
    $sth->bind_param( 2, $term_subsets,         SQL_VARCHAR );
    $sth->bind_param( 3, $term->{'accession'},  SQL_VARCHAR );
    $sth->bind_param( 4, $term->{'name'},       SQL_VARCHAR );
    $sth->bind_param( 5, $term->{'definition'}, SQL_VARCHAR );

    $sth->execute();

    if ( !defined($id) ) {
      $id = $dbh->last_insert_id( undef, undef, 'term', 'term_id' );
    } else {
      ++$id;
    }
    $term->{'id'} = $id;

    foreach my $syn ( @{ $term->{'synonyms'} } ) {
      $syn_sth->bind_param( 1, $id,  SQL_INTEGER );
      $syn_sth->bind_param( 2, $syn, SQL_VARCHAR );

      $syn_sth->execute();

      ++$syn_count;
    }

    ++$count;
  } ## end foreach my $accession ( sort...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("OPTIMIZE TABLE synonym");
  $dbh->do("UNLOCK TABLES");

  printf(
       "\tWrote %d entries into 'term' and %d entries into 'synonym'\n",
       $count, $syn_count );
} ## end sub write_term

#-----------------------------------------------------------------------

sub write_relation_type {
  my ( $dbh, $relation_types ) = @_;

  print("Writing to 'relation_type' table...\n");

  my $select_stmt =
    "SELECT relation_type_id FROM relation_type WHERE name = ?";
  my $select_sth = $dbh->prepare($select_stmt);

  my $insert_stmt = "INSERT INTO relation_type (name) VALUES (?)";
  my $insert_sth  = $dbh->prepare($insert_stmt);

  my $count = 0;

  local $SIG{ALRM} = sub {
    printf( "\t%d entries, %d to go...\n",
            $count, scalar( keys( %{$relation_types} ) ) - $count );
    alarm(10);
  };
  alarm(10);

  foreach my $relation_type ( sort( keys( %{$relation_types} ) ) ) {
    $select_sth->bind_param( 1, $relation_type, SQL_VARCHAR );
    $select_sth->execute();

    my $id;
    my $found = 0;

    $select_sth->bind_columns( \$id );

    while ( $select_sth->fetch() ) {
      $relation_types->{$relation_type} = { 'id' => $id };
      $found = 1;
    }

    if ( !$found ) {
      $insert_sth->bind_param( 1, $relation_type, SQL_VARCHAR );
      $insert_sth->execute();
      $relation_types->{$relation_type} = {
             'id' =>
               $dbh->last_insert_id( undef,           undef,
                                     'relation_type', 'relation_type_id'
               ) };
      ++$count;
    }
  } ## end foreach my $relation_type (...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE relation_type");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation_type

#-----------------------------------------------------------------------

sub write_relation {
  my ( $dbh, $terms, $relation_types ) = @_;

  print("Writing to 'relation' table...\n");

  $dbh->do("LOCK TABLE relation WRITE");

  my $statement =
      "INSERT INTO relation "
    . "(child_term_id, parent_term_id, relation_type_id) "
    . "VALUES (?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;

  local $SIG{ALRM} = sub {
    printf( "\t%d entries, %d to go...\n",
            $count, scalar( keys( %{$terms} ) ) - $count );
    alarm(10);
  };
  alarm(10);

  foreach my $child_term ( sort { $a->{'id'} <=> $b->{'id'} }
                           values( %{$terms} ) )
  {
    foreach my $relation_type (
                         sort( keys( %{ $child_term->{'parents'} } ) ) )
    {
      foreach my $parent_acc (
                 sort( @{ $child_term->{'parents'}{$relation_type} } ) )
      {
        if ( !defined( $terms->{$parent_acc} ) ) {
          printf( "WARNING: Parent accession '%s' does not exist!\n",
                  $parent_acc );
        } else {

	  if ( !( defined($child_term->{'id'}) && defined($terms->{$parent_acc}{'id'}) && defined($relation_types->{$relation_type}{'id'}) ) )
	  {
	      print "Null value encountered: child term id " . $child_term->{'id'} ." parent term id ". $terms->{$parent_acc}{'id'} . " relationship type " . $relation_types->{$relation_type}{'id'} . "\n";
	      exit;
	  }

          $sth->bind_param( 1, $child_term->{'id'}, SQL_INTEGER );
          $sth->bind_param( 2, $terms->{$parent_acc}{'id'},
                            SQL_INTEGER );
          $sth->bind_param( 3, $relation_types->{$relation_type}{'id'},
                            SQL_INTEGER );
          $sth->execute();
        }
      }
    }

    ++$count;
  } ## end foreach my $child_term ( sort...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE relation");
  $dbh->do("UNLOCK TABLES");

  printf( "\tWrote %d entries\n", $count );
} ## end sub write_relation

#-----------------------------------------------------------------------

my ( $dbhost, $dbport );
my ( $dbuser, $dbpass );
my ( $dbname, $obo_file_name );
my $ontology_name;

$dbport   = '3306';

if ( !GetOptions( 'dbhost|host|h=s' => \$dbhost,
                  'dbport|port|P=i' => \$dbport,
                  'dbuser|user|u=s' => \$dbuser,
                  'dbpass|pass|p=s' => \$dbpass,
                  'dbname|name|d=s' => \$dbname,
                  'file|f=s'        => \$obo_file_name,
                  'ontology|o=s'      => \$ontology_name,
		  'help|?'          => sub { usage(); exit } )
     || !defined($dbhost)
     || !defined($dbuser)
     || !defined($dbname)
     || !defined($obo_file_name)
     || !defined($ontology_name) )
{
  usage();
  exit;
}

#if parsing an EFO obo file delete xref lines - not compatible with OBO:Parser

my @returncode = `grep "default-namespace: efo" $obo_file_name`;
my $returncode = @returncode;
if ($returncode > 0) {
	open (OBO_FILE, "$obo_file_name") or die;
	my $terminator = $/;
	undef $/;
        my $buf = <OBO_FILE>;
	$buf =~ s/\cM\cJ/\n/g;
        $/ = $terminator;
        my @file_lines = split(/\n/, $buf);
	my $no_of_lines = @file_lines;
        my $new_obo_file_name = "new" . $obo_file_name;
        open (NEW_OBO_FILE, ">$new_obo_file_name");
	for (my $i = 0; $i < $no_of_lines; $i++) {
	   my $line = shift(@file_lines);
	   if ($line !~ /^xref/) {
		print NEW_OBO_FILE $line, "\n";
	   }
	}
	close OBO_FILE;
	close NEW_OBO_FILE;
	$obo_file_name = $new_obo_file_name;
 	
}



my $my_parser = OBO::Parser::OBOParser->new;

printf( "Reading OBO file '%s'...\n", $obo_file_name );

my $ontology = $my_parser->work($obo_file_name) or die;

my $dsn = sprintf( 'dbi:mysql:database=%s;host=%s;port=%s',
                   $dbname, $dbhost, $dbport );

my $dbh = DBI->connect( $dsn, $dbuser, $dbpass,
                        { 'RaiseError' => 1, 'PrintError' => 2 } );

my $statement = "SELECT name from ontology where name = ? group by name";

my $select_sth = $dbh->prepare($statement);

$select_sth->bind_param(1,$ontology_name,SQL_VARCHAR);
$select_sth->execute();

if ( $select_sth->fetch() ) {
     print("This ontology name already exists in the database.\nPlease run the program again with a new name.\n");
     $select_sth->finish();
     $dbh->disconnect();
     exit;
}

$statement =
  "SELECT meta_value FROM meta WHERE meta_key = 'OBO_file_date'";

my $stored_obo_file_date = $dbh->selectall_arrayref($statement)->[0][0];

my $obo_file_date = $ontology->date();

$obo_file_date = sprintf( "%s/%s", $obo_file_name, $obo_file_date );

if ( defined($stored_obo_file_date) ) {
    if ( $stored_obo_file_date eq $obo_file_date ){
         print("This OBO file has already been processed.\n");
         $dbh->disconnect();
         exit;
    } elsif ( index( $stored_obo_file_date, $obo_file_name ) != -1)
    {
      print <<EOT;
==> Trying to load a newer (?) OBO file that has already
==> been loaded.  Please clean the database manually of
==> data associated with this file and try again...
EOT
      $dbh->disconnect();
      exit;
    }
}

my ( %terms, %namespaces, %relation_types, %subsets );

my $default_namespace = $ontology->default_namespace();
my @subsets;
if ($OBO::Core::Ontology::VERSION <= 1.31) {
    my $set = $ontology->subset_def_set();
    @subsets = $set->get_set();
} else {
    @subsets = $ontology->subset_def_map()->values();
}

foreach my $subs (@subsets) {
        $subsets{$subs->name()}{'name'}  = $subs->name();
        $subsets{$subs->name()}{'definition'} = $subs->description();
}


# get all non obsolete terms
foreach my $t (@{$ontology->get_terms()}) {
     if (!($t->is_obsolete())) {
        
	my %term;

	my @t_namespace = $t->namespace();
        my $t_namespace_elm = @t_namespace;
	if ($t_namespace_elm > 0){
	  $term{'namespace'} = $t_namespace[0];
	}else {
	  $term{'namespace'} = $default_namespace;
        }
	$namespaces{ $term{'namespace'} } = $ontology_name;	
        
	$term{'accession'} = $t->id();
        $term{'name'} = $t->name();
        $term{'definition'}  = $t->def_as_string();

	my $rels = $ontology->get_relationships_by_source_term($t);
        foreach my $r (@{$rels}) {
		#get parent term
		my $pterm = $r->head();
		push( @{ $term{'parents'}{$r->type()} }, $pterm->id() );        
	}
	
	my @term_subsets = $t->subset();
	foreach my $term_subset (@term_subsets) {
		push( @{ $term{'subsets'} }, $term_subset );
	}      
	my @t_synonyms = $t->synonym_set();
	foreach my $t_synonym (@t_synonyms) {
		push( @{ $term{'synonyms'} }, $t_synonym->def_as_string() );
	}


	$terms{ $term{'accession'} } = {%term};

     }
}

#get all relationship types
foreach my $rel_type (@{$ontology->get_relationship_types()}) {
        my $rel_type_name = $rel_type->id();
	if ( !exists( $relation_types{$rel_type_name} ) ) {
          $relation_types{$rel_type_name} = 1;
        }
}

print("Finished reading OBO file, now writing to database...\n");

write_ontology( $dbh, \%namespaces );
write_subset( $dbh, \%subsets );
write_term( $dbh, \%terms, \%subsets, \%namespaces );
write_relation_type( $dbh, \%relation_types );
write_relation( $dbh, \%terms, \%relation_types );

print("Updating meta table...\n");

my $sth =
  $dbh->prepare(   "DELETE FROM meta "
                 . "WHERE meta_key = 'OBO_file_date' "
                 . "AND meta_value LIKE ?" );
$sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ), SQL_VARCHAR );
$sth->execute();
$sth->finish();

$sth = $dbh->prepare(   "INSERT INTO meta (meta_key, meta_value)"
                      . "VALUES ('OBO_file_date', ?)" );
$sth->bind_param( 1, $obo_file_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

my $obo_load_date =
  sprintf( "%s/%s", $obo_file_name, scalar( localtime() ) );

$sth =
  $dbh->prepare(   "DELETE FROM meta "
                 . "WHERE meta_key = 'OBO_load_date' "
                 . "AND meta_value LIKE ?" );
$sth->bind_param( 1, sprintf( "%s/%%", $obo_file_name ), SQL_VARCHAR );
$sth->execute();
$sth->finish();

$sth = $dbh->prepare(   "INSERT INTO meta (meta_key, meta_value)"
                      . "VALUES ('OBO_load_date', ?)" );
$sth->bind_param( 1, $obo_load_date, SQL_VARCHAR );
$sth->execute();
$sth->finish();

$dbh->disconnect();

print("Done.\n");


# $Id$
