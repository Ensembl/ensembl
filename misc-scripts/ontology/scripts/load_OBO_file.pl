#!/usr/bin/env perl

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
  my $pro = $0;
  my $off = q{ } x length($pro);
  print <<USAGE;
Usage:
  $pro  -h dbhost [-P dbport]
  $off  -u dbuser [-p dbpass]
  $off  -d dbname
  $off  -f file
  $off  -o ontology
  $off  [-?]

Arguments:
  -h/--host dbhost  Database server host name
  -P/--port dbport  Database server port (optional)
  -u/--user dbuser  Database user name
  -p/--pass dbpass  User password (optional)
  -d/--name dbname  Database name
  -f/--file file    The OBO file to parse
  -o/--ontology     Ontology name
  -?/--help         Displays this information
USAGE
}

#-----------------------------------------------------------------------

sub write_ontology {
  my ($dbh, $namespaces) = @_;

  print("Writing to 'ontology' table...\n");

  my $statement = "INSERT INTO ontology (name, namespace) VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $id;
  my $count = 0;

  local $SIG{ALRM} = sub {
    printf("\t%d entries, %d to go...\n", $count, scalar(keys(%{$namespaces})) - $count);
    alarm(10);
  };
  alarm(10);

  foreach my $namespace (sort(keys(%{$namespaces}))) {
    my $ontology = $namespaces->{$namespace};

    $sth->bind_param(1, $ontology,  SQL_VARCHAR);
    $sth->bind_param(2, $namespace, SQL_VARCHAR);

    $sth->execute();

    if (!defined($id)) {
      $id = $dbh->last_insert_id(undef, undef, 'ontology', 'ontology_id');
    }
    else {
      ++$id;
    }

    $namespaces->{$namespace} = {'id' => $id, 'name' => $ontology};

    ++$count;
  }
  alarm(0);

  #if unknown ontology not found, store it
  my ($unknown_onto_id) = $dbh->selectrow_array("SELECT ontology_id from ontology WHERE name = 'UNKNOWN'");

  if (!defined($unknown_onto_id)) {

    my $statement = "INSERT INTO ontology (name, namespace) VALUES ('UNKNOWN','unknown ontology')";
    my $sth       = $dbh->prepare($statement);
    $sth->execute();
    $unknown_onto_id = $dbh->last_insert_id(undef, undef, 'ontology', 'ontology_id');

  }
  $dbh->do("OPTIMIZE TABLE ontology");

  printf("\tWrote %d entries\n", $count);

  return $unknown_onto_id;
} ## end sub write_ontology

#-----------------------------------------------------------------------

sub write_subset {
  my ($dbh, $subsets) = @_;

  print("Writing to 'subset' table...\n");

  $dbh->do("LOCK TABLE subset WRITE");

  my $statement = "INSERT INTO subset " . "(name, definition) " . "VALUES (?,?)";

  my $sth = $dbh->prepare($statement);

  my $id;
  my $count = 0;

  local $SIG{ALRM} = sub {
    printf("\t%d entries, %d to go...\n", $count, scalar(keys(%{$subsets})) - $count);
    alarm(10);
  };
  alarm(10);

  foreach my $subset_name (sort(keys(%{$subsets}))) {
    my $subset = $subsets->{$subset_name};

    if (!(defined($subset->{'name'}) && defined($subset->{'definition'}))) {
      print "Null value encountered: subset name " . $subset->{'name'} . " subset definition " . $subset->{'definition'} . "\n";
      exit;
    }

    $sth->bind_param(1, $subset->{'name'},       SQL_VARCHAR);
    $sth->bind_param(2, $subset->{'definition'}, SQL_VARCHAR);

    $sth->execute();

    if (!defined($id)) {
      $id = $dbh->last_insert_id(undef, undef, 'subset', 'subset_id');
    }
    else {
      ++$id;
    }
    $subset->{'id'} = $id;

    ++$count;
  }
  alarm(0);

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("UNLOCK TABLES");

  printf("\tWrote %d entries\n", $count);
} ## end sub write_subset

#-----------------------------------------------------------------------

sub write_term {
  my ($dbh, $terms, $subsets, $namespaces, $unknown_onto_id) = @_;

  print("Writing to 'term' and 'synonym' tables...\n");

  $dbh->do("LOCK TABLES term WRITE, synonym WRITE");

  my $statement = "INSERT INTO term " . "(ontology_id, subsets, accession, name, definition) " . "VALUES (?,?,?,?,?)";

  my $syn_stmt = "INSERT INTO synonym (term_id, name) VALUES (?,?)";

  my $existing_term_st = "SELECT term_id, ontology_id FROM term WHERE accession = ?";

  my $update_stmt = "UPDATE term SET ontology_id = ?, subsets = ?, name = ?, definition = ? WHERE term_id = ?";

  my $sth               = $dbh->prepare($statement);
  my $update_sth        = $dbh->prepare($update_stmt);
  my $syn_sth           = $dbh->prepare($syn_stmt);
  my $existing_term_sth = $dbh->prepare($existing_term_st);

  my $id;
  my $count         = 0;
  my $updated_count = 0;
  my $syn_count     = 0;

  local $SIG{ALRM} = sub {
    printf("\t%d entries, %d to go...\n", $count, scalar(keys(%{$terms})) - $count);
    alarm(10);
  };
  alarm(10);

  foreach my $accession (sort(keys(%{$terms}))) {
    my $term = $terms->{$accession};

    my $term_subsets;

    if (exists($term->{'subsets'})) {
      $term_subsets = join(',', map { $subsets->{$_}{'name'} } @{$term->{'subsets'}});
    }

    if (!defined($term->{'name'})) {

      #check if we have this term in the db already
      $existing_term_sth->execute($term->{'accession'});
      my ($existing_term_id, $ontology_id) = $existing_term_sth->fetchrow_array;
      if (defined($existing_term_id)) {
        $term->{'id'} = $existing_term_id;
      }
      else {

        #if not link it to Unknown ontology

        $sth->bind_param(1, $unknown_onto_id,     SQL_INTEGER);
        $sth->bind_param(2, $term_subsets,        SQL_VARCHAR);
        $sth->bind_param(3, $term->{'accession'}, SQL_VARCHAR);
        $sth->bind_param(4, 'UNKNOWN NAME',       SQL_VARCHAR);
        $sth->bind_param(5, 'UNKNOWN DEFINITION', SQL_VARCHAR);

        $sth->execute();

        if (!defined($id)) {
          $id = $dbh->last_insert_id(undef, undef, 'term', 'term_id');
        }
        else {
          ++$id;
        }
        $term->{'id'} = $id;
        ++$count;
      }

    }
    else {

      #we have the term name, check if the term is linked to unknown ontology in the ontology db
      #if so, update it
      $existing_term_sth->execute($term->{'accession'});
      my ($existing_term_id, $ontology_id) = $existing_term_sth->fetchrow_array;
      if (defined($existing_term_id) && $ontology_id == $unknown_onto_id) {

        $update_sth->bind_param(1, $namespaces->{$term->{'namespace'}}{'id'}, SQL_INTEGER);
        $update_sth->bind_param(2, $term_subsets,                             SQL_VARCHAR);
        $update_sth->bind_param(3, $term->{'name'},                           SQL_VARCHAR);
        $update_sth->bind_param(4, $term->{'definition'},                     SQL_VARCHAR);
        $update_sth->bind_param(5, $existing_term_id,                         SQL_INTEGER);

        $update_sth->execute();

        $term->{'id'} = $existing_term_id;
        ++$updated_count;
      }
      else {

        $sth->bind_param(1, $namespaces->{$term->{'namespace'}}{'id'}, SQL_INTEGER);
        $sth->bind_param(2, $term_subsets,                             SQL_VARCHAR);
        $sth->bind_param(3, $term->{'accession'},                      SQL_VARCHAR);
        $sth->bind_param(4, $term->{'name'},                           SQL_VARCHAR);
        $sth->bind_param(5, $term->{'definition'},                     SQL_VARCHAR);

        $sth->execute();

        if (!defined($id)) {
          $id = $dbh->last_insert_id(undef, undef, 'term', 'term_id');
        }
        else {
          ++$id;
        }
        $term->{'id'} = $id;

        ++$count;
      }

      foreach my $syn (@{$term->{'synonyms'}}) {
        $syn_sth->bind_param(1, $id,  SQL_INTEGER);
        $syn_sth->bind_param(2, $syn, SQL_VARCHAR);

        $syn_sth->execute();

        ++$syn_count;
      }
    }
  } ## end foreach my $accession ( sort...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE term");
  $dbh->do("OPTIMIZE TABLE synonym");
  $dbh->do("UNLOCK TABLES");

  printf("\tWrote %d entries into 'term', updated %d entries in 'term' and wrote %d entries into 'synonym'\n", $count, $updated_count, $syn_count);
} ## end sub write_term

#-----------------------------------------------------------------------

sub write_relation_type {
  my ($dbh, $relation_types) = @_;

  print("Writing to 'relation_type' table...\n");

  my $select_stmt = "SELECT relation_type_id FROM relation_type WHERE name = ?";
  my $select_sth  = $dbh->prepare($select_stmt);

  my $insert_stmt = "INSERT INTO relation_type (name) VALUES (?)";
  my $insert_sth  = $dbh->prepare($insert_stmt);

  my $count = 0;

  local $SIG{ALRM} = sub {
    printf("\t%d entries, %d to go...\n", $count, scalar(keys(%{$relation_types})) - $count);
    alarm(10);
  };
  alarm(10);

  foreach my $relation_type (sort(keys(%{$relation_types}))) {
    $select_sth->bind_param(1, $relation_type, SQL_VARCHAR);
    $select_sth->execute();

    my $id;
    my $found = 0;

    $select_sth->bind_columns(\$id);

    while ($select_sth->fetch()) {
      $relation_types->{$relation_type} = {'id' => $id};
      $found = 1;
    }

    if (!$found) {
      $insert_sth->bind_param(1, $relation_type, SQL_VARCHAR);
      $insert_sth->execute();
      $relation_types->{$relation_type} = {'id' => $dbh->last_insert_id(undef, undef, 'relation_type', 'relation_type_id')};
      ++$count;
    }
  } ## end foreach my $relation_type (...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE relation_type");

  printf("\tWrote %d entries\n", $count);
} ## end sub write_relation_type

#-----------------------------------------------------------------------

sub write_relation {
  my ($dbh, $terms, $relation_types) = @_;

  print("Writing to 'relation' table...\n");

  $dbh->do("LOCK TABLE relation WRITE");

  my $statement = "INSERT INTO relation " . "(child_term_id, parent_term_id, relation_type_id, intersection_of) " . "VALUES (?,?,?,?)";

  my $sth = $dbh->prepare($statement);

  my $count = 0;

  local $SIG{ALRM} = sub {
    printf("\t%d entries, %d to go...\n", $count, scalar(keys(%{$terms})) - $count);
    alarm(10);
  };
  alarm(10);

  foreach my $child_term (sort { $a->{'id'} <=> $b->{'id'} } values(%{$terms})) {
    foreach my $relation_type (sort(keys(%{$child_term->{'parents'}}))) {
      foreach my $parent_acc (sort(@{$child_term->{'parents'}{$relation_type}})) {
        if (!defined($terms->{$parent_acc})) {
          printf("WARNING: Parent accession '%s' does not exist!\n", $parent_acc);
        }
        else {

          if (!(defined($child_term->{'id'}) && defined($terms->{$parent_acc}{'id'}) && defined($relation_types->{$relation_type}{'id'}))) {
            print "Null value encountered: child term id " . $child_term->{'id'} . " parent term id " . $terms->{$parent_acc}{'id'} . " relationship type " . $relation_types->{$relation_type}{'id'} . "\n";
            exit;
          }

          $sth->bind_param(1, $child_term->{'id'},                     SQL_INTEGER);
          $sth->bind_param(2, $terms->{$parent_acc}{'id'},             SQL_INTEGER);
          $sth->bind_param(3, $relation_types->{$relation_type}{'id'}, SQL_INTEGER);
          $sth->bind_param(4, 0,                                       SQL_INTEGER);
          $sth->execute();
        }
      }
    }
    foreach my $relation_type (sort(keys(%{$child_term->{'intersection_of_parents'}}))) {
      foreach my $parent_acc (sort(@{$child_term->{'intersection_of_parents'}{$relation_type}})) {
        if (!defined($terms->{$parent_acc})) {
          printf("WARNING: Parent accession '%s' does not exist!\n", $parent_acc);
        }
        else {

          if (!(defined($child_term->{'id'}) && defined($terms->{$parent_acc}{'id'}) && defined($relation_types->{$relation_type}{'id'}))) {
            print "Null value encountered: child term id " . $child_term->{'id'} . " parent term id " . $terms->{$parent_acc}{'id'} . " relationship type " . $relation_types->{$relation_type}{'id'} . "\n";
            exit;
          }

          $sth->bind_param(1, $child_term->{'id'},                     SQL_INTEGER);
          $sth->bind_param(2, $terms->{$parent_acc}{'id'},             SQL_INTEGER);
          $sth->bind_param(3, $relation_types->{$relation_type}{'id'}, SQL_INTEGER);
          $sth->bind_param(4, 1,                                       SQL_INTEGER);
          $sth->execute();
        }
      }
    }
    ++$count;
  } ## end foreach my $child_term ( sort...)
  alarm(0);

  $dbh->do("OPTIMIZE TABLE relation");
  $dbh->do("UNLOCK TABLES");

  printf("\tWrote %d entries\n", $count);
} ## end sub write_relation

#-----------------------------------------------------------------------

my ($dbhost, $dbport);
my ($dbuser, $dbpass);
my ($dbname, $obo_file_name);
my $ontology_name;
my $delete_unknown;

$dbport = '3306';

GetOptions(
  'dbhost|host|h=s' => \$dbhost,
  'dbport|port|P=i' => \$dbport,
  'dbuser|user|u=s' => \$dbuser,
  'dbpass|pass|p=s' => \$dbpass,
  'dbname|name|d=s' => \$dbname,
  'file|f=s'        => \$obo_file_name,
  'ontology|o=s'    => \$ontology_name,
  'delete_unknown'  => \$delete_unknown,
  'help|?'          => sub { usage(); exit }
);

if (!defined($dbhost) || !defined($dbuser) || !defined($dbname)) {
  usage();
  exit;
}

if (!$delete_unknown && (!defined($obo_file_name) || !defined($ontology_name))) {
  usage();
  exit;
}

my $dsn = sprintf('dbi:mysql:database=%s;host=%s;port=%s', $dbname, $dbhost, $dbport);

my $dbh = DBI->connect($dsn, $dbuser, $dbpass, {'RaiseError' => 1, 'PrintError' => 2});

#delete the 'UNKNOWN' ontology and exit
if ($delete_unknown) {

  my ($unknown_onto_id) = $dbh->selectrow_array("SELECT ontology_id from ontology WHERE name = 'UNKNOWN'");

  if (!defined($unknown_onto_id)) {
    print "'UNKNOWN' ontology doesn\'t exist - nothing to delete.\n";
  }
  else {
    my ($unknown_term_count) = $dbh->selectrow_array("select count(1) from term t join ontology o using(ontology_id) where o.name = 'UNKNOWN'");

    if ($unknown_term_count > 0) {
      print("Cannot delete ontology 'UNKNOWN' as $unknown_term_count terms are linked to this ontology.\n");

    }
    else {

      my $delete_unknown_sth = $dbh->prepare("delete from ontology where ontology_id = $unknown_onto_id");
      $delete_unknown_sth->execute();
      $delete_unknown_sth->finish();
      print("0 terms linked to 'UNKNOWN' ontology - ontology deleted.\n");
    }

  }

  $dbh->disconnect();
  exit;
}

#if parsing an EFO obo file delete xref lines - not compatible with OBO:Parser

my @returncode = `grep "default-namespace: efo" $obo_file_name`;
my $returncode = @returncode;
if ($returncode > 0) {
  open(OBO_FILE, "$obo_file_name") or die;
  my $terminator = $/;
  undef $/;
  my $buf = <OBO_FILE>;
  $buf =~ s/\cM\cJ/\n/g;
  $/ = $terminator;
  my @file_lines        = split(/\n/, $buf);
  my $no_of_lines       = @file_lines;
  my $new_obo_file_name = "new" . $obo_file_name;
  open(NEW_OBO_FILE, ">$new_obo_file_name");

  for (my $i = 0 ; $i < $no_of_lines ; $i++) {
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

printf("Reading OBO file '%s'...\n", $obo_file_name);

my $ontology = $my_parser->work($obo_file_name) or die;

my $statement = "SELECT name from ontology where name = ? group by name";

my $select_sth = $dbh->prepare($statement);

$select_sth->bind_param(1, $ontology_name, SQL_VARCHAR);
$select_sth->execute();

if ($select_sth->fetch()) {
  print("This ontology name already exists in the database.\nPlease run the program again with a new name.\n");
  $select_sth->finish();
  $dbh->disconnect();
  exit;
}

$statement = "SELECT meta_value FROM meta WHERE meta_key = 'OBO_file_date'";

my $stored_obo_file_date = $dbh->selectall_arrayref($statement)->[0][0];

my $obo_file_date = $ontology->date();

$obo_file_date = sprintf("%s/%s", $obo_file_name, $obo_file_date);

if (defined($stored_obo_file_date)) {
  if ($stored_obo_file_date eq $obo_file_date) {
    print("This OBO file has already been processed.\n");
    $dbh->disconnect();
    exit;
  }
  elsif (index($stored_obo_file_date, $obo_file_name) != -1) {
    print <<EOT;
==> Trying to load a newer (?) OBO file that has already
==> been loaded.  Please clean the database manually of
==> data associated with this file and try again...
EOT
    $dbh->disconnect();
    exit;
  }
}

my (%terms, %namespaces, %relation_types, %subsets);

my $default_namespace = $ontology->default_namespace();
my @subsets;
if ($OBO::Core::Ontology::VERSION <= 1.31) {
  my $set = $ontology->subset_def_set();
  @subsets = $set->get_set();
}
else {
  @subsets = $ontology->subset_def_map()->values();
}

foreach my $subs (@subsets) {
  $subsets{$subs->name()}{'name'}       = $subs->name();
  $subsets{$subs->name()}{'definition'} = $subs->description();
}

# get all non obsolete terms
foreach my $t (@{$ontology->get_terms()}) {
  if (!($t->is_obsolete())) {

    my %term;

    my @t_namespace     = $t->namespace();
    my $t_namespace_elm = @t_namespace;
    if ($t_namespace_elm > 0) {
      $term{'namespace'} = $t_namespace[0];
    }
    else {
      $term{'namespace'} = $default_namespace;
    }
    $namespaces{$term{'namespace'}} = $ontology_name;

    $term{'accession'}  = $t->id();
    $term{'name'}       = $t->name();
    $term{'definition'} = $t->def_as_string();

    my $rels = $ontology->get_relationships_by_source_term($t);
    foreach my $r (@{$rels}) {

      #get parent term
      my $pterm = $r->head();
      push(@{$term{'parents'}{$r->type()}}, $pterm->id());
    }

    my @intersection_of = $t->intersection_of();
    foreach my $r (@intersection_of) {

      if ($r->isa('OBO::Core::Relationship')) {

        #get parent term for relationship
        my $pterm = $r->head();
        my $type  = $r->type();
        if (!defined($type) || $type eq 'nil') {
          $type = 'is_a';
        }
        push(@{$term{'intersection_of_parents'}{$type}}, $pterm->id());
      }
      else {

        #this is a term
        push(@{$term{'intersection_of_parents'}{'is_a'}}, $r->id());
      }

    }

    my @term_subsets = $t->subset();
    foreach my $term_subset (@term_subsets) {
      push(@{$term{'subsets'}}, $term_subset);
    }
    my @t_synonyms = $t->synonym_set();
    foreach my $t_synonym (@t_synonyms) {
      push(@{$term{'synonyms'}}, $t_synonym->def_as_string());
    }

    $terms{$term{'accession'}} = {%term};

  }
}

#get all relationship types
foreach my $rel_type (@{$ontology->get_relationship_types()}) {
  my $rel_type_name = $rel_type->id();
  if (!exists($relation_types{$rel_type_name})) {
    $relation_types{$rel_type_name} = 1;
  }
}

print("Finished reading OBO file, now writing to database...\n");

my $unknown_onto_id = write_ontology($dbh, \%namespaces);
write_subset($dbh, \%subsets);
write_term($dbh, \%terms, \%subsets, \%namespaces, $unknown_onto_id);
write_relation_type($dbh, \%relation_types);
write_relation($dbh, \%terms, \%relation_types);

print("Updating meta table...\n");

my $sth = $dbh->prepare("DELETE FROM meta " . "WHERE meta_key = 'OBO_file_date' " . "AND meta_value LIKE ?");
$sth->bind_param(1, sprintf("%s/%%", $obo_file_name), SQL_VARCHAR);
$sth->execute();
$sth->finish();

$sth = $dbh->prepare("INSERT INTO meta (meta_key, meta_value)" . "VALUES ('OBO_file_date', ?)");
$sth->bind_param(1, $obo_file_date, SQL_VARCHAR);
$sth->execute();
$sth->finish();

my $obo_load_date = sprintf("%s/%s", $obo_file_name, scalar(localtime()));

$sth = $dbh->prepare("DELETE FROM meta " . "WHERE meta_key = 'OBO_load_date' " . "AND meta_value LIKE ?");
$sth->bind_param(1, sprintf("%s/%%", $obo_file_name), SQL_VARCHAR);
$sth->execute();
$sth->finish();

$sth = $dbh->prepare("INSERT INTO meta (meta_key, meta_value)" . "VALUES ('OBO_load_date', ?)");
$sth->bind_param(1, $obo_load_date, SQL_VARCHAR);
$sth->execute();
$sth->finish();

$dbh->disconnect();

print("Done.\n");

# $Id$
