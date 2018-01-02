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


=pod 

=head1 NAME

  compare_displays.pl

=head1 SYNOPSIS

  ./compare_displays.pl 
    -h host [-P port] -u user [-p password] -d database\\
    [-ph host -pP port -pu user -pu password -pd database] \\
    [-s species]
    [-v]

=head1 DESCRIPTION

  This script compares the display xrefs between two databases.
  By default, copmares databases on staging with databases on ens-livemirror


head1 OPTIONS

  -dbhost|h     User database server host
  -dbport|P     User database server port
  -dbname|d     User database name
  -dbuser|u     User username
  -dbpass|p     User password
  -species|s    Species to run the comparison on
  -dbphost|ph   Comparison database server host
  -dbpport|pP   Comparison database server port
  -dbpname|pd   Comparison database name
  -dbpuser|pu   Comparison username
  -dbppass|pp   Comparison password


=head1 EXAMPLES

  perl compare_displays.pl -dbhost my_host -dbuser ensro \
    -dbname my_database -species human 

=cut


use strict;
use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Utils::CliHelper;
use Bio::EnsEMBL::ApiVersion;

use Bio::EnsEMBL::Registry;
my $registry = "Bio::EnsEMBL::Registry";


my $dbhost = 'ens-staging';
my $dbport = 3306;
my $dbname;
my $dbuser = 'ensro';
my $dbpass;
my $species = 'homo_sapiens';
my $pdbhost = 'ens-livemirror';
my $pdbport = 3306;
my $pdbname;
my $pdbuser = 'ensro';
my $pdbpass;
my $verbose = 0;


if ( !GetOptions( 'dbhost|h=s'       => \$dbhost,
                  'dbport|P=i'       => \$dbport,
                  'dbuser|u=s'       => \$dbuser,
                  'dbpass|p=s'       => \$dbpass,
                  'dbname|d=s'       => \$dbname,
                  'species=s'        => \$species,
                  'pdbhost|ph=s'     => \$pdbhost,
                  'pdbport|pP=i'     => \$pdbport,
                  'pdbuser|pu=s'     => \$pdbuser,
                  'pdbpass|pp=s'     => \$pdbpass,
                  'pdbdatabase|pd=s' => \$pdbname,
                  'verbose|v!'       => \$verbose)
     || !(
           defined($dbhost)
        && defined($dbuser)
        && defined($pdbhost)
        && defined($pdbuser)
        && defined($species) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script compares the gene displays between two databases for a given species.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password]
  $indent -d database | --pattern pattern \\
  $indent -dp dumppath
  $indent [-dbf dbname override file]
  $indent [-mh host] [-mP port] \\
  $indent [-mu user] [-mp password] [-md database] \\
  $indent [-v]

  -h / --dbhost       User database server host
                      (default is ens-staging)
  -P / --dbport       User database server port
                      (optional, default is 3306)
  -u / --dbuser       User username
                      (no write-access required, default is 'ensro')
  -p / --dbpass       User password
                      (optional, default is undefined)
  -d / --dbname       User database name
                      (optional, can be inferred from the species name)

  -s / species        species for which to run the comparison
  
  -ph / --pdbhost     Previous database server host
                      (optional, default is 'ens-livemirror')
  -pP / --pdbport     Previous database server port
                      (optional, default is 3306)
  -pu / --pdbuser     Previous database username (no write-access required)
                      (optional, default is 'ensro')
  -pp / --pdbpass     Previous database password
                      (optional, default is undefined)
  -pd / --dbname      Previous database name
                      (optional, can be inferred from the species name)

  -v / --verbose      Be verbose, display every SQL statement as they
                      are executed (on standard error)

USAGE_END

  die(   "Need the following options: "
       . "-h -u -ph -pu and -s\n" );

}

my ($dba, $pdba, $previous_display, $previous_db, $new_display, $new_db);
my ($new_gene, $stable_sql, $display_sql, %sources);
my $lost_ids = 0;
my $changed_display = 0;
my $changed_db = 0;
my $changed = 0;

my $version = software_version() - 1;

if ($dbname) {
      $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
         '-host'     => $dbhost,
         '-user'     => $dbuser,
         '-pass'     => $dbpass,
         '-dbname'   => $dbname,
         '-species'  => $species,
         '-group'    => 'core',
       );
} else {
      $registry->load_registry_from_db(
        '-host'       => $dbhost,
        '-user'       => $dbuser,
        '-pass'       => $dbpass,
      );
      $dba = $registry->get_DBAdaptor($species, 'core');
}
Bio::EnsEMBL::Registry->clear();
print "Connected to database " . $dba->dbc()->dbname . " on " . $dba->dbc()->host() . "\n";

if ($pdbname) {
      $pdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
         '-host'     => $pdbhost,
         '-user'     => $pdbuser,
         '-pass'     => $pdbpass,
         '-dbname'   => $pdbname,
         '-species'  => $species,
         '-group'    => 'core',
       );
} else {
      Bio::EnsEMBL::Registry->load_registry_from_db(
        '-host'       => $pdbhost,
        '-user'       => $pdbuser,
        '-pass'       => $pdbpass,
        '-db_version' => $version,
      );
      $pdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
}
print "Connected to previous database " . $pdba->dbc()->dbname . " on " . $pdba->dbc->host . "\n";


my $previous_helper = $pdba->dbc()->sql_helper();
my $new_helper = $dba->dbc()->sql_helper();
my $gene_sql = 'SELECT stable_id FROM gene';
my @genes = @{ $previous_helper->execute_simple(-SQL => $gene_sql) };
my $nb_genes = scalar(@genes);

foreach my $gene (@genes) {
      $stable_sql = 'SELECT count(*) FROM gene WHERE stable_id = ?';
      $new_gene = $new_helper->execute_single_result(-SQL => $stable_sql, -PARAMS => [$gene]);
      if ($new_gene == 0) {
          print STDOUT "Stable id $gene does not exist in new database " . $dba->dbc()->dbname() . "\n";
          $lost_ids++;
          next;
      }

      $display_sql = 'SELECT display_label, db_name FROM gene LEFT JOIN xref ON display_xref_id = xref_id LEFT JOIN external_db USING(external_db_id) WHERE stable_id = ?';
      ($previous_display, $previous_db) = @{ $previous_helper->execute(-SQL => $display_sql, -PARAMS => [$gene])->[0] };
      ($new_display, $new_db) = @{ $new_helper->execute(-SQL => $display_sql, -PARAMS => [$gene])->[0] };

      if (!$new_display && $previous_display) {
          print STDOUT "New display for $gene: display $previous_display from $previous_db changed to null\n";
          $sources{$previous_db . " to null"}++;
          $changed++;
          next;
      }

      if (!$previous_display && $new_display) {
          print STDOUT "New display for $gene: display null changed to $new_display from $new_db\n";
          $sources{"null to $new_db"}++;
          $changed++;
          next;
      }

      if ($previous_db eq $new_db && $previous_display eq $new_display) {
          next;
      }
     
      if ($previous_db eq $new_db) {
          print STDOUT "Same external_db for $gene: external_db $previous_db changed from $previous_display to $new_display\n";
          $changed_display++;
      }

      elsif ($previous_display eq $new_display) {
         print STDOUT "Same display for $gene: display $previous_display changed from $previous_db to $new_db\n";
         $changed_db++;
      }

      else {
         print STDOUT "New display for $gene: display $previous_display from $previous_db changed to $new_display from $new_db\n" ;
         $sources{$previous_db . " to " . $new_db}++;
         $changed++;
      }

}

my $total = $changed_display + $changed_db + $changed;

print STDOUT "\n\nSUMMARY\n\n";

foreach my $k (keys %sources) {
    print $sources{$k} . " genes changed from $k\n" ;
}

print STDOUT "$lost_ids (" . sprintf("%.2f", $lost_ids/$nb_genes*100) . "%) stable ids were lost\n";
print STDOUT "$changed_display (" . sprintf("%.2f", $changed_display/$nb_genes*100) . "%) genes have the same external_db but different displays\n";
print STDOUT "$changed_db (" . sprintf("%.2f", $changed_db/$nb_genes*100) . "%) genes have the same display but different external_dbs\n";
print STDOUT "$changed (" . sprintf("%.2f", $changed/$nb_genes*100) . "%) genes have been assigned a fully new display\n";
print STDOUT "In total, $total (" . sprintf("%.2f", $total/$nb_genes*100) . "%) genes changed\n";




