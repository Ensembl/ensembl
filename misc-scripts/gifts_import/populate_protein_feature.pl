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

#The script populates the CORE db's protein_feature table with all the uniprot mappings found in GIFTS db's ensp_u_cigar table.


use strict;
use warnings;
use DBI qw( :sql_types );
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my ( $gfhost, $gfport, $gfuser, $gfpass, $gfdb );
my ( $pfhost, $pfport, $pfuser, $pfpass, $pfdb );

my $create_index = 0;
my $dryrun = 0;
#my $limit = 10;
my $limit = 0;


GetOptions(
  "gfhost=s"    => \$gfhost,
  "gfport=i"       => \$gfport,
  "gfuser=s"    => \$gfuser,
  "gfpass=s"    => \$gfpass,
  "gfdb=s"   => \$gfdb,

  "pfhost=s"      => \$pfhost,
  "pfport=i"      => \$pfport,
  "pfuser=s"      => \$pfuser,
  "pfpass=s"        => \$pfpass,
  "pfdb=s"   => \$pfdb,

  "dryrun!"         => \$dryrun,
  "help", \&usage,

);

usage() if ( !defined $gfhost || !defined $gfport || !defined $gfuser || !defined $gfdb );
usage() if ( !defined $pfhost  || !defined $pfport || !defined $pfuser || !defined $pfpass );


my ( $giftsDB, $protfeatureDB ) = init_db();
die unless defined($giftsDB) and defined($protfeatureDB);

print "giftsDB credentials:\n";
print Dumper($giftsDB);

print "protfeatureDB credentials:\n";
print Dumper($protfeatureDB);


my $registry = get_registry();
my $translation_adaptor = $registry->get_adaptor( 'human', 'core', 'translation' );


sub init_db {

# This is where the gifts info will be read from
  my $giftsDB = {
    user => $gfuser,
    pass => $gfpass,
    host => $gfhost,
    port => $gfport,
    dbname => $gfdb

  };

# This is where the gifts info will be written to
  my $protfeatureDB = {
    user   => $pfuser,
    pass   => $pfpass,
    host   => $pfhost,
    port   => $pfport,
    dbname => $pfdb
  };

  return ( $giftsDB, $protfeatureDB );
}


populate_protein_feature_db();

sub get_analysis_id {
  my ($logic_name) = @_;
  my $analysis_id;
  
  my $dba = get_dba($protfeatureDB);
  my $analysis_dba = $dba->get_AnalysisAdaptor();
  my $analysis = $analysis_dba->fetch_by_logic_name($logic_name);
  
  $analysis_id = $analysis->dbID() if defined $analysis;
  return $analysis_id;
}


sub populate_protein_feature_db {

  my $gifts_dbh = db_connect($giftsDB );
  my $protfeature_dbh = db_connect($protfeatureDB);
  
  if($gifts_dbh){
    print "Got connection to GiftsDB \n";
  }
  
  if($protfeature_dbh){
    print "Got connection to ProtfeatureDB \n";
  }
  
  my $insert_container = build_insert_sql($gifts_dbh, $protfeature_dbh);

  
  my $rows_inserted = 0;
  my @failed_stmts = [];
  unless($dryrun){

  foreach my $insert_stmt (@$insert_container) {
    eval{
    my $affected = $protfeature_dbh->do($insert_stmt);
    $rows_inserted += $affected;
    };
    if($@){
      print "$@\n";
      push(@failed_stmts, $@);
    }
    
  }
  print Dumper(\@failed_stmts);
  print "INSERT END\n";
  }

  print "Number of rows inserted into protein_feature $protfeatureDB->{'dbname'}.protein_feature ", $rows_inserted, "\n";
  return $rows_inserted;
  
  $gifts_dbh->disconnect();
  $protfeature_dbh->disconnect();
  
}

sub annotate_ref{
  my $ref = shift;
  my $ensp_id = $ref->[0];
  my $ensp = $translation_adaptor->fetch_by_stable_id( $ensp_id ) ;
  if($ensp){
  $ref->[0] = $ensp->dbID;
  $ref->[1] = 1;
  $ref->[2] = $ensp->length()-1;
  
  }
  
  return $ref;
}

sub build_insert_sql {
  my ($gifts_dbh, $protfeature_dbh) = @_;
  
  #get analysis id
  my $analysis_id = get_analysis_id("gifts_import");
  print("Analysis id $analysis_id \n");
  
  my $select_sql = qq{SELECT DISTINCT ensp_id, 0, 0, uniprot_acc, $analysis_id as analysis_id, mdz, 'mdtag' FROM ensp_u_cigar};
  $select_sql = $limit ? $select_sql . " limit 10" : $select_sql;
  
  my $insert_sql = qq{INSERT INTO $protfeatureDB->{'dbname'}.protein_feature(translation_id, seq_start, seq_end, hit_name, analysis_id, cigar_line, align_type) VALUES };


  my $import_sth = $gifts_dbh->prepare($select_sql);
  $import_sth->execute;
  my $rows = $import_sth->rows;

  my $start = 0;
  my $insert_values;
  my $max_rows_in_insert = 9999;    #Batch size of 10,000
  # my $max_rows_in_insert = 10;    #Batch size of 10
  my $rows_in_insert     = 0;

  my @insert_container = ();
  my $miss = 1;
  while ( my $ref = $import_sth->fetchrow_arrayref ) {
     my $not_ensp =1;
     if($ref->[0] !~ /^ENSP.*/){
      warn "\t\t$not_ensp ensp id not in right format for ", $ref->[0], " Uniprot id ", $ref->[3], "\n";
      $not_ensp++;
      next;
    }
    $ref = annotate_ref($ref);
    if($ref->[0] =~ /^ENSP.*/){
      warn "$miss No translation record found in core translation for ", $ref->[0], " Uniprot id ", $ref->[3], "\n";
      $miss++;
      next;
    }
    
    if ( $rows_in_insert < $max_rows_in_insert ) {
      $insert_values .= ',' if $start++;
      $insert_values .= '(' . ( join( ",", map { $gifts_dbh->quote($_) } @{$ref} ) ) . ')';
      $rows_in_insert++;
    }
    else {
      $insert_values .= ',' if $start++;
      $insert_values .= '(' . ( join( ",", map { $gifts_dbh->quote($_) } @{$ref} ) ) . ')';
      my $insert_sql_with_values = $insert_sql . $insert_values;
      push( @insert_container, $insert_sql_with_values );

      $insert_values  = undef;
      $rows_in_insert = 0;
      $start          = 0;
  }
}

  if ($insert_values) {
    my $insert_sql_with_values = $insert_sql . $insert_values;
    push( @insert_container, $insert_sql_with_values );
  }
  
  print "Number of rows fetched from GIFTS $rows \n";
  return \@insert_container;
 

}

#===========DB Connection routines=============
sub db_connect {
  my ($connectDB) = @_;
  my $dbname = $connectDB->{'dbname'};
  my $host = $connectDB->{'host'};
  my $user = $connectDB->{'user'};
  my $pass = $connectDB->{'pass'};
  my $port = $connectDB->{'port'};

  my $dsn = "DBI:mysql:host=$host;";
  if ($port) {
    $dsn .= "port=$port;";
  }
  $dsn .= "database=$dbname";

  my $dbh = DBI->connect( $dsn, $user, $pass,
  { 'PrintError' => 1, 'RaiseError' => 1 } );

  if ( !$dbh ) {
    die "ERROR: $DBI::errstr";
  }

  return $dbh;
}

sub get_dba{
  my ($connectDB) = @_;

  my $dbname = $connectDB->{'dbname'};
  my $host = $connectDB->{'host'};
  my $user = $connectDB->{'user'};
  my $pass = $connectDB->{'pass'};
  my $port = $connectDB->{'port'};

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host    => $host,
    -user    => $user,
    -dbname  => $dbname,
    -pass    => $pass,
    -port    => $port
  );
return $dba;
}

sub get_registry {

  my $registry = "Bio::EnsEMBL::Registry";

  $registry->load_registry_from_db(
  -host       => $protfeatureDB->{'host'},
  -user       => $protfeatureDB->{'user'},
  -pass       => $protfeatureDB->{'pass'},
  -port => $protfeatureDB->{'port'}
);

  return $registry;
}



sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);

The script populates the CORE db's protein_feature table with all the uniprot mappings found in GIFTS db's ensp_u_cigar table.


Options -gfhost -gfport -gfuser -gfpass -gfdb are mandatory and specify the credentials for the server on which a GIFTS database exists
Options -pfhost -pfport -pfuser -pfpass -pfdb are credentials of the core db where the GIFTS data will be copied into

To run the script cd into the directory where the script lives and run the script:
cd ensembl/misc-scripts/gifts_import/


eg: 
(with -dryrun flag won't insert)
perl populate_protein_feature.pl -gfhost localhost -gfport 3306 -gfuser prem -gfpass prem -gfdb carlos_ensembl_gifts -pfhost localhost -pfport 3306 -pfuser prem -pfpass prem -pfdb homo_sapiens_core_91_38 -dryrun

(without -dryrun flag creates a full version)
perl populate_protein_feature.pl -gfhost localhost -gfport 3306 -gfuser prem -gfpass prem -gfdb carlos_ensembl_gifts -pfhost localhost -pfport 3306 -pfuser prem -pfpass prem -pfdb homo_sapiens_core_91_38


Usage:

  $0 -gfhost giftsdb_host -gfport giftsdb_port -gfuser giftsdb_user -gfpass giftsdb_pass -gfdb giftsdb_dbname -pfhost protfeatdb_host -pfport protfeatdb_port -pfuser protfeatdb_user -pfpass protfeatdb_pass -pfdb protfeatdb_dbname 
  $indent [-help]  
  
  -gfhost              GIFTS database host

  -gfport              GIFTS database port

  -gfuser              GIFTS database user

  -gfpass              GIFTS database password

  -gfdb                GIFTS database name

  -pfhost              core database host

  -pfport              core database port

  -pfuser              core database user
  
  -pfpass              core database pass

  -pfdb                core database name

  -dryrun              dryrun, no insert

  -help                This message


EOF

}
