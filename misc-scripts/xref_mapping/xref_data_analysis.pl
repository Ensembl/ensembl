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

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long qw(:config pass_through);

my $species;
my $version;
my $host;
my $user;
my $dbname;
my $display_name;
my $source;
my $show_links;

my $ret = Getopt::Long::GetOptions ('dbname=s'        => \$dbname,
                                    'species=s'       => \$species,
                                    'host=s'          => \$host,
                                    'user=s'          => \$user,
                                    'version=s'       => \$version,
                                    'name=s'          => \$display_name,
                                    'source=s'        => \$source,
                                    'show_links'      => \$show_links,
                                    'help'            => sub {usage(); exit(0); } );


if(!defined $host){
  $host = "ensembldb.ensembl.org";
  $user = "anonymous";
}

if(defined $dbname){
  if(!defined $species){
    $species = "no_eye_deer";
  }
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
			     -user    => $user,
			     -dbname  => $dbname,
			     -host    => $host,
			     -species => $species,
			     -group   => "core");
}
elsif(defined $version){
 Bio::EnsEMBL::Registry->load_registry_from_db( -host => $host,  -user => $user, -db_version => $version);
}
else{
 Bio::EnsEMBL::Registry->load_registry_from_db( -host => $host,  -user => $user);
}

if(defined $display_name){
  specific_example($display_name, $source);
}
#elsif(defined $display_name){
#  die "must specify source if using -name option\n";
#}

my $adap = Bio::EnsEMBL::Registry->get_adaptor($species,"core","gene");

if(!defined $adap){
  die "Could not get adaptor for $species\n";
}

my $dbi = $adap->dbc();

# General
my $gen_sql =(<<'GEN');
SELECT  e.db_name, x.info_type, count(*) 
  FROM xref x, external_db e 
    WHERE x.external_db_id = e.external_db_id 
      GROUP BY  e.db_name, x.info_type
GEN

my %big_hash;

my $gen_sth = $dbi->prepare($gen_sql);
my ($name, $type, $count);
$gen_sth->execute();
$gen_sth->bind_columns(\$name, \$type, \$count);
while($gen_sth->fetch){
  $type ||= "GeneBuilder";
  $big_hash{$name}{$type} = $count;
}
$gen_sth->finish;


# dependent xrefs
my $dependent_sql =(<<'DEP');
SELECT master.db_name, dependent.db_name, count(1)
  FROM xref as xd, xref as xm,
       external_db as master, external_db as dependent,
       dependent_xref dx
    WHERE dx.dependent_xref_id = xd.xref_id AND
          dx.master_xref_id = xm.xref_id AND
          xm.external_db_id = master.external_db_id AND
          xd.external_db_id = dependent.external_db_id
     GROUP BY  master.db_name, dependent.db_name
DEP

my $dep_sth = $dbi->prepare($dependent_sql);
$dep_sth->execute();
my ($master, $dependent);
$dep_sth->bind_columns(\$master, \$dependent, \$count);

my %dependents;
while($dep_sth->fetch()){
  $dependents{$dependent}{$master} = $count;
}
$dep_sth->finish;

foreach my $dbname (sort keys %big_hash){
  foreach my $info_type (keys %{$big_hash{$dbname}}){
    print "$dbname\t$info_type\t".$big_hash{$dbname}{$info_type}."\n";
    if($info_type eq "DEPENDENT"){
      foreach my $master_db (keys %{$dependents{$dbname}} ){
	print "\tvia $master_db\t".$dependents{$dbname}{$master_db}."\n";
      }
    }
  }
}

sub specific_example {

  my ($search_name, $source_name) = @_;

  my $db_adap = Bio::EnsEMBL::Registry->get_adaptor($species,"core","dbentry");
  my $gene_adap = Bio::EnsEMBL::Registry->get_adaptor($species,"core","gene");

  my $dbentrys = $db_adap->fetch_all_by_name($search_name, $source_name);

  my $found = 0;
  foreach my $dbentry (@$dbentrys){
    $found = 1;
    print "##############################\ndbname :".$dbentry->db_display_name."\n";
    print "accession is: ".$dbentry->primary_id."\n";
    print "display id: ".$dbentry->display_id."\n";
    print "info type: ".$dbentry->info_type."\n";
    if($dbentry->info_type =~ /DEPENDENT/ms){
      my @masters = @{$dbentry->get_all_masters()};
      foreach my $entry (@masters){
	print "\tMaster: ".$entry->db_display_name." : ".$entry->primary_id."\n";
      }
    }
    if(defined $dbentry->info_text){
      print "info text: ".$dbentry->info_text."\n";
    }
    if($show_links && !($dbentry->info_type =~ /UNMAPPED/ms)){
      print "Linked to the following genes: ";
      foreach my $gene (@{$gene_adap->fetch_all_by_external_name($dbentry->primary_id, $dbentry->dbname)}){
	print $gene->stable_id." ";
      }	
      print "\n";
    }
    print "##############################\n";
  }
  if(!$found){
    if(defined $source_name){
      print "$search_name NOT FOUND for $source_name\n";
    }
    else{
      print "$search_name NOT FOUND\n";
    }
  }
  exit(0);
}

sub usage {

  print << 'EOF';
  This script will report where external database references come from and give numbers
  for each of these. If a name is given, information about only those ones are given.

  Options

  -species    The species you want to look at. This MUST be set unless you use dbname

  -dbname     The database you want data on.

  -user       The user to use to read the core database in ro mode (default anonymous)

  -host       The mysql server (default ensembldb.ensembl.org)

  -version    Can be used to change the version of the database to be examined.

  -name       Get data for this display_label/name.

  -source     the source name for the accession (i.e. HGNC, MGI, RefSeq_mRNA)

  -show_links If name used then show which ensembl objects is is linked to.

   A typical run would be to see what xrefs are there for human :-

      perl xref_data_analysis.pl -species human


  this would generate a list of the number of xrefs for the human database on 
  ensembldb.ensembl.org that matches the API version you are running.


  If you want to look at the 64 version and you are not using the 64 API :-

      perl xref_data_analysis.pl -species human -version 64


  If you want to look at a test database

     perl xref_data_analysis.pl -dbname ianl_human_core_65 -host ens-research -user ro


  To find how a partcular xref was mapped use -name to specify this and -source if more than
  one xref may have the same accession

  e.g. how was accession BRCA2 in HGNC  mapped to ensembl and to which genes

       perl xref_data_analysis.pl -species human -name BRCA2 -source HGNC -show_links


EOF
  return;
}
