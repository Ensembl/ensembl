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
use Getopt::Long;
my $reg = "Bio::EnsEMBL::Registry";



my ($user, $pass, $host, $port, $dbname, $type, $stable_id);
my ($species, $db_version, $ignore);

GetOptions('user=s'        => \$user,
           'pass=s'        => \$pass,
           'host=s'        => \$host,
           'port=i'        => \$port,
           'dbname=s'      => \$dbname,
           'type=s'        => \$type,
           'species=s'     => \$species,
           'version=s'  => \$db_version,
           'id=s'          => \$stable_id,
           'ignore=s'      => \$ignore,
           'help'          => sub {usage(); exit(0); } );



if(defined $species and defined $dbname){
  die "You can only specify a species OR a dbname\n";
}

my %ignore_sources;
if(defined $ignore){
  foreach my $item (split(',',$ignore)){
    $ignore_sources{$item} = 1;
    print "ignoring $item sources\n";
  }
}
if(defined $species){
   $reg->load_registry_from_db( -host => $host || 'ensembldb.ensembl.org',
                                -user => $user || 'anonymous',
                                -db_version => $db_version || undef)
}
else{
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
             -species => "test",
	     -group   => 'core',
             -host    => $host||'ensembldb.ensembl.org',
	     -port    => $port||3306,
             -user    => $user||'anonymous',
	     -pass    => $pass||undef,
             -dbname  => $dbname);
  my $species = "test";
}

my $gene_adap =  $reg->get_adaptor($species,"core","gene");
if(!defined($gene_adap)){
  die "unable to connect to database and get adaptor\n";
}
my @transcripts ;
if(lc($type) eq "gene"){
  my $gene = $gene_adap->fetch_by_stable_id($stable_id);
  if(!defined($gene)){
    die "unable to find stable_id $stable_id\n";
  }	
  my $status = $gene->status;
  my $display = $gene->display_xref;
  if(defined($display)){
    $display = $display->display_id;
  }
  else{
    $display = "NO_DISPLAY_XREF";
  }	
  print "$status Gene $stable_id display = $display status = $status\n";
  print "position\t".$gene->slice->seq_region_name."\n";
  if(defined($gene->description)){
    print $gene->description."\n";
  }
  my @db_entries = @{$gene->get_all_DBEntries()}; 
  foreach my $entry (@db_entries){
    if(!defined $ignore_sources{$entry->dbname}){
      print "\t".$entry->dbname."\t".$entry->primary_id."\t".$entry->display_id."\t".$entry->info_text."\n";
    }
  }
  @transcripts = @{$gene->get_all_Transcripts()};
}	
elsif(lc($type) eq "name"){
  my @genes = @{$gene_adap->fetch_all_by_external_name($stable_id)};

  foreach my $gene (@genes){
    my $status = $gene->status;
    my $display = $gene->display_xref;
    if(defined($display)){
      $display = $display->display_id;
    }
    else{
      $display = "NO_DISPLAY_XREF";
    }	
    print "$status Gene $stable_id display = $display status = $status\n";
    print "position\t".$gene->slice->seq_region_name."\n";
    my @db_entries = @{$gene->get_all_DBEntries()}; 
    foreach my $entry (@db_entries){
      if(!defined $ignore_sources{$entry->dbname}){
	print "\t".$entry->dbname."\t".$entry->primary_id."\t".$entry->display_id."\t".$entry->info_text."\n";
      }
    }
    if(defined($gene->description)){
      print $gene->description."\n";
    }
    push  @transcripts, @{$gene->get_all_Transcripts()};  
  }
}
else{  
  my $tran_adap =  $reg->get_adaptor($species,"core","transcript");
  my $tran = $tran_adap->fetch_by_stable_id($stable_id);
  push @transcripts, $tran; 
}

foreach my $tran (@transcripts){
  my $display = $tran->display_xref;
  if(defined($display)){
    $display = $display->display_id;
  }
  else{
    $display = "NO_DISPLAY_XREF";
  }
  print $tran->stable_id."\t".$display."\tstatus=".($tran->status||"NOT SET")."\n";
  my @db_entries = @{$tran->get_all_DBLinks()};
  foreach my $dbe (@db_entries){
    if(!defined $ignore_sources{$dbe->dbname}){
      print "\t*".$dbe->dbname."*\t".$dbe->primary_id."\t".$dbe->display_id."\t".$dbe->info_text."\n";
    }
  }
}


sub usage {

print << 'EOF';

  This script will report the external database references for a particular
  gene or transcript. The gene or transcript is specified by setting the type
  (gene, transcript or name) and also the identifier being looked up. i.e.
  -type gene id ENSG00000053918
  -type transcript ENST00000155840
  -type name BRCA2


  -species    The species you want to look at. This MUST be set unless you
               use dbname

  -type       Options are gene, transcript or name.(must be set)

  -id         Identifier to look up. must be of the type specified above.

  -ignore     external database source to ignore (seperated by comas)

  -dbname     The database to look for the xrefs. (if species not set)

  -user       The user to use to read the core database in ro mode 
              (default anonymous)

  -pass       The password for the user.

  -host       The mysql server (default ensembldb.ensembl.org)

  -version    Can be used to change the version of the database to be examined.
              (only used when species is set)


  Examples.

  1) Find the external databases references for the human gene ENSG00000053918.

     perl list_xrefs.pl -species human -type gene -id ENSG00000053918

  2) Find the xrefs for the mouse gene that is named Cntnap1

     perl list_xrefs.pl -species mouse -type name -id cntnap1

  3) as above but we are not interested in GO or goslim_goa references

     perl list_xrefs.pl -species mouse -type name -id cntnap1
                        -ignore GO,goslim_goa

EOF
}


