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
use Bio::SeqIO;
my $reg = "Bio::EnsEMBL::Registry";

my $species = shift;

if($species eq "mouse"){
  
  $reg->load_registry_from_db(
			      -host => 'ens-staging2',
			      -user => 'ensro');
}
else{
  $reg->load_registry_from_db(
			      -host => 'ens-staging1',
			      -user => "ensro");
}

my $core_sa = $reg->get_adaptor($species,"core","slice");
my $vega_sa = $reg->get_adaptor($species,"vega","slice");

if(!defined($core_sa)){
  die "Could not get core slice adaptor for $species???";
}

if(!defined($vega_sa)){
  die "Could not get vega slice adaptor for $species???";
}


my $sql = 'select t.stable_id, x.display_label from xref x, object_xref ox , transcript t, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and t.transcript_id = ox.ensembl_id and e.db_name like ?';



my $sth = $vega_sa->dbc->prepare($sql) || die "Could not prepare $sql for vega";


my %vega_ott_to_enst;
my %core_ott_to_enst;

$sth->execute("ENST_CDS") or croak( $vega_sa->dbc->errstr() );
while ( my @row = $sth->fetchrow_array() ) {
  $vega_ott_to_enst{$row[0]} = $row[1];
}

$sth->execute("ENST_ident") or croak( $vega_sa->dbc->errstr() );
while ( my @row = $sth->fetchrow_array() ) {
  $vega_ott_to_enst{$row[0]} = $row[1];
}

print "We have from the vega database ".scalar(%vega_ott_to_enst)." ott to enst entries\n ";



$sth = $core_sa->dbc->prepare($sql) || die "Could not prepare $sql for core ";
$sth->execute("OTTT") or croak( $core_sa->dbc->errstr() );
my $core_extra_count = 0;
while ( my @row = $sth->fetchrow_array() ) {
  if(!defined($vega_ott_to_enst{$row[1]})){
    $core_extra_count++;
    if($core_extra_count < 10){
      print "core extra ".$row[1]." ". $row[0]."\n";
    }
  }
  $core_ott_to_enst{$row[1]} = $row[0];
}
print "Core extra tags (OTTT) -> $core_extra_count \n";
$core_extra_count = 0;

$sth = $core_sa->dbc->prepare($sql) || die "Could not prepare $sql for core ";
$sth->execute("shares_CDS_and_UTR_with_OTTT") or croak( $core_sa->dbc->errstr() );
my $core_extra_count = 0;
while ( my @row = $sth->fetchrow_array() ) {
  if(!defined($vega_ott_to_enst{$row[1]})){
    $core_extra_count++;
    if($core_extra_count < 10){
      print "core extra ".$row[1]." ". $row[0]."\n";
    }
  }
  $core_ott_to_enst{$row[1]} = $row[0];
}
print "Core extra tags (shares_CDS_and_UTR_with_OTTT) -> $core_extra_count \n";
$core_extra_count = 0;

#$sth = $core_sa->dbc->prepare($sql) || die "Could not prepare $sql for core ";
#$sth->execute("shares_CDS_with_ENST") or croak( $core_sa->dbc->errstr() );
#my $core_extra_count = 0;
#while ( my @row = $sth->fetchrow_array() ) {
#  if(!defined($vega_ott_to_enst{$row[1]})){
#    $core_extra_count++;
#    if($core_extra_count < 10){
#      print "core extra ".$row[1]." ". $row[0]."\n";
#    }
#  }
#  $core_ott_to_enst{$row[1]} = $row[0];
#}
#print "Core extra tags (shares_CDS_with_ENST) -> $core_extra_count \n";
#$core_extra_count = 0;



#$sth->execute("Vega_transcript") or croak( $core_sa->dbc->errstr() );
#while ( my @row = $sth->fetchrow_array() ) {
#  if(!($row[1] =~ /^OTT/)){
#    next;
#  }
#  if(!defined($vega_ott_to_enst{$row[1]})){
#    $core_extra_count++;
#    if($core_extra_count < 10){
#      print "core extra ".$row[1]." ". $row[0]."\n";
#    }
#  }
#  $core_ott_to_enst{$row[1]} = $row[0];
#}
#print "Core extra tags (Vega_transcript) -> $core_extra_count \n";
$core_extra_count = 0;

$sth->execute("shares_CDS_with_OTTT") or croak( $core_sa->dbc->errstr() );
while ( my @row = $sth->fetchrow_array() ) {
  if(!($row[1] =~ /^OTT/)){
    next;
  }
  if(!defined($vega_ott_to_enst{$row[1]})){
    $core_extra_count++;
    if($core_extra_count < 10){
      print "core extra ".$row[1]." ". $row[0]."\n";
    }
  }
  $core_ott_to_enst{$row[1]} = $row[0];
}

print "Core extra tags (shares_CDS_with_OTTT) -> $core_extra_count \n";
$core_extra_count = 0;
print "We have from the core database ".scalar(%core_ott_to_enst)." ott to enst entries\n ";



my $vega_extra_count = 0;
foreach my $key (keys %vega_ott_to_enst){
  if(!defined($core_ott_to_enst{$key})){
    $vega_extra_count++;
    if($vega_extra_count < 10){
      print "vega extra ".$key." ". $vega_ott_to_enst{$key}."\n";
    }
  }
}

print "Vega extra tags -> $vega_extra_count \n";

