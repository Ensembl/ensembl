=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::HGNC_CCDSParser;

use strict;
use warnings;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port;
  my $dbname;
  my $pass;
  my $wget = "";

  if($file =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($file =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($file =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($file =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }
  if($file =~ /wget[=][>](\S+?)[,]/){
    $wget = $1;
  }


  my $ua = LWP::UserAgent->new();
  $ua->timeout(10);
  $ua->env_proxy();

  my %ccds_to_hgnc;

  my $response = $ua->get($wget);

  if ( !$response->is_success() ) {
    die $response->status_line;
  }
  else{
    my @lines = split(/\n/,$response->content);
    foreach my $line (@lines){
      my($hgnc, $junk, $ccds) = split(/\t/,$line);
      my @ccds_list = split(/, /,$ccds);
      foreach my $c (@ccds_list){
	$ccds_to_hgnc{$c} = $hgnc;
      }
    }

  }

  my $ccds_db =  XrefParser::Database->new({ host   => $host,
					     port   => $port,
					     user   => $user,
					     dbname => $dbname,
					     pass   => $pass});
  my $dbi2 = $ccds_db->dbi();

  if(!defined($dbi2)){
    return 1;
  }



  my $sql =(<<'SQL');
SELECT ox.ensembl_id, x.dbprimary_acc
  FROM object_xref ox, xref x, external_db e
    WHERE x.xref_id = ox.xref_id and
          x.external_db_id = e.external_db_id AND
          e.db_name like "Ens_%_transcript" AND
          x.dbprimary_acc like "ENST%"
SQL

  my %trans_id_to_stable_id;
  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $trans_id_to_stable_id{$row[0]} = $row[1];
  }
  $sth->finish;  


  
  $sql = 'select ox.ensembl_id, x.dbprimary_acc from object_xref ox, xref x, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and ox.ensembl_object_type = "Transcript" and e.db_name like "CCDS"'; 

  my %ccds_to_stable_id;
  $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    if(defined($trans_id_to_stable_id{$row[0]})){
      push @{$ccds_to_stable_id{$row[1]}}, $trans_id_to_stable_id{$row[0]};
    }
    else{
      print "NO transcript_stable_id for  for ".$row[0]."\n";
    }
  }
  $sth->finish;
  


  # becouse the direct mapping have no descriptions etc
  # we have to steal these fromt he previous HGNC parser.
  # This is why the order states this is after the other one.
  # maybe 1091,1092 is not right maybe should use name = HGNC and priority = 30r4 ??

  my %label;
  my %version;
  my %description;

  my $sql_syn = "insert ignore into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql_syn);
  
  my $syn_hash = $self->get_ext_synonyms("HGNC", $dbi);

  $sql = 'select source_id, priority_description from source where name like "HGNC"';
  $sth = $dbi->prepare($sql);
  
  $sth->execute();
  my ($hgnc_source_id, $desc);
  $sth->bind_columns(\$hgnc_source_id, \$desc);
  my @arr;
  while($sth->fetch()){
    push @arr, $hgnc_source_id;
  }
  $sth->finish;
  
  $sql = "select accession, label, version,  description from xref where source_id in (".join(", ",@arr).")";

  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    if(defined($desc)){
      $label{$acc} = $lab;
      $version{$acc} = $ver;
      $description{$acc} = $desc;
    }
  }
  $sth->finish;

  my $xref_count = 0;
  my $no_ccds_to_hgnc = 0;
  my $direct_count = 0;
  foreach my $ccds (keys %ccds_to_stable_id){
    if(defined($ccds_to_hgnc{$ccds})){
      my $hgnc = $ccds_to_hgnc{$ccds};
      my $xref_id = $self->add_xref({ acc        => $hgnc,
				      version    => $version{$hgnc} ,
				      label      => $label{$hgnc}||$hgnc ,
				      desc       => $description{$hgnc},
				      source_id  => $source_id,
				      species_id => $species_id,
                                      dbi        => $dbi,
				      info_type  => "DIRECT"} );

      foreach my $stable_id (@{$ccds_to_stable_id{$ccds}}){
	$self->add_direct_xref($xref_id, $stable_id, "Transcript", "", $dbi);
	$direct_count++;
      }	

      $xref_count++;

      if(defined($syn_hash->{$hgnc})){
	foreach my $syn (@{$syn_hash->{$hgnc}}){
	  $add_syn_sth->execute($xref_id, $syn);
	}
      }
      
    }
    else{
         $no_ccds_to_hgnc++;
#      print "no ccds to hgnc for $ccds\n";

    }
  }
  $add_syn_sth->finish;
  print "$no_ccds_to_hgnc missed as no hgnc for the ccds. Added $xref_count HGNC xrefs via CCDS and $direct_count direct xrefs\n" if($verbose);
  return 0;
}

1;
