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

package XrefParser::curated_transcriptParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

use XrefParser::Database;
use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";


sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my ($type, $my_args) = split(/:/,$file);
  
  my ($host, $source_prefix);
  my $user = "ensro";

  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }
  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if ($my_args =~ /source[=][>](\S+?)[,]/){
    $source_prefix = $1;
  }
  my %id2name = $self->species_id2name;
  my $species_name = $id2name{$species_id}[0];
  if (!$source_prefix){
    die "Species is $species_name and is not homo_sapiens, mus_musculus, danio_rerio or sus_scrofa, the only four valid species\n";
  }

  my $vuser = "ensro";
  my $vhost;
  my $vport = 3306;
  my $vdbname;
  my $vpass;
 
  my $cuser  ="ensro";
  my $chost;
  my $cport = 3306;
  my $cdbname;
  my $cpass;

  if($my_args =~ /chost[=][>](\S+?)[,]/){
    $chost = $1;
  }
  if($my_args =~ /cport[=][>](\S+?)[,]/){
    $cport =  $1;
  }
  if($my_args =~ /cdbname[=][>](\S+?)[,]/){
    $cdbname = $1;
  }
  if($my_args =~ /cpass[=][>](\S+?)[,]/){
    $cpass = $1;
  }
  if($my_args =~ /cuser[=][>](\S+?)[,]/){
    $cuser = $1;
  }
  if($my_args =~ /vhost[=][>](\S+?)[,]/){
    $vhost = $1;
  }
  if($my_args =~ /vport[=][>](\S+?)[,]/){
    $vport =  $1;
  }
  if($my_args =~ /vdbname[=][>](\S+?)[,]/){
    $vdbname = $1;
  }
  if($my_args =~ /vpass[=][>](\S+?)[,]/){
    $vpass = $1;
  }
  if($my_args =~ /vuser[=][>](\S+?)[,]/){
    $vuser = $1;
  }

  my $vega_dbc;
  my $core_dbc;
  if(defined($vdbname)){
    print "Using $vhost $vdbname for Vega\n";

    my $vega_db =  XrefParser::Database->new({ host   => $vhost,
					       port   => $vport,
					       user   => $vuser,
					       dbname => $vdbname,
					       pass   => $vpass});

    $vega_dbc = $vega_db->dbi;
    if(!defined($vega_dbc)){
      print "Problem could not open connection to $vhost, $vport, $vuser, $vdbname, $vpass\n";
      return 1;
    }
  } else {
    $reg->load_registry_from_db(
                                -host => $host,
                                -user => $user,
                                -port => 4519,
                                -species => $species_name);

    $vega_dbc = $reg->get_adaptor($species_name,"vega","slice");
    if(!defined($vega_dbc)){
      print "Could not connect to $species_name vega database using load_registry_from_db $host $user\n";
      return 1;
    }
    $vega_dbc = $vega_dbc->dbc;
  }

   if (defined($cdbname)){
    my $core_db =  XrefParser::Database->new({ host   => $chost,
					       port   => $cport,
					       user   => $cuser,
					       dbname => $cdbname,
					       pass   => $cpass});
    $core_dbc = $core_db->dbi;
    if(!defined($core_dbc)){
      print "Problem could not open connectipn to $chost, $cport, $cuser, $cdbname, $cpass\n";
      return 1;
    }
  } else{
    $reg->load_registry_from_db(
                                -host => $host,
                                -user => $user,
                                -port => 4519,
			        -species => $species_name);

    $core_dbc = $reg->get_adaptor($species_name,"core","slice");
    if(!defined($core_dbc)){
      print "Could not connect to $species_name core database using load_registry_from_db $host $user\n";
      return 1;
    }
    $core_dbc= $core_dbc->dbc;
  }


  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');
  my $curated_source_id =
    $self->get_source_id_for_source_name($source_prefix."_curated_transcript_notransfer");
 
 print "source id is $source_id, curated_source_id is $curated_source_id\n";

  my $sql = 'select distinct t.stable_id, x.display_label, t.status from analysis a, xref x, object_xref ox , external_db e, transcript t where t.analysis_id = a.analysis_id and a.logic_name like "%havana%" and e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and t.transcript_id = ox.ensembl_id and ox.ensembl_object_type = "Transcript" and e.db_name like ?';

  my $sql_vega = 'select t.stable_id, x.display_label, t.status from xref x, object_xref ox , external_db e, transcript t where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and t.transcript_id = ox.ensembl_id and t.stable_id <> x.display_label and ox.ensembl_object_type = "Transcript" and e.db_name like ?';


  my %ott_to_vega_name;
  my %ott_to_enst;


  my $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";

  my $count_ott_to_enst;
  foreach my $external_db (qw(Vega_transcript shares_CDS_with_OTTT shares_CDS_and_UTR_with_OTTT OTTT)){
    $sth->execute($external_db) or croak( $core_dbc->errstr());
    while ( my @row = $sth->fetchrow_array() ) {
      push (@{$ott_to_enst{$row[1]}},$row[0]);
      $count_ott_to_enst++;
    }
  }

  print "We have $count_ott_to_enst ott to enst entries\n " if($verbose);


  my $dbi = $self->dbi();

  my $status_insert_sth = $dbi->prepare("INSERT IGNORE INTO havana_status (stable_id, status) values(?, ?)")
    || die "Could not prepare status_insert_sth";

  my %ott_to_status;
  $sth = $vega_dbc->prepare($sql_vega);   # funny number instead of stable id ?????
  $sth->execute("Vega_transcript") or croak( $vega_dbc->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_vega_name{$row[0]} = $row[1];
    $ott_to_status{$row[0]} = $row[2];
  }
  $sth->finish;

  my $xref_count = 0;

  foreach my $ott (keys %ott_to_enst){
    if(defined($ott_to_vega_name{$ott})){
      my $id = $curated_source_id;
      my $name  = $ott_to_vega_name{$ott};
      $name =~ s/WU://;
      if($name =~ /[.]/){
	$id = $clone_source_id;
# number is no longer the clone version but the gene number so we need to keep it now.
#        $name =~ s/[.]\d+//;    #remove .number  #
      }
      my $xref_id = $self->add_xref({ acc        => $name,
				      label      => $name,
				      source_id  => $id,
				      species_id => $species_id,
				      info_type  => "DIRECT"} );
      $xref_count++;
      
      foreach my $ensembl_stable_id (@{$ott_to_enst{$ott}}) {
	  $self->add_direct_xref($xref_id, $ensembl_stable_id, "transcript", "");
      }
    }
    if(defined($ott_to_status{$ott})){
      foreach my $ensembl_stable_id (@{$ott_to_enst{$ott}}) {
        $status_insert_sth->execute($ensembl_stable_id, $ott_to_status{$ott});
      }
    }
    
  }
 

  # need to add gene info to havana_status table
  $sql = 'select g.stable_id, x.display_label from xref x, object_xref ox , external_db e, gene g where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and g.gene_id = ox.ensembl_id and e.db_name like "OTTG" and ox.ensembl_object_type = "Gene" ';

  $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";
  $sth->execute() or croak( $core_dbc->errstr());
  my %ottg_to_ensg;
  while ( my @row = $sth->fetchrow_array() ) {
    $ottg_to_ensg{$row[1]} = $row[0];
  }

  $sth = $vega_dbc->prepare("select stable_id, status from gene");
  $sth->execute() or croak( $core_dbc->errstr());
  while ( my @row = $sth->fetchrow_array() ) {
    if(defined($ottg_to_ensg{$row[0]}) and defined($row[1])){
      $status_insert_sth->execute($ottg_to_ensg{$row[0]}, $row[1]);
    }
  }

  print "$xref_count direct xrefs succesfully parsed\n" if($verbose);
  return 0;
}





1;
