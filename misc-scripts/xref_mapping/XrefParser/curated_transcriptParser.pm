package XrefParser::curated_transcriptParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use Bio::EnsEMBL::Registry;
my $reg = "Bio::EnsEMBL::Registry";


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: curated_transcriptParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run_script {
  my $self = shift if (defined(caller(1)));

  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;
  my $verbose = shift;

  my ($type, $my_args) = split(/:/,$file);
  
  my $host = "ens-staging1";
  my $user = "ensro";

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }
  my %id2name = $self->species_id2name;
  my $species_name = $id2name{$species_id}[0];
  my $source_prefix;
  if($species_name eq "homo_sapiens" ){
    $source_prefix = "HGNC";
  }
  elsif($species_name eq "mus_musculus" ){
    $source_prefix = "MGI";
  }
  else{
    die "Species is $species_name and is not homo_sapines or mus_musculus the only two valid species\n";
  }

  my $vuser  ="ensro";
  my $vhost;
  my $vport;
  my $vdbname;
  my $vpass;
 
  my $cuser  ="ensro";
  my $chost;
  my $cport;
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
    print "Using $host $vdbname for Vega and cdbname for Core\n";
    $vega_dbc = $self->dbi2($vhost, $vport, $vuser, $vdbname, $vpass);
    if(!defined($vega_dbc)){
      print "Problem could not open connectipn to $vhost, $vport, $vuser, $vdbname, $vpass\n";
      return 1;
    }
    $core_dbc = $self->dbi2($chost, $cport, $cuser, $cdbname, $cpass);
    if(!defined($core_dbc)){
      print "Problem could not open connectipn to $chost, $cport, $cuser, $cdbname, $cpass\n";
      return 1;
    }

  }
  else{
    my @dbas = @{$reg->get_all_adaptors()};
    foreach my $db (@dbas){
      print "DB: $db\n";
    }
    $reg->load_registry_from_db(
                                -host => $host,
                                -user => $user);

    @dbas = @{$reg->get_all_adaptors(-species => $species_name)};
#    foreach my $db (@dbas){
#      print "DB2: ".$db->db->group."\t$db\n";
#    }
    $vega_dbc = $reg->get_adaptor($species_name,"vega","slice");
    if(!defined($vega_dbc)){
      print "Could not connect to $species_name vega database using load_registry_from_db $host $user\n";
      return 1;
    }
    $vega_dbc = $vega_dbc->dbc;
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
    $self->get_source_id_for_source_name($source_prefix."_curated_transcript");
 
 print "source id is $source_id, curated_source_id is $curated_source_id\n";

  my $sql = 'select tsi.stable_id, x.display_label, t.status from xref x, object_xref ox , transcript_stable_id tsi, external_db e, transcript t where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and t.transcript_id = tsi.transcript_id and e.db_name like ?';


  my %ott_to_vega_name;
  my %ott_to_enst;


  my $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";

  foreach my $external_db (qw(Vega_transcript OTTT shares_CDS_and_UTR_with_OTTT shares_CDS_with_OTTT)){
    $sth->execute($external_db) or croak( $core_dbc->errstr());
    while ( my @row = $sth->fetchrow_array() ) {
      $ott_to_enst{$row[1]} = $row[0];
    }
  }

  print "We have ".scalar(%ott_to_enst)." ott to enst entries\n " if($verbose);

  my $sth = $vega_dbc->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("Vega_transcript") or croak( $vega_dbc->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_vega_name{$row[0]} = $row[1];
  }
  $sth->finish;

  my $xref_count = 0;

  foreach my $ott (keys %ott_to_enst){
    if(defined($ott_to_vega_name{$ott})){
      my $id = $curated_source_id;
      my $name  = $ott_to_vega_name{$ott};
      if($name =~ /[.]/){
	$id = $clone_source_id;
# number is no longer the clone version but the gene number so we need to keep it now.
#        $name =~ s/[.]\d+//;    #remove .number  #
      }
      my $xref_id = $self->add_xref($name, "" , $name , "", $id, $species_id, "DIRECT");
      $xref_count++;
      
      $self->add_direct_xref($xref_id, $ott_to_enst{$ott}, "transcript", "");
    }
  }

  print "$xref_count direct xrefs succesfully parsed\n" if($verbose);
  return 0;
}





1;

