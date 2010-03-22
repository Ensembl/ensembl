package XrefParser::MGI_curated_transcriptParser;

use strict;
use File::Basename;
use Bio::EnsEMBL::Registry;
use base qw( XrefParser::BaseParser );
my $reg = "Bio::EnsEMBL::Registry";


if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: MGI_curated_transcriptParser.pm file <source_id> <species_id>\n\n";
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

  my $host = "ens-staging2";
  my $user = "ensro";

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
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
    $reg->load_registry_from_db(
				-host => $host,
				-user => $user);

    $vega_dbc = $reg->get_adaptor("mouse","vega","slice");
    if(!defined($vega_dbc)){
      print "Could not connect to mouse vega database using load_registry_from_db $host $user\n";
      return 1;
    } 
    $vega_dbc = $vega_dbc->dbc;
    $core_dbc = $reg->get_adaptor("mouse","core","slice");
    if(!defined($core_dbc)){
      print "Could not connect to mouse core database using load_registry_from_db $host $user\n";
      return 1;
    }
    $core_dbc= $core_dbc->dbc;
  }

  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');
  my $curated_source_id =
    $self->get_source_id_for_source_name('MGI_curated_transcript');
 
  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my %ott_to_vega_name;
  my %ott_to_enst;


  my $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";

  foreach my $external_db (qw(OTTT shares_CDS_and_UTR_with_OTTT shares_CDS_with_ENST)){
    $sth->execute($external_db) or croak( $core_dbc->errstr());
    while ( my @row = $sth->fetchrow_array() ) {
      $ott_to_enst{$row[1]} = $row[0];
    }
  }

  
  $sth = $vega_dbc->prepare($sql) || die "Could not prepare for vega $sql\n";
  $sth->execute("Vega_transcript") or croak( $vega_dbc->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
     next if ($row[1] =~ /^OTT/); 
    $ott_to_vega_name{$row[0]} = $row[1];
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
# keep it apparently       $name =~ s/[.]\d+//;    #remove .number
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

