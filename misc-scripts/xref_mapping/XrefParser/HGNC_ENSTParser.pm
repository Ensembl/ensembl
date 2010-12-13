package XrefParser::HGNC_ENSTParser;

use strict;

use DBI;

use Bio::EnsEMBL::Registry;
use base qw( XrefParser::BaseParser );
my $reg = "Bio::EnsEMBL::Registry";

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes



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

  my $vuser;
  my $vhost;
  my $vport;
  my $vdbname;
  my $vpass;

  if($my_args =~ /vhost[=][>](\S+?)[,]/){
    $vhost = $1;
  }
  if($my_args =~ /vuser[=][>](\S+?)[,]/){
    $vuser = $1;
  }
  if($my_args =~ /vport[=][>](\S+?)[,]/){
    $vport =  $1;
  }
  if($my_args =~ /vdbname[=][>](\S+?)[,]/){
    $vdbname = $1;
  }
  else{
    print "No vdbname??? $my_args??\n";
  }	
  if($my_args =~ /vpass[=][>](\S+?)[,]/){
    $vpass = $1;
  }

  my $cuser;
  my $chost;
  my $cport;
  my $cdbname;
  my $cpass;

  if($my_args =~ /chost[=][>](\S+?)[,]/){
    $chost = $1;
  }
  if($my_args =~ /cuser[=][>](\S+?)[,]/){
    $cuser = $1;
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
  my $vega_dbc;
  my $core_dbc;
  if(defined($vdbname)){
    print "Using $host $vdbname for Vega and $cdbname for Core\n";
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

    $vega_dbc = $reg->get_adaptor("human","vega","slice");
    if(!defined($vega_dbc)){
      print "Could not connect to human vega database using load_registry_from_db $host $user\n";
      return 1;
    }
    $vega_dbc = $vega_dbc->dbc;
    $core_dbc = $reg->get_adaptor("human","core","slice");
    if(!defined($core_dbc)){
      print "Could not connect to human core database using load_registry_from_db $host $user\n";
      return 1;
    }
    $core_dbc= $core_dbc->dbc;
  }



  $source_id = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","vega");


  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my $hgnc_sql = 'select tsi.stable_id, x.dbprimary_acc from xref x, object_xref ox , transcript_stable_id tsi, gene g, external_db e, transcript t  where t.gene_id = g.gene_id and g.gene_id = ox.ensembl_id and tsi.transcript_id = t.transcript_id and e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and ox.ensembl_object_type = "Gene" and e.db_name like "HGNC"';


  my %ott_to_hgnc;
  my %ott_to_enst;

  my $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";

  foreach my $external_db (qw(OTTT shares_CDS_and_UTR_with_OTTT shares_CDS_with_OTTT)){
    $sth->execute($external_db) or croak( $core_dbc->errstr());
    while ( my @row = $sth->fetchrow_array() ) {
      $ott_to_enst{$row[1]} = $row[0];
    }
  }

  print "We have ".scalar(%ott_to_enst)." ott to enst entries\n " if($verbose);

  $sth = $vega_dbc->prepare($hgnc_sql);
  $sth->execute() or croak( $vega_dbc->errstr() );

  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_hgnc{$row[0]} = $row[1];
  }
  $sth->finish;
  print "We have ".scalar(%ott_to_hgnc)." ott to hgnc entries\n" if($verbose);


  my $line_count = 0;
  my $xref_count = 0;

  # becouse the direct mapping have no descriptions etc
  # we have to steal these fromt he previous HGNC parser.
  # This is why the order states this is after the other one.
  # maybe 1091,1092 is not right maybe should use name = HGNC and priority = 30r4 ??

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  


  my $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    
  
  my $syn_hash = $self->get_hgnc_synonyms();
  


  #get the source ids for HGNC refseq, entrezgene and unitprot
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
#  print "$sql\n";;
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver);
  my $hgnc_loaded_count = 0;
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
    $hgnc_loaded_count++;
  }
  $sth->finish;
  if($hgnc_loaded_count == 0){
    die "No point continuing no hgncs there\n";
  }


  my $ignore_count = 0;
  my $ignore_examples ="";
  my %acc;

  foreach my $key ( keys %ott_to_hgnc){
    if(defined($ott_to_enst{$key} )){
      
      my $hgnc = $ott_to_hgnc{$key};
      $hgnc =~ s/HGNC://;
      my $stable_id = $ott_to_enst{$key};

      if(!defined($label{$hgnc})){
	$ignore_count++;
	if($ignore_count < 10){
	  $ignore_examples .= " ".$hgnc;
	}
	next;
      }
      if(!defined($acc{$hgnc})){
	$acc{$hgnc} = 1;
	my $version ="";
	$line_count++;
	
	my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , $description{$hgnc}, $source_id, $species_id, "DIRECT");
	$xref_count++;
	
	
	$self->add_direct_xref($xref_id, $stable_id, "transcript", "");

	if(defined($syn_hash->{$hgnc})){
	  foreach my $syn (@{$syn_hash->{$hgnc}}){
	    $add_syn_sth->execute($xref_id, $syn);
	  }
	}

      }
    }
  }
  $add_syn_sth->finish;
  print "Parsed $line_count HGNC identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n" if($verbose);
  if($ignore_count){
    print $ignore_count." ignoreed due to numbers no identifiers being no longer valid :- $ignore_examples\n" if($verbose);
  }
  
  
  return 0;
}


1;
