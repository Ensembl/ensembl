package XrefParser::VegaOfficialNameParser;

use strict;
use warnings;
use Carp;
use DBI;
use Bio::EnsEMBL::Registry;
use base qw( XrefParser::BaseParser );
my $reg = "Bio::EnsEMBL::Registry";

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  my ($type, $my_args) = split(/:/,$file);

  my $host;
  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }

  my %id2name = $self->species_id2name;
  my $species_name = $id2name{$species_id}[0];

  my $source_name;
  my $prepend = 1;


  if($species_name eq "homo_sapiens" ){
    $source_name = "HGNC";
    $host = "ens-staging1";
  }
  elsif($species_name eq "mus_musculus" ){
    $source_name = "MGI";
    $prepend = 0;
    $host = "ens-staging2";
  }
  elsif($species_name eq "danio_rerio" ){
    $source_name = "ZFIN_ID";
    $host = "ens-staging1";
  }
  else{
    die "Species is $species_name and is not homo_sapines, mus_musculus or danio_rerio the only three valid species\n";
  }

  my $user = 'ensro';
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

#populate species hash - source name is the hash key
  my %species_names = (
      HGNC => "human",
      ZFIN_ID => "zebrafish",
      MGI => "mouse"
  );

  my $vega_dbc;
  my $core_dbc;

  if(defined($vdbname)){
    print "Using $host $vdbname for Vega and $cdbname for Core\n";
    my $vega_db =  XrefParser::Database->new({ host   => $vhost,
					       port   => $vport,
					       user   => $vuser,
					       dbname => $vdbname,
					       pass   => $vpass});

    $vega_dbc = $vega_db->dbi();
    if(!defined($vega_dbc)){
      print "Problem could not open connectipn to $vhost, $vport, $vuser, $vdbname, $vpass\n";
      return 1;
    }
    my $core_db =  XrefParser::Database->new({ host   => $chost,
					       port   => $cport,
					       user   => $cuser,
					       dbname => $cdbname,
					       pass   => $cpass});

    my $core_dbc = $core_db->dbi();
    if(!defined($core_dbc)){
      print "Problem could not open connectipn to $chost, $cport, $cuser, $cdbname, $cpass\n";
      return 1;
    }

  }
  else{
    $reg->load_registry_from_db(
                                -host => $host,
                                -user => $user);

    $vega_dbc = $reg->get_adaptor($species_names{$source_name},"vega","slice");
    if(!defined($vega_dbc)){
      print "Could not connect to $species_names{$source_name} vega database using load_registry_from_db $host $user\n";
      return 1;
    }
    $vega_dbc = $vega_dbc->dbc;
    $core_dbc = $reg->get_adaptor($species_names{$source_name},"core","slice");
    if(!defined($core_dbc)){
      print "Could not connect to $species_names{$source_name} core database using load_registry_from_db $host $user\n";
      return 1;
    }
    $core_dbc= $core_dbc->dbc;
  }



  $source_id = $self->get_source_id_for_source_name($source_name,"vega");


  my $sql =(<<SQL);
SELECT tsi.stable_id, x.display_label 
  FROM analysis a, xref x, object_xref ox , transcript_stable_id tsi, external_db e , transcript t
   WHERE a.analysis_id = t.analysis_id and 
         t.transcript_id = tsi.transcript_id and
         e.external_db_id = x.external_db_id and 
         x.xref_id = ox.xref_id and 
         tsi.transcript_id = ox.ensembl_id and
         a.logic_name like "%havana%" and 
         e.db_name like ?
SQL

  my $ext_sql =(<<EXT);
SELECT tsi.stable_id, x.dbprimary_acc 
  FROM xref x, object_xref ox , transcript_stable_id tsi, gene g, external_db e, transcript t  
    WHERE t.gene_id = g.gene_id and 
          g.gene_id = ox.ensembl_id and 
          tsi.transcript_id = t.transcript_id and 
          e.external_db_id = x.external_db_id and 
          x.xref_id = ox.xref_id and 
          ox.ensembl_object_type = "Gene" and 
          e.db_name like '$source_name'
EXT

  my %vega_to_ext;
  my %ext_to_core;

  my $sth = $core_dbc->prepare($sql) || die "Could not prepare for core $sql\n";


  foreach my $external_db (qw(Vega_transcript shares_CDS_with_OTTT shares_CDS_and_UTR_with_OTTT OTTT)){
    $sth->execute($external_db) or croak( $core_dbc->errstr());
    while ( my @row = $sth->fetchrow_array() ) {
      $ext_to_core{$row[1]} = $row[0];
    }
  }

  print "We have ".scalar(%ext_to_core)." vega to external source entries\n " if($verbose);

  $sth = $vega_dbc->prepare($ext_sql) || die "could not prepare vega sql: $ext_sql";
  $sth->execute() or croak( $vega_dbc->errstr() );

  while ( my @row = $sth->fetchrow_array() ) {
    $vega_to_ext{$row[0]} = $row[1];
  }
  $sth->finish;
  print "We have ".scalar(%vega_to_ext)." vega to external source entries\n" if($verbose);


  my $line_count = 0;
  my $xref_count = 0;

  # becouse the direct mapping have no descriptions etc
  # we have to steal these fromt he previous external source parser.
  # This is why the order states this is after the other one.

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  


  my $syn_sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($syn_sql);    
  
  my $syn_hash = $self->get_ext_synonyms($source_name);

  $sql = 'select source_id, priority_description from source where name like "' . $source_name . '"';
  $sth = $dbi->prepare($sql);
  
  $sth->execute();

  my ($ext_source_id, $desc);
  $sth->bind_columns(\$ext_source_id, \$desc);
  my @arr;
  while($sth->fetch()){
    push @arr, $ext_source_id;
  }
  $sth->finish;
  
  $sql = "select accession, label, version,  description from xref where source_id in (".join(", ",@arr).")";
#  print "$sql\n";;
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver);
  my $ext_loaded_count = 0;
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while ($sth->fetch()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc if(defined($desc));
    $ext_loaded_count++;
  }
  $sth->finish;
  if($ext_loaded_count == 0){
    warn "No point continuing no external references there\n";
    return 1;
  }

  my $ignore_count = 0;
  my $ignore_examples ="";
  my %acc;

  my $at_least_1_xref_loaded = 0;

  foreach my $key ( keys %vega_to_ext){
    if(defined($ext_to_core{$key} )){
      
      my $ext = $vega_to_ext{$key};
      if($prepend){
	my $regex = $source_name . ':';
	$ext =~ s/$regex//;
      }
      my $stable_id = $ext_to_core{$key};

      if(!defined($label{$ext})){
	$ignore_count++;
	if($ignore_count < 10){
	  $ignore_examples .= " ".$ext." (". $vega_to_ext{$key} ." )";
	}
	next;
      }

      
      my $version ="";
      $line_count++;
      if(!defined($acc{$ext})){
	my $xref_id = $self->add_xref({ acc        => $ext,
					version    => $version{$ext} ,
					label      => $label{$ext}||$ext ,
					desc       => $description{$ext},
					source_id  => $source_id,
					species_id => $species_id,
					info_type  => "DIRECT"} );

	$acc{$ext} = $xref_id;
	$xref_count++;
	if ((defined($xref_id)) and ($at_least_1_xref_loaded == 0)){
	   $at_least_1_xref_loaded = 1;
	}
      }
	
	
      $self->add_direct_xref($acc{$ext}, $stable_id, "transcript", "");

      if(defined($syn_hash->{$ext})){
	foreach my $syn (@{$syn_hash->{$ext}}){
	  $add_syn_sth->execute($acc{$ext}, $syn);
	}
      }
      
    }
  }
  $add_syn_sth->finish;

  print "Parsed $line_count $source_name identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n" if($verbose);
  if($ignore_count){
    print $ignore_count." ignoreed due to numbers no identifiers being no longer valid :- $ignore_examples\n" if($verbose);
  }
  
  if ($at_least_1_xref_loaded == 0) {    
      print "None of the external references were loaded to the core database.\n";
      return 1;
  }
  return 0;
}

1;
