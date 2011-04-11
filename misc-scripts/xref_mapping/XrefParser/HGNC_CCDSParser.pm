package XrefParser::HGNC_CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes

sub run_script {

  my ($self, $file, $source_id, $species_id, $verbose) = @_;

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
#      print "ccds:$ccds\n";
      my @ccds_list = split(/, /,$ccds);
      foreach my $c (@ccds_list){
#	print $c."\t".$hgnc."\n";
	$ccds_to_hgnc{$c} = $hgnc;
      }
    }
    
  }
  
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }



  my $sql = 'select ox.ensembl_id, x.dbprimary_acc from object_xref ox, xref x, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like "Ens_%_transcript" and x.dbprimary_acc like "ENST%"'; 


  my %trans_id_to_stable_id;
  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $trans_id_to_stable_id{$row[0]} = $row[1];
  }
  $sth->finish;  


  
  $sql = 'select ox.ensembl_id, x.display_label from object_xref ox, xref x, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like "CCDS"'; 

  my %ccds_to_stable_id;
  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    if(defined($trans_id_to_stable_id{$row[0]})){
      $ccds_to_stable_id{$row[1]} = $trans_id_to_stable_id{$row[0]};
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

  my $dbi = $self->dbi();  

  my $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $dbi->prepare($sql);    
  
  my $syn_hash = $self->get_ext_synonyms("HGNC");

  $sql = 'select source_id, priority_description from source where name like "HGNC"';
  my $sth = $dbi->prepare($sql);
  
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
  my ($acc, $lab, $ver, $desc);
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
  foreach my $ccds (keys %ccds_to_stable_id){
    if(defined($ccds_to_hgnc{$ccds})){
      my $hgnc = $ccds_to_hgnc{$ccds};
      $hgnc =~ s/HGNC://;
      my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , 
				      $description{$hgnc}, $source_id, $species_id, "DIRECT");
      $self->add_direct_xref($xref_id, $ccds_to_stable_id{$ccds}, "Transcript", "");
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
  print "$no_ccds_to_hgnc missed as no hgnc for the ccds. Added $xref_count HGNC xrefs via CCDS\n" if($verbose);
  return 0;
}

1;
