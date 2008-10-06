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
  

  my %refseq_to_hgnc;
  my %refseq_to_ccds;

  my $response = $ua->get($wget);
  
  if ( !$response->is_success() ) {
    die $response->status_line;
  }
  else{
    my @lines = split(/\n/,$response->content);
    foreach my $line (@lines){
      my($hgnc, $refseq) = split(/\s+/,$line);
      if(defined($refseq) and $refseq ne ""){
	$refseq_to_hgnc{$refseq} = $hgnc;
      }
    }
    
  }

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }

  my $sql = "select  cu.ccds_uid, a.nuc_acc from Accessions a, Accessions_GroupVersions agv, GroupVersions  gv, CcdsUids cu where a.accession_uid = agv.accession_uid and a.organization_uid=1 and agv.group_version_uid=gv.group_version_uid and gv.ccds_status_val_uid in (3) and cu.group_uid=gv.group_uid  order by gv.ccds_status_val_uid, cu.ccds_uid";


  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $refseq_to_ccds{$row[1]} = $row[0];
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

  my $sql = 'select source_id, priority_description from source where name like "HGNC"';
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
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
  }
  $sth->finish;
 
  $sql = 'select x.accession, d.ensembl_stable_id, d.type 
            from xref x, direct_xref d, source s 
             where s.source_id = x.source_id and 
                   x.xref_id = d.general_xref_id and s.name like "CCDS"'; 
 
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($access, $stable_id, $type);
  $sth->bind_columns(\$access, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  while (my @row = $sth->fetchrow_array()) {
    $ensembl_stable_id{$access} = $stable_id;
    $ensembl_type{$access} = $type;
  }
  $sth->finish;
  
  my $line_count = 0;
  my $xref_count = 0;
  my %seen;
  my $ignore_count = 0;
  my $ignore_examples ="";

  foreach my $refseq (keys %refseq_to_ccds){
    if(defined($refseq_to_hgnc{$refseq})){
      my $ccds = $refseq_to_ccds{$refseq};
      my $hgnc = $refseq_to_hgnc{$refseq};
      my $key = "CCDS".$ccds;
      if(defined($ensembl_stable_id{$key})){
	my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , 
				      $description{$hgnc}, $source_id, $species_id);
	$self->add_direct_xref($xref_id, $ensembl_stable_id{$key}, $ensembl_type{$key}, "");
	$xref_count++;
      }
    }
  }


}

1;
