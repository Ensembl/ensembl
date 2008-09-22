package XrefParser::HGNC_ENSTParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of HGNC records and assign direct xrefs
# All assumed to be linked to genes



sub run_script {

  my ($self, $file, $source_id, $species_id, $verbose) = @_;

  my ($type, $my_args) = split(/:/,$file);

  my $user = "ensro";
  my $host ="ens-staging";
  my $port = "3306";
  my $dbname = "homo_sapiens_vega_51_36m";
  my $pass;

  if($my_args =~ /host[=][>](\S+?)[,]/){
    $host = $1;
  }
  if($my_args =~ /port[=][>](\S+?)[,]/){
    $port =  $1;
  }
  if($my_args =~ /dbname[=][>](\S+?)[,]/){
    $dbname = $1;
  }
  if($my_args =~ /pass[=][>](\S+?)[,]/){
    $pass = $1;
  }

  $source_id = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","havana");


  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my $hgnc_sql = 'select tsi.stable_id, x.dbprimary_acc from xref x, object_xref ox , transcript_stable_id tsi, gene g, external_db e, transcript t  where t.gene_id = g.gene_id and g.gene_id = ox.ensembl_id and tsi.transcript_id = t.transcript_id and e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and ox.ensembl_object_type = "Gene" and e.db_name like "HGNC"';


  my %ott_to_hgnc;
  my %ott_to_enst;

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);


  my $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("ENST_CDS") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;

  $sth = $dbi2->prepare($sql);
  $sth->execute("ENST_ident") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[0]} = $row[1];
  }
  $sth->finish;
  print "We have ".scalar(%ott_to_enst)." ott to enst entries\n " if($verbose);

  $sth = $dbi2->prepare($hgnc_sql);
  $sth->execute() or croak( $dbi2->errstr() );

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

  #get the source ids for HGNC refseq, entrezgene and unitprot
  my $sql = 'select source_id, priority_description from source where name like "HGNC"';
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
	
	my $xref_id = $self->add_xref($hgnc, $version{$hgnc} , $label{$hgnc}||$hgnc , $description{$hgnc}, $source_id, $species_id);
	$xref_count++;
	
	
	$self->add_direct_xref($xref_id, $stable_id, "transcript", "");
      }
    }
  }
  print "Parsed $line_count HGNC identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n" if($verbose);
  if($ignore_count){
    print $ignore_count." ignoreed due to numbers no identifiers being no longer valid :- $ignore_examples\n" if($verbose);
  }
  
  
  return 0;
}


1;
