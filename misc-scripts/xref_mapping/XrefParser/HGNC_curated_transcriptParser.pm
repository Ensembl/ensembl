package XrefParser::HGNC_curated_transcriptParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: HGNC_curated_transcriptParser.pm file <source_id> <species_id>\n\n";
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

  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');
  my $curated_source_id =
    $self->get_source_id_for_source_name('HGNC_curated_transcript');
 
  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my %ott_to_vega_name;
  my %ott_to_enst;
  
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  

  my $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
  $sth->execute("Vega_transcript") or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_vega_name{$row[0]} = $row[1];
  }
  $sth->finish;

  $sth = $dbi2->prepare($sql);   # funny number instead of stable id ?????
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

  my $xref_count = 0;
  foreach my $ott (keys %ott_to_enst){
    if(defined($ott_to_vega_name{$ott})){
      my $id = $curated_source_id;
      my $name  = $ott_to_vega_name{$ott};
      my $primary_acc = $name;
      if($name =~ /[.]/){
	$id = $clone_source_id;
      }
      my $xref_id = $self->add_xref($name, "" , $name , "", $id, $species_id);
      $xref_count++;
      
      
      $self->add_direct_xref($xref_id, $ott_to_enst{$ott}, "transcript", "");
      
    }
  }

  print "xref_count direct xrefs succesfully parsed\n" if($verbose);
  return 0;
}





1;

