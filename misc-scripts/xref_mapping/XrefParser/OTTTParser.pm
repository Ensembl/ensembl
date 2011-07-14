package XrefParser::OTTTParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;


sub run_script {
  my $self = shift;
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;

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

 
  my $sql = 'select tsi.stable_id, x.display_label from xref x, object_xref ox , transcript_stable_id tsi, external_db e where e.external_db_id = x.external_db_id and x.xref_id = ox.xref_id and tsi.transcript_id = ox.ensembl_id and e.db_name like ?';


  my %ott_to_enst;
  
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }
  
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

  my $xref_count = 0;
  foreach my $ott (keys %ott_to_enst){
  
    my $xref_id = $self->add_xref($ott, "" , $ott , "", $source_id, $species_id, "DIRECT");
    $xref_count++;
    
    
    $self->add_direct_xref($xref_id, $ott_to_enst{$ott}, "transcript", "");
    
  }
}




#sub dbi2{

#  my $self = shift;
#  my ($host, $port, $user, $dbname, $pass) = @_;

#  if ( !defined $dbi2 || !$dbi2->ping() ) {
#    my $connect_string =
#      sprintf( "dbi:mysql:host=%s;port=%s;database=%s",
#	       $host, $port, $dbname );

#    $dbi2 =
#      DBI->connect( $connect_string, $user, $pass,
#		    {
#		     'RaiseError' => 1 } )
#	or croak( "Can't connect to database: " . $DBI::errstr );
#    $dbi2->{'mysql_auto_reconnect'} = 1; # Reconnect on timeout
#  }
    
#  return $dbi2;
#}

1;

