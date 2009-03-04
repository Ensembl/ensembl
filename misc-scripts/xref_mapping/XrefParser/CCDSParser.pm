package XrefParser::CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts
# The same CCDS may be linked to more than one transcript, but need to only
# add the xref once, so check if it already exists before adding it.

sub run_script {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id  = shift;
  my $species_id = shift;
  my $verbose    = shift;

  my $user = "ensro";
  my $host;
  my $port = 3306;
  my $dbname;
  my $pass;
  my $tran_name;


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
  if($file =~ /tran_name[=][>](\S+?)[,]/){
    $tran_name = $1;
  }

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);

  if(!defined($dbi2)){
    return 1;
  }


  my $line_count = 0;
  my $xref_count = 0;

 my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND version=? AND source_id=$source_id AND species_id=$species_id");


#
# Need to get the stable_id via the ENST xrefs!!!!!
#

#
# 
#

  my $sql = 'select ox.ensembl_id, x.dbprimary_acc from object_xref ox, xref x, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like "ENST" and x.dbprimary_acc like "'.$tran_name.'%"'; 


  my %trans_id_to_stable_id;
  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $trans_id_to_stable_id{$row[0]} = $row[1];
#    print $row[0]."\t".$row[1]."\n";
  }
  $sth->finish;

  $sql = 'select ox.ensembl_id, x.dbprimary_acc, x.display_label, x.version from object_xref ox, xref x, external_db e where x.xref_id = ox.xref_id and x.external_db_id = e.external_db_id and e.db_name like "CCDS"';
  
  $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  my $xref_count = 0;
  my $direct_count=0;
  while ( my @row = $sth->fetchrow_array() ) {
#    print "Processing ".$row[0]."\n";
    if(defined($trans_id_to_stable_id{$row[0]})){
      my $acc = $row[1];
      my $version = $row[3];
      my $display_label = $row[2];
      my $tran_id = $row[0];
      my $stable_id = $trans_id_to_stable_id{$tran_id};

      # check if an xref already exists
      $xref_sth->execute($acc, $version);
      my $xref_id = ($xref_sth->fetchrow_array())[0];
      if (!$xref_id) {
	$xref_id = $self->add_xref($acc, $version, $display_label, "", $source_id, $species_id, "DIRECT");
	$xref_count++;
      }
      
      $self->add_direct_xref($xref_id, $stable_id, "Transcript", "");
      $direct_count++;
    }
    else{
      print "Could not find trans_id_to_stable_id for ".$row[0]."\n";
    }
  }

  print "Parsed CCDS identifiers from $file, added $xref_count xrefs and $direct_count direct_xrefs\n" if($verbose);

  return 0;
}

1;
