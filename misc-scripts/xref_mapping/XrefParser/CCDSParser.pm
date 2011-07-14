package XrefParser::CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts
# The same CCDS may be linked to more than one transcript, but need to only
# add the xref once, so check if it already exists before adding it.

sub run_script {

  my $self = shift;
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

  my $sql = 'select tsi.stable_id, x.dbprimary_acc from xref x, object_xref ox, transcript_stable_id tsi, external_db e where x.xref_id=ox.xref_id and  ox.ensembl_object_type = "Transcript" and ox.ensembl_id = tsi.transcript_id and e.external_db_id = x.external_db_id and e.db_name like "Ens_%_transcript"';


  my %seen;

  my $sth = $dbi2->prepare($sql) or "Could not prepare sql $sql\n";;
  $sth->execute() or die "Could not execute $sql\n";;
  my $xref_count = 0;
  my $direct_count=0;
  my ($stable_id, $display_label);
  $sth->bind_columns( \$display_label,\$stable_id);
  while ( $sth->fetch ) {

    my ($acc, $version) = split (/\./,$display_label);

    my $xref_id;
    if (!defined($seen{$display_label})) {
      $xref_id = $self->add_xref($acc, $version, $display_label, "", $source_id, $species_id, "DIRECT");
      $xref_count++;
      $seen{$display_label} = $xref_id;
    }
    else{
      $xref_id = $seen{$display_label};
    }

    $self->add_direct_xref($xref_id, $stable_id, "Transcript", "");
    $direct_count++;
  }

  print "Parsed CCDS identifiers from $file, added $xref_count xrefs and $direct_count direct_xrefs\n" if($verbose);

  return 0;
}

1;
