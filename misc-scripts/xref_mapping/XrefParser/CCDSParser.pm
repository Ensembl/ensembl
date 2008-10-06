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
  my $port;
  my $dbname;
  my $pass;


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

  my $line_count = 0;
  my $xref_count = 0;

  my $xref_sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND version=? AND source_id=$source_id AND species_id=$species_id");

  my $sql = "select display_label, stable_id from transcript_stable_id, translation, object_xref, xref where transcript_stable_id.transcript_id = translation.transcript_id && translation.translation_id = object_xref.ensembl_id && object_xref.ensembl_object_type = 'Translation' && object_xref.xref_id = xref.xref_id && external_db_id = 3800";

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);

  if(!defined($dbi2)){
    return 1;
  }

  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $stable_id = $row[1];
    my $ccds = $row[0];

    my ($acc, $version) = split (/\./, $ccds);
    $line_count++;

    # check if an xref already exists
    $xref_sth->execute($acc, $version);
    my $xref_id = ($xref_sth->fetchrow_array())[0];
    if (!$xref_id) {
      $xref_id = $self->add_xref($acc, $version, $ccds, "", $source_id, $species_id);
      $xref_count++;
    }

    $self->add_direct_xref($xref_id, $stable_id, "transcript", "");

  }

  print "Parsed $line_count CCDS identifiers from $file, added $xref_count xrefs and $line_count direct_xrefs\n" if($verbose);

  return 0;
}

1;
