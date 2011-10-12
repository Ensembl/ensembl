package XrefParser::CCDSParser;

use strict;
use Carp;
use DBI;

use base qw( XrefParser::BaseParser );
use XrefParser::Database;
# Parse file of CCDS records and assign direct xrefs
# All assumed to be linked to transcripts
# The same CCDS may be linked to more than one transcript, but need to only
# add the xref once, so check if it already exists before adding it.

sub run_script {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file  as pairs";
  }
  $verbose |=0;

  my $user = "ensro";
  my $host;
  my $port = 3306;
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

  my $ccds_db =  XrefParser::Database->new({ host   => $host,
					     port   => $port,
					     user   => $user,
					     dbname => $dbname,
					     pass   => $pass});

  my $dbi2 = $ccds_db->dbi();

  if(!defined($dbi2)){
    return 1;
  }


  my $line_count = 0;
  my $xref_count = 0;

  my $sql =(<<'SCD');
SELECT t.stable_id, x.dbprimary_acc 
  FROM xref x, object_xref ox, transcript t, external_db e
    WHERE x.xref_id=ox.xref_id AND
          ox.ensembl_object_type = "Transcript" AND
          ox.ensembl_id = t.transcript_id AND
          e.external_db_id = x.external_db_id AND
          e.db_name like "Ens_%_transcript"
SCD

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
      $xref_id = $self->add_xref({ acc        => $acc,
				   version    => $version,
				   label      => $display_label,
				   source_id  => $source_id,
				   species_id => $species_id,
				   info_type  => "DIRECT"} );
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
