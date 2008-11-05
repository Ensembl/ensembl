package XrefParser::MGI_curated_transcriptParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: MGI_curated_transcriptParser.pm file <source_id> <species_id>\n\n";
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
  
  my $cuser = "ensro";
  my $chost ="ens-staging";
  my $cport = "3306";
  my $cdbname = "";
  my $cpass;

  my $vuser = "ensro";
  my $vhost ="ens-staging";
  my $vport = "3306";
  my $vdbname = "mus_musculus_vega_51_37d";
  my $vpass;

  if($my_args =~ /chost[=][>](\S+?)[,]/){
    $chost = $1;
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

  if($my_args =~ /vhost[=][>](\S+?)[,]/){
    $vhost = $1;
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


  my $xref_count = 0;


  my $clone_source_id =
    $self->get_source_id_for_source_name('Clone_based_vega_transcript');
  my $curated_source_id =
    $self->get_source_id_for_source_name('MGI_curated_transcript');
 
  my %name_to_mgi_number = %{ $self->get_label_to_acc( "MGI", $species_id ) };
  my %name_to_mgi_desc = %{ $self->get_label_to_desc( "MGI", $species_id, "descriptions") };

  my $core_sql = 'select tsi.stable_id, x.dbprimary_acc from transcript_stable_id tsi, transcript t, object_xref ox, xref x, external_db e where tsi.transcript_id = t.transcript_id and ox.ensembl_id = t.transcript_id and ox.ensembl_object_type = "Transcript" and ox.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name like "%OTTT"';

  my $vega_sql = 'select x.dbprimary_acc, x.display_label from xref x, external_db e where  x.external_db_id = e.external_db_id and e.db_name like "Vega_transcript" and display_label not like "OTT%"';
  
  my %ott_to_vega_name;
  my %ott_to_enst;
  
  my $dbi2 = $self->dbi2($chost, $cport, $cuser, $cdbname, $cpass);
  if(!defined($dbi2)){
    return 1;
  }
  

  my $sth = $dbi2->prepare($core_sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[1]} = $row[0];
  }
  $sth->finish;
  
  my $dbi3 = $self->dbi2($vhost, $vport, $vuser, $vdbname, $vpass);
  if(!defined($dbi3)){
    return 1;
  }
  
  
  $sth = $dbi3->prepare($vega_sql); 
  $sth->execute() or croak( $dbi3->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    # [0] OTTMUST...,   [1] Mrpl15-001 etc 
    my $desc= "";
    if(defined($ott_to_enst{$row[0]})){
      my $id = $curated_source_id;
      my $name = $row[1];
      $name =~ s/^WU://;
      my $prim_acc = $name;
      my $stable_id = $ott_to_enst{$row[0]};
      if($name =~ /[.]/){
        $id = $clone_source_id;
      }
      else{
	my $mgi_name = $name;
	# find MGI name 
	my($mgi_bit, $num) = split(/-\d\d\d/,$name);
	if(defined($name_to_mgi_number{$mgi_bit})){
	  $prim_acc = $name_to_mgi_number{$mgi_bit};
	}
	if(defined($name_to_mgi_desc{$mgi_bit})){
	  $desc = $name_to_mgi_desc{$mgi_bit};
	}
	else{
	  print "$mgi_bit has no description\n" if($verbose);
	}
      }
      my $xref_id = $self->add_xref($prim_acc, "" , $name , $desc, $id, $species_id);
      $xref_count++;
      
      $self->add_direct_xref($xref_id, $ott_to_enst{$row[0]}, "Transcript", "");

      
    }
  }

  print "$xref_count direct xrefs succesfully parsed\n" if($verbose);

  #Finally add the synonyms:-
  my $synonym_sql = (<<SYNO);
SELECT x2.xref_id, s.synonym
  FROM synonym s
    INNER JOIN xref x1 ON  x1.xref_id = s.xref_id
    INNER JOIN xref x2 ON  x2.accession = x1.accession 
    INNER JOIN source s1 ON s1.source_id = x1.source_id
    INNER JOIN source s2 ON s2.source_id = x2.source_id    
      WHERE x2.xref_id != x1.xref_id
	AND s2.name = "MGI_curated_transcript" 
        AND s1.name = "MGI" 
        AND s1.priority_description = "descriptions"
SYNO
 

  $sth = $self->dbi()->prepare($synonym_sql);
  
  $sth->execute() or croak( $self->dbi()->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $self->add_synonym($row[0], $row[1]);
  }
  $sth->finish;  


  return 0;
}





1;

