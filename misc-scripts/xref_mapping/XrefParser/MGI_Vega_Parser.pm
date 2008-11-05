package XrefParser::MGI_Vega_Parser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: MGI_Vega_Parser.pm file <source_id> <species_id>\n\n";
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
 

  #
  # need to get label and derscriptions fro primary acc.
  #

  my %mgi_to_label;
  my %mgi_to_desc;
  my %mgi_syn;

  my $sth = $self->dbi()->prepare("SELECT x.accession, x.label, x.description from xref x, source s where x.source_id = s.source_id and s.name like 'MGI' and s.priority_description like 'descriptions'");
  
  $sth->execute() or croak( $self->dbi()->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $mgi_to_label{$row[0]} = $row[1];
    $mgi_to_desc{$row[0]} = $row[2];
  }
  $sth->finish;

  #
  # Also add synonyms
  #

  $sth = $self->dbi()->prepare("SELECT sy.synonym, x.accession from xref x, source s, synonym sy where sy.xref_id = x.xref_id and x.source_id = s.source_id and s.name like 'MGI' and s.priority_description like 'descriptions'");
  
  $sth->execute() or croak( $self->dbi()->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $mgi_syn{$row[0]} = $row[1]; 
  }
  $sth->finish;




  my $core_sql = 'select tsi.stable_id, x.dbprimary_acc from transcript_stable_id tsi, transcript t, object_xref ox, xref x, external_db e where tsi.transcript_id = t.transcript_id and ox.ensembl_id = t.transcript_id and ox.ensembl_object_type = "Transcript" and ox.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name like "%OTTT"';


  my %ott_to_enst;
  
  my $dbi2 = $self->dbi2($chost, $cport, $cuser, $cdbname, $cpass);
  if(!defined($dbi2)){
    return 1;
  }
  
  
  $sth = $dbi2->prepare($core_sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $ott_to_enst{$row[1]} = $row[0];
  }
  $sth->finish;
  

  #
  # get the enst->ensg mappings.
  #
  my %enst_to_ensg;
  $sth = $dbi2->prepare("select gsi.stable_id, tsi.stable_id from transcript t, gene_stable_id gsi, transcript_stable_id tsi where tsi.transcript_id = t.transcript_id and t.gene_id = gsi.gene_id"); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    $enst_to_ensg{$row[1]} = $row[0];
  }
  $sth->finish;


  #
  # Get the ott -> mgi mappings
  #


  my $vega_sql = (<<VSQL);
  SELECT DISTINCT(tsi.stable_id) , x.dbprimary_acc, x.display_label    
    FROM        transcript_stable_id tsi      
     INNER JOIN transcript t              ON tsi.transcript_id = t.transcript_id      
     INNER JOIN gene g                    ON g.gene_id = t.gene_id      
     INNER JOIN object_xref ox            ON ox.ensembl_id = g.gene_id      
     INNER JOIN xref x                    ON x.xref_id = ox.xref_id      
     INNER JOIN external_db e             ON e.external_db_id = x.external_db_id      
     WHERE ox.ensembl_object_type = "Gene"         
       AND e.db_name like "MGI"
VSQL

  my $dbi3 = $self->dbi2($vhost, $vport, $vuser, $vdbname, $vpass);
  if(!defined($dbi3)){
    return 1;
  }


  my %seen;
  $sth = $dbi3->prepare($vega_sql); 
  $sth->execute() or croak( $dbi3->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    # [0] OTTMUST...,   [1] MGI:123456 [2]  Asx15 etc 
    my $desc= "";
    my $prim_acc = $row[1];
    if(defined($ott_to_enst{$row[0]})){
      my $tran_stable_id = $ott_to_enst{$row[0]};
      my $name = $prim_acc;
      my $desc = "";
      my $label = "";
      if(defined($mgi_to_desc{$name})){
	$desc = $mgi_to_desc{$name};
	$label  = $mgi_to_label{$name};
      }	
      elsif( defined( $mgi_syn{$row[2]} ) and defined( $mgi_to_desc{$mgi_syn{$row[2]}})){ # synonym
        $prim_acc = $mgi_syn{$row[2]};
	$desc = $mgi_to_desc{$prim_acc};
        $label = $mgi_to_label{$name};
      }
      else{
	print "VEGA: $name [".$row[2]."} has no description\n" if($verbose);
      }
      my $xref_id = $self->add_xref($prim_acc, "" , $label , $desc, $source_id, $species_id);
      my $ensg = $enst_to_ensg{$tran_stable_id};
      if(!defined($seen{$xref_id.$ensg})){
	$xref_count++;
	$self->add_direct_xref($xref_id, $ensg , "Gene", "");
	$seen{$xref_id.$ensg} = 1;
      }	
    }
  }

  print "$xref_count direct xrefs succesfully parsed\n" if($verbose);
 

# Done in the mapper
#  #
#  # Finally addd the synonyms
#  #

#  my $synonym_sql = (<<SYNO);
#SELECT x2.xref_id, s.synonym
#  FROM synonym s
#    INNER JOIN xref x1 ON  x1.xref_id = s.xref_id
#    INNER JOIN xref x2 ON  x2.accession = x1.accession 
#    INNER JOIN source s1 ON s1.source_id = x1.source_id
#    INNER JOIN source s2 ON s2.source_id = x2.source_id    
#      WHERE x2.xref_id != x1.xref_id
#	AND s2.name = "MGI" 
#	AND s2.priority_description = "vega" 
#        AND s1.name = "MGI" 
#        AND s1.priority_description = "descriptions"
#SYNO

#  $sth = $self->dbi()->prepare($synonym_sql);
  
#  $sth->execute() or croak( $self->dbi()->errstr() );
#  while ( my @row = $sth->fetchrow_array() ) {
#    $self->add_synonym($row[0], $row[1]);
#  }
#  $sth->finish;  


return 0;
}





1;

