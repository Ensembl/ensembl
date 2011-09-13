package XrefParser::EG_DBParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub run_script {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $file         = $ref_arg->{file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $file) ){
    croak "Need to pass source_id, species_id and file as pairs";
  }
  $verbose |=0;

  print STDERR "parsing EG_DB_Xrefs...\n";

  my ($type, $my_args) = split(/:/,$file);

  my $user  ="ensro";
  my $host;
  my $port;
  my $dbname;
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
  if($my_args =~ /user[=][>](\S+?)[,]/){
    $user = $1;
  }

  print STDERR "species_id, $species_id\n";

  my %source;

  # Todo: get the whole list of sources for All EG set
  
  my $sources_aref = ['BROAD_U_maydis','CADRE','CADRE_Afum_A1163','AspGD','ENA_GENE','BROAD_F_oxysporum','BROAD_G_moniliformis','BROAD_G_zeae','GeneDB','BROAD_P_infestans','phyra_jgi_v1.1','physo1_jgi_v1.1','PGD_GENE','phatr_jgi_v2_bd','phatr_jgi_v2','thaps_jgi_v2','thaps_jgi_v2_bd'];
  foreach my $source_name (@$sources_aref) {
      $source{$source_name}     = $self->get_source_id_for_source_name($source_name)    || die "Could not get source_id for $source_name\n";
  }
 
  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
      print STDERR "failed to connect to EG_Xrefs database!\n";
    return 1;
  }

  my $sql=(<<"SQL");
  SELECT gene_stable_id, transcript_stable_id, gene_dbid, transcript_dbid, 
         source, xref_name, xref_primary_id, xref_description
    FROM EG_Xref
      WHERE taxonomy_id = $species_id
SQL

  my ($gene_stable_id, $transcript_stable_id, $gene_dbid, $transcript_dbid, $source_name, $label, $acc, $desc);
  my $sth = $dbi2->prepare($sql);
  $sth->execute() or croak( $dbi2->errstr() );
  $sth->bind_columns(\$gene_stable_id, \$transcript_stable_id, \$gene_dbid, \$transcript_dbid, \$source_name, \$label, \$acc, \$desc);
  my $added = 0;
  while ( $sth->fetch() ) {

    my ($description,$junk) = split("[[]Source:",$desc);

    if(!defined($source{$source_name})){
      print STDERR "Could not find source_id for source $source_name for xref $acc\n";
      next;
    }
    my $xref_id = $self->get_xref($acc,$source{$source_name}, $species_id);
    if(!defined($xref_id)){
      $xref_id = $self->add_xref($acc,"",$label,$description,$source{$source_name}, $species_id,"DIRECT");
      $added++;
    }
    my $transcript_id = $transcript_dbid;
    if(defined($transcript_stable_id) and $transcript_stable_id ne ""){
      $transcript_id = $transcript_stable_id;
    }
    my $gene_id = $gene_dbid;
    if(defined($gene_stable_id) and $gene_stable_id ne ""){
	$gene_id = $gene_stable_id;
    }
    
    #$self->add_direct_xref($xref_id, $transcript_id, "Transcript", "") if (defined($transcript_id));    
    $self->add_direct_xref($xref_id, $gene_id, "Gene", "") if (defined($gene_id));  
  }
  $sth->finish;

  print "Added $added Xrefs for EGs\n" if($verbose);
  return 0;
}

1;
