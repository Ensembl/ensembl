package XrefParser::EG_DBParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: EG_DBParser.pm file <source_id> <species_id> <verbose>\n\n";
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

  $source{"BROAD_U_maydis"}     = $self->get_source_id_for_source_name('BROAD_U_maydis')    || die "Could not get source_id for BROAD_U_maydis\n";
 
  print STDERR "source_id for Broad Umaydis: " . $source{"BROAD_U_maydis"} . "\n";

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
      print STDERR "failed to connect to EG_Xrefs database!\n";
    return 1;
  }

  my $sql=(<<SQL);
  SELECT gene_stable_id, transcript_stable_id, gene_dbid, transcript_dbid, source, xref_name, xref_primary_id, xref_description
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
