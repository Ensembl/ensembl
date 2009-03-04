package XrefParser::RefSeq_CCDSParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Parse file of Refseq records and assign direct xrefs

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

  my $dna_pred = XrefParser::BaseParser->get_source_id_for_source_name("RefSeq_dna_predicted");

  # becouse the direct mapping have no descriptions etc
  # we have to steal these from the previous Refseq parser.

  my %label;
  my %version;
  my %description;

  my $dbi = $self->dbi();  
  my $sql = "select xref.accession, xref.label, xref.version,  xref.description from xref, source where xref.source_id = source.source_id and source.name = 'RefSeq_dna'";
  my $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($acc, $lab, $ver, $desc);
  $sth->bind_columns(\$acc, \$lab, \$ver, \$desc);
  while (my @row = $sth->fetchrow_array()) {
    $label{$acc} = $lab;
    $version{$acc} = $ver;
    $description{$acc} = $desc;
  }
  $sth->finish;
 
  $sql = 'select x.accession, x.xref_id, d.ensembl_stable_id, "Transcript"
            from xref x, transcript_direct_xref d, source s 
             where s.source_id = x.source_id and 
                   x.xref_id = d.general_xref_id and s.name like "CCDS"'; 
 
  $sth = $dbi->prepare($sql);
  $sth->execute();
  my ($access, $old_xref_id, $stable_id, $type);
  $sth->bind_columns(\$access, \$old_xref_id, \$stable_id, \$type);
  my %ensembl_stable_id;
  my %ensembl_type;
  my %old_xref;
  while (my @row = $sth->fetchrow_array()) {
      $ensembl_stable_id{$access} = $stable_id;
      $ensembl_type{$access} = $type;
      $old_xref{$access} = $old_xref_id; 
  }
  $sth->finish;
  
 

  my $line_count = 0;
  my $xref_count = 0;
  my %seen;
  my %old_to_new;

  my $dbi2 = $self->dbi2($host, $port, $user, $dbname, $pass);
  if(!defined($dbi2)){
    return 1;
  }
  
  my $sql = "select  cu.ccds_uid, a.nuc_acc from Accessions a, Accessions_GroupVersions agv, GroupVersions  gv, CcdsUids cu where a.accession_uid = agv.accession_uid and a.organization_uid=1 and agv.group_version_uid=gv.group_version_uid and gv.ccds_status_val_uid in (3) and cu.group_uid=gv.group_uid  order by gv.ccds_status_val_uid, cu.ccds_uid";
  
  
  my $sth = $dbi2->prepare($sql); 
  $sth->execute() or croak( $dbi2->errstr() );
  while ( my @row = $sth->fetchrow_array() ) {
    my $ccds = $row[0];
    my $refseq = $row[1];
    
    $line_count++;
    if(!defined($seen{$refseq})){
      $seen{$refseq} = 1;
      my $key = "CCDS".$ccds;
      if(defined($ensembl_stable_id{$key})){
	my $new_source_id = $source_id;
	if($refseq =~ /^XM/){
	  $new_source_id = $dna_pred;
	}
	my $xref_id = $self->add_xref($refseq, $version{$refseq} , $label{$refseq}||$refseq , 
				      $description{$refseq}, $new_source_id, $species_id, "DIRECT");
	$self->add_direct_xref($xref_id, $ensembl_stable_id{$key}, $ensembl_type{$key}, "");
	$old_to_new{$old_xref{$refseq}} = $xref_id;
	$xref_count++;
      }
    }
  }
  
  #for each one seen get all its dependent xrefs and load them fro the new one too;

  my $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref VALUES(?,?,?,?)");
  my $get_dependent_xref_sth = $dbi->prepare("SELECT dependent_xref_id, linkage_annotation "
					    .  "FROM  dependent_xref where master_xref_id = ?");

  foreach my $old_xref (keys %old_to_new){
      my $linkage;
      my $dependent_id;
      $get_dependent_xref_sth->execute($old_xref);
      $get_dependent_xref_sth->bind_columns(\$dependent_id, \$linkage);
      while(my @row = $get_dependent_xref_sth->fetchrow_array()){
	  $add_dependent_xref_sth->execute($old_to_new{$old_xref}, $dependent_id, $linkage, $source_id); 
      }   
  }


  print "Parsed $line_count RefSeq_dna identifiers from $file, added $xref_count xrefs and $xref_count direct_xrefs  from $line_count lines.\n" if ($verbose);


  return 0;

}

1;
