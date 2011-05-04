package XrefMapper::BasicMapper;

use strict;
use warnings;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;


=head2 new

  Description: Constructor for BasicMapper.
  Returntype : BasicMapper
  Exceptions : none
  Caller     : general

=cut

sub new{
  my($class, @args) = @_;

  my $self ={};
  bless $self,$class;
  return $self;
}



=head2 xref

  Arg [1]    : (optional)
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core.
               info for the xref database.
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub xref{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_xref} = $arg );
  return $self->{_xref};
}

=head2 farm_queue

  Arg [1]    : (optional)
  Example    : $mapper->farm_queue("long");
  Description: Getter / Setter for the farm queue.
  Returntype : string
  Exceptions : none

=cut

sub farm_queue{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_queue} = $arg );
  return $self->{_queue};
}

=head2 exonerate

  Arg [1]    : (optional)
  Example    : $mapper->exonerate("/usr/local/exonerate1.1.1");
  Description: Getter / Setter for the exonerate executable with full path.
  Returntype : string
  Exceptions : none

=cut

sub exonerate{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_exonerate} = $arg );
  return $self->{_exonerate};
}

=head2 core

  Arg [1]    : (optional)
  Example    : $mapper->core($new_core);
  Description: Getter / Setter for the core.
               info for the ensembl core database.
  Returntype : XrefMapper::db
  Exceptions : none

=cut

sub core{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_core} = $arg );
  return $self->{_core};
}


sub add_meta_pair {

  my ($self, $key, $value) = @_;

  my $sth = $self->xref->dbc->prepare('insert into meta (meta_key, meta_value, date) values("'.$key.'", "'.$value.'", now())');
  $sth->execute;
  $sth->finish;

}


sub xref_latest_status { 
  my $self = shift;
  my $verbose = shift || 0;

  my $sth = $self->xref->dbc->prepare("select id, status, date from process_status order by id");
  
  $sth->execute();
  my ($xref_id, $acc);
  my ($id, $status, $date);
  $sth->bind_columns(\$id, \$status,\$date);
  while($sth->fetch){
    print "$status\t$date\n" if($verbose and $self->verbose);
  }
  return $status;

}

sub get_meta_value {
  my ($self, $key) = @_;

  my $sth = $self->xref->dbc->prepare('select meta_value from meta where meta_key like "'.$key.'" order by meta_id');
  
  $sth->execute();
  my $value;
  $sth->bind_columns(\$value);
  while($sth->fetch){   # get the last one
  }
  $sth->finish;

  return $value;  
}

sub process_file {
  my $self = shift;
  my $file = shift;
  my $verbose = shift;
  my $no_xref = shift;

  open(FILE, $file) or die ("\nCannot open input file '$file':\n $!\n");
  

  my $xref=undef;
  my $ensembl=undef;
  my $type;
  
  my %xref_hash=();
  my %species_hash=();
  my %farm_hash=();
  
  while( my $line = <FILE> ) {
    
    chomp($line);
    next if $line =~ /^#/;
    next if !$line;
    
    #  print $line."\n";
    my ($key, $value) = split("=",$line);
    if($key eq "species"){
      $type = "species";
      $species_hash{'species'} = $value;
    }
    elsif($key eq "xref"){
      $type = "xref";
    }
    elsif($key eq "farm"){
      $type = "farm";
    }
    elsif($type eq "species"){ # processing species data
      $species_hash{lc($key)} = $value;
    }
    elsif($type eq "xref"){    # processing xref data
      $xref_hash{lc($key)} = $value;
    }
    elsif($type eq "farm"){
      $farm_hash{lc($key)} = $value;
    }
  }
  

  my $value = $species_hash{'species'};
  my $taxon = $species_hash{'taxon'};

  if ($value !~ /_/) {
    print STDERR "\'$value\' is not a recognised species - please use full species name (e.g. homo_sapiens) in $file\n";
    exit(1);
  }
    
  my $mapper;
  my $module;
  my $class = "XrefMapper/$value.pm";
  eval {
    require $class;
  };
  if($@) {
    if ($@ =~ /Can\'t locate $class/) {
      $class = "XrefMapper/$taxon.pm";
      eval {
         require $class;
      };
      if($@) {
         if ($@ =~ /Can\'t locate $class/) {
             warn("Did not find a specific mapping module XrefMapper::$value or XrefMapper::$taxon - using XrefMapper::BasicMapper instead\n") if(defined($verbose) and $verbose);
             require XrefMapper::BasicMapper;
	     $module = "BasicMapper";
	 } else { die "$@"; }
       } else {
         $module = $taxon; 
       }
    } else {
      die "$@";
    }
      
  } else{
    $module = $value;
  }
  
  
  $mapper = "XrefMapper::$module"->new();

  if(defined($farm_hash{'queue'})){
    $mapper->farm_queue($farm_hash{'queue'});
  }
  if(defined($farm_hash{'exonerate'})){
    $mapper->exonerate($farm_hash{'exonerate'});
  }
  

  if(defined($xref_hash{host}) and !defined($no_xref)){
    my ($host, $user, $dbname, $pass, $port);
    $host = $xref_hash{'host'};
    $user = $xref_hash{'user'};
    $dbname = $xref_hash{'dbname'};
    if(defined($xref_hash{'password'})){
      $pass = $xref_hash{'password'};
    }
    else{
      $pass = '';
    }
    if(defined($xref_hash{'port'})){
      $port = $xref_hash{'port'};
    }
    else{
      $port = 3306;
    }
    
    $xref = new XrefMapper::db(-host => $host,
			       -port => $port,
			       -user => $user,
			       -pass => $pass,
			       -group   => 'core',
			       -dbname => $dbname);
    
    $mapper->xref($xref);
    $mapper->add_meta_pair("xref", $host.":".$dbname);
    if(defined($xref_hash{'dir'})){
      $xref->dir($xref_hash{'dir'});
      if(!-d $xref_hash{'dir'}){
	die "directory ".$xref_hash{'dir'}." does not exist please create this\n";
      }
    }
    else{
      die "No directory specified for the xref fasta files\n";
    }	
    
  }
  elsif(!defined($no_xref)){
    die "No host name given for xref database\n";
  }
  else{
    print "No xref database is too be used\n" if ($verbose)
  }
  
  
  if(defined($species_hash{'species'})){

    my ($host, $port, $user, $dbname, $pass);
    $host = $species_hash{'host'};
    $user = $species_hash{'user'};
    $dbname = $species_hash{'dbname'};
    if(defined($species_hash{'password'})){
      $pass = $species_hash{'password'};
    }
    else{
      $pass = '';
    }
    if(defined($species_hash{'port'})){
      $port = $species_hash{'port'};
    }
    else{
      $port = '';
    }

    my $core = new XrefMapper::db(-host => $host,
				  -port => $port,
				  -user => $user,
				  -pass => $pass,
				  -group   => 'core',
				  -dbname => $dbname);
    
    $mapper->core($core);
    if(!defined($no_xref)){
      $mapper->add_meta_pair("species", $host.":".$dbname);

      if(defined($species_hash{'dir'})){
	$core->dir($species_hash{'dir'});
	if(!-d $species_hash{'dir'}){
	  die "directory ".$species_hash{'dir'}." does not exist please create this\n";
	}
      }    
      else{
	die "No directory specified for the ensembl fasta files\n";
      }	
    }
    $core->species($value);
  }

  return $mapper;
}


=head2 dumpcheck

  Arg [1]    : (optional)
  Example    : $mapper->dumpcheck("yes");
  Description: Getter / Setter for dumpcheck.
               If set the mapper will not dump fasta files
               if they exist already.
  Returntype : scalar
  Exceptions : none

=cut

sub dumpcheck {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_dumpcheck} = $arg );
  return $self->{_dumpcheck};
}

sub nofarm {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_nofarm} = $arg );
  return $self->{_nofarm};
}

sub verbose {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_verbose} = $arg );
  return $self->{_verbose};
}

sub species_id {
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_species_id} = $arg );
  return $self->{_species_id};
}

sub get_id_from_species_name {
 my ($self, $species_name) = @_;

 my $sql = "select species_id from species where name = '".$species_name."'";
 my $sth = $self->xref->dbc->prepare($sql);
 $sth->execute();
 my @row = $sth->fetchrow_array();
 my $species_id;
 if (@row) {
   $species_id = $row[0];
 } else {
   print STDERR "Couldn't get ID for species ".$species_name."\n";
   print STDERR "It must be one of :-\n";
   $sql = "select name from species";
   $sth = $self->dbc->prepare($sql);
   $sth->execute();
   while(my @row = $sth->fetchrow_array()){
     print STDERR $row[0]."\n";
   }
   die("Please try again :-)\n");
 }
 $sth->finish();
 
 return $species_id;
 

}

#
# Alt alleles
#

sub get_alt_alleles {
  my $self =  shift;

  my $gene_id;
  my $alt_id;
  my $sth = $self->core->dbc->prepare("select alt_allele_id, gene_id from alt_allele");
  $sth->execute;
  $sth->bind_columns(\$alt_id,\$gene_id);
  my $count = 0 ;
  my %alt_id_to_gene_id;
  my %gene_id_to_alt_id;
  my $max_alt_id = 0;
  while($sth->fetch){
    $count++;
    push @{$alt_id_to_gene_id{$alt_id}}, $gene_id;
    $gene_id_to_alt_id{$gene_id} = $alt_id;
      if($alt_id > $max_alt_id){
	$max_alt_id = $alt_id;
      }
  }

  my $insert_sth = $self->xref->dbc->prepare("insert into alt_allele (alt_allele_id, gene_id, is_reference) values (?, ?,?)");



  if($count){
    my %non_reference;

    my $sql = (<<SEQ);
SELECT g.gene_id
  FROM gene g, seq_region_attrib sra, attrib_type at
    WHERE g.seq_region_id = sra.seq_region_id AND
          at.attrib_type_id = sra.attrib_type_id AND
          at.code = 'non_ref'
SEQ

    $sth = $self->core->dbc->prepare($sql);
    $sth->execute;
    $sth->bind_columns(\$gene_id);
    while($sth->fetch()){
      $non_reference{$gene_id} = 1;
    }
    
    $sth = $self->xref->dbc->prepare("delete from alt_allele");
    $sth->execute;


    my $alt_added = 0;
    my $num_of_genes = 0;
    my $alt_failed = 0;
    foreach my $alt_id (keys %alt_id_to_gene_id){

      # make sure one and only one is on the reference
      my $ref_count = 0;
      foreach my $gene (@{$alt_id_to_gene_id{$alt_id}}){
	if(!defined($non_reference{$gene})){
	  $ref_count++;
	}
      }
      if($ref_count == 1){
	$alt_added++;
	foreach my $gene (@{$alt_id_to_gene_id{$alt_id}}){
	  $num_of_genes++;
	  my $ref =0 ;
	  if(!defined($non_reference{$gene})){
	    $ref = 1;
	  }
	  $insert_sth->execute($alt_id, $gene, $ref);
	}		
      }
      else{
	$alt_failed++;
      }
    }
    print "$alt_added alleles found containing $num_of_genes genes\n";
  }
  else{
    print "No alt_alleles found for this species.\n" ;
  }


  ### LRGs added as alt_alleles in the XREF system but never added to core.

  #
  # Use $max_alt_id for new ones.
  #
  
  my $sql =(<<LRG);
SELECT  ox.ensembl_id, gsi.gene_id
  FROM xref x, object_xref ox, external_db e, gene_stable_id gsi
    WHERE x.xref_id = ox.xref_id AND
          e.external_db_id = x.external_db_id AND
          e.db_name like "Ens_Hs_gene" AND
          ox.ensembl_object_type = "Gene" AND
           x.display_label = gsi.stable_id
LRG
  
  $sth = $self->core->dbc->prepare($sql);
  my ($core_gene_id, $lrg_gene_id);
  $sth->execute();
  $sth->bind_columns(\$lrg_gene_id, \$core_gene_id);
  
  $count =0;
  
  my $old_count = 0;
  my $new_count = 0;
  my $lrg_count = 0;
  #
  # If the core gene is already in an alt_allele set then use that alt_id for the LRG gene only.
  # Else use a new one and add both core and LRG.
  #
  

  while ($sth->fetch()){
    if(defined($gene_id_to_alt_id{$core_gene_id})){
      $insert_sth->execute($gene_id_to_alt_id{$core_gene_id}, $lrg_gene_id, 0);
      $old_count++;
    }
    elsif(defined($gene_id_to_alt_id{$lrg_gene_id})){
      $insert_sth->execute($gene_id_to_alt_id{$lrg_gene_id}, $core_gene_id, 1);
      print "LRG perculiarity\t$core_gene_id\t$lrg_gene_id\n";
      $lrg_count++;
    }
    else{ # new one.
      $max_alt_id++;
      $insert_sth->execute($max_alt_id, $lrg_gene_id, 0);
      $insert_sth->execute($max_alt_id, $core_gene_id, 1);
      $new_count++;
    }
    $count++;
  }

  if($count){
    print "Added $count alt_allels for the lrgs. $old_count added to previous alt_alleles and $new_count new ones\n";
    print "LRG problem count = $lrg_count\n";
  }



  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('alt_alleles_added',now())");
  $sth_stat->execute();
  $sth_stat->finish;
  return;
  
}


#
# official naming
#

sub get_official_name {
  return undef;
}



#
# NOW need to set display_xrefs and gene descriptions too to make it easier.
# But set them in the xxx_stable_id table for now otherwise it is too complicated
# to dump the xrefs etc
#

sub official_naming{
  my $self = shift;
  my $dbname = $self->get_official_name();
  
  if(!defined($dbname)){
    my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('official_naming_done',now())");
    $sth_stat->execute();
    $sth_stat->finish;    
    return;
  }
  if($dbname eq "MGI"){ # first Copy MGI to Genes
    $self->biomart_fix("MGI","Translation","Gene");
    $self->biomart_fix("MGI","Transcript","Gene");
  }   
  if($dbname eq "ZFIN_ID"){ # first Copy MGI to Genes
    $self->biomart_fix("ZFIN_ID","Translation","Gene");
    $self->biomart_fix("ZFIN_ID","Transcript","Gene");
  } 
  #  print "Official naming started. Copy $dbname from gene to canonical transcript\n" if($self->verbose);
  my ($max_object_xref_id, $max_xref_id,$max_object_xref_id2);


  $self->species_id($self->get_id_from_species_name($self->core->species));
			 
  my $sth = $self->xref->dbc->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;

  $sth = $self->xref->dbc->prepare("SELECT MAX(object_xref_id) FROM identity_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id2);
  $sth->fetch;

  
  
  $sth = $self->xref->dbc->prepare("SELECT MAX(xref_id) FROM xref");
  $sth->execute();
  $sth->bind_columns(\$max_xref_id);
  $sth->fetch;


  print "MAX xref_id = $max_xref_id MAX object_xref_id = $max_object_xref_id, max_object_xref from identity_xref = $max_object_xref_id2\n";

  #get the xref synonyms store as a hash of arrays.

#  my $syn_hash = $self->get_official_synonyms();
  my $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $self->xref->dbc->prepare($sql);    
   

  #
  # get the OFDN (HGNC/MGI) xrefs in the xrefs and add the synonyms, plus create hash to get id from the display name and desc
  # This is becouse not all MGI.HGNC's are loaded as these are priority xrefs.
  #

  my %display_label_to_id;
  my %display_label_to_desc;
  

  $sql = 'select x.accession, sy.synonym, x.description from synonym sy, xref x, source so where x.xref_id = sy.xref_id and so.source_id = x.source_id and so.name like "'.$dbname.'"';

  
  $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($display_label, $acc, $syn, $desc);
  $sth->bind_columns(\$acc,\$display_label, \$desc);
  while($sth->fetch){
    $display_label_to_id{$display_label} = $acc;
    $display_label_to_desc{$display_label} = $desc;
  }
  $sth->finish;

  # get label to id from xref database to start with.
  $sql = 'select x.accession, x.label, x.description from xref x, source s where s.source_id = x.source_id and s.name like "'.$dbname.'"';

  
  $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  $sth->bind_columns(\$acc,\$display_label, \$desc);
  while($sth->fetch){
    $display_label_to_id{$display_label} = $acc;
    if(!defined($desc)){
      print "undef desc for $display_label\n";
    }
    else{
      $display_label_to_desc{$display_label} = $desc;
    }
  }	
  $sth->finish;



  my %synonym;
  $sth = $self->xref->dbc->prepare('select es.synonym, x.label from synonym es, xref x, source s where x.xref_id = es.xref_id and x.source_id = s.source_id and
 s.name = "'.$dbname.'"' );
  $sth->execute();
  my ($name);
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;

#  $sth = $self->xref->dbc->prepare('select es.synonym, x.label from synonym es, xref x, source s where x.xref_id = es.xref_id and x.source_id = s.source_id and
# s.name = "EntrezGene"' );
#  $sth->execute();
#  $sth->bind_columns(\$syn,\$name);
#  while($sth->fetch){
#    $synonym{$syn} = $name;
#  }
#  $sth->finish;





#######################
#Do the naming bit now.
#######################

# get the vega external_sources


  my (              $clone_based_vega_gene_id, $clone_based_ensembl_gene_id);
  my ($odn_tran_id, $clone_based_vega_tran_id, $clone_based_ensembl_tran_id);
  my ($mirbase_gene_id, $rfam_gene_id);
  my ($mirbase_tran_id, $rfam_tran_id);

  $sth = $self->xref->dbc->prepare("select source_id from source where name like ?");

  $sth->execute("Clone_based_vega_gene");
  $sth->bind_columns(\$clone_based_vega_gene_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_gene");
  $sth->bind_columns(\$clone_based_ensembl_gene_id);
  $sth->fetch;
  
  $sth->execute("RFAM_gene_name");
  $sth->bind_columns(\$rfam_gene_id);
  $sth->fetch;

  $sth->execute("miRBase_gene_name");
  $sth->bind_columns(\$mirbase_gene_id);
  $sth->fetch;

  if(!defined($clone_based_vega_gene_id)){
    die "Could not find external database name Clone_based_vega_gene\n";
  }
  if(!defined($clone_based_ensembl_gene_id)){
    die "Could not find external database name Clone_based_ensembl_gene\n";
  }
  if(!defined($rfam_gene_id)){
    die "Could not find external database name RFAM_gene_name\n";
  }
  if(!defined($mirbase_gene_id)){
    die "Could not find external database name miRBase_gene_name\n";
  }


  $sth->execute($dbname."_transcript_name");
  $sth->bind_columns(\$odn_tran_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_vega_transcript");
  $sth->bind_columns(\$clone_based_vega_tran_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_transcript");
  $sth->bind_columns(\$clone_based_ensembl_tran_id);
  $sth->fetch;
  
  $sth->execute("RFAM_transcript_name");
  $sth->bind_columns(\$rfam_tran_id);
  $sth->fetch;

  $sth->execute("miRBase_transcript_name");
  $sth->bind_columns(\$mirbase_tran_id);
  $sth->fetch;

  if(!defined($odn_tran_id)){
    die "Could not find external database name ".$dbname."_transcript_name\n";
  }
  if(!defined($clone_based_vega_tran_id)){
    die "Could not find external database name Clone_based_vega_transcript\n";
  }
  if(!defined($clone_based_ensembl_tran_id)){
    die "Could not find external database name Clone_based_ensembl_transcript\n";
  }
  if(!defined($rfam_tran_id)){
    die "Could not find external database name RFAM_transcript_name\n";
  }
  if(!defined($mirbase_tran_id)){
    die "Could not find external database name miRBase_transcript_name\n";
  }


  #
  # Set up quick hashes for getting source id from database_name
  #
  my %dbname_tran_source;
  $dbname_tran_source{$dbname}   = $odn_tran_id;
  $dbname_tran_source{"miRBase"} = $mirbase_tran_id;
  $dbname_tran_source{"RFAM"}    = $rfam_tran_id;


  ###########################
  # Delete the old ones.
  ###########################

  # delete the synonyms first


### AHHH gene ones are not new!!!! or may not be!! need a way to differentiate them

  my $list =  "$odn_tran_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id, $clone_based_ensembl_tran_id, $rfam_tran_id, $rfam_gene_id, $mirbase_tran_id, $mirbase_gene_id";

  $sql =(<<DE1);
DELETE s
  FROM synonym s, xref x 
    WHERE s.xref_id = x.xref_id AND 
          x.source_id in ( $list );
DE1

  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();

 
  my $del_identity_sql =(<<DE2);
DELETE i 
  FROM object_xref o, xref x, identity_xref i
    WHERE i.object_xref_id = o.object_xref_id AND
           x.xref_id = o.xref_id AND
            x.source_id in ( $list )
DE2
  $sth = $self->xref->dbc->prepare($del_identity_sql);
  $sth->execute();
 
  my $del_ox_sql = (<<DE3);
DELETE o 
  FROM object_xref o, xref x 
    WHERE x.xref_id = o.xref_id AND
           x.source_id in ( $list )
DE3
  $sth = $self->xref->dbc->prepare($del_ox_sql);
  $sth->execute();
 
  my $del_x_sql = "delete x from xref x where x.source_id in ( $list )";

  $sth = $self->xref->dbc->prepare($del_x_sql);
  $sth->execute();


  my $del_synonym_sql = "delete s from xref x, synonym s where s.xref_id = x.xref_id and x.source_id = $clone_based_vega_tran_id and x.info_type = 'MISC'"; # original ones added have info type of "DIRECT"

  $sth = $self->xref->dbc->prepare($del_synonym_sql);
  $sth->execute();

  $del_x_sql = "delete x from xref x where x.source_id = $clone_based_vega_tran_id and x.info_type = 'MISC'"; # original ones added have info type of "DIRECT"
  $sth = $self->xref->dbc->prepare($del_x_sql);
  $sth->execute();


  $sth =  $self->xref->dbc->prepare("update transcript_stable_id set display_xref_id = null");
  $sth->execute;

  $sth =  $self->xref->dbc->prepare("update gene_stable_id set display_xref_id = null");
  $sth->execute;


  ######################################################
  # Get the current max values for xref and object_xref
  ######################################################


  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $self->core->dbc);
  my $ga = $db->get_GeneAdaptor();


  ###########################
  # Process each Gene
  ###########################

  

  my %gene_to_transcripts;
  my %gene_id_to_stable_id;
  my %tran_id_to_stable_id;


  $sql = "UPDATE gene_stable_id SET display_xref_id = null, desc_set =0";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute;

  $sql = "UPDATE transcript_stable_id SET display_xref_id = null";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute;

  

  $sql =(<<SQ0);
SELECT gtt.gene_id, gtt.transcript_id, gsi.stable_id, tsi.stable_id 
  FROM gene_transcript_translation gtt, gene_stable_id gsi, transcript_stable_id tsi 
    WHERE gtt.gene_id = gsi.internal_id AND
          gtt.transcript_id = tsi.internal_id
    ORDER BY gsi.stable_id, tsi.stable_id
SQ0

  $sth = $self->xref->dbc->prepare($sql);

  $sth->execute;
  my $gene_id;
  my $tran_id;
  my $gsi;
  my $tsi;
  $sth->bind_columns(\$gene_id, \$tran_id, \$gsi, \$tsi);
  my @sorted_gene_ids;
  while ($sth->fetch){
    if(!defined($gene_to_transcripts{$gene_id})){
      push @sorted_gene_ids, $gene_id;
    }
    push @{$gene_to_transcripts{$gene_id}}, $tran_id;
    $gene_id_to_stable_id{$gene_id} = $gsi; 
    $tran_id_to_stable_id{$tran_id} = $tsi; 
  }
  


  $sql =(<<SQ1);
SELECT x.label, x.xref_id, ox.object_xref_id, s.priority 
  FROM xref x, object_xref ox, source s
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND
          s.name = ? AND
          ox.ox_status = 'DUMP_OUT' AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = ?
SQ1
  my $dbentrie_sth = $self->xref->dbc->prepare($sql);


  $sql=(<<SQ2);
SELECT x.label, x.xref_id, ox.object_xref_id, s.priority 
  FROM xref x, object_xref ox, source s 
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND 
          s.name = ? AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = ?
SQ2
  my $lrg_find_sth = $self->xref->dbc->prepare($sql);

  my $lrg_set_status_sth = $self->xref->dbc->prepare("update object_xref set ox_status = 'NO_DISPLAY' where object_xref_id = ?");

  $sql=(<<SQ4);
SELECT x.xref_id, s.priority 
  FROM xref x,source s, object_xref ox
    WHERE x.xref_id = ox.xref_id AND
          x.source_id = s.source_id AND
          x.label = ? AND
          s.name = ? AND
          ox.ox_status = 'DUMP_OUT'
    ORDER BY s.priority
SQ4
  my $lrg_to_hgnc_sth  = $self->xref->dbc->prepare($sql);



  $sql = "insert into xref (xref_id, source_id, accession, label, version, species_id, info_type, info_text, description) values (?, ?, ?, ?,  0, ".$self->species_id.", 'MISC', ?, ? )";
  my $ins_xref_sth = $self->xref->dbc->prepare($sql); 

  $sql = "insert into xref (xref_id, source_id, accession, label, version, species_id, info_type, info_text, description) values (?, ?, ?, ?,  ?, ".$self->species_id.", 'MISC', ?, ? )";
  my $ins_xref_ver_sth = $self->xref->dbc->prepare($sql);

  my $ins_dep_ix_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)");

  my $delete_odn_sth = $self->xref->dbc->prepare('UPDATE object_xref SET ox_status = "MULTI_DELETE" where object_xref_id = ?');

  $sql=(<<SQ5);
SELECT x.xref_id 
  FROM xref x, source s, object_xref ox 
    WHERE ox.xref_id = x.xref_id AND
           x.source_id = s.source_id AND
           x.label = ? AND 
           s.name like '$dbname' AND
           s.priority_description like "vega" AND
           ox.ox_status ="DUMP_OUT"
SQ5


  my $find_hgnc_sth = $self->xref->dbc->prepare($sql);

  my $set_gene_display_xref_sth =  $self->xref->dbc->prepare('UPDATE gene_stable_id SET display_xref_id =? where internal_id = ?');

  my $set_tran_display_xref_sth =  $self->xref->dbc->prepare('UPDATE transcript_stable_id SET display_xref_id =? where internal_id = ?');

  # important or will crash and burn!!!
  my %xref_added; # store those added already $xref_added{$accession:$source_id} = $xref_id;

  my $get_xref_info_sth =  $self->xref->dbc->prepare("select x.label, x.accession, s.priority_description, x.description  from xref x, source s where xref_id = ? and s.source_id = x.source_id");

  #
  # Okay we assign unused_priority to be the number of time a vega transcript is attached to make sure
  # we get the BEST name for the gene (i.e. the one that appears the most)
  #
  my $ins_object_xref_sth =  $self->xref->dbc->prepare("insert into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_type, ox_status, unused_priority) values (?, ?, ?, ?, 'MISC', 'DUMP_OUT', ?)");
  foreach my $gene_id (@sorted_gene_ids){
    
    my @ODN=();
    my $xref_id;
    my $object_xref_id;
    my $display;
    my $level;
    my $tran_source = $dbname;

    my $best_info = undef;

    $dbentrie_sth->execute($dbname, $gene_id, "Gene");
    $dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
    my $best_level=999;

    my $count = 0;
    my @list=();
    my @list_ox=();
    my %xref_id_to_display;
    my %best_list;
    while($dbentrie_sth->fetch){

      push @list, $xref_id;
      push @list_ox, $object_xref_id;
      $count++;
      $xref_id_to_display{$xref_id} = $display;
      if($level < $best_level){
	@ODN = ();
	push @ODN, $xref_id;
	$best_level = $level;
      }
      elsif($level == $best_level){
	push @ODN, $xref_id;
      }
    }

    foreach my $x (@ODN){
      $best_list{  $xref_id_to_display{$x} } = 1;
    }

    my $gene_symbol = undef;
    my $gene_symbol_xref_id;
    my $other_symbol = undef;
    my $other_xref_id = undef;
    my $other_source = undef;
    my $clone_name = undef;
    my $vega_clone_name = undef;

    if($count > 1){
      print "For gene ".$gene_id_to_stable_id{$gene_id}." we have mutiple ".$dbname."'s\n";
      if(scalar(@ODN) == 1){ # found one that is "best"
#	foreach my $x (@list){
	my $i=0;
	while ($i < scalar(@list)){
	  my $x = $list[$i];
	  if($x != $ODN[0]){
	    print "\tremoving ".$xref_id_to_display{$x}." from gene\n";
	    #remove object xref....
	    $delete_odn_sth->execute($list_ox[$i])|| die "Could not set staus to MULTI_DELETE for object_xref ".$list_ox[$i]."\n";
	  }
	  else{
	    print "\tKeeping the best one ".$xref_id_to_display{$x}."\n";
	    $gene_symbol = $xref_id_to_display{$x};
	    $gene_symbol_xref_id = $x;
	  }
	  $i++;
	}
      }
    }

    my %name_count;
    my %tran_to_vega_ext;
    foreach my $tran_id ( @{$gene_to_transcripts{$gene_id}} ){
      $dbentrie_sth->execute($dbname."_curated_transcript_notransfer", $tran_id, "Transcript");
      $dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
      while($dbentrie_sth->fetch){
	my $symbol_bit;
	my $num;
	if($display =~ /(.+)-(\d+)$/){
	  $symbol_bit = $1;
	  $num = $2;
	}
	else{
	  print STDERR "$display does not have the usual expected regex\n";
	  next;
	}	
	if(!defined($num) or !$num or $num eq ""){
	  print "Problem finding number for $display\n";
	}

	$tran_to_vega_ext{$tran_id} = $num;

	if(defined($display_label_to_desc{$symbol_bit})){
	}
	elsif(defined($synonym{$symbol_bit})){
	  $symbol_bit = $synonym{$symbol_bit};
	}
	else{
	  # -ps added as MGI have added -ps to the pseudo genes but vega has not caught
	  # up with this yet so check for this.
	
	  if(!defined($display_label_to_desc{$symbol_bit."-ps"})){
	    print STDERR "Warning Could not find id for $symbol_bit came from $display for $dbname\n";
	    next;
	  }
	}
	if($best_list{$symbol_bit}){
	  $name_count{$symbol_bit}++;
	}
      }
    }

    if(scalar(@ODN) == 1){  # one hgnc to this gene - perfect case :-)
      # $ODN[0] is an xref_id i need the display name
      $gene_symbol = $xref_id_to_display{$ODN[0]};
      $gene_symbol_xref_id = $ODN[0];
    }
    elsif(scalar(@ODN) > 1){ # try to use vega to find the most common one
      print "Multiple best ".$dbname."'s using vega to find the most common for ".$gene_id_to_stable_id{$gene_id}."\n";


      #############################################################################################
      # need some way to get the xref_id for the hgnc from the curated??
      if(scalar(%name_count)){

	my $top =0;
	foreach my $vn ( keys %name_count){
	  if($name_count{$vn} > $top){
	    $top = $name_count{$vn};
	    $gene_symbol = $vn;
	    foreach my $y (@ODN){
	      if($vn eq $xref_id_to_display{$y}){
		$gene_symbol_xref_id = $y;
	      }
	    }
	  }
	}
	print "\t$gene_symbol chosen from vega\n";
      }
      else{  # take the first one ??
	my $i = 0;
	foreach my $x (@ODN){
	  print "\t".$xref_id_to_display{$x};
	  if(!$i){
	    print "  (chosen as first)\n";
	    $gene_symbol =  $xref_id_to_display{$x};
	    $gene_symbol_xref_id =  $x;
	  }
	  else{
	    print "  (left as $dbname reference but not gene symbol)\n";
	  }
	  $i++;
	}
      }
    }
    else{# look for LRG
      # look for LRG_HGNC_notransfer, if found then find HGNC equiv and set to this
      print "LRG FOUND with no HGNC, should have gotten this via the alt allele table?? gene_id = $gene_id\n";
      $lrg_find_sth->execute("LRG_HGNC_notransfer", $gene_id, "Gene");
      $lrg_find_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
      while($lrg_find_sth->fetch){
	$lrg_set_status_sth->execute($object_xref_id); # set oc_status to no _display as we do not want this transferred, just the equivalent hgnc
	my $new_xref_id  = undef;
	my $pp;
	$lrg_to_hgnc_sth->execute($display,"HGNC");
	$lrg_to_hgnc_sth->bind_columns(\$new_xref_id,\$pp);
	$lrg_to_hgnc_sth->fetch;
	if(defined($new_xref_id)){
	  $gene_symbol = $display;
	  $gene_symbol_xref_id = $new_xref_id;
	}
      }
    }
    if(!defined($gene_symbol)){ # try ther database source (should be RFAm and mirbase only)
      #set $other_symbol if RfAM or miRBase found
      $other_symbol = undef;
      $other_xref_id = undef;
      foreach my $ext_db_name (qw(miRBase RFAM)){
	$dbentrie_sth->execute($ext_db_name, $gene_id, "Gene");
	$dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
	while($dbentrie_sth->fetch){
	  $gene_symbol = $display;
	  $gene_symbol_xref_id = $xref_id;
	  $tran_source = $ext_db_name;
	  next;
	}
      }	
    }

    foreach my $tran_id  (@{$gene_to_transcripts{$gene_id}}){
      $dbentrie_sth->execute("Clone_based_vega_transcript", $tran_id, "Transcript");
      $dbentrie_sth->bind_columns(\$display, \$xref_id, \$object_xref_id, \$level);
      while($dbentrie_sth->fetch){
	$display =~ /(.+)-(\d\d\d)$/;
	my $acc_bit =$1;
	my $num = $2;
	$tran_to_vega_ext{$tran_id} = $num;
	$xref_added{$display.":".$clone_based_vega_tran_id} = $xref_id;
	$vega_clone_name = $acc_bit;
      }
    }
    if(!defined($gene_symbol) and !defined($other_symbol) ){   # No HGNC or other so look for vega clone names

      if(!defined($vega_clone_name)){ #if no vega clone name use the ensembl clone name


        my $gene= $ga->fetch_by_dbID($gene_id);
	my $slice = $gene->slice->sub_Slice($gene->start,$gene->end,$gene->strand);
	my $len = 0;
	if($dbname ne "ZFIN_ID"){
	  my $clone_projection = $slice->project('clone'); 
	  foreach my $seg (@$clone_projection) {
	    my $clone = $seg->to_Slice();
	    if($clone->length() > $len){
	      $clone_name = $clone->seq_region_name;
	      $len = $clone->length;
	    }
	   }
	}
	if(!defined($clone_name)){
	  # try going via contig
	  my $super_projection = $slice->project('contig');
	  my $super;
	  $len = 0;
	  foreach my $seg (@$super_projection) {
	    $super = $seg->to_Slice();
	    if($super->length() > $len){
	      $clone_name = $super->seq_region_name;
	      $len = $super->length;
	    }
	  }
	  $len = 0;
	  if($dbname ne "ZFIN_ID"){
	    my $clone_projection = $super->project('clone');
	    foreach my $seg (@$clone_projection) {
	      my $clone = $seg->to_Slice();
	      if($clone->length() > $len){
		$clone_name = $clone->seq_region_name;
		$len = $clone->length;
	      }
	    }
	  }	
	  if(!defined($clone_name)){
	    print STDERR "PROJECT failed for ".$gene->stable_id."\n";
	    next;
	  }

	}

	if(defined($clone_name)){
	  $clone_name =~ s/[.]\d+//;    #remove .number	
	}
      } 
    }


    #
    # Set the names now that we know which to use.
    #
    if( !(defined($clone_name) or defined($vega_clone_name)) and !defined($gene_symbol) and !defined($other_symbol)){
      print STDERR "Problem gene ".$gene_id_to_stable_id{$gene_id}." could not get a clone name or ".$dbname." symbol\n";
      next;
    }
    if(defined($gene_symbol)){
      #gene symbol already set as it is HGNC or MGI so do not need add a new xref anything for the gene;

      my $desc = $display_label_to_desc{$name};


      if(!defined($gene_symbol_xref_id)){
	$find_hgnc_sth->execute($gene_symbol);
	$find_hgnc_sth->bind_columns(\$gene_symbol_xref_id);
	$find_hgnc_sth->fetch();
	if(!defined($gene_symbol_xref_id)){  # remember mouse has -ps newly added and vega may not be uptodate
	  $find_hgnc_sth->execute($gene_symbol."-ps");
	  $find_hgnc_sth->bind_columns(\$gene_symbol_xref_id);
	  $find_hgnc_sth->fetch();
	}	
	if(!defined($gene_symbol_xref_id)){
	  print "BLOOMING NORA could not find $gene_symbol in $dbname\n";
	  next;
	}
      }
      $set_gene_display_xref_sth->execute($gene_symbol_xref_id,$gene_id);


      my $no_vega_ext = 201;
      foreach my $tran_id ( @{$gene_to_transcripts{$gene_id}} ){
	my $ext;
	my $source_id = $dbname_tran_source{$tran_source};
	if(defined($tran_to_vega_ext{$tran_id})){
	  $ext = $tran_to_vega_ext{$tran_id};
	}
	else{
	  $ext = $no_vega_ext;
	  $no_vega_ext++;
	}
	my $id = $gene_symbol."-".$ext;
	if(!defined($xref_added{$id.":".$source_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $source_id, $id, $id, "", $desc);
	  $xref_added{$id.":".$source_id} = $max_xref_id;
	}
	$set_tran_display_xref_sth->execute($xref_added{$id.":".$source_id}, $tran_id);
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $tran_id, 'Transcript', $xref_added{$id.":".$source_id},undef);
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
      }
    }
    else{ # use clone name
      my $t_source_id;
      my $g_source_id;
      my $desc;
      if(defined($vega_clone_name)){
	$name = $vega_clone_name;
	$t_source_id = $clone_based_vega_tran_id;
	$g_source_id = $clone_based_vega_gene_id;
	$desc = "via havana clone name";
      }
      else{
	if(defined($clone_name)){
	  $name = $clone_name;
	  $t_source_id = $clone_based_ensembl_tran_id;
	  $g_source_id = $clone_based_ensembl_gene_id;
	  $desc = "via ensembl clone name";
	}
	elsif(defined($other_symbol)){
	  die "SHOULD NOT GET HERE NOW??\n";
	  $name = $other_symbol;
	  if($other_source =~ /RFAM/){
	    $t_source_id = $rfam_tran_id;
	    $g_source_id = $rfam_gene_id;
	  }
	  else{
	    $t_source_id = $mirbase_tran_id;
	    $g_source_id = $mirbase_gene_id;
	  }
	}
	else{
	  die "No name";
	}
	my $num = 1;
	my $unique_name = $name.".".$num;
	while(defined($xref_added{$unique_name.":".$g_source_id})){
	  $num++;
	  $unique_name = $name.".".$num;
	}
	$name = $unique_name;
      }
      
      # first add the gene xref and set display_xref_id
      # store the data
      if(!defined($xref_added{$name.":".$g_source_id})){
	$max_xref_id++;
	$ins_xref_sth->execute($max_xref_id, $g_source_id, $name, $name, $desc, undef);
	$xref_added{$name.":".$g_source_id} = $max_xref_id;
      }	
      $set_gene_display_xref_sth->execute($xref_added{$name.":".$g_source_id},$gene_id);
	
      $max_object_xref_id++;
      $ins_object_xref_sth->execute($max_object_xref_id, $gene_id, 'Gene', $xref_added{$name.":".$g_source_id}, undef);	
      $ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
      
      # then add transcript names.
      my $ext = 201;
      foreach my $tran_id ( sort @{$gene_to_transcripts{$gene_id}}){
	my $id =  $name."-".$ext;
	if(defined($tran_to_vega_ext{$tran_id})){
	  $id = $name."-".$tran_to_vega_ext{$tran_id};
	  if(!defined($xref_added{$id.":".$t_source_id})){
	    $max_xref_id++;
	    $ins_xref_sth->execute($max_xref_id, $t_source_id, $id, $id, $desc, undef);
	    $xref_added{$id.":".$t_source_id} = $max_xref_id;

	    $max_object_xref_id++;
	    $ins_object_xref_sth->execute($max_object_xref_id, $tran_id, 'Transcript', $xref_added{$id.":".$t_source_id}, undef);	
	    $ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
	  }
	  $set_tran_display_xref_sth->execute($xref_added{$id.":".$t_source_id}, $tran_id);
	}
	else{
	  while(defined($xref_added{$id.":".$t_source_id})){
	    $ext++;
	    $id = $name."-".$ext;
	  }	
	  $ext++;
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $t_source_id, $id, $id, $desc, undef);
	  $xref_added{$id.":".$t_source_id} = $max_xref_id;

	  $max_object_xref_id++;
#	  print "$id\t$t_source_id\t".$xref_added{$id.":".$t_source_id}."\n";
	  $ins_object_xref_sth->execute($max_object_xref_id, $tran_id, 'Transcript', $xref_added{$id.":".$t_source_id}, undef);	
	  $ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
	  $set_tran_display_xref_sth->execute($max_xref_id, $tran_id);
	}

#	print $id."\ts=".$t_source_id."\n";

	
      }
      
    }
  } # for each gene



  ########################################################
  # Copy $dbname from gene to the canonical transcripts. #
  ########################################################

  # remove the ignore later on after testing
  my $sth_add_ox = $self->xref->dbc->prepare("insert ignore into object_xref (object_xref_id, xref_id, ensembl_id, ensembl_object_type, linkage_type, ox_status, master_xref_id) values(?, ?, ?, 'Transcript', ?, ?, ?)");
  
  #  my $object_xref_id = $max_object_xref_id + 1;
  
  $max_object_xref_id++;
  
  if($max_object_xref_id == 1){
    die "max_object_xref_id should not be 1\n";
  }
  
  my $object_sql = (<<FSQL);
select x.xref_id, o.ensembl_id, o.linkage_type, o.ox_status, o.master_xref_id, ix.query_identity, ix.target_identity
  from xref x, source s, object_xref o, identity_xref ix
    where x.source_id = s.source_id and 
      ix.object_xref_id  = o.object_xref_id and
      o.ox_status = "DUMP_OUT" and
      s.name like "$dbname" and 
      o.xref_id = x.xref_id  and
      o.ensembl_object_type = "Gene";
FSQL
  
  $sql = "select gene_id, canonical_transcript_id from gene";
  $sth = $self->core->dbc->prepare($sql);
  
  $sth->execute();

  $sth->bind_columns(\$gene_id,\$tran_id);
  my %gene_to_tran_canonical;
  while ($sth->fetch){
    $gene_to_tran_canonical{$gene_id} = $tran_id;
  }
  $sth->finish;

  $sth = $self->xref->dbc->prepare($object_sql);

  $sth->execute();
  my ($xref_id, $linkage_type, $ox_status, $q_id, $t_id, $master_id);
  $sth->bind_columns(\$xref_id, \$gene_id, \$linkage_type, \$ox_status, \$q_id, \$t_id, \$master_id);


  while ($sth->fetch){
    if(defined($gene_to_tran_canonical{$gene_id})){
      $max_object_xref_id++;
      $sth_add_ox->execute($max_object_xref_id, $xref_id, $gene_to_tran_canonical{$gene_id}, $linkage_type, $ox_status, $master_id) || print STDERR "(Gene id - $gene_id) Could not add  $max_object_xref_id, .".$gene_to_tran_canonical{$gene_id}.", $xref_id, $linkage_type, $ox_status to object_xref, master_xref_id to $master_id\n";
      $ins_dep_ix_sth->execute($max_object_xref_id, $q_id, $t_id);
    }
    else{
      print STDERR "Could not find canonical for gene $gene_id\n";
    }
  }
  $sth->finish;

  print "Copied all $dbname from gene to canonical transcripts\n" if($self->verbose);
  
  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('official_naming_done',now())");
  $sth_stat->execute();
  $sth_stat->finish;
  return;
}



sub biomart_fix{
  my ($self, $db_name, $type1, $type2) = @_;
  my $xref_dbc = $self->xref->dbc;

  print "$db_name is associated with both $type1 and $type2 object types\n" if($self->verbose);
  
  my $to;
  my $from;
  my $to_id;
  my $from_id;
  if($type1 eq "Gene" or $type2 eq "Gene"){
    $to = "Gene";
    $to_id = "gene_id";
    if($type1 eq "Translation" or $type2 eq "Translation"){
      $from = "Translation";
      $from_id = "translation_id"
    }
    else{
      $from = "Transcript";
      $from_id = "transcript_id";
    }
  }
  else{
    $to = "Transcript";
    $to_id = "transcript_id";
    $from = "Translation";
    $from_id = "translation_id";
  }
  
  print "Therefore moving all associations from $from to ".$to."\n" if($self->verbose);
  

  my $sql =(<<EOF);
  UPDATE IGNORE object_xref, gene_transcript_translation, xref, source
    SET object_xref.ensembl_object_type = "$to",
      object_xref.ensembl_id = gene_transcript_translation.$to_id 
	WHERE object_xref.ensembl_object_type = "$from" AND
	  object_xref.ensembl_id = gene_transcript_translation.$from_id AND
	    xref.xref_id = object_xref.xref_id AND
	      xref.source_id = source.source_id AND
                object_xref.ox_status = "DUMP_OUT"  AND
		  source.name = "$db_name";
EOF
  my $result =  $xref_dbc->do($sql) ;
#  print "\n$sql\n";

  if($db_name eq "GO"){
    $sql =(<<EOF2);
  DELETE object_xref, identity_xref, go_xref
    FROM object_xref, xref, source, identity_xref, go_xref
      WHERE object_xref.ensembl_object_type = "$from" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
	xref.xref_id = object_xref.xref_id AND
          go_xref.object_xref_id = object_xref.object_xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF2
    
  $result = $xref_dbc->do($sql);  
  }
  else{
    $sql =(<<EOF3);
  DELETE object_xref, identity_xref
    FROM object_xref, xref, source, identity_xref
      WHERE object_xref.ensembl_object_type = "$from" AND
        identity_xref.object_xref_id = object_xref.object_xref_id AND
	xref.xref_id = object_xref.xref_id AND
	  xref.source_id = source.source_id AND
            object_xref.ox_status = "DUMP_OUT"  AND
	      source.name = "$db_name";
EOF3
    
  $result = $xref_dbc->do($sql);
  }
#  print "\n$sql\n";

  #delete dependent_xref 
  $sql =(<<EOF4);
  DELETE FROM dependent_xref WHERE object_xref_id NOT IN 
   (SELECT object_xref_id FROM object_xref);
EOF4
      
}

sub biomart_testing{
  my ($self) = @_;

  my $sql = 'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';


  my $again = 1;
  while ($again){
    $again = 0;
    
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    my ($type, $count, $name);
    my ($last_type, $last_count, $last_name);
    $sth->bind_columns(\$type,\$count,\$name);
    $last_name = "DEFAULT";
    while (!$again and $sth->fetch){
      if($last_name eq $name){
	$again  = 1;
	$self->biomart_fix($name,$last_type, $type);
#	$again = 0; # remove this line after testing
      }
      $last_name = $name;
      $last_type= $type;
      $last_count = $count;
    }
    $sth->finish;  
  }

  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('biomart_test_finished',now())");
  $sth_stat->execute();
  $sth_stat->finish;
}
  


sub biomart_test{
  my ($self) = @_;

  my $sql = 'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';


  my $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($type, $count, $name);
  my ($last_type, $last_count, $last_name);
  $sth->bind_columns(\$type,\$count,\$name);
  $last_name = "NOTKNOWN";
  my $first = 1;
  while ($sth->fetch){
    if($last_name eq $name){
      if($first){
	print STDERR "\nProblem Biomart test fails\n";
	$first=0;
      }
      print STDERR "$last_name\t$last_count\t$last_type\n";
      print STDERR "$name\t$count\t$type\n";
    }
    $last_name = $name;
    $last_type= $type;
    $last_count = $count;
  }
  $sth->finish;

}

# remove a list of patterns from a string
sub filter_by_regexp {

  my ($self, $str, $regexps) = @_;

  foreach my $regexp (@$regexps) {
    $str =~ s/$regexp//ig;
  }

  return $str;

}



sub get_official_synonyms{
  my $self = shift;
  my %hgnc_syns;
  my %seen;          # can be in more than one for each type of hgnc.

  my $dbname = $self->get_official_name();

  my $sql = (<<SYN);
SELECT  x.accession, x.label, sy.synonym 
  FROM xref x, source so, synonym sy
    WHERE x.xref_id = sy.xref_id
      AND so.source_id = x.source_id
      AND so.name like "$dbname"
SYN

  my $sth = $self->xref->dbc->prepare($sql);    

  $sth->execute;
  my ($acc, $label, $syn);
  $sth->bind_columns(\$acc, \$label, \$syn);

  my $count = 0;
  while($sth->fetch){
    if(!defined($seen{$acc.":".$syn})){
      push @{$hgnc_syns{$acc}}, $syn;
      push @{$hgnc_syns{$label}}, $syn;
      $count++;
    }
    $seen{$acc.":".$syn} = 1;
  }
  $sth->finish;

  return \%hgnc_syns;

}

sub get_species_id_from_species_name{
  my ($self,$species) = @_;


  my $sql = "select species_id from species where name = '".$species."'";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute();
  my @row = $sth->fetchrow_array();
  my $species_id;
  if (@row) {
    $species_id = $row[0];
  } else {
    print STDERR "Couldn't get ID for species ".$species."\n";
    print STDERR "It must be one of :-\n";
    $sql = "select name from species";
    $sth = $self->dbc->prepare($sql);
    $sth->execute();
    while(my @row = $sth->fetchrow_array()){
      print STDERR $row[0]."\n";
    }
    die("Please try again :-)\n");
  }
  $sth->finish();

  return $species_id;
}


sub clean_up{
  my $self = shift;
  my $stats = shift;
  

  # remove all object_xref, identity_xref  entries

  my $sql = "DELETE from object_xref";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from go_xref";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from identity_xref";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 
 
  # remove all xrefs after PARSED_xref_id
  # set dumped to NULL fro all xrefs.

  my $max_xref_id = $self->get_meta_value("PARSED_xref_id");

  if($max_xref_id){
    $sql = "DELETE from xref where xref_id > $max_xref_id";
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute(); 
  }

  $sql = "UPDATE xref set dumped = null";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 


  
  # remove all from core_info tables
  #        gene_transcript_translation
  #        [gene/transcript/translation]_stable_id
  #
  $sql = "DELETE from gene_transcript_translation";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from gene_stable_id";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 
 
  $sql = "DELETE from transcript_stable_id";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 
 
  $sql = "DELETE from translation_stable_id";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 
 

}

sub remove_mapping_data{
  my $self = shift;

  my $sql = "DELETE from mapping_jobs";
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  $sql = "DELETE from mapping";
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

}


sub revert_to_parsing_finished{
  my $self = shift;


  $self->clean_up();
  $self->remove_mapping_data();
  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('parsing_finished',now())");
  $sth_stat->execute();
  $sth_stat->finish;    
}


sub revert_to_mapping_finished{
  my $self = shift;

  $self->clean_up();

  # set mapping jobs to SUBMITTED
  my $sql = 'UPDATE mapping_jobs set status = "SUBMITTED"';;
  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute(); 

  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('mapping_finished',now())");
  $sth_stat->execute();
  $sth_stat->finish;    

}

#
# In case we have alt alleles with xefs, these will be direct ones
# we need to move all xrefs on to the reference
#


sub get_alt_allele_hashes{
  my $self= shift;

  my %alt_to_ref;
  my %ref_to_alts;

  my $sql = "select gene_id, is_reference from alt_allele order by alt_allele_id, is_reference DESC";

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
  my ($gene_id, $is_ref);
  $sth->bind_columns(\$gene_id, \$is_ref);
  my $ref_gene;
  while($sth->fetch()){
    if($is_ref){
      $ref_gene = $gene_id;
    }
    else{
      $alt_to_ref{$gene_id} = $ref_gene;
      push @{$ref_to_alts{$ref_gene}}, $gene_id;
    }
  }
  $sth->finish;

  return \%alt_to_ref, \%ref_to_alts;
}


#sub move_xrefs_from_alt_allele_to_reference {
sub process_alt_alleles{
  my $self = shift;

  # ALL are on the Gene level now. This may change but for now it is okay.
  my ($alt_to_ref, $ref_to_alts) = $self->get_alt_allele_hashes();

  my $tester = XrefMapper::TestMappings->new($self);
  if($tester->unlinked_entries){
    die "Problems found before process_alt_alleles\n";
  }
  #
  # Move the xrefs on to the reference Gene.
  # NOTE: Igonore used as the xref might already be on this Gene already and we do not want it to crash
  #
  my $move_sql =(<<MOVE);
UPDATE IGNORE object_xref ox, xref x, source s 
  SET ox.ensembl_id = ? 
    WHERE x.source_id = s.source_id AND 
          ox.xref_id = x.xref_id AND
          ox.ensembl_id = ? AND
          ox.ensembl_object_type = 'Gene' AND
          ox.ox_status = 'DUMP_OUT' AND 
          s.name in (
MOVE
$move_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

print "MOVE SQL\n$move_sql\n";

  #
  # Now where it was already on the Gene the ignore will have stopped the move
  # so we now want to just remove those ones as they already exist.
  #
  my $del_ix_sql =(<<DIX);
DELETE ix 
  FROM identity_xref ix, object_xref ox, xref x, source s 
    WHERE x.source_id = s.source_id AND
          ox.object_xref_id = ix.object_xref_id AND
          ox.xref_id = x.xref_id AND 
          ox.ensembl_id = ? AND 
          ox.ensembl_object_type = 'Gene' AND 
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DIX
$del_ix_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $del_sql =(<<DEL);
DELETE ox 
  FROM object_xref ox, xref x, source s 
    WHERE x.source_id = s.source_id AND
          ox.xref_id = x.xref_id AND 
          ox.ensembl_id = ? AND 
          ox.ensembl_object_type = 'Gene' AND 
          ox.ox_status = 'DUMP_OUT' AND
           s.name in (
DEL
$del_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $move_sth = $self->xref->dbc->prepare($move_sql)  || die "$move_sql cannot be prepared";
  my $del_ix_sth = $self->xref->dbc->prepare($del_ix_sql)    || die "$del_ix_sql cannot be prepared";
  my $del_sth = $self->xref->dbc->prepare($del_sql)    || die "$del_sql cannot be prepared";

  my $move_count = 0;
  my $del_ix_count = 0;
  my $del_ox_count = 0;
  foreach my $key (keys %$alt_to_ref){
    $move_sth->execute($alt_to_ref->{$key}, $key);
    $move_count += $move_sth->rows;

    $del_ix_sth->execute($key);
    $del_ix_count += $del_ix_sth->rows;

    $del_sth->execute($key);
    $del_ox_count += $del_sth->rows;
  }
  $move_sth->finish;
  $del_sth->finish;
  $del_ix_sth->finish;

  print "Number of rows:- moved = $move_count, identitys deleted = $del_ix_count, object_xrefs deleted = $del_ox_count\n";
  if($tester->unlinked_entries){
    die "Problems found mid process_alt_alleles\n";
  }
  #
  # Now we have all the data on the reference Gene we want to copy all the data
  # onto the alt alleles.
  #


  my $get_data_sql=(<<GET);
SELECT ox.object_xref_id, ox.ensembl_object_type, ox.xref_id, ox.linkage_annotation, 
       ox.linkage_type, ox.ox_status, ox.unused_priority, ox.master_xref_id,
       ix.query_identity, ix.target_identity, ix.hit_start, ix.hit_end,
       ix.translation_start, ix.translation_end, ix.cigar_line, ix.score, ix.evalue
  FROM xref x, source s, object_xref ox 
    LEFT JOIN identity_xref ix ON ox.object_xref_id =ix.object_xref_id 
      WHERE  x.source_id = s.source_id AND
             ox.xref_id = x.xref_id AND
             ox.ensembl_id = ? AND
             ox.ox_status = 'DUMP_OUT' AND
             ox.ensembl_object_type = 'Gene' AND
              s.name in (
GET

  $get_data_sql .= "'".join("', '",$self->get_gene_specific_list()) . "')";

  my $get_data_sth = $self->xref->dbc->prepare($get_data_sql) || die "Could not prepare $get_data_sql";



  my $insert_object_xref_sql =(<<INO);
INSERT INTO object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_annotation, 
            linkage_type, ox_status, unused_priority, master_xref_id) 
       VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
INO

  my $insert_ox_sth = $self->xref->dbc->prepare($insert_object_xref_sql) || die "Could not prepare $insert_object_xref_sql";


  my $insert_identity_xref_sql = (<<INI);
INSERT INTO identity_xref (object_xref_id, query_identity, target_identity, hit_start, hit_end,
            translation_start, translation_end, cigar_line, score, evalue ) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
INI

  my $insert_ix_sth = $self->xref->dbc->prepare($insert_identity_xref_sql) || die "Could not prepare $insert_identity_xref_sql";



  my $max_object_xref_id;

  my $sth = $self->xref->dbc->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;
  if(!defined($max_object_xref_id) or !$max_object_xref_id){
    die "Problem getting max object_xref_id";
  }
  $max_object_xref_id++;

  my $added_count = 0;
  my $ignored = 0;
  foreach my $key (keys %$ref_to_alts){
    $get_data_sth->execute($key);
    my ($object_xref_id, $ensembl_object_type, $xref_id, $linkage_annotation,
	$linkage_type, $ox_status, $unused_priority, $master_xref_id,
	$query_identity, $target_identity, $hit_start, $hit_end,
	$translation_start, $translation_end, $cigar_line, $score, $evalue);

    $get_data_sth->bind_columns(\$object_xref_id, \$ensembl_object_type, \$xref_id, \$linkage_annotation,
				\$linkage_type, \$ox_status, \$unused_priority, \$master_xref_id,
				\$query_identity, \$target_identity, \$hit_start, \$hit_end,
				\$translation_start, \$translation_end, \$cigar_line, \$score, \$evalue);

    while( $get_data_sth->fetch()){
      foreach my $alt (@{$ref_to_alts->{$key}}){
	$max_object_xref_id++;
        $insert_ox_sth->execute($max_object_xref_id, $alt, $ensembl_object_type, $xref_id, $linkage_annotation,
				$linkage_type, $ox_status, $unused_priority, $master_xref_id) || die "Could not insert object_xref data";

#ONLY add identity xref if object_xref was added successfully.
	if( $insert_ox_sth->rows){
	  $added_count++;
	  $insert_ix_sth->execute($max_object_xref_id, $query_identity, $target_identity, $hit_start, $hit_end,
				$translation_start, $translation_end, $cigar_line, $score, $evalue) ||  die "Could not insert identity_xref data";
	}
	else{
	  $ignored++;
	}
      }
    }	
  }
  print "Added $added_count new mapping but ignored $ignored\n";
  
  if($tester->unlinked_entries){
    die "Problems found after process_alt_alleles\n";
  }
  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('alt_alleles_processed',now())");
  $sth_stat->execute();
  $sth_stat->finish;

}



sub get_gene_specific_list {
  my $self = shift;

  my @list = qw(DBASS3 DBASS5 EntrezGene miRBase RFAM UniGene Uniprot_genename WikiGene MIM_GENE MIM_MORBID HGNC);

  return @list;
}

sub source_defined_move{
  my $self = shift;

  my $tester = XrefMapper::TestMappings->new($self);
  if($tester->unlinked_entries){
    die "Problems found before source_defined_move\n";
  }
  foreach my $source ($self->get_gene_specific_list()){
    $self->biomart_fix($source,"Translation","Gene");
    $self->biomart_fix($source,"Transcript","Gene");
  }
  if($tester->unlinked_entries){
    die "Problems found after source_defined_move\n";
  }
  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('source_level_move_finished',now())");
  $sth_stat->execute();
  $sth_stat->finish;
}



#sub process_alt_alleles {
#  my $self = shift;
#
#
#
#
#
#  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('alt_alleles_processed',now())");
#  $sth_stat->execute();
#  $sth_stat->finish;
#}


1;
