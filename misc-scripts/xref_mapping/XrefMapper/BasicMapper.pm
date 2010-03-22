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
      warn("Did not find a specific mapping module XrefMapper::$value - using XrefMapper::BasicMapper instead\n") if(defined($verbose) and $verbose);
      require XrefMapper::BasicMapper;
	$module = "BasicMapper";
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
# official naming
#

sub get_official_name {
  return undef;
}


sub official_naming {
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
  print "Official naming started. Copy $dbname from gene to canonical transcript\n" if($self->verbose);
  my ($max_object_xref_id, $max_xref_id);


  $self->species_id($self->get_id_from_species_name($self->core->species));
			 
  my $sth = $self->xref->dbc->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;

  
  
  $sth = $self->xref->dbc->prepare("SELECT MAX(xref_id) FROM xref");
  $sth->execute();
  $sth->bind_columns(\$max_xref_id);
  $sth->fetch;



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

  my $sql = "select gene_id, canonical_transcript_id from gene";
  $sth = $self->core->dbc->prepare($sql);
  
  $sth->execute();
  my ($gene_id, $tran_id);
  $sth->bind_columns(\$gene_id,\$tran_id);
  my %gene_to_tran_canonical;
  while ($sth->fetch){
    $gene_to_tran_canonical{$gene_id} = $tran_id;
  }
  $sth->finish;

  print "Pre copy all $dbname from gene to canonical transcripts\n" if($self->verbose);

  $self->biomart_test();

  my $ins_dep_ix_sth = $self->xref->dbc->prepare("insert into identity_xref (object_xref_id, query_identity, target_identity) values(?, ?, ?)");

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

  print "Post copy all $dbname from gene to canonical transcripts\n" if($self->verbose);

  $self->biomart_test();


##########################################################################################
# HGNC/MGI synonyms are special as they are only on some prioritys so copy all the synonyms
# from the xref database to the xref database. Granted some will already be there but use
# update ignore just in case.
###########################################################################################


  #get the xref synonyms store as a hash of arrays.

  my $syn_hash = $self->get_official_synonyms();
  $sql = "insert into synonym (xref_id, synonym) values (?, ?)";
  my $add_syn_sth = $self->xref->dbc->prepare($sql);    
   

  #
  # get the OFDN (HGNC/MGI) xrefs in the xrefs and add the synonyms, plus create hash to get id from the display name
  # This is becouse not all MGI.HGNC's are loaded as these are priority xrefs.
  #

  my %display_label_to_id;
  

  $sql = 'select x.accession, sy.synonym from synonym sy, xref x, source so where x.xref_id = sy.xref_id and so.source_id = x.source_id and so.name like "'.$dbname.'"';

  
  $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($display_label, $acc, $syn,);
  $sth->bind_columns(\$acc,\$display_label);
  while($sth->fetch){
    $display_label_to_id{$display_label} = $acc;
  }
  $sth->finish;

  # get label to id from xref database to start with.
  $sql = 'select x.accession, x.label from xref x, source s where s.source_id = x.source_id and s.name like "'.$dbname.'"';

  
  $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  $sth->bind_columns(\$acc,\$display_label);
  while($sth->fetch){
    $display_label_to_id{$display_label} = $acc;
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

  $sth = $self->xref->dbc->prepare('select es.synonym, x.label from synonym es, xref x, source s where x.xref_id = es.xref_id and x.source_id = s.source_id and
 s.name = "EntrezGene"' );
  $sth->execute();
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;





#######################
#Do the naming bit now.
#######################

# get the vega external_sources


  my ($odn_curated_gene_id, $odn_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id);
  my ($odn_curated_tran_id, $odn_automatic_tran_id, $clone_based_vega_tran_id, $clone_based_ensembl_tran_id);
  my $odn_id;

  $sth = $self->xref->dbc->prepare("select source_id from source where name like ?");
  

  $sth->execute($dbname."_curated_gene");
  $sth->bind_columns(\$odn_curated_gene_id);
  $sth->fetch;
  
  $sth->execute($dbname."_automatic_gene");
  $sth->bind_columns(\$odn_automatic_gene_id);
  $sth->fetch;

  $sth->execute("Clone_based_vega_gene");
  $sth->bind_columns(\$clone_based_vega_gene_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_gene");
  $sth->bind_columns(\$clone_based_ensembl_gene_id);
  $sth->fetch;
  
  if(!defined($odn_curated_gene_id)){
    die "Could not find external database name ".$dbname."_curated_gene\n";
  }
  if(!defined($odn_automatic_gene_id)){
    die "Could not find external database name ".$dbname."_automatic_gene\n";
  }
  if(!defined($clone_based_vega_gene_id)){
    die "Could not find external database name Clone_based_vega_gene\n";
  }
  if(!defined($clone_based_ensembl_gene_id)){
    die "Could not find external database name Clone_based_ensembl_gene\n";
  }


  $sth->execute($dbname."_curated_transcript");
  $sth->bind_columns(\$odn_curated_tran_id);
  $sth->fetch;
  
  $sth->execute($dbname."_automatic_transcript");
  $sth->bind_columns(\$odn_automatic_tran_id);
  $sth->fetch;

  $sth->execute("Clone_based_vega_transcript");
  $sth->bind_columns(\$clone_based_vega_tran_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_transcript");
  $sth->bind_columns(\$clone_based_ensembl_tran_id);
  $sth->fetch;
  
  if(!defined($odn_curated_tran_id)){
    die "Could not find external database name ".$dbname."_curated_transcript\n";
  }
  if(!defined($odn_automatic_tran_id)){
    die "Could not find external database name ".$dbname."_automatic_transcript\n";
  }
  if(!defined($clone_based_vega_tran_id)){
    die "Could not find external database name Clone_based_vega_transcript\n";
  }
  if(!defined($clone_based_ensembl_tran_id)){
    die "Could not find external database name Clone_based_ensembl_transcript\n";
  }



  ###########################
  # Delete the old ones.
  ###########################

  # delete the synonyms first

  my $del_synonym_sql = "delete s from synonym s, xref x where s.xref_id = x.xref_id and x.source_id in ( $odn_curated_gene_id, $odn_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$odn_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->xref->dbc->prepare($del_synonym_sql);
  $sth->execute();

 
  my $del_identity_sql = "delete i from object_xref o, xref x, identity_xref i where i.object_xref_id = o.object_xref_id and x.xref_id = o.xref_id and x.source_id in ( $odn_curated_gene_id, $odn_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$odn_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->xref->dbc->prepare($del_identity_sql);
  $sth->execute();
 
  my $del_vega_sql = "delete o from object_xref o, xref x where x.xref_id = o.xref_id and x.source_id in ( $odn_curated_gene_id, $odn_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$odn_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->xref->dbc->prepare($del_vega_sql);
  $sth->execute();
 
  $del_vega_sql = "delete x from xref x where x.source_id in ($odn_curated_gene_id, $odn_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$odn_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->xref->dbc->prepare($del_vega_sql);
  $sth->execute();





  $del_vega_sql = "delete s from xref x, synonym s where s.xref_id = x.xref_id and x.source_id = $clone_based_vega_tran_id and x.info_type = 'MISC'"; # original ones added have info type of "DIRECT"

  $sth = $self->xref->dbc->prepare($del_vega_sql);
  $sth->execute();

  $del_vega_sql = "delete x from xref x where x.source_id = $clone_based_vega_tran_id and x.info_type = 'MISC'"; # original ones added have info type of "DIRECT"

  $sth = $self->xref->dbc->prepare($del_vega_sql);
  $sth->execute();



  ######################################################
  # Get the current max values for xref and object_xref
  ######################################################


  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $self->core->dbc);
  my $ga = $db->get_GeneAdaptor();

#  if(!defined($max_xref_id) or $max_xref_id == 0){
#    die "Sorry i dont believe there are no xrefs  max xref = 0\n";
#  }
#  if(!defined($max_object_xref_id) or $max_object_xref_id == 0){
#    die "Sorry i dont believe there are no object_xrefs  max object_xref = 0\n";
#  }




  ###########################
  # Process each Gene
  ###########################

  

  my %gene_to_transcripts;
  my %gene_id_to_stable_id;
  my %tran_id_to_stable_id;

  $sth = $self->xref->dbc->prepare("select gtt.gene_id, gtt.transcript_id, gsi.stable_id, tsi.stable_id 
                                        from gene_transcript_translation gtt, gene_stable_id gsi, transcript_stable_id tsi 
                                          where gtt.gene_id = gsi.internal_id and gtt.transcript_id = tsi.internal_id");

  $sth->execute;
  my $gsi;
  my $tsi;
  $sth->bind_columns(\$gene_id, \$tran_id, \$gsi, \$tsi);
  while ($sth->fetch){
    push @{$gene_to_transcripts{$gene_id}}, $tran_id;
    $gene_id_to_stable_id{$gene_id} = $gsi; 
    $tran_id_to_stable_id{$tran_id} = $tsi; 
  }
  

  my $total_gene_vega = 0;
  my $total_gene = 0;
  my $total_clone_name = 0;
  my $odn_count = 0;

  my $dbentrie_sth = $self->xref->dbc->prepare("select x.label, x.xref_id, s.priority from xref x, object_xref ox, source s where x.xref_id = ox.xref_id and 
                                                x.source_id = s.source_id and s.name = ? and ox.ox_status = 'DUMP_OUT' and ox.ensembl_id = ? and ox.ensembl_object_type = ? ");



  $sql = "insert into xref (xref_id, source_id, accession, label, version, species_id, info_type, info_text) values (?, ?, ?, ?,  0, ".$self->species_id.", 'MISC', ? )";
  my $ins_xref_sth = $self->xref->dbc->prepare($sql); 

  $sql = "insert into xref (xref_id, source_id, accession, label, version, species_id, info_type, info_text) values (?, ?, ?, ?,  ?, ".$self->species_id.", 'MISC', ? )";
  my $ins_xref_ver_sth = $self->xref->dbc->prepare($sql);


  # important or will crash and burn!!!
  my %xref_added; # store those added already $xref_added{$accession:$source_id} = $xref_id;

  my $get_xref_info_sth =  $self->xref->dbc->prepare("select x.label, x.accession, s.priority_description  from xref x, source s where xref_id = ? and s.source_id = x.source_id");

  #
  # Okay we assign unused_priority to be the number of time a vega transcript is attached to make sure
  # we get the BEST name for the gene (i.e. the one that appears the most)
  #
  my $ins_object_xref_sth =  $self->xref->dbc->prepare("insert into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id, linkage_type, ox_status, unused_priority) values (?, ?, ?, ?, 'MISC', 'DUMP_OUT', ?)");
  my %gene_clone_name_count;  
  foreach my $gene_id (keys %gene_to_transcripts){
    
    my @ODN=();
    my @VEGA_NAME=();
    my $CLONE_NAME = undef;
    my $xref_id;
    my $display;
    my $level;

    my $best_info = undef;

    $dbentrie_sth->execute($dbname, $gene_id, "Gene");
    $dbentrie_sth->bind_columns(\$display, \$xref_id,\$level);
    my $best_level=999;
    while($dbentrie_sth->fetch){
      if($level < $best_level){
	@ODN = ();
	push @ODN, $xref_id;
      }
      elsif($level == $best_level){
	push @ODN, $xref_id;
      }
    }
    
    if(scalar(@ODN)){
      $odn_count++;
    }
    
    my %no_vega; # hash now as we want to sort by stable id and hence need a key value pair
    my %vega_clone;
    my $vega_count = 0;
    my %name_count;
    
    
    my $count = 0;
    foreach my $tran_id ( @{$gene_to_transcripts{$gene_id}} ){
      my $VEGA = undef;
      $count = 0 ;

      $dbentrie_sth->execute($dbname."_curated_transcript", $tran_id, "Transcript");
      $dbentrie_sth->bind_columns(\$display, \$xref_id, \$level);
      while($dbentrie_sth->fetch){
	my($hgnc_bit, $num) = split(/-\d\d\d/,$display);
	$VEGA = $hgnc_bit;
        $name_count{$VEGA}++;
	$vega_count++;
	my $multiple = 1;
	if(scalar(@VEGA_NAME)){
	  foreach my $name (@VEGA_NAME){
	    if($name ne $VEGA){
	      if(defined($synonym{$name}) and uc($synonym{$name}) eq $VEGA){
		$multiple = 0;
	      }
	      if(defined($synonym{$VEGA}) and uc($synonym{$VEGA}) eq uc($name)){
		$multiple = 0;
	      }
	    }
	    else{
	      $multiple = 0;
	    }
	  }
	  if($multiple){
	    push @VEGA_NAME , $VEGA;
	  }
	}
	else{
	  push @VEGA_NAME , $VEGA;
	}
	$count++;
      }
       
      $dbentrie_sth->execute("Clone_based_vega_transcript", $tran_id, "Transcript");
      $dbentrie_sth->bind_columns(\$display, \$xref_id, \$level);
      while($dbentrie_sth->fetch){
	my($hgnc_bit, $num) = split(/-0\d\d/,$display);
	$CLONE_NAME = $hgnc_bit;
        $vega_clone{$tran_id_to_stable_id{$tran_id}} = $tran_id;
      }

      if($count == 0){
	$no_vega{$tran_id_to_stable_id{$tran_id}} = $tran_id;
       }
      if($count > 1){
	print STDERR "Problem: ".$tran_id_to_stable_id{$tran_id}." has more than one vega_transcript\n";
       }
      if($count == 1){
	my $name ="noidea";
 	if(defined($VEGA) and scalar(@ODN)){
	  my $found = 0;
	  foreach my $xref_id (@ODN){
	    my ($display, $acc, $text);
	    $get_xref_info_sth->execute($xref_id);
	    $get_xref_info_sth->bind_columns(\$display, \$acc, \$text);
	    $get_xref_info_sth->fetch();
	    if(uc($VEGA) eq uc($display)){
	      $found = 1;
	    }
	    elsif(defined($synonym{$VEGA}) and uc($synonym{$VEGA}) eq uc($display)){
	      $found = 1;
	    }
	    else{
	      $name = $display;
	    }
	  }
	  if(!$found){
	    print STDERR "Problem: ".$gene_id_to_stable_id{$gene_id}." linked to ".$dbname." (".$name.")   BUT ".$tran_id_to_stable_id{$tran_id}." linked to vega_transcript $VEGA????\n";	
	  }
	}
      }
    } # end for each transcript
    ####################################################################################
    # if there is at least one transcript
    # set vega_gene
    # loop through no_vega array and set vega_transcript_like for each starting at 201
    ####################################################################################
    if(scalar(@VEGA_NAME) > 1){
      print STDERR "Warning: gene ".$gene_id_to_stable_id{$gene_id}." has more than one vega_transcript these are (".join(', ',@VEGA_NAME).")\n";
    }	
    if($vega_count){

#
# Find the most common one
#

      my $v_name = $VEGA_NAME[0];
      my $top =0;
      foreach my $vn ( keys %name_count){
	if($name_count{$vn} > $top){
	  $top = $name_count{$vn};
	  $v_name = $vn;
	}
      }


      foreach my $name (keys %name_count){
	my $id = $display_label_to_id{$name};
	if(!defined($id)){
	  $id = $name;
	  print STDERR "Warning Could not find id for $name\n";
	}
	if(!defined($xref_added{$id.":".$odn_curated_gene_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $odn_curated_gene_id, $id, $name, "via havana");

	  $xref_added{$id.":".$odn_curated_gene_id} = $max_xref_id;

 	  if(defined($syn_hash->{$name})){
	    foreach my $syn (@{$syn_hash->{$name}}){
	      $add_syn_sth->execute($max_xref_id, $syn);
	    }
	  }
	  
	}
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $gene_id, 'Gene', $xref_added{$id.":".$odn_curated_gene_id}, $name_count{$name});
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
      }

      $name = $v_name;
      my $tran_name_ext = 201;
      foreach my $tran (sort keys %no_vega){
	my $id = $name."-".$tran_name_ext;
	if(!defined($xref_added{$id.":".$odn_automatic_tran_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $odn_automatic_tran_id, $id, $id, "via havana");
	  $xref_added{$id.":".$odn_automatic_tran_id} = $max_xref_id;
	}	
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $no_vega{$tran}, 'Transcript', $xref_added{$id.":".$odn_automatic_tran_id},undef);
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);

	$tran_name_ext++;
      }
    }
    
    ####################################################################################
    # if no vega_transcript but hgnc/mgi
    # set vega_gene_like to hgnc/mgi
    # loop through both arrays and set vega_transcript_like to hgnc-101 etc
    ####################################################################################
    elsif(scalar(@ODN)){
      foreach my $xref_id (@ODN){
	# now store xref_id so get xref details from this.

	# $get_xref_info_sth =  $self->xref->dbc->prepare("select x.label, x.accession, x.info_text  from xref x where xref_id = ?");
	my ($display, $acc, $text);
	$get_xref_info_sth->execute($xref_id);
	$get_xref_info_sth->bind_columns(\$display, \$acc, \$text);
	$get_xref_info_sth->fetch();
	my $id = $display_label_to_id{$display};
	if(!defined($display)){
	  print STDERR "Warning Could not find id for xref  $xref_id\n";
	}

	if(!defined($xref_added{$acc.":".$odn_automatic_gene_id})){
	  $max_xref_id++;
	  #   "insert into xref (xref_id, source_id, accession, label, version, species_id          , info_type, info_text) 
          #              values (?      , ?        , ?        , ?    ,  0    , ".$self->species_id.", 'MISC'   , ? )";
	  $ins_xref_sth->execute($max_xref_id, $odn_automatic_gene_id, $acc, $display, "via ".$text);
	  $xref_added{$acc.":".$odn_automatic_gene_id} = $max_xref_id;
	  if(defined($syn_hash->{$name})){
	    foreach my $syn (@{$syn_hash->{$name}}){
	      $add_syn_sth->execute($max_xref_id, $syn);
	    }
	  }
	}	
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $gene_id, 'Gene', $xref_added{$acc.":".$odn_automatic_gene_id}, undef);
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);


      }

      my $xref_id = $ODN[0];
      my ($display, $acc, $text);
      $get_xref_info_sth->execute($xref_id);
      $get_xref_info_sth->bind_columns(\$display, \$acc, \$text);
      $get_xref_info_sth->fetch();

      my $name = $display;


      my $tran_name_ext = 201;
      foreach my $tran (sort keys %no_vega){
	my $id = $name."-".$tran_name_ext;
	while(defined($xref_added{$id.":".$odn_automatic_tran_id})){
	  $tran_name_ext++;
	  $id = $name."-".$tran_name_ext;
	}
	if(!defined($xref_added{$id.":".$odn_automatic_tran_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $odn_automatic_tran_id, $id, $id, "via ".$text);
	  $xref_added{$id.":".$odn_automatic_tran_id} = $max_xref_id;
	}
	else{
	  print "ERROR: should not get here $id already defined?\n";
	}
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $no_vega{$tran}, 'Transcript', $xref_added{$id.":".$odn_automatic_tran_id}, undef);
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
	$tran_name_ext++;
      }
    }	
    
    
    ####################################################################################
    # if no vega_transcript and no hgnc/mgi use clone name
    # set vega_gene_like to clone name
    # loop through both arrays and set vega_transcript_like to clone_name-101 etc
    ####################################################################################
    else{
      if(defined($CLONE_NAME)){
	if(defined($gene_clone_name_count{$CLONE_NAME})){
	  $gene_clone_name_count{$CLONE_NAME}++;
	}
	else{
	  $gene_clone_name_count{$CLONE_NAME} = 1;
	}
#	$CLONE_NAME .= ".".$gene_clone_name_count{$CLONE_NAME};
	my $id = $CLONE_NAME.".".$gene_clone_name_count{$CLONE_NAME};
	if(!defined($xref_added{$id.":".$clone_based_vega_gene_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $clone_based_vega_gene_id, $id, $id, "via Vega clonename");
	  $xref_added{$id.":".$clone_based_vega_gene_id} = $max_xref_id;
	}	

	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $gene_id, 'Gene', $xref_added{$id.":".$clone_based_vega_gene_id}, undef);	
	$ins_dep_ix_sth->execute($max_object_xref_id, 100 , 100);


	my $tran_name_ext = 201;
	foreach my $tran (sort keys %no_vega){
	  next if(defined($vega_clone{$tran}));
	  my $id = $CLONE_NAME."-".$tran_name_ext;
	  if(!defined($xref_added{$id.":".$clone_based_vega_tran_id})){
	    $max_xref_id++;
	    $ins_xref_sth->execute($max_xref_id, $clone_based_vega_tran_id, $id, $id, "via Vega clonename");
	    $xref_added{$id.":".$clone_based_vega_tran_id} = $max_xref_id;
	  }	
	  
	  $max_object_xref_id++;
	  $ins_object_xref_sth->execute($max_object_xref_id, $no_vega{$tran}, 'Transcript', $xref_added{$id.":".$clone_based_vega_tran_id}, undef);	
	  $ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
	  

	  $tran_name_ext++;
	}
	
      }
      else{
	#get the clone name
	my $gene= $ga->fetch_by_dbID($gene_id);
	

#	my $new_gene = $gene->transform('clone');
	my $new_clone_name = undef;
#	if(defined($new_gene)){
#	  $new_clone_name = $new_gene->slice->seq_region_name;
#	}
#	else{
	  # on more than one clone?? try  
	  my $slice = $gene->slice->sub_Slice($gene->start,$gene->end,$gene->strand);
	  my $clone_projection = $slice->project('clone');
	  foreach my $seg (@$clone_projection) {
	    my $clone = $seg->to_Slice();
	    $new_clone_name = $clone->seq_region_name;
	  }
	  if(!defined($new_clone_name)){
	    # try going via contig
	    my $super_projection = $slice->project('contig');
	    my $super;
	    foreach my $seg (@$super_projection) {
	      $super = $seg->to_Slice();
	      $new_clone_name = $super->seq_region_name;
	    }
	    my $clone_projection = $super->project('clone');
	    foreach my $seg (@$clone_projection) {
	      my $clone = $seg->to_Slice();
	      $new_clone_name = $clone->seq_region_name;
	    }
	    if(!defined($new_clone_name)){
	      print STDERR "PROJECT failed for ".$gene->stable_id."\n";
	      next;
	    }
	  }
#	}

	if(defined($new_clone_name)){
	  $new_clone_name =~ s/[.]\d+//;    #remove .number	
	}
	else{
	  print "ERROR: No clone name can be found for ".$gene->stable_id."\n";
	}

	if(defined($gene_clone_name_count{$new_clone_name})){
	  $gene_clone_name_count{$new_clone_name}++;
	}
	else{
	  $gene_clone_name_count{$new_clone_name} = 1;
	}
	my $gene_name = $new_clone_name . ".".$gene_clone_name_count{$new_clone_name};

	# store the data
	if(!defined($xref_added{$gene_name.":".$clone_based_ensembl_gene_id})){
	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $clone_based_ensembl_gene_id, $gene_name, $gene_name, "via clonename");
	  $xref_added{$gene_name.":".$clone_based_ensembl_gene_id} = $max_xref_id;
	}	
	
	$max_object_xref_id++;
	$ins_object_xref_sth->execute($max_object_xref_id, $gene->dbID, 'Gene', $xref_added{$gene_name.":".$clone_based_ensembl_gene_id}, undef);	
	$ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);
	  	
	my $tran_name_ext = 201;
	foreach my $tran (sort keys %no_vega){
	  my $id = $new_clone_name."-".$tran_name_ext;
	  while(defined($xref_added{$id.":".$clone_based_ensembl_tran_id})){
	    $tran_name_ext++;
	    $id = $new_clone_name."-".$tran_name_ext;
	  }

	  $max_xref_id++;
	  $ins_xref_sth->execute($max_xref_id, $clone_based_ensembl_tran_id, $id, $id, "via clonename");
	  $xref_added{$id.":".$clone_based_ensembl_tran_id} = $max_xref_id;

	  $max_object_xref_id++;
	  $ins_object_xref_sth->execute($max_object_xref_id, $no_vega{$tran}, 'Transcript', $xref_added{$id.":".$clone_based_ensembl_tran_id}, undef);	
	  $ins_dep_ix_sth->execute($max_object_xref_id, 100, 100);



	  $tran_name_ext++;
	}
	

	
      }


    }
    
    if($vega_count){
      $total_gene_vega++;
    }
    $total_gene++;
    #  print "Finished Gene ".$gene->stable_id."\t$ODN\t $vega_count\n"
  }
  $add_syn_sth->finish;


  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('official_naming_done',now())");
  $sth_stat->execute();
  $sth_stat->finish;

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
}
  


sub biomart_test{
  my ($self) = @_;

  my $sql = 'SELECT ox.ensembl_object_type, COUNT(*), s.name  FROM xref x, object_xref ox, source s  WHERE x.xref_id = ox.xref_id AND s.source_id = x.source_id  and ox.ox_status = "DUMP_OUT" GROUP BY s.name, ox.ensembl_object_type';


  my $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($type, $count, $name);
  my ($last_type, $last_count, $last_name);
  $sth->bind_columns(\$type,\$count,\$name);
  $last_name = "ARSE";
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

1;
