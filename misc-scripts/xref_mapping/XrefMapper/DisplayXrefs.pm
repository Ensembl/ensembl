package XrefMapper::DisplayXrefs;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

use strict;
use warnings;
use XrefMapper::BasicMapper;

use Cwd;
use DBI;
use File::Basename;
use IPC::Open3;

my %genes_to_transcripts;
my %translation_to_transcript;
my %transcript_to_translation;
my %transcript_length;


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  return $self;
}


sub mapper{
  my ($self, $arg) = @_;

  (defined $arg) &&
    ($self->{_mapper} = $arg );
  return $self->{_mapper};
}



sub genes_and_transcripts_attributes_set{
  # Runs build_transcript_and_gene_display_xrefs,
  # new_build_gene_descriptions, build_gene_transcript_status and
  # build_meta_timestamp, and, if "-upload" is set, uses the SQL files
  # produced to update the core database.

  my ($self) = @_;

  my $status = $self->mapper->xref_latest_status();
  if($status ne "display_xref_done"){
    $self->build_transcript_and_gene_display_xrefs();
    
    my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('display_xref_done',now())");
    $sth_stat->execute();
    $sth_stat->finish;
    
  }

  $self->new_build_gene_descriptions();
  $self->build_gene_transcript_status();
  $self->build_meta_timestamp;

  my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('gene_description_done',now())");
  $sth_stat->execute();
  $sth_stat->finish;
  
  return 1;
}




sub build_transcript_and_gene_display_xrefs {
  my ($self) = @_;

  print "Building Transcript and Gene display_xrefs\n";


  my %external_name_to_id;
  my $sql1 = "SELECT external_db_id, db_name, status from external_db";
  
  my $sth1 = $self->core->dbc->prepare($sql1) || die "prepare failed for $sql1\n";
  $sth1->execute() || die "execute failed";
  my ($db_id, $name, $status);
  $sth1->bind_columns(\$db_id, \$name, \$status);
  while($sth1->fetch()){
    $external_name_to_id{$name}  = $db_id;
  }
  $sth1->finish;


  #
  # create a prioritised list of sources to use
  # and also a list of those xrefs to ignore 
  # where the source name is the key and the value is the string to test for 
  # 
  my ($presedence, $ignore) = @{$self->transcript_display_xref_sources()};
  my $i=0;
  my %level;
  print "precedense in reverse order:-\n";
  foreach my $name (reverse (@$presedence)){
    $i++;
    if(!defined($external_name_to_id{$name})){
      print STDERR "unknown external database name *$name* being used\n";
    }
    $level{$external_name_to_id{$name}} = $i;
    print "\t".$name."\t$i\n";
  }

  $self->build_genes_to_transcripts();

  $self->load_translation_to_transcript();
  

  # Gets object and identity xref data for 'SEQUENCE_MATCH' xrefs
  # for a given db_id of a given object type (gene/transcript/translation) 
  my $sql = (<<ESQL);
  SELECT ox.xref_id, ix.xref_identity, ix.ensembl_identity, 
         x.external_db_id, x.display_label, 
         e.db_name, ox.linkage_annotation
    FROM (object_xref ox, xref x, external_db e) 
      LEFT JOIN identity_xref ix 
        ON (ox.object_xref_id = ix.object_xref_id) 
   WHERE x.xref_id = ox.xref_id 
     AND ox.ensembl_object_type = ? 
     AND ox.ensembl_id = ? 
     AND x.info_type = 'SEQUENCE_MATCH'
     AND e.external_db_id = x.external_db_id
ESQL
  my $primary_sth = $self->core->dbc->prepare($sql)
    || die "prepare failed for $sql\n";


  # Gets object and identity xref data for 'DEPENDENT' xrefs
  # for a given db_id of a given object type (gene/transcript/translation) 
  $sql = (<<ZSQL);
  SELECT ox.xref_id, ix.xref_identity, ix.ensembl_identity, 
         x.external_db_id, x.display_label, 
         e.db_name, ox.linkage_annotation
    FROM (object_xref ox, xref x, external_db e, dependent_xref dx) 
      LEFT JOIN identity_xref ix 
        ON (ox.object_xref_id = ix.object_xref_id) 
   WHERE dx.master_xref_id = ox.xref_id
     AND dx.dependent_xref_id = x.xref_id
     AND ox.ensembl_object_type = ? 
     AND ox.ensembl_id = ? and x.info_type = 'DEPENDENT'
     AND e.external_db_id = x.external_db_id
ZSQL
  my $dependent_sth = $self->core->dbc->prepare($sql) 
    || die "prepare failed for $sql\n";


  # Gets object xref data for 'DIRECT' xrefs
  # for a given db_id of a given object type (gene/transcript/translation) 
  $sql = (<<QSQL);
  SELECT  x.xref_id, x.external_db_id, x.display_label, 
          e.db_name, o.linkage_annotation
     FROM object_xref o, xref x, external_db e 
    WHERE x.xref_id = o.xref_id 
      AND o.ensembl_object_type = ?
      AND o.ensembl_id = ? 
      AND (x.info_type = 'DIRECT' or x.info_type = 'MISC')
      AND e.external_db_id = x.external_db_id
QSQL
                             
  my $direct_sth = $self->core->dbc->prepare($sql) 
    || die "prepare failed for $sql\n";

  # Gets object xrefs data for any type of xref
  # for a given gene id.
  $sql = (<<GSQL);
  SELECT x.xref_id, x.external_db_id, e.db_name, o.linkage_annotation
    FROM object_xref o, xref x, external_db e   
   WHERE x.xref_id = o.xref_id 
     AND o.ensembl_object_type = 'Gene'
     AND o.ensembl_id = ?
     AND e.external_db_id = x.external_db_id
GSQL
       


  my $gene_sth = $self->core->dbc->prepare($sql) 
    || die "prepare failed for $sql\n";


  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $self->core->dbc->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  my $update_gene_sth = $self->core->dbc->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_tran_sth = $self->core->dbc->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");


  my $count =0;
  
  my ($xref_id, $qid, $tid, $ex_db_id, $display_label, $external_db_name, $linkage_annotation);
      
  # Loop for each gene_id
  foreach my $gene_id (keys %genes_to_transcripts) {
    my %percent_id;
    my %level_db;
    my %parent;
    my %percent_id_via_acc;
    my @gene_xrefs = ();
    
    # Query for object_xrefs attached directly o gene
    $gene_sth->execute($gene_id) || die "execute failed";
    $gene_sth->bind_columns
        (\$xref_id, \$ex_db_id, \$external_db_name, \$linkage_annotation);
    
    my $best_gene_xref  = 0;    # store xref
    my $best_gene_level = 0;    # store level
    my $best_gene_percent = 0;  # additoon of precentage ids

    while($gene_sth->fetch()){

      # Skip certain hard-coded external_db.names. 
      if(defined($$ignore{$external_db_name})){
	if($linkage_annotation =~ /$$ignore{$external_db_name}/){
	  #print( "Ignoring $xref_id as linkage_annotation has ",
          #       $$ignore{$external_db_name}." in it." ); # DEBUG

	  next;
	}
      }

      # 
      if(defined($level{$ex_db_id})){
	if($level{$ex_db_id} > $best_gene_level){
	  $best_gene_xref = $xref_id;
	  $best_gene_level = $level{$ex_db_id};
	}
      }
    }
    

    my @transcripts = @{$genes_to_transcripts{$gene_id}};
    foreach my $transcript_id (@transcripts) {

      my @transcript_xrefs = ();
      
      foreach my $type ("Transcript", "Translation"){
	my $ens_id;
	if($type eq "Transcript"){
	  $ens_id = $transcript_id;
	}
	else{
	  if(defined($transcript_to_translation{$transcript_id})){
	    $ens_id=$transcript_to_translation{$transcript_id};
	  }
	  else{
	    next;
	  }
	}
	$primary_sth->execute($type, $ens_id ) || die "execute failed";
	$primary_sth->bind_columns(\$xref_id, \$qid, \$tid, \$ex_db_id, 
				   \$display_label, \$external_db_name, 
				   \$linkage_annotation);
	while($primary_sth->fetch()){
	  if($level{$ex_db_id}  and $display_label =~ /\D+/ ){ #correct level and label is not just a number 	
	    if(defined($$ignore{$external_db_name})){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
#		print "Ignoring $xref_id as linkage_annotation has ".$$ignore{$external_db_name}." in it. DELETE THIS MESSAGE AFTER TESTING\n";
		next;
	      }
	    }

	    push @transcript_xrefs, $xref_id;
	    if(!defined($qid) || !defined($tid)){
	      print "PRIMARY $xref_id\n";
	      $percent_id{$xref_id} = 0;
	    }
	    else{
	      $percent_id{$xref_id}  = $qid + $tid;
	    }
	  
	    $level_db{$xref_id}  = $level{$ex_db_id};
	  }  
	}
	
	$dependent_sth->execute($type, $ens_id ) || die "execute failed";
	$dependent_sth->bind_columns(\$xref_id, \$qid, \$tid, \$ex_db_id, 
				     \$display_label, \$external_db_name, 
				     \$linkage_annotation);
	while($dependent_sth->fetch()){
	  if($level{$ex_db_id}  and $display_label =~ /\D+/){
	    if( defined($$ignore{$external_db_name}) and defined($linkage_annotation) ){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
#		print "Ignoring $xref_id as linkage_annotation has ".$$ignore{$external_db_name}." in it. DELETE THIS MESSAGE AFTER TESTING\n";
		next;
	      }
	    }
	    push @transcript_xrefs, $xref_id;
	    if(!defined($qid) || !defined($tid)){
#	      print "DEPENDENT $xref_id\n" if($ex_db_id != 1100); #HGNC has added one with no %ids.
	      $percent_id{$xref_id} = 0;
	    }
	    else{
	      $percent_id{$xref_id}  = $qid + $tid;
	    }
	    $level_db{$xref_id}  = $level{$ex_db_id};	    
	  }  
	}
	
	$direct_sth->execute($type, $ens_id ) || die "execute failed";
	$direct_sth->bind_columns(\$xref_id, \$ex_db_id, \$display_label,
				  \$external_db_name, \$linkage_annotation);
	while($direct_sth->fetch()){
	  if($level{$ex_db_id}  and $display_label =~ /\D+/){ 	
	    if(defined($$ignore{$external_db_name})){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
#		print "Ignoring $xref_id as linkage_annotation has ".$$ignore{$external_db_name}." in it. DELETE THIS MESSAGE AFTER TESTING\n";
		next;
	      }
	    }
	    push @transcript_xrefs, $xref_id;
	    $percent_id{$xref_id} = 0;
	    $level_db{$xref_id}  = $level{$ex_db_id};
	  }  
	}
      
      }      
      
      my $best_tran_xref  = 0; # store xref
      my $best_tran_level = 0; # store level
      my $best_tran_percent = 0; # store best %id total

      foreach my $t_xref_id (@transcript_xrefs) {
	if(defined($level_db{$t_xref_id}) and $level_db{$t_xref_id}){
	  if($level_db{$t_xref_id} < $best_tran_level){
	    next;
	  }

	  if($level_db{$t_xref_id} == $best_tran_level){
	    if($percent_id{$t_xref_id} < $best_tran_percent){
	      next;
	    }
	  }
	  $best_tran_percent = $percent_id{$t_xref_id};
	  $best_tran_level = $level_db{$t_xref_id};
	  $best_tran_xref  = $t_xref_id;
	}
      }       
      
      if($best_tran_xref){
	$update_tran_sth->execute($best_tran_xref, $transcript_id);
#        print TRANSCRIPT_DX "UPDATE transcript SET display_xref_id=" .$best_tran_xref. 
#            " WHERE transcript_id=" . $transcript_id . ";\n";
      }

      if($best_tran_level < $best_gene_level){
         next;
      }
      if($best_tran_level == $best_gene_level){
        if($best_tran_percent < $best_gene_percent){
          next;
        }
      }

      $best_gene_percent = $best_tran_percent;
      $best_gene_level   = $best_tran_level;
      $best_gene_xref    = $best_tran_xref;
    }
  
    if($best_gene_xref){
      $update_gene_sth->execute($best_gene_xref, $gene_id);
#      print GENE_DX "UPDATE gene g SET g.display_xref_id=" . $best_gene_xref . 
#	" WHERE g.gene_id=" . $gene_id . ";\n";
    }
  }

  # Done
  return;

}


sub load_translation_to_transcript{
  my ($self) = @_;

  my $sth = $self->core->dbc->prepare("SELECT translation_id, transcript_id FROM translation");
  $sth->execute();
  
  my ($translation_id, $transcript_id);
  $sth->bind_columns(\$translation_id, \$transcript_id);
  
  while ($sth->fetch()) {
    $translation_to_transcript{$translation_id} = $transcript_id;
    $transcript_to_translation{$transcript_id} = $translation_id if ($translation_id);
  }
}


sub build_genes_to_transcripts {

  my ($self) = @_;

  my $sql = "SELECT gene_id, transcript_id, seq_region_start, seq_region_end FROM transcript";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();

  my ($gene_id, $transcript_id, $start, $end);
  $sth->bind_columns(\$gene_id, \$transcript_id, \$start, \$end);

  # Note %genes_to_transcripts is global
  while ($sth->fetch()) {
    push @{$genes_to_transcripts{$gene_id}}, $transcript_id;
    $transcript_length{$transcript_id} = $end- $start;
  }

  $sth->finish
}



sub new_build_gene_descriptions{
  my ($self) = @_;
  

  print "Building gene Descriptions\n";

  my $update_gene_desc_sth =  $self->core->dbc->prepare("UPDATE gene SET description = ? where gene_id = ?");

  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET description = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  my @regexps = $self->gene_description_filter_regexps();

  if(scalar(@regexps) == 0){
    warn "no reg exps\n";
  }
  my @presedence = $self->gene_description_sources();



 if(!scalar(keys %translation_to_transcript)){
   $self->load_translation_to_transcript();
 }

  my %external_name_to_id;  
  my %ex_db_id_to_status;
  my %ex_db_id_to_display_name;
  my $sql1 = "SELECT external_db_id, db_name, status, db_display_name from external_db";
  
  my $sth1 = $self->core->dbc->prepare($sql1) 
    || die "prepare failed for $sql1\n";
  $sth1->execute() || die "execute failed";
  my ($db_id, $name, $status, $display_name);
  $sth1->bind_columns(\$db_id, \$name, \$status, \$display_name);
  while($sth1->fetch()){
    $external_name_to_id{$name} = $db_id;
    $ex_db_id_to_status{$db_id} = $status;
    $ex_db_id_to_display_name{$db_id}      = $display_name;
  }
  $sth1->finish;

  my $i=0;
  my %level;
  
  foreach my $ord (reverse (@presedence)){
    $i++;
    if(!defined($external_name_to_id{$ord})){
      print STDERR "Unknown external database name  *$ord* being used\n";
    }
    $level{$external_name_to_id{$ord}} = $i;

  }

  if(!scalar(keys %genes_to_transcripts)){
    $self->build_genes_to_transcripts();
  }


  my $sql = (<<ESQL);
  SELECT ox.xref_id, ix.xref_identity, ix.ensembl_identity,  x.external_db_id, x.description, x.dbprimary_acc
    FROM (object_xref ox, xref x) 
      LEFT JOIN identity_xref ix ON (ox.object_xref_id = ix.object_xref_id) 
	WHERE x.xref_id = ox.xref_id and ox.ensembl_object_type = ? 
              and ox.ensembl_id = ? and x.info_type = 'SEQUENCE_MATCH'
ESQL

  my $primary_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";

  $sql = (<<ZSQL);
  SELECT ox.xref_id, ix.xref_identity, ix.ensembl_identity, 
         x.external_db_id, x.description, x.dbprimary_acc 
    FROM (object_xref ox, xref x, external_db e, dependent_xref dx) 
      LEFT JOIN identity_xref ix 
        ON (ox.object_xref_id = ix.object_xref_id) 
   WHERE dx.master_xref_id = ox.xref_id
     AND dx.dependent_xref_id = x.xref_id
     AND ox.ensembl_object_type = ? 
     AND ox.ensembl_id = ? and x.info_type = 'DEPENDENT'
     AND e.external_db_id = x.external_db_id
ZSQL
#  $sql = (<<ZSQL);
#  SELECT ox.xref_id, ix.query_identity, ix.target_identity, x.external_db_id, x.description, x.dbprimary_acc
#    FROM (object_xref ox, xref x) 
#      LEFT JOIN identity_xref_temp ix ON (ox.object_xref_id = ix.object_xref_id) 
#	WHERE x.xref_id = ox.xref_id and ox.ensembl_object_type = ? 
#              and ox.ensembl_id = ? and x.info_type = 'DEPENDENT'
#ZSQL
 
  my $dependent_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";


  $sql = (<<QSQL);
  SELECT x.xref_id, x.external_db_id, x.description, x.dbprimary_acc
   FROM object_xref o, xref x  
    WHERE x.xref_id = o.xref_id 
        and o.ensembl_object_type = ? and o.ensembl_id = ? and x.info_type = 'DIRECT'
QSQL
                             
  my $direct_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";

  $sql = (<<GSQL);
  SELECT x.xref_id, x.external_db_id, x.description, x.dbprimary_acc
   FROM object_xref o, xref x  
    WHERE x.xref_id = o.xref_id 
        and o.ensembl_object_type = 'Gene' and o.ensembl_id = ?
GSQL

  my $gene_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";
 
  my $count =0;
  
  my ($xref_id, $qid, $tid, $ex_db_id, $description, $acc);
  

  # Do not use the %id here just use the Best one using the @words to do this
  my @words = qw(unknown hypothetical putative novel probable [0-9]{3} kDa fragment cdna protein);
  my $trembl_id = $external_name_to_id{"Uniprot/SPTREMBL"};


  my $checked = 0;
  my $added   = 0;
  my $removed = 0; 

  foreach my $gene_id (keys %genes_to_transcripts) {
    
    my %percent_id;
    my %level_db;
    my %parent;
    my %ex_db;
    my %xref_descriptions;
    my %xref_accessions;
    my @gene_xrefs = ();
    my @transcript_xrefs = ();

    my $best_gene_xref  = 0;    # store xref
    my $best_gene_level = 0;    # store level
    my $best_gene_percent = 0;  # additoon of precentage ids
    my $best_gene_length  = 0;  # best transcript for the genes length

    $gene_sth->execute($gene_id) || die "execute failed";
    $gene_sth->bind_columns(\$xref_id, \$ex_db_id, \$description, \$acc);
    
    while($gene_sth->fetch()){
      $checked++;
      if ($description and defined($level{$ex_db_id})) {
	my $filtered_description = $self->filter_by_regexp($description, \@regexps);
	if ($filtered_description ne "") {
	  $xref_descriptions{$xref_id} = $description;
	  $xref_accessions{$xref_id} = $acc;
	  if($level{$ex_db_id} > $best_gene_level){
	    $best_gene_xref = $xref_id;
	    $best_gene_level = $level{$ex_db_id};
	  }
	  $added++;
	} else {
          $removed++;
	}
      }
    }
    
    
    my @transcripts = @{$genes_to_transcripts{$gene_id}};
    foreach my $transcript_id (@transcripts) {
      foreach my $type("Transcript", "Translation"){
	my $ens_id;
	if($type eq "Transcript"){
	  $ens_id = $transcript_id;
	}
	else{
	  if(defined($transcript_to_translation{$transcript_id})){
	    $ens_id=$transcript_to_translation{$transcript_id};
	  }
	  else{
	    next;
	  }
	}
	
	$primary_sth->execute($type, $ens_id) || die "execute failed";
	$primary_sth->bind_columns(\$xref_id, \$qid, \$tid, \$ex_db_id, \$description, \$acc);
	while($primary_sth->fetch()){
	
	  if($level{$ex_db_id}){
	    $checked++;
	    if ($description) {
	      my $filtered_description = $self->filter_by_regexp($description, \@regexps);
	      if ($filtered_description ne "") {
		$xref_descriptions{$xref_id} = $description;
		$xref_accessions{$xref_id} = $acc;
		if(defined($qid)){
		  $percent_id{$xref_id}  = $qid + $tid;
		}
		else{
		  print "WARN: xref_id $xref_id PRIMARY can't find percrnt id\n";
		  $percent_id{$xref_id} = 50;
		}
		push @transcript_xrefs, $xref_id;  # added ?
		$percent_id{$xref_id}  = $qid + $tid;
		$ex_db{$xref_id} = $ex_db_id;
		$level_db{$xref_id}  = $level{$ex_db_id};
		$added++;
	      } else {
		$removed++;
	      }
	    }
	  }  
	}
      
	$dependent_sth->execute($type, $ens_id) || die "execute failed";
	$dependent_sth->bind_columns(\$xref_id, \$qid, \$tid, \$ex_db_id, \$description, \$acc);
	while($dependent_sth->fetch()){
	  if($level{$ex_db_id}){
	    $checked++;
	    if ($description) {
	      my $filtered_description = $self->filter_by_regexp($description, \@regexps);
	      if ($filtered_description ne "") {
		$xref_descriptions{$xref_id} = $description;
		$xref_accessions{$xref_id} = $acc;
		push @transcript_xrefs, $xref_id;
		if(defined($qid)){
		  $percent_id{$xref_id}  = $qid + $tid;
		}
		else{
#		  print "WARN: xref_id $xref_id DEPEND can't find percrnt id. Type = $type, ens_id == $ens_id, acc =$acc\n";
		  $percent_id{$xref_id} = 50;
		}
		$ex_db{$xref_id} = $ex_db_id;	
		$level_db{$xref_id}  = $level{$ex_db_id};
		$added++;
	      } else {
		$removed++;
	      }
	    }
	  }  
	}	
	
	$direct_sth->execute($type, $ens_id) || die "execute failed";
	$direct_sth->bind_columns(\$xref_id, \$ex_db_id, \$description, \$acc);
	while($direct_sth->fetch()){
	  if($level{$ex_db_id}){
	    $checked++;
	    if ($description) {
	      my $filtered_description = $self->filter_by_regexp($description, \@regexps);
	      if ($filtered_description ne "") {
		$xref_descriptions{$xref_id} = $description;
		$xref_accessions{$xref_id} = $acc;
		push @transcript_xrefs, $xref_id;
		$percent_id{$xref_id}  = 200;
		$ex_db{$xref_id} = $ex_db_id;
		$level_db{$xref_id}  = $level{$ex_db_id};
		$added++;
	      } else {
		$removed++;
	      }
	    }
	  }  
	}
      
      }
      
      my $best_tran_xref  = 0;   # store xref
      my $best_tran_level = 0;   # store level
      my $best_tran_percent = 0; # store best %id total
      my $best_tran_length =0 ;  # store length of the best
      
      foreach my $xref_id (@transcript_xrefs) {
	if(defined($xref_descriptions{$xref_id})){
	  if(defined($level_db{$xref_id} and $level_db{$xref_id})){
	    if($level_db{$xref_id} < $best_tran_level){
	      next;
	    }
	    if($level_db{$xref_id} == $best_tran_level){
	      if($percent_id{$xref_id} < $best_tran_percent){
		next;
	      }
	    }
	  
	    $best_tran_percent = $percent_id{$xref_id};
	    $best_tran_level = $level_db{$xref_id};
	    $best_tran_xref  = $xref_id;
	    $best_tran_length =  $transcript_length{$transcript_id};
	  }      
	}  
      }       
      
      if($best_tran_level < $best_gene_level){
	next;
      }
      if($best_tran_level == $best_gene_level){
	if($best_tran_percent < $best_gene_percent){
	  next;
	}
	elsif($best_tran_percent == $best_gene_percent){
	if($transcript_length{$transcript_id} < $best_gene_length){
	  next;
	}
      } 
      
      }
    
      $best_gene_percent = $best_tran_percent;
      $best_gene_level   = $best_tran_level;
      $best_gene_xref    = $best_tran_xref;
      $best_gene_length  = $transcript_length{$transcript_id};
    }
    
    if($best_gene_xref){
      my $description = $xref_descriptions{$best_gene_xref};
      my $acc = $xref_accessions{$best_gene_xref};
      
      $description =~ s/\"//ig; # remove " as they will cause problems in .sql files
      
      my $desc = $description . " [Source:".$ex_db_id_to_display_name{$ex_db{$best_gene_xref}}.";Acc:$acc]";
 
      $update_gene_desc_sth->execute($desc, $gene_id) if($description);
      
    }
  }

  return;
}

sub build_gene_transcript_status{
  # Creates the files that contain the SQL needed to (re)set the
  # gene.status and transcript.status values
  my $self = shift;
  
  my $reset_sth = $self->core->dbc->prepare('UPDATE gene SET status = "NOVEL"');
  $reset_sth->execute();
  $reset_sth->finish;

  $reset_sth = $self->core->dbc->prepare('UPDATE transcript SET status = "NOVEL"');
  $reset_sth->execute();
  $reset_sth->finish;
  
  my $update_gene_sth = $self->core->dbc->prepare('UPDATE gene SET status = "KNOWN" where gene_id = ?');
  my $update_tran_sth = $self->core->dbc->prepare('UPDATE transcript SET status = "KNOWN" where transcript_id = ?');

  #create a hash known which ONLY has databases names of those that are KNOWN and KNOWNXREF
  my %known;
  my $sth = $self->core->dbc->prepare("select db_name from external_db where status in ('KNOWNXREF','KNOWN')");
  $sth->execute();
  my ($name);
  $sth->bind_columns(\$name);
  while($sth->fetch){
    $known{$name} = 1;
  }
  $sth->finish;
  
  
  # loop throught the gene and all transcript until you find KNOWN/KNOWNXREF as a status
  my $ensembl = $self->core;
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $ensembl->dbc);
  my $gene_adaptor = $db->get_GeneAdaptor();

  my @genes = @{$gene_adaptor->fetch_all()};

  while (my $gene = shift @genes){
    my $gene_found = 0;
    my @dbentries = @{$gene->get_all_DBEntries()};
    foreach my $dbe (@dbentries){
      if(defined($known{$dbe->dbname})){
	$gene_found =1;
      }
    }
    my $one_tran_found = 0;
    foreach my $tr (@{$gene->get_all_Transcripts}){
      my $tran_found = 0;
      foreach my $dbe (@{$tr->get_all_DBLinks}){
	if(defined($known{$dbe->dbname})){
	  $tran_found = 1;
	  $one_tran_found = 1;
	}
      }
      if($tran_found or $gene_found){
	$update_tran_sth->execute($tr->dbID);
      }
    }
    if($gene_found or $one_tran_found){
      $update_gene_sth->execute($gene->dbID);
    }
  }

  return;
}

sub build_meta_timestamp{
  # Creates a file that contains the SQL needed to (re)set the 
  # 'xref.timestamp' key of the meta table.
  my $self = shift;


  my $sth = $self->core->dbc->prepare("DELETE FROM meta WHERE meta_key='xref.timestamp'");
  $sth->execute();
  $sth->finish;

  $sth = $self->core->dbc->prepare("INSERT INTO meta (meta_key,meta_value) VALUES ('xref.timestamp', NOW())");
  $sth->execute();
  $sth->finish;

  return;
}


1;
