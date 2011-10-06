package XrefMapper::DisplayXrefs;
use strict;

use vars '@ISA';
@ISA = qw{ XrefMapper::BasicMapper };

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


#
# ignore should be some sql to return object_xref_ids that should be ignored. FOR full mode METHOD 
# ignore should return regexp and source name as key for update METHODS
#

sub gene_description_sources {

  return ("RFAM",
          "RNAMMER",
          "TRNASCAN_SE",
	  "miRBase",
          "HGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
	  "RefSeq_dna",
	  "Uniprot/Varsplic",
	  "Uniprot/SPTREMBL");

}

sub gene_description_filter_regexps {

  return ();

}

sub transcript_display_xref_sources {
  my $self     = shift;
  my $fullmode = shift;

  my @list = qw(HGNC
                MGI
                Clone_based_vega_gene
                Clone_based_ensembl_gene
                HGNC_transcript_name
                MGI_transcript_name
                Clone_based_vega_transcript
                Clone_based_ensembl_transcript
                miRBase
                RFAM
                IMGT/GENE_DB
                SGD
                flybase_symbol
                Anopheles_symbol
                Genoscope_annotated_gene
                Uniprot/SWISSPROT
                Uniprot/Varsplic
                Uniprot/SPTREMBL
                EntrezGene);

  my %ignore;
  

  # Both methods

  if(!$fullmode){
    $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  }
  else{
    $ignore{"EntrezGene"} =(<<AIGN);

SELECT ox.object_xref_id 
    FROM object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 
     WHERE ox.object_xref_id = dx.object_xref_id AND dx.dependent_xref_id = x1.xref_id 
     AND x1.source_id = s1.source_id and s1.name = 'EntrezGene' 
     AND x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id 
     AND (s2.name like 'Refseq_dna_predicted' or s2.name like 'RefSeq_peptide_predicted') 
     AND ox.ox_status = 'DUMP_OUT'
AIGN

    $ignore{"Uniprot/SPTREMBL"} =(<<BIGN);
SELECT object_xref_id
    FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
     WHERE ox_status = 'DUMP_OUT' AND name = 'Uniprot/SPTREMBL' 
      AND priority_description = 'protein_evidence_gt_3'
BIGN

  }

  return [\@list,\%ignore];

}


sub new {
  my($class, $mapper) = @_;

  my $self ={};
  bless $self,$class;
  $self->core($mapper->core);
  $self->xref($mapper->xref);
  $self->mapper($mapper);
  $self->verbose($mapper->verbose);
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

  my ($self, $fullmode, $noxref_database) = @_;

  my $status;
  if(defined($noxref_database)){
    $status = "none";
  }
  else{
    $status = $self->mapper->xref_latest_status();
  }
     
  if(!$fullmode){
    
    if($status ne "display_xref_done"){
      $self->build_transcript_and_gene_display_xrefs();
      
      if(!defined($noxref_database)){
	my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('display_xref_done',now())");
	$sth_stat->execute();
	$sth_stat->finish;
      }
      
    }
    
    $self->new_build_gene_descriptions();
    $self->build_gene_transcript_status();
  }
  else{
    if($self->mapper->can("set_display_xrefs")){
      $self->mapper->set_display_xrefs();
    }
    else{
      $self->set_display_xrefs();
    }	
    my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('display_xref_done',now())");
    $sth_stat->execute();
    $sth_stat->finish;
    if($self->mapper->can("set_gene_descriptions")){
      $self->mapper->set_gene_descriptions();
    }
    else{
      $self->set_gene_descriptions();
    }	
    $self->set_status(); # set KNOWN,NOVEL etc 
  }

  $self->build_meta_timestamp;

  # Special removal of LRG transcript display xref, xref and object_xrefs;

  my $sth_lrg =  $self->core->dbc->prepare('DELETE ox, x  FROM object_xref ox, xref x, transcript t, transcript_stable_id tsi WHERE ox.xref_id = x.xref_id and tsi.transcript_id = t.transcript_id and t.display_xref_id = x.xref_id and tsi.stable_id like "LRG%"');
  $sth_lrg->execute;

  $sth_lrg = $self->core->dbc->prepare('UPDATE transcript t, transcript_stable_id tsi SET t.display_xref_id = null WHERE  tsi.stable_id like "LRG%" and tsi.transcript_id = t.transcript_id');
  $sth_lrg->execute;

  #End of Special

  $sth_lrg = $self->core->dbc->prepare("UPDATE xref SET info_text=null WHERE info_text=''");
  $sth_lrg->execute;


  if(!defined($noxref_database)){
    my $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('gene_description_done',now())");
    $sth_stat->execute();
    $sth_stat->finish;
  }

  return 1;
}

sub set_gene_descriptions_from_display_xref{
  my $self = shift;
  
  $self->set_gene_descriptions(1);
#  print "NO NEED to write sub as done when display_xref was set :-)\n";
}




sub set_display_xrefs_from_stable_table{
  my $self = shift;
  print "Setting Transcript and Gene display_xrefs from xref database into core and setting the desc\n" if ($self->verbose);

  my $xref_offset = $self->get_meta_value("xref_offset");

  print "Using xref_off set of $xref_offset\n" if($self->verbose);

  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $self->core->dbc->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  $reset_sth = $self->core->dbc->prepare("UPDATE gene SET description = null");
  $reset_sth->execute();
  $reset_sth->finish;


  my %name_to_external_name;
  my $sql = "select external_db_id, db_name, db_display_name from external_db";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  my ($id, $name, $display_name);
  $sth->bind_columns(\$id, \$name, \$display_name);
  while($sth->fetch()){
    $name_to_external_name{$name} = $display_name;
   }
  $sth->finish;

  my %source_id_to_external_name;

  $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$id, \$name);

  while($sth->fetch()){
     if(defined($name_to_external_name{$name})){
      $source_id_to_external_name{$id} = $name_to_external_name{$name};
    }
  }
  $sth->finish;


  my $update_gene_sth = $self->core->dbc->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_gene_desc_sth = $self->core->dbc->prepare("UPDATE gene g SET g.description= ? WHERE g.gene_id=?");

  my $update_tran_sth = $self->core->dbc->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

  my $get_gene_display_xref = $self->xref->dbc->prepare("SELECT gsi.internal_id, gsi.display_xref_id, x.description ,x.source_id, x.accession
                                                              FROM gene_stable_id gsi, xref x 
                                                                 WHERE gsi.display_xref_id = x.xref_id");

  my $get_tran_display_xref = $self->xref->dbc->prepare("SELECT gsi.internal_id, gsi.display_xref_id from transcript_stable_id gsi");

  $reset_sth = $self->xref->dbc->prepare("UPDATE gene_stable_id gsi SET gsi.desc_set=0");
  $reset_sth->execute();

  my $set_desc_done_sth = $self->xref->dbc->prepare("UPDATE gene_stable_id gsi SET gsi.desc_set=1 WHERE gsi.internal_id=?");

  $get_gene_display_xref->execute();
  my $xref_id;
  my $desc;
  my $gene_id;
  my $source_id;
  my $label;
  $get_gene_display_xref->bind_columns(\$gene_id, \$xref_id, \$desc, \$source_id, \$label);
  my $gene_count =0;
  while($get_gene_display_xref->fetch()){

    $update_gene_sth->execute($xref_id+$xref_offset, $gene_id);

    if (defined($desc) and $desc ne "") {
      $desc .= " [Source:".$source_id_to_external_name{$source_id}.";Acc:".$label."]";
      $update_gene_desc_sth->execute($desc,$gene_id);
      $set_desc_done_sth->execute($gene_id);
      $gene_count++;
    }
#    else{
#      $update_gene_desc_set($gene_id);
#    }
  }

  $update_gene_desc_sth->finish;
  $update_gene_sth->finish;

  print "$gene_count gene descriptions added\n" if($self->verbose);

  $get_tran_display_xref->execute();
  my $tran_id;
  $get_tran_display_xref->bind_columns(\$tran_id, \$xref_id);

  while($get_tran_display_xref->fetch()){
    if(defined($xref_id)){
      $update_tran_sth->execute($xref_id+$xref_offset, $tran_id);
      if(!defined($tran_id) || !defined($xref_id) || !defined($xref_offset)){
	print "PROB: tran_id = $tran_id\nxref_id = $xref_id\n$xref_offset = $xref_offset\n";
      }
    }
  }	
}



sub set_status{
  my $self = shift;

# set all genes to NOVEL

  
  my $reset_sth = $self->core->dbc->prepare('UPDATE gene SET status = "NOVEL"');
  $reset_sth->execute();
  $reset_sth->finish;

  $reset_sth = $self->core->dbc->prepare('UPDATE transcript SET status = "NOVEL"');
  $reset_sth->execute();
  $reset_sth->finish;
  
  my $update_gene_sth = $self->core->dbc->prepare('UPDATE gene SET status = ? where gene_id = ?');
  my $update_tran_sth = $self->core->dbc->prepare('UPDATE transcript SET status = ? where transcript_id = ?');

  
my $known_xref_sql =(<<DXS);
select  distinct 
        IF (ox.ensembl_object_type = 'Gene',        gtt_gene.gene_id,
        IF (ox.ensembl_object_type = 'Transcript',  gtt_transcript.gene_id,
                                                    gtt_translation.gene_id)) AS gene_id,

        IF (ox.ensembl_object_type = 'Gene',        gtt_gene.transcript_id,
        IF (ox.ensembl_object_type = 'Transcript',  gtt_transcript.transcript_id,
                                                    gtt_translation.transcript_id)) AS transcript_id
from    (   source s
      join    (   xref x
        join      (   object_xref ox
                  ) using (xref_id)
              ) using (source_id)
          )
  left join gene_transcript_translation gtt_gene
    on (gtt_gene.gene_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_transcript
    on (gtt_transcript.transcript_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_translation
    on (gtt_translation.translation_id = ox.ensembl_id)
where   ox.ox_status = 'DUMP_OUT'
        AND s.status like "KNOWN%"
        ORDER BY gene_id DESC, transcript_id DESC
DXS


  my $last_gene = 0;

  my $known_xref_sth = $self->xref->dbc->prepare($known_xref_sql);

  $known_xref_sth->execute();
  my ($gene_id, $transcript_id);  # remove labvel after testig it is not needed
  $known_xref_sth->bind_columns(\$gene_id, \$transcript_id);
  while($known_xref_sth->fetch()){
    if($gene_id != $last_gene){
      $update_gene_sth->execute("KNOWN",$gene_id);
      $last_gene = $gene_id;
    } 
    $update_tran_sth->execute("KNOWN",$transcript_id);
  }


  # 1) load list of stable_gene_id from xref database and covert to internal id in
  #    new core database table.
  #    Use this table to reset havana gene/transcript status.

  if(!scalar(keys %genes_to_transcripts)){
    $self->build_genes_to_transcripts();
  }

  #
  # Reset status for those from vega
  #

#  my %gene_id_to_status;
  my $gene_status_sth = $self->xref->dbc->prepare("SELECT gsi.internal_id, hs.status FROM gene_stable_id gsi, havana_status hs WHERE hs.stable_id = gsi.stable_id") 
    || die "Could not prepare gene_status_sth";

  $gene_status_sth->execute();
  my ($internal_id, $status);
  $gene_status_sth->bind_columns(\$internal_id,\$status);
  while($gene_status_sth->fetch()){
#    $gene_id_to_status{$internal_id} = $status;
    $update_gene_sth->execute($status, $internal_id);
  }
  $gene_status_sth->finish();

  #
  # need to create a transcript_id to status hash
  #

#  my %transcript_id_to_status;
  my $transcript_status_sth = $self->xref->dbc->prepare("SELECT tsi.internal_id, hs.status FROM transcript_stable_id tsi, havana_status hs WHERE hs.stable_id = tsi.stable_id") 
    || die "Could not prepare transcript_status_sth";

  $transcript_status_sth->execute();
  $transcript_status_sth->bind_columns(\$internal_id,\$status);
  while($transcript_status_sth->fetch()){
    #    $transcript_id_to_status{$internal_id} = $status;
    $update_tran_sth->execute($status,$internal_id);  
  }
  $transcript_status_sth->finish();


#  #
#  # Get some stats
#  #
#  my %count;

#  # Loop for each gene_id
#  foreach my $gene_id (keys %genes_to_transcripts) {
#    # if gene havana status is set
#    if(defined($gene_id_to_status{$gene_id})){
#      #if set check each transcript
#      my $missing = 0;
#      foreach my $tran_id (@{$genes_to_transcripts{$gene_id}}){
#	if(!defined($transcript_id_to_status{$tran_id})){
#	  $missing++;
#	}	
#      }
#      #if all transcript have havana status
#      if(!$missing){
#	#   set the status for each transcript and the gene
#	foreach my $tran_id (@{$genes_to_transcripts{$gene_id}}){
#	  $update_tran_sth->execute($transcript_id_to_status{$tran_id},$tran_id);
#	}
#	$update_gene_sth->execute($gene_id_to_status{$gene_id},$gene_id);
#	$count{"Setting for all transcripts and gene"}++;
#      }
#      else{
#	$count{"One or more transcripts failed test"}++;
#      }
#    }
#    else{
#      $count{"No havana gene status"}++;
#    }
#  }
      
#  print "\n";
#  foreach my $key (keys %count){
#    print "$key\t".$count{$key}."\n";
#  }

  $known_xref_sth->finish;
  $update_gene_sth->finish;
  $update_tran_sth->finish;
  

}


sub build_transcript_and_gene_display_xrefs {
  my ($self) = @_;

  print "Building Transcript and Gene display_xrefs Using OLD methods accessing core database alone.\n" if ($self->verbose);

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

#if mapper can transcript_display_xref_sources # from species.pm file
#do it else.....
  my $presedence;
  my $ignore; 
  if( $self->mapper->can("transcript_display_xref_sources") ){
    ($presedence, $ignore) = @{$self->mapper->transcript_display_xref_sources(0)}; # UPDATE MODE
  }
  else{
    ($presedence, $ignore) = @{$self->transcript_display_xref_sources(0)}; # UPDATE MODE
  }
  my $i=0;
  my %level;

  my $last_name = "";
  print "precedense in reverse order:-\n" if($self->verbose);
  foreach my $name (reverse (@$presedence)){
    $i++;
    if(!defined($external_name_to_id{$name})){
      print STDERR "unknown external database name *$name* being used\n";
    }
    $level{$external_name_to_id{$name}} = $i;
    if($name ne $last_name){
      print "\t".$name."\t$i\n" if($self->verbose);
    }
    $last_name = $name;
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
	  if(defined($level{$ex_db_id})  and $display_label =~ /\D+/ ){ #correct level and label is not just a number 	
	    if(defined($$ignore{$external_db_name})){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
		next;
	      }
	    }

	    push @transcript_xrefs, $xref_id;
	    if(!defined($qid) || !defined($tid)){
	      print STDERR "PROBLEM:: PRIMARY $xref_id\n";
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
		next;
	      }
	    }
	    push @transcript_xrefs, $xref_id;
	    if(!defined($qid) || !defined($tid)){
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
	  if(defined($level{$ex_db_id})  and $display_label =~ /\D+/){ 	
	    if(defined($$ignore{$external_db_name})){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
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
  

  print "Building gene Descriptions the OLD way from core database\n" if ($self->verbose);


  my $update_gene_desc_sth =  $self->core->dbc->prepare("UPDATE gene SET description = ? where gene_id = ?");

  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET description = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  my @presedence;
  my @regexps;
  if( $self->mapper->can("gene_description_sources") ){
    @presedence = $self->mapper->gene_description_sources();
  }
  else{
    @presedence = $self->gene_description_sources();
  }

  if( $self->mapper->can("gene_description_filter_regexps") ){
    @regexps = $self->mapper->gene_description_filter_regexps();
  }
  else{
    @regexps = $self->gene_description_filter_regexps();
  }

  if(scalar(@regexps) == 0){
    warn "no reg exps\n" if($self->verbose);
  }



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

  my %no_source_name_in_desc;
  if( $self->mapper->can("no_source_label_list") ){
    foreach my $name (@{$self->mapper->no_source_label_list()}){
      my $id = $external_name_to_id{$name};
      print "$name will have no [Source:..] info\n";
      $no_source_name_in_desc{$id} = 1;
    }
  }

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
	  $ex_db{$xref_id} = $ex_db_id;
	  if($level{$ex_db_id} > $best_gene_level){
	    $best_gene_xref = $xref_id;
	    $best_gene_level = $level{$ex_db_id};
	  }
          # if same level then use prioritys to work out best
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
		  print STDERR "WARN: xref_id $xref_id PRIMARY can't find percent id\n";
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
      
      if(!defined($ex_db_id_to_display_name{$ex_db{$best_gene_xref}})){
	print STDERR "Could not find display name for gene $best_gene_xref for external db ".$ex_db{$best_gene_xref}."\n";
      }
      my $desc = $description;
      if(!defined($no_source_name_in_desc{$ex_db{$best_gene_xref}})){
	$desc .= " [Source:".$ex_db_id_to_display_name{$ex_db{$best_gene_xref}}.";Acc:$acc]";
      }

      $update_gene_desc_sth->execute($desc, $gene_id) if($description);
      
    }
  }
  # remove until dependent_xref added to core database

  # Now is part of the xref system.

#  my $sth = $self->core->dbc->prepare("drop table dependent_xref");
#  $sth->execute || die "Could not drop temp table dependent_xref\n";
#  $sth->finish;  

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



sub set_display_xrefs{
  my $self = shift;


  print "Building Transcript and Gene display_xrefs using xref database\n" if ($self->verbose);

  my $xref_offset = $self->get_meta_value("xref_offset");

  print "Using xref_off set of $xref_offset\n" if($self->verbose);

  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $self->core->dbc->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  my $update_gene_sth = $self->core->dbc->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_tran_sth = $self->core->dbc->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

 #get hash for sources in hash
  #get priority description

my $sql =(<<SQL); 
  CREATE TABLE display_xref_prioritys(
    source_id INT NOT NULL,
    priority       INT NOT NULL,
    PRIMARY KEY (source_id)
  ) COLLATE=latin1_swedish_ci TYPE=InnoDB
SQL

  my $sth = $self->xref->dbc->prepare($sql);
  $sth->execute;
  $sth->finish;

  my $presedence;
  my $ignore; 
  if( $self->mapper->can("transcript_display_xref_sources") ){
    ($presedence, $ignore) = @{$self->mapper->transcript_display_xref_sources(1)}; # FULL update mode pass 1
  }
  else{
    ($presedence, $ignore) = @{$self->transcript_display_xref_sources(1)}; # FULL update mode pass 1
  }
#  my ($presedence, $ignore) = @{$self->transcript_display_xref_sources(1)};  # FULL update mode pass 1
  my $i=0;
  
  my $ins_p_sth = $self->xref->dbc->prepare("INSERT into display_xref_prioritys (source_id, priority) values(?, ?)");
  my $get_source_id_sth = $self->xref->dbc->prepare("select source_id from source where name like ? order by priority desc");

#
# So the higher the number the better then 
#


  my $last_name = "";
  print "Presedence for the display xrefs\n" if($self->verbose);
  foreach my $name (reverse (@$presedence)){
    $i++;
    $get_source_id_sth->execute($name);
    my $source_id;
    $get_source_id_sth->bind_columns(\$source_id);
    while($get_source_id_sth->fetch){
      $ins_p_sth->execute($source_id, $i);
      if($name ne $last_name){
	print "\t$name\t$i\n" if ($self->verbose);
      }	
      $last_name = $name;
    }
  }
  $ins_p_sth->finish;
  $get_source_id_sth->finish;


#
# Set status to 'NO_DISPLAY' for those that match the ignore REGEXP in object_xref
# Xrefs have already been dump to core etc so no damage done.
#

  my $update_ignore_sth = $self->xref->dbc->prepare('UPDATE object_xref SET ox_status = "NO_DISPLAY" where object_xref_id = ?');

  foreach my $ignore_sql (values %$ignore){
    print "IGNORE SQL: $ignore_sql\n" if($self->verbose);
    my $ignore_sth = $self->xref->dbc->prepare($ignore_sql);

    my $gene_count = 0;
    $ignore_sth->execute();
    my ($object_xref_id); 
    $ignore_sth->bind_columns(\$object_xref_id);
    while($ignore_sth->fetch()){    
      $update_ignore_sth->execute($object_xref_id);
    }
    $ignore_sth->finish;
  }
  $update_ignore_sth->finish;

#
# Do a similar thing for those with a display_label that is just numeric;
#

  $update_ignore_sth = $self->xref->dbc->prepare('UPDATE object_xref ox, source s, xref x SET ox_status = "NO_DISPLAY" where ox_status like "DUMP_OUT" and s.source_id = x.source_id and x.label REGEXP "^[0-9]+$" and ox.xref_id = x.xref_id');

  $update_ignore_sth->execute();
  $update_ignore_sth->finish;


#######################################################################

my $display_xref_sql =(<<DXS);
select  IF (ox.ensembl_object_type = 'Gene',        gtt_gene.gene_id,
        IF (ox.ensembl_object_type = 'Transcript',  gtt_transcript.gene_id,
          gtt_translation.gene_id)) AS gene_id,
        IF (ox.ensembl_object_type = 'Gene',        gtt_gene.transcript_id,
        IF (ox.ensembl_object_type = 'Transcript',  gtt_transcript.transcript_id,
          gtt_translation.transcript_id)) AS transcript_id,
        p.priority as priority,
        x.xref_id, 
        ox.ensembl_object_type as object_type,
        x.label  as label
from    (   display_xref_prioritys p
    join  (   source s
      join    (   xref x
        join      (   object_xref ox
          join        (   identity_xref ix
                      ) using (object_xref_id)
                  ) using (xref_id)
              ) using (source_id)
          ) using (source_id)
        )
  left join gene_transcript_translation gtt_gene
    on (gtt_gene.gene_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_transcript
    on (gtt_transcript.transcript_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_translation
    on (gtt_translation.translation_id = ox.ensembl_id)
where   ox.ox_status = 'DUMP_OUT'
order by    gene_id DESC, p.priority DESC, (ix.target_identity+ix.query_identity) DESC, ox.unused_priority DESC

DXS
#  SELECT gtt.gene_id, gtt.transcript_id, p.priority, x.xref_id, ox.ensembl_object_type, x.label  
#    FROM source s, xref x, object_xref ox, identity_xref ix, gene_transcript_translation gtt, display_xref_prioritys p
#     WHERE  x.source_id = s.source_id 
#       AND x.xref_id = ox.xref_id 
#       AND ox.ox_status = "DUMP_OUT"   
#       AND (          (ox.ensembl_object_type = "Transcript" and gtt.transcript_id = ox.ensembl_id)       
#                 OR   (ox.ensembl_object_type = "Translation" and gtt.translation_id = ox.ensembl_id)
#                 OR   (ox.ensembl_object_type = "Gene" and gtt.gene_id = ox.ensembl_id)     
#           )
#       AND ox.object_xref_id = ix.object_xref_id 
#       AND p.source_id = s.source_id
#  ORDER BY gtt.gene_id DESC, p.priority DESC, (ix.target_identity+ix.query_identity) DESC
#DXS

########################################################################

  my %seen_transcript; # first time we see it is the best due to ordering :-)
                         # so either write data to database or store

  
#  my $gene_sth = $self->core->dbc->prepare("select x.display_label from gene g, xref x where g.display_xref_id = x.xref_id and g.gene_id = ?"); 
#  my $tran_sth = $self->core->dbc->prepare("select x.display_label from transcript t, xref x where t.display_xref_id = x.xref_id and t.transcript_id = ?"); 


  my $last_gene = 0;

  my $display_xref_sth = $self->xref->dbc->prepare($display_xref_sql);

  my $gene_count = 0;
  $display_xref_sth->execute();
  my ($gene_id, $transcript_id, $p, $xref_id, $type, $label);  # remove labvel after testig it is not needed
  $display_xref_sth->bind_columns(\$gene_id, \$transcript_id, \$p, \$xref_id, \$type, \$label);
  while($display_xref_sth->fetch()){
    if($gene_id != $last_gene){
      $update_gene_sth->execute($xref_id+$xref_offset, $gene_id);
      $last_gene = $gene_id;
      $gene_count++;
    } 
    if($type ne "Gene"){
      if(!defined($seen_transcript{$transcript_id})){ # not seen yet so its the best
	$update_tran_sth->execute($xref_id+$xref_offset, $transcript_id);
      }
      $seen_transcript{$transcript_id} = $xref_id+$xref_offset;
      
    }
  }
  $display_xref_sth->finish;
  $update_gene_sth->finish;
  $update_tran_sth->finish;

  #
  # reset the status to DUMP_OUT fro thise that where ignored for the display_xref;
  #

  my $reset_status_sth = $self->xref->dbc->prepare('UPDATE object_xref SET ox_status = "DUMP_OUT" where ox_status = "NO_DISPLAY"');
  $reset_status_sth->execute();
  $reset_status_sth->finish;

  $sth = $self->xref->dbc->prepare("drop table display_xref_prioritys");
  $sth->execute || die "Could not drop temp table display_xref_prioritys\n";
  $sth->finish;  


  print "Updated $gene_count display_xrefs for genes\n" if($self->verbose);
}


# Remove after sure everything is cool
sub check_label{
  my $self  = shift;
  my $id    = shift;
  my $label = shift;
  my $sth   = shift;
  my $type  = shift;

  $sth->execute($id);
  my $old_label;
  $sth->bind_columns(\$old_label);
  $sth->fetch;

  if($old_label ne $label){
    print "ERROR: $type ($id) has different display_xrefs ???  old:$old_label   new:$label\n";
  }
}



sub set_source_id_to_external_name {
    
    my $self = shift;
    my $name_to_external_name_href = shift;

    my $source_id_to_external_name_href = {};
    my $name_to_source_id_href = {};
    
    my $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
    
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute();
    my ($id, $name);
    $sth->bind_columns(\$id, \$name);
    while($sth->fetch()){
	if(defined($name_to_external_name_href->{$name})){
	    $source_id_to_external_name_href->{$id} = $name_to_external_name_href->{$name};
	    $name_to_source_id_href->{$name} = $id;
	}
	elsif($name =~ /notransfer$/){
	}
	else{
	    die "ERROR: Could not find $name in external_db table please add this too continue";
	}
    }
    
    $sth->finish;
    
    return ($source_id_to_external_name_href, $name_to_source_id_href);
}



sub set_gene_descriptions{
  my $self = shift;
  my $only_those_not_set = shift || 0;

  my $update_gene_desc_sth =  $self->core->dbc->prepare("UPDATE gene SET description = ? where gene_id = ?");

  if(!$only_those_not_set){
    my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET description = null");
    $reset_sth->execute();
    $reset_sth->finish;
  }

  my %ignore;
  if($only_those_not_set){
    print "Only setting those not already set\n";
    my $sql = "select internal_id from gene_stable_id where desc_set = 1";
    my $sql_sth = $self->xref->dbc->prepare($sql);
    $sql_sth->execute;
    my $id;
    $sql_sth->bind_columns(\$id);
    while($sql_sth->fetch){
      $ignore{$id} = 1;
    }
    $sql_sth->finish;
  }	

  ##########################################
  # Get source_id to external_disaply_name #
  ##########################################

  my %name_to_external_name;
  my $sql = "select external_db_id, db_name, db_display_name from external_db";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();
  my ($id, $name, $display_name);
  $sth->bind_columns(\$id, \$name, \$display_name);
  while($sth->fetch()){
    $name_to_external_name{$name} = $display_name;
   }
  $sth->finish;
  
  my ($source_id_to_external_name_href, $name_to_source_id_href);
  if( $self->mapper->can("set_source_id_to_external_name") ){
      ($source_id_to_external_name_href, $name_to_source_id_href) = $self->mapper->set_source_id_to_external_name (\%name_to_external_name);
  }
  else{
      ($source_id_to_external_name_href, $name_to_source_id_href) = $self->set_source_id_to_external_name (\%name_to_external_name);
  }

  my %source_id_to_external_name = %$source_id_to_external_name_href;
  my %name_to_source_id = %$name_to_source_id_href;
  
  $sql =(<<SQL); 
  CREATE TABLE gene_desc_prioritys(
    source_id INT NOT NULL,
    priority       INT NOT NULL,
    PRIMARY KEY (source_id)
  ) COLLATE=latin1_swedish_ci TYPE=InnoDB
SQL

  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute;
  $sth->finish;

  my @presedence;
  my @regexps;
  if( $self->mapper->can("gene_description_sources") ){
    @presedence = $self->mapper->gene_description_sources();
  }
  else{
    @presedence = $self->gene_description_sources();
  }

  if( $self->mapper->can("gene_description_filter_regexps") ){
    @regexps = $self->mapper->gene_description_filter_regexps();
  }
  else{
    @regexps = $self->gene_description_filter_regexps();
  }


  my $i=0;
  
  my $ins_p_sth = $self->xref->dbc->prepare("INSERT into gene_desc_prioritys (source_id, priority) values(?, ?)");
  my $get_source_id_sth = $self->xref->dbc->prepare("select source_id from source where name like ?");

#
# So the higher the number the better then 
#


  print "Presedence for Gene Descriptions\n" if($self->verbose);
  my $last_name = "";
  foreach my $name (reverse (@presedence)){
    $i++;
    $get_source_id_sth->execute($name);
    my $source_id;
    $get_source_id_sth->bind_columns(\$source_id);
    while($get_source_id_sth->fetch){
      $ins_p_sth->execute($source_id, $i);
      if($last_name ne $name){
	print "\t$name\t$i\n" if ($self->verbose);
      }
      $last_name = $name;
    }
  }
  $ins_p_sth->finish;
  $get_source_id_sth->finish;

#  print "REgular expressions to be uased:=\n";
#  foreach my $reg (@regexps){
#    print "\t".$reg."\n";
#  }

#######################################################################
my $gene_desc_sql =(<<DXS);
select  IF (ox.ensembl_object_type = 'Gene',        gtt_gene.gene_id,
        IF (ox.ensembl_object_type = 'Transcript',  gtt_transcript.gene_id,
          gtt_translation.gene_id)) AS gene_id,
        x.description AS description,
        s.source_id AS source_id,
        x.accession AS accession
from    (   gene_desc_prioritys p
    join  (   source s
      join    (   xref x
        join      (   object_xref ox
          join        (   identity_xref ix
                      ) using (object_xref_id)
                  ) using (xref_id)
              ) using (source_id)
          ) using (source_id)
        )
  left join gene_transcript_translation gtt_gene
    on (gtt_gene.gene_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_transcript
    on (gtt_transcript.transcript_id = ox.ensembl_id)
  left join gene_transcript_translation gtt_translation
    on (gtt_translation.translation_id = ox.ensembl_id)
where   ox.ox_status = 'DUMP_OUT'
order by    gene_id desc,
            p.priority desc,
            (ix.target_identity+ix.query_identity) desc
#  SELECT gtt.gene_id, x.description, s.source_id, x.accession
#    FROM source s, xref x, object_xref ox, identity_xref ix, gene_transcript_translation gtt, gene_desc_prioritys p
#     WHERE  x.source_id = s.source_id 
#       AND s.source_id = p.source_id
#       AND x.xref_id = ox.xref_id 
#       AND ox.ox_status = "DUMP_OUT"   
#       AND (          (ox.ensembl_object_type = "Transcript" and gtt.transcript_id = ox.ensembl_id)       
#                 OR   (ox.ensembl_object_type = "Translation" and gtt.translation_id = ox.ensembl_id)
#                 OR   (ox.ensembl_object_type = "Gene" and gtt.gene_id = ox.ensembl_id)     
#           )
#       AND ox.object_xref_id = ix.object_xref_id 
#  ORDER BY gtt.gene_id DESC, p.priority DESC, (ix.target_identity+ix.query_identity) DESC
DXS

########################################################################

#  my $get_desc_sql = "select description from xref where xref_id = ?";
  
  
  my $gene_sth = $self->core->dbc->prepare("select g.description from gene g where g.gene_id = ?"); 


  my $last_gene = 0;

  my %no_source_name_in_desc;
  if( $self->mapper->can("no_source_label_list") ){
    foreach my $name (@{$self->mapper->no_source_label_list()}){
      my $id = $name_to_source_id{$name};
      print "$name will not have [Source:...] info in desc\n";
      $no_source_name_in_desc{$id} = 1;
    }
  }

  my $gene_desc_sth = $self->xref->dbc->prepare($gene_desc_sql);

  $gene_desc_sth->execute();
  my ($gene_id, $desc,$source_id,$label);  # remove labvel after testig it is not needed
  $gene_desc_sth->bind_columns(\$gene_id, \$desc, \$source_id, \$label);
  
  my $gene_count = 0;
  while($gene_desc_sth->fetch()){
    #    print "$gene_id, $transcript_id, $p, $xref_id, $type, $label\n";
    
    next if(defined($ignore{$gene_id}));
    
    if($gene_id != $last_gene and defined($desc) ){
      my $filtered_description = $self->filter_by_regexp($desc, \@regexps);
      if ($filtered_description ne "") {
	if(!defined($no_source_name_in_desc{$source_id})){
	  $desc .= " [Source:".$source_id_to_external_name{$source_id}.";Acc:".$label."]";
	}
	$update_gene_desc_sth->execute($desc,$gene_id);
        $gene_count++;
	$last_gene = $gene_id;
      }
    }
  }
  $update_gene_desc_sth->finish;
  $gene_desc_sth->finish;
  print "$gene_count gene descriptions added\n";# if($self->verbose);
  


  $sth = $self->xref->dbc->prepare("drop table gene_desc_prioritys");
  $sth->execute || die "Could not drop temp table gene_desc_prioritys\n";
  $sth->finish;  
}

sub filter_by_regexp {

  my ($self, $str, $regexps) = @_;

  foreach my $regexp (@$regexps) {
    $str =~ s/$regexp//ig;
  }

  return $str;

}

sub check_desc{
  my $self  = shift;
  my $id    = shift;
  my $desc = shift;
  my $sth   = shift;
  my $type  = shift;

  $sth->execute($id);
  my $old_desc;
  $sth->bind_columns(\$old_desc);
  $sth->fetch;

  if($old_desc ne $desc){
    print "ERROR: $type ($id) has different descriptions ???  \n\told:$old_desc \n\tnew:$desc\n";
  }
}


#sub no_source_label_list{
#  my $self = shift;
#  my @list;

## this needs to be added to the species .pm file 
##  @list = qw(Uniprot/SWISSPROT Refseq_dna);
##	
#  return \@list;
#}

1;
