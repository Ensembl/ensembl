package XrefMapper::drosophila;
use strict;

use  XrefMapper::BasicMapper;
use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

use XrefMapper::BasicMapper qw(%stable_id_to_internal_id %object_xref_mappings %xref_to_source %xref_accessions %source_to_external_db);

my %genes_to_transcripts;
my %transcript_to_translation;
my %translation_to_transcript;
my %transcript_length;

sub gene_description_filter_regexps {

  return ();

}

# Special logic for drosophila display_xrefs:
#
# gene: flybase_name if present, else gadfly_gene_cgid
#
# transcript: flybase_name if present, else gadfly_transcript_cgid
sub xref_offset{
  my ($self, $val) = @_;

  if(defined($val)){
    $self->{'_xref_offset'} = $val;
  }
  return $self->{'_xref_offset'};
}


sub gene_display_xref_sources {

  my @list = qw(
								FlyBaseName_gene
                FlyBaseCGID_gene
                flybase_gene_id
                );

  my %ignore;
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}


sub build_transcript_and_gene_display_xrefs {
  my ($self) = @_;
  my $dir = $self->core->dir();

  my %external_name_to_id;
  my %ex_db_id_to_status;
  my $sql1 = "SELECT external_db_id, db_name, status from external_db";
  
  my $sth1 = $self->core->dbc->prepare($sql1) || die "prepare failed for $sql1\n";
  $sth1->execute() || die "execute failed";
  my ($db_id, $name, $status);
  $sth1->bind_columns(\$db_id, \$name, \$status);
  while($sth1->fetch()){
    $external_name_to_id{$name}  = $db_id;
    $ex_db_id_to_status{$db_id} = $status;
  }
  $sth1->finish;


  #############################
  #create the tempory table
  #############################

  my $sth = $self->core->dbc->prepare("create table identity_xref_temp like identity_xref");
  print "creating table identity_xref_temp\n";
  $sth->execute() || die "Could not \ncreate table identity_xref_temp like identity_xref\n";


  #############################
  #populate the tempory table
  #############################
  my $file = $dir."/identity_xref_temp.txt";
  
  if(-s $file){
    my $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$file\' IGNORE INTO TABLE identity_xref_temp");
    print "Uploading data in $file to identity_xref_temp\n";
    $sth->execute();
  }
  else{
    print "NO file or zero size file, so not able to load file $file to identity_xref_temp\n";
  }

  #
  # get a list of sources to use
  # and also a list of those xrefs to ignore 
  # where the source name is the key and the value is the string to test for 
  # 
  my ($genepresedence, $geneignore) = @{$self->gene_display_xref_sources()};
  my ($presedence, $ignore) = @{$self->transcript_display_xref_sources()};

  my $i=0;
  my %level;

  foreach my $ord (reverse (@$presedence)){
    $i++;
    if(!defined($external_name_to_id{$ord})){
      print STDERR "unknown external database name *$ord* being used\n";
    }
    $level{$external_name_to_id{$ord}} = $i;
  }

  foreach my $ord (reverse (@$genepresedence)){
    $i++;
    if(!defined($external_name_to_id{$ord})){
      print STDERR "unknown external database name *$ord* being used\n";
    }
    $level{$external_name_to_id{$ord}} = $i;
  }

  if(!scalar(keys %genes_to_transcripts)){
    $self->build_genes_to_transcripts();
  }

  if(!scalar(keys %translation_to_transcript)){
    $self->load_translation_to_transcript();
  }



  my $sql = (<<ESQL);
  SELECT ox.xref_id, ix.query_identity, ix.target_identity,  x.external_db_id, x.display_label, e.db_name, ox.linkage_annotation
    FROM (object_xref ox, xref x, external_db e) 
      LEFT JOIN identity_xref ix ON (ox.object_xref_id = ix.object_xref_id) 
	WHERE x.xref_id = ox.xref_id AND ox.ensembl_object_type = ? 
              AND ox.ensembl_id = ? AND x.info_type = 'SEQUENCE_MATCH'
              AND e.external_db_id = x.external_db_id
ESQL

  my $primary_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";



  $sql = (<<ZSQL);
  SELECT ox.xref_id, ix.query_identity, ix.target_identity, x.external_db_id, x.display_label, e.db_name, ox.linkage_annotation
    FROM (object_xref ox, xref x, external_db e) 
      LEFT JOIN identity_xref_temp ix ON (ox.object_xref_id = ix.object_xref_id) 
	WHERE x.xref_id = ox.xref_id and ox.ensembl_object_type = ? 
              and ox.ensembl_id = ? and x.info_type = 'DEPENDENT'
              AND e.external_db_id = x.external_db_id
ZSQL
 
  my $dependent_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";



  $sql = (<<QSQL);
  SELECT  x.xref_id, x.external_db_id, x.display_label, e.db_name, o.linkage_annotation
   FROM object_xref o, xref x, external_db e 
    WHERE x.xref_id = o.xref_id 
        and o.ensembl_object_type = ? and o.ensembl_id = ? and x.info_type = 'DIRECT'
              AND e.external_db_id = x.external_db_id
QSQL
                             
  my $direct_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";



# get xrefs connect directly to the gene.

  $sql = (<<GSQL);
  SELECT x.xref_id, x.external_db_id, e.db_name, o.linkage_annotation
   FROM object_xref o, xref x, external_db e   
    WHERE x.xref_id = o.xref_id 
        and o.ensembl_object_type = 'Gene' and o.ensembl_id = ?
              AND e.external_db_id = x.external_db_id
GSQL
                             
  my $gene_sth = $self->core->dbc->prepare($sql) || die "prepare failed for $sql\n";

  my $count =0;
  
  my ($xref_id, $qid, $tid, $ex_db_id, $display_label, $external_db_name, $linkage_annotation);
  
  

  # Open file handles to recieve SQL and text data used to set 
  # display_xrefs 
  my $gene_dx_file       = "$dir/gene_display_xref.sql"; 
  my $tran_dx_file       = "$dir/transcript_display_xref.sql";
  my $unset_gene_dx_file = "$dir/gene_unset_display_xref.sql";
  my $unset_tran_dx_file = "$dir/transcript_unset_display_xref.sql";

  open (my $GENE_DX, ">", $gene_dx_file)
      or die( "Could not open $gene_dx_file: $!" );
  open (my $TRANSCRIPT_DX, ">", $tran_dx_file) 
      or die( "Could not open $tran_dx_file: $!" );
  open (my $GENE_DX_UNSET, ">", $unset_gene_dx_file)
      or die( "Could not open $unset_gene_dx_file: $!" );
  open (my $TRAN_DX_UNSET, ">", $unset_tran_dx_file) 
      or die( "Could not open $unset_tran_dx_file: $!" );
  open (my $GENE_DX_TXT, ">", "$dir/gene_display_xref.txt");
  open (my $TRANSCRIPT_DX_TXT, ">", "$dir/transcript_display_xref.txt");

  # These are the files that this method will return
  my @files = ($unset_gene_dx_file,$gene_dx_file, 
               $unset_tran_dx_file,$tran_dx_file);

  # Write the 'unset' sql to the files, and cose them
  print $GENE_DX_UNSET qq(UPDATE gene       SET display_xref_id=NULL;\n);
  print $TRAN_DX_UNSET qq(UPDATE transcript SET display_xref_id=NULL;\n);
  close( $GENE_DX_UNSET );
  close( $TRAN_DX_UNSET );
  
  foreach my $gene_id (keys %genes_to_transcripts) {
    my %percent_id;
    my %level_db;
    my %parent;
    my %percent_id_via_acc;
    my @gene_xrefs = ();
    
    $gene_sth->execute($gene_id) || die "execute failed";
    $gene_sth->bind_columns(\$xref_id, \$ex_db_id, \$external_db_name, \$linkage_annotation);
    
    
    my $best_gene_xref  = 0;    # store xref
    my $best_gene_level = 0;    # store level
    my $best_gene_percent = 0;  # additoon of precentage ids

    while($gene_sth->fetch()){
      if(defined($$ignore{$external_db_name})){
	if($linkage_annotation =~ /$$ignore{$external_db_name}/){
#	  print "Ignoring $xref_id as linkage_annotation has ".$$ignore{$external_db_name}." in it. DELETE THIS MESSAGE AFTER TESTING\n";
	  next;
	}
      }
      if($level{$ex_db_id} > $best_gene_level){
	$best_gene_xref = $xref_id;
	$best_gene_level = $level{$ex_db_id};
      }
      if($best_gene_xref){
        print $GENE_DX "UPDATE gene g SET g.display_xref_id=" . $best_gene_xref .
          " WHERE g.gene_id=" . $gene_id . ";\n";
        print $GENE_DX_TXT $best_gene_xref . "\t" . $gene_id ."\n";
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
	    if(defined($$ignore{$external_db_name})){
	      if($linkage_annotation =~ /$$ignore{$external_db_name}/){
#		print "Ignoring $xref_id as linkage_annotation has ".$$ignore{$external_db_name}." in it. DELETE THIS MESSAGE AFTER TESTING\n";
		next;
	      }
	    }
	    push @transcript_xrefs, $xref_id;
	    if(!defined($qid) || !defined($tid)){
	      print "DEPENDENT $xref_id\n" if($ex_db_id != 1100); #HGNC has added one with no %ids.
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

      foreach my $xref_id (@transcript_xrefs) {
	if(defined($level_db{$xref_id}) and $level_db{$xref_id}){
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
	}
      }       
      
      if($best_tran_xref){
        print $TRANSCRIPT_DX "UPDATE transcript SET display_xref_id=" .$best_tran_xref. 
            " WHERE transcript_id=" . $transcript_id . ";\n";
        print $TRANSCRIPT_DX_TXT  $best_tran_xref. "\t" . $transcript_id . "\n";
      }

      if($best_tran_level < $best_gene_level){
         next;
      }
      if($best_tran_level == $best_gene_level){
        if($best_tran_percent < $best_gene_percent){
          next;
        }
      }

    }
  
  }
  close $TRANSCRIPT_DX;
  close $TRANSCRIPT_DX_TXT;
  close $GENE_DX;
  close $GENE_DX_TXT;
  

  return @files;
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

sub gene_description_sources {
  return (
          "FlyBaseName_gene",
#          "gadfly_gene_cgid",
          "FlyBaseCGID_gene",
         );
}

sub transcript_display_xref_sources {

  my @list = qw(FlyBaseName_transcript FlyBaseCGID_transcript);

 # gadfly_transcript_cgid flybase_annotation_id

  my %ignore;
  $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';

  return [\@list,\%ignore];

}

1;
