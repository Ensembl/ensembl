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

sub gene_description_sources {
  return (
          "FlyBaseName_gene",
          "FlyBaseCGID_gene",
         );
}

sub transcript_display_xref_sources {
	  my $self     = shift;
		my $fullmode = shift;
	
		my @list = qw(FlyBaseName_transcript 
                  FlyBaseCGID_transcript);

 # gadfly_transcript_cgid flybase_annotation_id

		print "Get FlyBase Transcript display_xrefs_sources with fullmode=$fullmode\n" if ($self->verbose);

  my %ignore;

  return [\@list,\%ignore];

}

sub gene_display_xref_sources {
		my $self     = shift;
		my $fullmode = shift;

  my @list = qw(
								FlyBaseName_gene
                FlyBaseCGID_gene
                flybase_gene_id
                );

  my %ignore;
  
		print "Get FlyBase Gene display_xrefs_sources with fullmode=$fullmode\n" if ($self->verbose);

  # Both methods

  if(!$fullmode){
    $ignore{"EntrezGene"}= 'FROM:RefSeq_[pd][en][pa].*_predicted';
  }
  else{
    $ignore{"EntrezGene"} = 'select ox.object_xref_id from object_xref ox, dependent_xref dx, source s1, xref x1, source s2, xref x2 where ox.object_xref_id = dx.object_xref_id and dx.dependent_xref_id = x1.xref_id and x1.source_id = s1.source_id and s1.name = "EntrezGene" and x2.xref_id = dx.master_xref_id and x2.source_id = s2.source_id and (s2.name like "Refseq_dna_predicted" or s2.name like "RefSeq_peptide_predicted") and ox.ox_status = "DUMP_OUT"';

  }


  return [\@list,\%ignore];

}


sub set_display_xrefs{
  my $self = shift;

  print "Building FlyBase Transcript and Gene display_xrefs using xref database\n" if ($self->verbose);

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

	##############################################
	# Careful: The table might have been created #
	# already if the pipeline was run previously #
	# delete the table first (if not here)       #
	##############################################

	my $sth = $self->xref->dbc->prepare("DROP TABLE IF EXISTS display_xref_prioritys");
  $sth->execute;
  $sth->finish;

my $sql =(<<SQL); 
  CREATE TABLE display_xref_prioritys(
    source_id INT NOT NULL,
    priority       INT NOT NULL,
    PRIMARY KEY (source_id)
  ) COLLATE=latin1_swedish_ci TYPE=InnoDB
SQL

  $sth = $self->xref->dbc->prepare($sql);
  $sth->execute;
  $sth->finish;

	############################################
	# We also create an extra table for the    #
	# gene priority                            #
	############################################

	my $sth = $self->xref->dbc->prepare("DROP TABLE IF EXISTS gene_display_xref_prioritys");
  $sth->execute;
  $sth->finish;

my $g_sql =(<<GSQL); 
  CREATE TABLE gene_display_xref_prioritys(
    source_id INT NOT NULL,
    priority       INT NOT NULL,
    PRIMARY KEY (source_id)
  ) COLLATE=latin1_swedish_ci TYPE=InnoDB
GSQL

  $sth = $self->xref->dbc->prepare($g_sql);
  $sth->execute;
  $sth->finish;

  ############################################
	# OK comes the specific code for FlyBase:  #
	# Separate Transcript display xrefs from   #
	# gene display xrefs.                      #
	# get a list of sources to use.            #
	# Gene and Trans. have different levels:   #
  # transcript presedence should never be    #
	# mixed with gene presedence               #
	# for this reason, we have to create a     #
	# distinct gene_level / transcript_level   #
	############################################

  # in we are here, it means that fullmode = 1;
	my $fullmode = 1;
  my ($gene_presedence, $gene_ignore) = @{$self->gene_display_xref_sources($fullmode)};
  my ($transcript_presedence, $transcript_ignore) = @{$self->transcript_display_xref_sources($fullmode)};

  my $i=0;
  my $j=0;

  my $ins_p_sth = $self->xref->dbc->prepare("INSERT into display_xref_prioritys (source_id, priority) values(?, ?)");
	my $ins_g_p_sth = $self->xref->dbc->prepare("INSERT into gene_display_xref_prioritys (source_id, priority) values(?, ?)");
  my $get_source_id_sth = $self->xref->dbc->prepare("select source_id from source where name like ? order by priority desc");

  ############################################
  # So the higher the number the better then #
	# Do it for transcripts and then for genes #
  ############################################

  my $last_name = "";
  print "Transcript presedence for the display xrefs\n" if($self->verbose);
  foreach my $name (reverse (@$transcript_presedence)){
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

	$last_name = "";
  print "Gene presedence for the display xrefs\n" if($self->verbose);
  foreach my $name (reverse (@$gene_presedence)){
			$j++;
			$get_source_id_sth->execute($name);
			my $source_id;
			$get_source_id_sth->bind_columns(\$source_id);
			while($get_source_id_sth->fetch){
					$ins_g_p_sth->execute($source_id, $j);
					if($name ne $last_name){
							print "\t$name\t$j\n" if ($self->verbose);
					}	
					$last_name = $name;
			}
  }
  $ins_p_sth->finish;

  $get_source_id_sth->finish;

  ############################################
  # Set status to 'NO_DISPLAY' for those that#
  # match the ignore REGEXP in object_xref   #
  # Xrefs have already been dump to core etc #
	# so no damage done.                       #
  ############################################

  my $update_ignore_sth = $self->xref->dbc->prepare('UPDATE object_xref SET ox_status = "NO_DISPLAY" where object_xref_id = ?');

  ############################################
  # Gene and transcript ignore are the same  #
  ############################################
	
	foreach my $ignore_sql (values %$transcript_ignore){
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

my $gene_display_xref_sql =(<<G_DXS);
select  gtt_gene.gene_id AS gene_id,
        p.priority as priority,
        x.xref_id, 
        ox.ensembl_object_type as object_type,
        x.label  as label
from    (   gene_display_xref_prioritys p
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
where   ox.ox_status = 'DUMP_OUT' 
and     ox.ensembl_object_type = 'Gene' 
order by    gene_id DESC, p.priority DESC, (ix.target_identity+ix.query_identity) DESC, ox.unused_priority DESC

G_DXS
  

  my $gene_display_xref_sth = $self->xref->dbc->prepare($gene_display_xref_sql);

	my $gene_count = 0;
	my $last_gene = 0;
	my $transcript_count = 0;

	$gene_display_xref_sth->execute();

	my ($gene_id, $transcript_id, $p, $xref_id, $type, $label);  # remove labvel after testig it is not needed
  $gene_display_xref_sth->bind_columns(\$gene_id, \$p, \$xref_id, \$type, \$label);
  while($gene_display_xref_sth->fetch()){
			if($gene_id != $last_gene){
					$update_gene_sth->execute($xref_id+$xref_offset, $gene_id);
					$last_gene = $gene_id;
					$gene_count++;
			} 
  }

  $gene_display_xref_sth->finish;

  print "Updated $gene_count display_xrefs for FlyBase genes\n" if($self->verbose);

  my %seen_transcript; # first time we see it is the best due to ordering :-)
                         # so either write data to database or store

  $last_gene = 0;

  my $display_xref_sth = $self->xref->dbc->prepare($display_xref_sql);


  $display_xref_sth->execute();
  $display_xref_sth->bind_columns(\$gene_id, \$transcript_id, \$p, \$xref_id, \$type, \$label);
  while($display_xref_sth->fetch()){
			if($gene_id != $last_gene){
					$last_gene = $gene_id;
			} 
			if($type ne "Gene"){
					if(!defined($seen_transcript{$transcript_id})){ # not seen yet so its the best
							$update_tran_sth->execute($xref_id+$xref_offset, $transcript_id);
					}
					$seen_transcript{$transcript_id} = $xref_id+$xref_offset;
					$transcript_count++;
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
	
	# also drop the gene_display_xref_prioritys
	$sth = $self->xref->dbc->prepare("drop table gene_display_xref_prioritys");
  $sth->execute || die "Could not drop temp table display_xref_prioritys\n";
  $sth->finish;  

	print "Updated $transcript_count display_xrefs for FlyBase transcripts\n" if($self->verbose);
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

1;
