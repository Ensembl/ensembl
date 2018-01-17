=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
          "VGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
	  "RefSeq_mRNA");

}

sub gene_description_filter_regexps {

  return ();

}

sub transcript_display_xref_sources {
  my $self     = shift;

  my @list = qw(RFAM
                miRBase
                Uniprot/SWISSPROT
                Uniprot/Varsplic
               
  );

  my %ignore;

  return [\@list,\%ignore];

}


sub gene_display_xref_sources {
  my $self     = shift;
	
  my @list = qw(VGNC
                RFAM
                miRBase
                Uniprot_gn
                EntrezGene);

  my %ignore;

  #don't use EntrezGene labels dependent on predicted RefSeqs

$ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" AND
          ox.master_xref_id = dx.master_xref_id 
IEG

  #don't use labels starting with LOC

$ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

  return [\@list,\%ignore];

}

sub remove_source_priorities {
    my $self = shift;

    my $sql = "DELETE from display_xref_priority";
    my $sth = $self->xref->dbc->prepare($sql);
    $sth->execute(); 

    $sql = "DELETE from gene_desc_priority";
    $sth = $self->xref->dbc->prepare($sql);
    $sth->execute(); 

    return;
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
  # Runs build_transcript_and_gene_display_xrefs and
  # build_meta_timestamp, and, if "-upload" is set, uses the SQL files
  # produced to update the core database.

  my ($self) = @_;

  my $status = $self->mapper->xref_latest_status();
    
  if($self->mapper->can("set_display_xrefs")){
    $self->mapper->set_display_xrefs();
  }
  else{
    $self->set_display_xrefs();
  }	
  if ($self->mapper->can("transcript_names_from_gene")) {
    $self->mapper->transcript_names_from_gene();
  } else {
    $self->transcript_names_from_gene();
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

  $self->build_meta_timestamp;

  $sth_stat = $self->xref->dbc->prepare("insert into process_status (status, date) values('gene_description_done',now())");
  $sth_stat->execute();
  $sth_stat->finish;
  

  return 1;
}

sub set_gene_descriptions_from_display_xref{
  my $self = shift;
  
  $self->set_gene_descriptions(1);
}




sub set_display_xrefs_from_stable_table{
  my $self = shift;
  print "Setting Transcript and Gene display_xrefs from xref database into core and setting the desc\n" if ($self->verbose);

  my $xref_offset = $self->get_meta_value("xref_offset");
  my $core_dbi = $self->core->dbc;
  my $xref_dbi = $self->xref->dbc;

  print "Using xref_off set of $xref_offset\n" if($self->verbose);

  my $reset_sth = $core_dbi->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $core_dbi->prepare("UPDATE transcript SET display_xref_id = null WHERE biotype NOT IN ('LRG_gene')");
  $reset_sth->execute();
  $reset_sth->finish;

  $reset_sth = $core_dbi->prepare("UPDATE gene SET description = null");
  $reset_sth->execute();
  $reset_sth->finish;


  my %name_to_external_name;
  my $sql = "select external_db_id, db_name, db_display_name from external_db";
  my $sth = $core_dbi->prepare($sql);
  $sth->execute();
  my ($id, $name, $display_name);
  $sth->bind_columns(\$id, \$name, \$display_name);
  while($sth->fetch()){
    $name_to_external_name{$name} = $display_name;
   }
  $sth->finish;

  my %source_id_to_external_name;

  $sql = 'select s.source_id, s.name from source s, xref x where x.source_id = s.source_id group by s.source_id'; # only get those of interest
  $sth = $xref_dbi->prepare($sql);
  $sth->execute();
  $sth->bind_columns(\$id, \$name);

  while($sth->fetch()){
     if(defined($name_to_external_name{$name})){
      $source_id_to_external_name{$id} = $name_to_external_name{$name};
    }
  }
  $sth->finish;


  my $update_gene_sth = $core_dbi->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_gene_desc_sth = $core_dbi->prepare("UPDATE gene g SET g.description= ? WHERE g.gene_id=?");

  my $update_tran_sth = $core_dbi->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

  my $get_gene_display_xref = $xref_dbi->prepare("SELECT gsi.internal_id, gsi.display_xref_id, x.description ,x.source_id, x.accession
                                                              FROM gene_stable_id gsi, xref x 
                                                                 WHERE gsi.display_xref_id = x.xref_id");

  my $get_tran_display_xref = $xref_dbi->prepare("SELECT gsi.internal_id, gsi.display_xref_id from transcript_stable_id gsi");

  $reset_sth = $xref_dbi->prepare("UPDATE gene_stable_id gsi SET gsi.desc_set=0");
  $reset_sth->execute();

  my $set_desc_done_sth = $xref_dbi->prepare("UPDATE gene_stable_id gsi SET gsi.desc_set=1 WHERE gsi.internal_id=?");

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

  #
  # Clean up synonyms linked to xrefs which are not the display xref
  # Synonyms are only used as alternative gene names, so should be synonyms of the gene symbol chosen
  #

  my $syn_clean_sth = $core_dbi->prepare("DELETE es FROM external_synonym es, xref x LEFT JOIN gene g ON g.display_xref_id = x.xref_id WHERE es.xref_id = x.xref_id AND isnull(g.display_xref_id)");
  $syn_clean_sth->execute();
  $syn_clean_sth->finish();

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


  print "Building Transcript and Gene display_xrefs\n" if ($self->verbose);

  my $xref_offset = $self->get_meta_value("xref_offset");
  my $core_dbi = $self->core->dbc();
  my $xref_dbi = $self->xref->dbc();

  print "Using xref_off set of $xref_offset\n" if($self->verbose);

  my $reset_sth = $core_dbi->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $core_dbi->prepare("UPDATE transcript SET display_xref_id = null WHERE biotype NOT IN ('LRG_gene')");
  $reset_sth->execute();
  $reset_sth->finish;

  my $update_gene_sth = $core_dbi->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_tran_sth = $core_dbi->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");



  # Set status to 'NO_DISPLAY' for object_xrefs with a display_label that is just numeric;

  my $update_ignore_sth = $xref_dbi->prepare('UPDATE object_xref ox, source s, xref x SET ox_status = "NO_DISPLAY" where ox_status like "DUMP_OUT" and s.source_id = x.source_id and x.label REGEXP "^[0-9]+$" and ox.xref_id = x.xref_id');

  $update_ignore_sth->execute();
  $update_ignore_sth->finish;


  my $ins_p_sth = $xref_dbi->prepare("INSERT ignore into display_xref_priority (ensembl_object_type,source_id, priority) values(?, ?, ?)");
  my $get_source_id_sth = $xref_dbi->prepare("select source_id from source where name like ? order by priority");
  my $list_sources_sth = $xref_dbi->prepare("select distinct name from display_xref_priority d join source using(source_id) where ensembl_object_type = ? order by d.priority");

  $update_ignore_sth = $xref_dbi->prepare('UPDATE object_xref SET ox_status = "NO_DISPLAY" where object_xref_id = ?');


  my %object_types = ('gene' => 'Gene', 'transcript' => 'Transcript'); 
  
  foreach my $object_type (keys %object_types) {

      my $precedence;
      my $ignore; 
      my $method = $object_type . '_display_xref_sources';
      if( $self->mapper->can($method) ){
	  ($precedence, $ignore) = @{$self->mapper->$method()};
      }
      else{
	  ($precedence, $ignore) = @{$self->$method()};
      }

      # The lower the priority number the better then 
      my $i=0;
      foreach my $name (@$precedence){
	  $i++;
	  $get_source_id_sth->execute($name);
	  my $source_id;
	  $get_source_id_sth->bind_columns(\$source_id);
	  while($get_source_id_sth->fetch){
	      $ins_p_sth->execute($object_types{$object_type},$source_id, $i);
	  }
      }
      $ins_p_sth->finish;
      $get_source_id_sth->finish;

      $i = 0;
      if ($self->verbose) {
	  print "Precedence for $object_type display xrefs (1- best name)\n";
	  $list_sources_sth->execute($object_types{$object_type});
	  my $source_name;
	  $list_sources_sth->bind_columns(\$source_name);
	  while ($list_sources_sth->fetch() ) {
	      $i++;
	      print "\t$i\t$source_name\n";
	  }
 
      }

      # Set status to 'NO_DISPLAY' for those that match the ignore REGEXP in object_xref
      # Xrefs have already been dump to core etc so no damage done.

      foreach my $ignore_sql (values %$ignore){
	  print "IGNORE SQL: $ignore_sql\n" if($self->verbose);
	  my $ignore_sth = $xref_dbi->prepare($ignore_sql);
	  $ignore_sth->execute();
	  my ($object_xref_id); 
	  $ignore_sth->bind_columns(\$object_xref_id);
	  while($ignore_sth->fetch()){    
	      $update_ignore_sth->execute($object_xref_id);
	  }
	  $ignore_sth->finish;
      }
      $update_ignore_sth->finish;


#look at sources of display xrefs which are relevant for this object type
#(listed in gene_display_xref_sources() or transcript_display_xref_sources() )
#but get xrefs for all levels Gene, its Transcripts and Translations 
#######################################################################
my $display_xref_sql =(<<DXS);
select  CASE ox.ensembl_object_type
           WHEN 'Gene' THEN gtt_gene.gene_id
	   WHEN 'Transcript' THEN gtt_transcript.gene_id
	   WHEN 'Translation' THEN gtt_translation.gene_id
	END AS d_gene_id,
        CASE ox.ensembl_object_type
           WHEN 'Gene' THEN gtt_gene.transcript_id
	   WHEN 'Transcript' THEN gtt_transcript.transcript_id
	   WHEN 'Translation' THEN gtt_translation.transcript_id
        END AS d_transcript_id,
        p.priority as priority,
        x.xref_id
from    (   display_xref_priority p
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
        and p.ensembl_object_type = ?
order by d_gene_id, ox.ensembl_object_type, 
	p.priority, (ix.target_identity + ix.query_identity) DESC, unused_priority DESC;

DXS


########################################################################

      my %object_seen;
 
      my $display_xref_sth = $xref_dbi->prepare($display_xref_sql);

      my $display_xref_count = 0;
      $display_xref_sth->execute($object_type);
      my ($gene_id, $transcript_id, $priority, $xref_id);  
      $display_xref_sth->bind_columns(\$gene_id, \$transcript_id, \$priority, \$xref_id);
      while($display_xref_sth->fetch()){
	  my $object_id;
	  if ($object_type eq 'gene') {
	      $object_id = $gene_id;
	  } elsif ($object_type eq 'transcript') {
	      $object_id = $transcript_id;
	  }

	  if (!exists($object_seen{$object_id}) ) {
	      if ($object_type eq 'gene') {
		  $update_gene_sth->execute($xref_id+$xref_offset, $object_id);
	      } elsif ($object_type eq 'transcript') {
		  $update_tran_sth->execute($xref_id+$xref_offset, $object_id);
	      }    
	      $display_xref_count++;
	      $object_seen{$object_id} = 1;
	  }
      } 
  
      $display_xref_sth->finish;
      $update_gene_sth->finish;
      $update_tran_sth->finish;

      print "Updated $display_xref_count $object_type display_xrefs\n" if($self->verbose);

  }

  #
  # reset the status to DUMP_OUT fro object_xrefs that where ignored for the display_xref;
  #

  my $reset_status_sth = $xref_dbi->prepare('UPDATE object_xref SET ox_status = "DUMP_OUT" where ox_status = "NO_DISPLAY"');
  $reset_status_sth->execute();
  $reset_status_sth->finish;

  #
  # Clean up synonyms linked to xrefs which are not the display xref
  # Synonyms are only used as alternative gene names, so should be synonyms of the gene symbol chosen
  #

  my $syn_clean_sth = $core_dbi->prepare("DELETE es FROM external_synonym es, xref x LEFT JOIN gene g ON g.display_xref_id = x.xref_id WHERE es.xref_id = x.xref_id AND isnull(g.display_xref_id)");
  $syn_clean_sth->execute();
  $syn_clean_sth->finish();


}


sub transcript_names_from_gene {
  my $self = shift;
  my $core_dbi = $self->core->dbc;
  my $xref_dbi = $self->xref->dbc;

  print "Assigning transcript names from gene names\n" if ($self->verbose);

  my $reset_sth = $core_dbi->prepare("UPDATE transcript SET display_xref_id = null WHERE biotype NOT IN ('LRG_gene')");
  $reset_sth->execute();
  $reset_sth->finish;

  my $xref_id_sth = $core_dbi->prepare("SELECT max(xref_id) FROM xref");
  my $ox_id_sth = $core_dbi->prepare("SELECT max(object_xref_id) FROM object_xref");
  my $del_xref_sth = $core_dbi->prepare("DELETE x FROM xref x, object_xref ox WHERE x.xref_id = ox.xref_id AND ensembl_object_type = 'Transcript' AND display_label REGEXP '-2[0-9]{2}\$'");
  my $reuse_xref_sth = $core_dbi->prepare("SELECT xref_id FROM xref x WHERE external_db_id = ? AND display_label = ? AND description = ? AND info_type = 'MISC'");
  my $del_ox_sth = $core_dbi->prepare("DELETE ox FROM object_xref ox LEFT JOIN xref x ON x.xref_id = ox.xref_id WHERE isnull(x.xref_id)");
  my $ins_xref_sth = $core_dbi->prepare("INSERT IGNORE into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values(?, ?, ?, ?, 0, ?, 'MISC', ?)");
  my $ins_ox_sth = $core_dbi->prepare("INSERT into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id) values(?, ?, 'Transcript', ?)");
  my $update_tran_sth = $core_dbi->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

  my $get_genes = $core_dbi->prepare("SELECT g.gene_id, e.db_name, x.dbprimary_acc, x.display_label, x.description FROM gene g, xref x, external_db e where g.display_xref_id = x.xref_id and e.external_db_id = x.external_db_id");
  my $get_transcripts = $core_dbi->prepare("SELECT transcript_id FROM transcript WHERE gene_id = ? ORDER BY seq_region_start, seq_region_end");
  my $get_source_id = $core_dbi->prepare("SELECT external_db_id FROM external_db WHERE db_name like ?");

  $get_genes->execute();
  my ($gene_id, $external_db, $external_db_id, $acc, $label, $description, $transcript_id, $xref_id, $ox_id, $ext, $reuse_xref_id, $info_text);
  $get_genes->bind_columns(\$gene_id, \$external_db, \$acc, \$label, \$description);
  $xref_id_sth->execute();
  $xref_id_sth->bind_columns(\$xref_id);
  $xref_id_sth->fetch();
  $ox_id_sth->execute();
  $ox_id_sth->bind_columns(\$ox_id);
  $ox_id_sth->fetch();
  $del_xref_sth->execute();
  while ($get_genes->fetch()) {
    $ext = '201';
    $get_source_id->execute($external_db . "_trans_name");
    $get_source_id->bind_columns(\$external_db_id);
    $get_source_id->fetch();
    $get_transcripts->execute($gene_id);
    $get_transcripts->bind_columns(\$transcript_id);
    while ($get_transcripts->fetch) {
      $xref_id++;
      $ox_id++;
      $reuse_xref_sth->execute($external_db_id, $label . '-' . $ext, $description);
      $reuse_xref_sth->bind_columns(\$reuse_xref_id);
      if ($reuse_xref_sth->fetch()) {
        $ins_ox_sth->execute($ox_id, $transcript_id, $reuse_xref_id);
        $update_tran_sth->execute($reuse_xref_id, $transcript_id);
      } else {
        $info_text = 'via gene ' . $acc;
        $ins_xref_sth->execute($xref_id, $external_db_id, $label. "-" . $ext, $label . "-" . $ext, $description, $info_text);
        $ins_ox_sth->execute($ox_id, $transcript_id, $xref_id); 
        $update_tran_sth->execute($xref_id, $transcript_id);
      }
      $ext++;
    }
  }

  $del_xref_sth->finish();
  $del_ox_sth->execute();
  $del_ox_sth->finish();
  $reuse_xref_sth->finish();
  $xref_id_sth->finish();
  $ox_id_sth->finish();
  $get_genes->finish();
  $get_source_id->finish();
  $get_transcripts->finish();
  $ins_xref_sth->finish();
  $ins_ox_sth->finish();
  $update_tran_sth->finish();
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
  my $sql;
  my $core_dbi = $self->core->dbc;
  my $xref_dbi = $self->xref->dbc;

  my $update_gene_desc_sth =  $core_dbi->prepare("UPDATE gene SET description = ? where gene_id = ?");

  if(!$only_those_not_set){
    my $reset_sth = $core_dbi->prepare("UPDATE gene SET description = null");
    $reset_sth->execute();
    $reset_sth->finish;
  }

  my %ignore;
  if($only_those_not_set){
    print "Only setting those not already set\n";
    $sql = "select internal_id from gene_stable_id where desc_set = 1";
    my $sql_sth = $xref_dbi->prepare($sql);
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
  $sql = "select external_db_id, db_name, db_display_name from external_db";
  my $sth = $core_dbi->prepare($sql);
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
  

  my @precedence;
  my @regexps;
  if( $self->mapper->can("gene_description_sources") ){
    @precedence = $self->mapper->gene_description_sources();
  }
  else{
    @precedence = $self->gene_description_sources();
  }

  if( $self->mapper->can("gene_description_filter_regexps") ){
    @regexps = $self->mapper->gene_description_filter_regexps();
  }
  else{
    @regexps = $self->gene_description_filter_regexps();
  }
  
  my $ins_p_sth = $xref_dbi->prepare("INSERT ignore into gene_desc_priority (source_id, priority) values(?, ?)");
  my $get_source_id_sth = $xref_dbi->prepare("select source_id from source where name like ?");
  my $list_sources_sth = $xref_dbi->prepare("select distinct name from gene_desc_priority d join source using(source_id) order by d.priority");

  # The lower the priority number the better then 
  my $i=0;
  foreach my $name (@precedence){
     $i++;
     $get_source_id_sth->execute($name);
     my $source_id;
     $get_source_id_sth->bind_columns(\$source_id);
     while($get_source_id_sth->fetch){
	  $ins_p_sth->execute($source_id, $i);
     }
  }
  $ins_p_sth->finish;
  $get_source_id_sth->finish;


  $i = 0;
  if ($self->verbose) {
      print "Precedence for gene descriptions (1- best description)\n";
      $list_sources_sth->execute();
      my $source_name;
      $list_sources_sth->bind_columns(\$source_name);
      while ($list_sources_sth->fetch() ) {
	  $i++;
	  print "\t$i\t$source_name\n";
      }
 
  }



#######################################################################
my $gene_desc_sql =(<<DXS);
select  CASE ox.ensembl_object_type
           WHEN 'Gene' THEN gtt_gene.gene_id
	   WHEN 'Transcript' THEN gtt_transcript.gene_id
	   WHEN 'Translation' THEN gtt_translation.gene_id
	END AS d_gene_id,
        x.description AS description,
        s.source_id AS source_id,
        x.accession AS accession
from    (   gene_desc_priority p
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
order by    d_gene_id,
            ox.ensembl_object_type,
            p.priority,
            (ix.target_identity+ix.query_identity) desc

DXS

######################################################################## 
  
  my $gene_sth = $core_dbi->prepare("select g.description from gene g where g.gene_id = ?"); 

  my %no_source_name_in_desc;
  if( $self->mapper->can("no_source_label_list") ){
    foreach my $name (@{$self->mapper->no_source_label_list()}){
      my $id = $name_to_source_id{$name};
      print "$name will not have [Source:...] info in desc\n";
      $no_source_name_in_desc{$id} = 1;
    }
  }

  my $gene_desc_sth = $xref_dbi->prepare($gene_desc_sql);

  $gene_desc_sth->execute();
  my ($gene_id, $desc,$source_id,$label);  
  $gene_desc_sth->bind_columns(\$gene_id, \$desc, \$source_id,\$label);
  
  my %gene_desc_updated;

  while($gene_desc_sth->fetch()){
     
    next if(exists($ignore{$gene_id}) || exists($gene_desc_updated{$gene_id}));
    
    if(defined($desc) ){
      my $filtered_desc = $self->filter_by_regexp($desc, \@regexps);
      if ($filtered_desc ne "") {
	if(!defined($no_source_name_in_desc{$source_id})){
	  $filtered_desc .= " [Source:".$source_id_to_external_name{$source_id}.";Acc:".$label."]";
	}
	$update_gene_desc_sth->execute($filtered_desc,$gene_id);
	$gene_desc_updated{$gene_id} = 1;
      }
    }
  }
  $update_gene_desc_sth->finish;
  $gene_desc_sth->finish;
  print scalar(keys %gene_desc_updated) ." gene descriptions added\n";
  

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


1;
