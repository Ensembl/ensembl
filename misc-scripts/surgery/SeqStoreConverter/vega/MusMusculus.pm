use strict;
use warnings;

use SeqStoreConverter::MusMusculus;

package SeqStoreConverter::vega::MusMusculus;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::MusMusculus);

# create_coord_systems, create_seq_regions, create assembly and clone_to_seq_region modified to ignore supercontigs.


sub create_coord_systems {
  my $self = shift;

  $self->debug("MusMusculus Specific: loading assembly data");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords = 
    (["chromosome" , $ass_def, "default_version"          ,1    ],
     ['clone'      , undef   , 'default_version'          ,3     ],
     ["contig"     , undef   , "default_version,sequence_level",4]);

  my @assembly_mappings =  ("chromosome:$ass_def|contig",
                            "clone|contig",
                            "chromosome:$ass_def|contig|clone",
			   );
  my %cs = (gene                  => 'chromosome',
            transcript            => 'chromosome',
            exon                  => 'chromosome',
            dna_align_feature     => 'contig',
            protein_align_feature => 'contig',
            marker_feature        => 'contig',
            simple_feature        => 'contig',
            repeat_feature        => 'contig',
            qtl_feature           => 'chromosome',
            misc_feature          => 'chromosome',
            prediction_transcript => 'contig',
            prediction_exon       => 'contig',
            karyotype             => 'chromosome');

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib, rank) VALUES (?,?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
  }
  $sth->finish();

  $self->debug("Building meta_coord table");
  $sth = $dbh->prepare("INSERT INTO $target.meta_coord VALUES (?, ?)");
  foreach my $val (keys %cs) {
    $sth->execute($val, $coord_system_ids{$cs{$val}});
  }
  $sth->finish();

  $self->debug("Adding assembly.mapping entries to meta table");

  $sth = $dbh->prepare("INSERT INTO $target.meta(meta_key, meta_value) " .
                       "VALUES ('assembly.mapping', ?)");

  foreach my $mapping (@assembly_mappings) {
    $sth->execute($mapping);
  }
  
  $sth->finish();

  return;
}



sub create_seq_regions {
  my $self = shift;

  $self->debug("MusMusculus Specific: creating contig, " .
               "clone, and chromosome seq_regions");

  $self->contig_to_seq_region();
  $self->clone_to_seq_region();
  $self->chromosome_to_seq_region();
}

sub create_assembly {
  my $self = shift;

  $self->debug("MusMusculus Specific: loading assembly data");

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();
}




sub clone_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  # target coord_system will have a different ID
  $target_cs_name ||= "clone";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("MusMusculus Specific:Transforming clones " .
               "into $target_cs_name seq_regions");

  #
  # We don't want to make clones out of the WGS contigs, only out of
  # the actual BACs with proper embl accessions
  #
  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id,
             CONCAT(cl.embl_acc, '.', cl.embl_version),
             MAX(ctg.embl_offset+ctg.length-1)
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id
     AND    cl.embl_acc not like 'C%'
     GROUP BY ctg.clone_id");
 
  $select_sth->execute();

  my ($clone_id, $embl_acc, $length);
  $select_sth->bind_columns(\$clone_id, \$embl_acc, \$length);

  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES(?,?,?)");

  my $tmp_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?, ?)");

  while ($select_sth->fetch()) {
    $insert_sth->execute("$embl_acc", $cs_id, $length);

    #store mapping of old -> new ids in temp table
    $tmp_insert_sth->execute($clone_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $insert_sth->finish();
  $tmp_insert_sth->finish();

  return;
}





sub copy_other_tables {
  my $self = shift;

  #xref tables
  $self->copy_tables("xref",
                     "go_xref",
                     "identity_xref",
                     "object_xref",
                     "external_db",
                     "external_synonym",
  #marker/qtl related tables
                     "map",
                     "marker",
                     "marker_synonym",
                     "qtl",
                     "qtl_synonym",
  #misc other tables
		     "supporting_feature",
		     "analysis",
		     "exon_transcript",
		     "interpro",
		     "gene_description",
		     "protein_feature",
  #vega tables
		     "gene_synonym",
		     "transcript_info",
		     "current_gene_info",
		     "current_transcript_info",
		     "author",
		     "gene_name",
		     "transcript_class",
		     "gene_remark",
		     "gene_info",
		     "evidence",
		     "transcript_remark",
		     "clone_remark",
		     "clone_info",
		     "clone_info_keyword",
		     "clone_lock",
		     "current_clone_info",
		     "keyword",
		     "job",
		     "job_status",
		     "input_id_analysis");
}

sub update_clone_info {
  my $self = shift;
  return;
}


1;
