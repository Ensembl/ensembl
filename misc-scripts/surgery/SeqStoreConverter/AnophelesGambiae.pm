use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::AnophelesGambiae;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);


sub create_coord_systems {
  my $self = shift;

  $self->debug("AnophelesGambiae Specific: creating scaffold, chunk and, " .
	       "chromosome coord systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords =
    (['chunk',         undef, 'default_version,sequence_level'],
     ['chromosome', $ass_def, 'top_level,default_version'],
     ["scaffold" ,     undef, "default_version,sequence_level"]);

  my @assembly_mappings = ("chromosome:$ass_def|chunk",
                           "scaffold|chunk");

  my %cs = (gene                  => 'chromosome',
            transcript            => 'chromosome',
            exon                  => 'chromosome',
            dna_align_feature     => 'scaffold',
            protein_align_feature => 'scaffold',
            marker_feature        => 'scaffold',
            simple_feature        => 'scaffold',
            repeat_feature        => 'scaffold',
            qtl_feature           => 'chromosome',
            misc_feature          => 'chromosome',
            prediction_transcript => 'chromosome',
            karyotype             => 'chromosome');

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib) VALUES (?,?,?)");

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

  $self->debug("AnophelesGambiae Specific: creating seq_regions");

  $self->contig_to_seq_region('chunk');
  $self->clone_to_seq_region('scaffold');
  $self->chromosome_to_seq_region();
}



sub create_assembly {
  my $self = shift;

  $self->debug("AnophelesGambiae Specific: loading assembly table");

  $self->assembly_contig_chromosome();
  $self->assembly_contig_clone();

  return;
}


sub transfer_prediction_transcripts {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();
  
  $self->debug("AnophelesGambiae Specific: building prediction_exon table");

  #
  # In Anopheles the predicion transcripts were computed in chromosomal
  # coords, so convert them to chromosomal coords and merge any adjacent
  # exons
  #

  my $sql = 
    "SELECT pt.prediction_transcript_id, tcm.new_id as seq_region_id, " .
    "       IF(a.contig_ori=1,(pt.contig_start+a.chr_start-a.contig_start),".
    "         (a.chr_start+a.contig_end-pt.contig_end)) as start, " .
    "       IF(a.contig_ori=1,(pt.contig_end+a.chr_start-a.contig_start)," .
    "         (a.chr_start+a.contig_end-pt.contig_start)) as end, " .
    "       a.contig_ori * pt.contig_strand as strand, " .
    "       pt.start_phase, pt.score, pt.p_value " .
    "FROM $source.assembly a, $target.tmp_chr_map tcm, " .
    "     $source.prediction_transcript pt " . 
    "WHERE pt.contig_id = a.contig_id " .
    "AND   a.chromosome_id = tcm.old_id " .
    "ORDER BY pt.prediction_transcript_id, exon_rank";

  my $sth = $dbh->prepare($sql);
  $sth->execute();

  my $prev_end   = undef;
  my $prev_start = undef;
  my $prev_id    = undef;
  my $rank       = 0;

  my %prev_exon = ();

  while(my $row = $sth->fetchrow_arrayref()) {
    my ($pt_id, $sr_id, $sr_start, $sr_end, $sr_strand, $start_phase,
        $score, $p_value) = @$row;
       
    if(defined($prev_id) && ($prev_id == $pt_id)) {
      #still in the same transcript

      if($sr_strand == 1 && 
         defined($prev_end) && $prev_end == $sr_start-1) {
        $self->debug("merged exon $rank in prediction_transcript $pt_id\n");
        #adjacent exons forward strand - merge them
        $prev_exon{'seq_region_end'} = $sr_end;
        $prev_end = $sr_end;
      } elsif($sr_strand == -1 && 
              defined($prev_start) && $prev_start == $sr_end+1) {
        $self->debug("merged exon $rank in prediction_transcript $pt_id\n");
        #adjacent exons negative strand - merge them
        $prev_exon{'seq_region_start'} = $sr_start;
        $prev_start = $sr_start;
      } else {
        #non-adjacent exons in the same transcript - no merge
        $rank++;

        #store the previous exon
        $self->store_pexon(\%prev_exon);

        #make current exon the previous exon
        %prev_exon = ('prediction_transcript_id' => $pt_id,
                      'seq_region_id'            => $sr_id,
                      'seq_region_start'         => $sr_start,
                      'seq_region_end'           => $sr_end,
                      'seq_region_strand'        => $sr_strand,
                      'start_phase'              => $start_phase,
                      'score'                    => $score,
                      'p_value'                  => $p_value,
                      'rank'                     => $rank);
      } 
    } else {
      #store previous exon
      $self->store_pexon(\%prev_exon) if(%prev_exon);

      #new ptranscript
      $rank      = 0;
      $prev_id   = $pt_id;
      $prev_end  = $sr_end;
      $prev_start = $sr_start;
      %prev_exon = ('prediction_transcript_id' => $pt_id,
                    'seq_region_id'            => $sr_id,
                    'seq_region_start'         => $sr_start,
                    'seq_region_end'           => $sr_end,
                    'seq_region_strand'        => $sr_strand,
                    'start_phase'              => $start_phase,
                    'score'                    => $score,
                    'p_value'                  => $p_value,
                    'rank'                     => $rank);
    }
  }

  #store the very last exon in the table
  $self->store_pexon(\%prev_exon) if(%prev_exon);

  $sth->finish();


  $self->debug("AnophelesGambiae Specific: building prediction_transcript " .
               "table");

  $dbh->do
   ("INSERT INTO $target.prediction_transcript (prediction_transcript_id, " .
    "       seq_region_id, seq_region_start, seq_region_end, " .
    "       seq_region_strand, analysis_id ) " .
    "SELECT pt.prediction_transcript_id, tcm.new_id as seq_region_id, " .
    "    MIN(IF(a.contig_ori=1,(pt.contig_start+a.chr_start-a.contig_start),".
    "          (a.chr_start+a.contig_end-pt.contig_end))) as start, " .
    "    MAX(IF(a.contig_ori=1,(pt.contig_end+a.chr_start-a.contig_start)," .
    "          (a.chr_start+a.contig_end-pt.contig_start))) as end, " .
    "       a.contig_ori * pt.contig_strand as strand, " .
    "       pt.analysis_id " .
    "FROM $source.assembly a, $target.tmp_chr_map tcm, " .
    "     $source.prediction_transcript pt " . 
    "WHERE pt.contig_id = a.contig_id " .
    "AND   a.chromosome_id = tcm.old_id " .
    "GROUP BY prediction_transcript_id");  

  return;
}


#
# helper function to store prediction exon
#
sub store_pexon {
  my $self = shift;
  
  my $pexon = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();
  
  my $store_sth = $dbh->prepare
    ("INSERT INTO $target.prediction_exon (prediction_transcript_id, " .
     "       exon_rank, seq_region_id, seq_region_start, seq_region_end, " .
     "       seq_region_strand, start_phase, score, p_value) " .
     "VALUES (?,?,?,?,?,?,?,?,?)");

  $store_sth->execute($pexon->{'prediction_transcript_id'},
                      $pexon->{'rank'},
                      $pexon->{'seq_region_id'},
                      $pexon->{'seq_region_start'},
                      $pexon->{'seq_region_end'},
                      $pexon->{'seq_region_strand'},
                      $pexon->{'start_phase'},
                      $pexon->{'score'},
                      $pexon->{'p_value'});
  $store_sth->finish();

  return;
}


#
# override the clone to seq region so that the version is not
# tacked on the end
#
sub clone_to_seq_region {
  my $self = shift;
  my $target_cs_name = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  # target coord_system will have a different ID
  $target_cs_name ||= "clone";
  my $cs_id = $self->get_coord_system_id($target_cs_name);

  $self->debug("AnophelesGambiae Specific: Transforming clones into " .
               "$target_cs_name seq_regions");

  my $select_sth = $dbh->prepare
    ("SELECT cl.clone_id, cl.embl_acc,
             MAX(ctg.embl_offset)+ctg.length-1
     FROM   $source.clone cl, $source.contig ctg
		 WHERE  cl.clone_id = ctg.clone_id GROUP BY ctg.clone_id");
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



1;
