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
    (['chunk',         undef, 'default_version,sequence_level', 3],
     ['chromosome', $ass_def, 'default_version', 1],
     ["scaffold" ,     undef, "default_version", 2]);

  my @assembly_mappings = ("chromosome:$ass_def|chunk",
                           "chromosome:$ass_def|scaffold",
                           "scaffold|chromosome:$ass_def|chunk");

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib,rank) VALUES (?,?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
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
  $self->supercontig_to_seq_region('scaffold');
  $self->chromosome_to_seq_region();
}



sub create_assembly {
  my $self = shift;

  $self->debug("AnophelesGambiae Specific: loading assembly table");

  $self->assembly_contig_chromosome();
  $self->assembly_supercontig_chromosome();

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
  my $rank       = undef;

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
      $rank      = 1;
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



1;
