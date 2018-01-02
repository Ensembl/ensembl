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

package Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor;

=pod


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor

=head1 SYNOPSIS

  my $isea = $dba->get_IntronSupportingEvidenceAdaptor();
  my $ise = $isea->fetch_by_dbID(1);
  my $ise_array = $dfa->fetch_all();

=head1 METHODS

=cut

use strict;
use warnings;
use base qw/Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor/;

use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::IntronSupportingEvidence;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw/throw/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

=head2 list_linked_transcript_ids

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to search with
  Example     : my $transcript_ids = @{$isea->list_linked_transcript_ids($ise)};
  Description : Uses the given IntronSupportingEvidence to find all linked
                transcript ids 
  Returntype  : ArrayRef[Integer] of transcript_id
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut

sub list_linked_transcript_ids {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  my $query = <<'SQL';
select transcript_id from transcript_intron_supporting_evidence 
where intron_supporting_evidence_id =? 
SQL
  return $self->dbc()->sql_helper()->execute_simple(-SQL => $query, -PARAMS => [$sf->dbID()]);
}

=head2 fetch_all_by_Transcript

  Arg[1]      : Bio::EnsEMBL::Transcript Transcript to search with
  Example     : my $ises = $isea->fetch_all_by_Transcript($transcript);
  Description : Uses the given Transcript to search for all instances of
                IntronSupportingEvidence linked to the transcript in the 
                database 
  Returntype  : ArrayRef of IntronSupportingEvidence objects
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut

sub fetch_all_by_Transcript {
  my ($self, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  my $query = <<'SQL';
select intron_supporting_evidence_id from transcript_intron_supporting_evidence where transcript_id =?
SQL
  my $ids = $self->dbc()->sql_helper()->execute_simple(-SQL => $query, -PARAMS => [$transcript->dbID()]);
  return $self->fetch_all_by_dbID_list($ids);
}

=head2 fetch_flanking_exon_ids

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to search with
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to search with
  Example     : my ($prev_id, $next_id) = @{$isea->fetch_flanking_exon_ids($ise, $transcript)};
  Description : Uses the given IntronSupportingEvidence and Transcript to search
                for the recorded previous and next exon database ids 
  Returntype  : ArrayRef 1 row long but with 2 columns representing previous 
                and next IDs respectivly
  Exceptions  : Thrown if arguments are not as stated and for DB errors 

=cut

sub fetch_flanking_exon_ids {
  my ($self, $sf, $transcript) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  my $query = <<'SQL';
select previous_exon_id, next_exon_id 
from transcript_intron_supporting_evidence 
where transcript_id =? and intron_supporting_evidence_id =?
SQL
  my $ids = $self->dbc()->sql_helper()->execute(-SQL => $query, -PARAMS => [$transcript->dbID(), $sf->dbID()]);
  return unless @{$ids}; 
  return @{$ids->[0]};
}

sub _tables {
  return ( [ 'intron_supporting_evidence', 'ise' ] );
}

sub _columns {
  return qw/
    ise.intron_supporting_evidence_id 
    ise.analysis_id 
    ise.seq_region_id ise.seq_region_start ise.seq_region_end ise.seq_region_strand
    ise.hit_name ise.score ise.score_type
    ise.is_splice_canonical
  /;
}

# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Intron supporting features
#  Returntype : listref of Bio::EnsEMBL::IntronSupportingFeatures in target coordinate system
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #
  
  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();

  my @features;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;
  
  my(
    $id,               $analysis_id,    $seq_region_id,
    $seq_region_start, $seq_region_end, $seq_region_strand,
    $hit_name,         $score,          $score_type,
    $splice_canonical);

  $sth->bind_columns(\(
                     $id,               $analysis_id,    $seq_region_id,
                     $seq_region_start, $seq_region_end, $seq_region_strand,
                     $hit_name,         $score,          $score_type,
                     $splice_canonical));
  
  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_cs;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  my $asma;

  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
    $dest_slice_length  = $dest_slice->length();
    $dest_slice_cs      = $dest_slice->coord_system();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id   = $dest_slice->get_seq_region_id();
    $asma               = $self->db->get_AssemblyMapperAdaptor();
  }

  FEATURE: while($sth->fetch()) {

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
    $analysis_hash{$analysis_id} = $analysis;

    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if (!$slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }

    #obtain a mapper if none was defined, but a dest_seq_region was
    if(!$mapper && $dest_slice && !$dest_slice_cs->equals($slice->coord_system)) {
      $mapper = $asma->fetch_by_CoordSystems($dest_slice_cs, $slice->coord_system);
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #

    if ($mapper) {

      if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper') ) {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->map($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs, 1, $dest_slice);

      } else {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
      }

      #skip features that map to gaps or coord system boundaries
      next FEATURE if (!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
      $slice = $slice_hash{"ID:".$seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
    }

    #
    # If a destination slice was provided convert the coords.
    #
    if (defined($dest_slice)) {
      my $seq_region_len = $dest_slice->seq_region_length();

      if ( $dest_slice_strand == 1 ) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end - $dest_slice_start + 1;

        if ( $dest_slice->is_circular ) {
        # Handle circular chromosomes.

          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromosome origin.

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the chromosome
              $seq_region_start -= $seq_region_len;
            }
            if ( $seq_region_end < 0 ) {
              $seq_region_end += $seq_region_len;
            }
          } else {
            if ($dest_slice_start > $dest_slice_end && $seq_region_end < 0) {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $seq_region_len;
              $seq_region_end   += $seq_region_len;
            }
          }
        }
      } else {

        my $start = $dest_slice_end - $seq_region_end + 1;
        my $end = $dest_slice_end - $seq_region_start + 1;

        if ($dest_slice->is_circular()) {

          if ($dest_slice_start > $dest_slice_end) {
            # slice spans origin or replication

            if ($seq_region_start >= $dest_slice_start) {
              $end += $seq_region_len;
              $start += $seq_region_len if $seq_region_end > $dest_slice_start;

            } elsif ($seq_region_start <= $dest_slice_end) {
              # do nothing
            } elsif ($seq_region_end >= $dest_slice_start) {
              $start += $seq_region_len;
              $end += $seq_region_len;

            } elsif ($seq_region_end <= $dest_slice_end) {
              $end += $seq_region_len if $end < 0;

            } elsif ($seq_region_start > $seq_region_end) {
              $end += $seq_region_len;
            }

          } else {

            if ($seq_region_start <= $dest_slice_end and $seq_region_end >= $dest_slice_start) {
              # do nothing
            } elsif ($seq_region_start > $seq_region_end) {
              if ($seq_region_start <= $dest_slice_end) {
                $start -= $seq_region_len;
              } elsif ($seq_region_end >= $dest_slice_start) {
                $end += $seq_region_len;
              }
            }
          }
        }

        $seq_region_start = $start;
        $seq_region_end = $end;
        $seq_region_strand *= -1;

      } ## end else [ if ( $dest_slice_strand...)]

      # Throw away features off the end of the requested slice or on
      # different seq_region.
      if ($seq_region_end < 1
          || $seq_region_start > $dest_slice_length
          || ($dest_slice_sr_id != $seq_region_id)) {
        next FEATURE;
      }
      $slice = $dest_slice;
    }
    
    push( @features,
          $self->_create_feature_fast(
                                    'Bio::EnsEMBL::IntronSupportingEvidence', {
                                      'start'       => $seq_region_start,
                                      'end'         => $seq_region_end,
                                      'strand'      => $seq_region_strand,
                                      'slice'       => $slice,
                                      'analysis'    => $analysis,
                                      'adaptor'     => $self,
                                      'dbID'        => $id,
                                      'hit_name'    => $hit_name,
                                      'score'       => $score,
                                      'score_type'  => $score_type,
                                      'is_splice_canonical' => $splice_canonical,
                                    } ) );
  }
  
  return \@features;
}

####### STORAGE

=head2 store

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to store
  Example     : $isea->store($ise);
  Description : Stores the IntronSupportingEvidence in the database. Duplicates
                are ignored.
  Returntype  : Integer The assigned database identifier
  Exceptions  : Thrown if the given object is not a IntronSupportingEvidence, 
                and for any DB exception. 

=cut

sub store {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_feature');
  
  my $db = $self->db();
  
  if($sf->is_stored($db)) {
    return $sf->dbID();
  }
  
  my $analysis = $sf->analysis();
  my $analysis_id = $analysis->is_stored($db) ? $analysis->dbID() : $db->get_AnalysisAdaptor()->store($analysis);
  
  my $seq_region_id;
  ($sf, $seq_region_id) = $self->_pre_store($sf);

  my $insert_ignore = $self->insert_ignore_clause();
  my $sql = qq{
${insert_ignore} into intron_supporting_evidence
(analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, score, score_type, is_splice_canonical)
values (?,?,?,?,?,?,?,?,?)
  };

  #Used later on for duplicate entry retrieval
  my $query_params = [
    [$analysis_id, SQL_INTEGER],
    [$seq_region_id, SQL_INTEGER],
    [$sf->start(), SQL_INTEGER],
    [$sf->end(), SQL_INTEGER],
    [$sf->strand(), SQL_INTEGER],
    [$sf->hit_name(), SQL_VARCHAR],
  ];

  my $params = [
    @{$query_params},
    [$sf->score(), SQL_FLOAT],
    [$sf->score_type(), SQL_VARCHAR],
    [$sf->is_splice_canonical(), SQL_INTEGER],
  ];
    
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params, -CALLBACK => sub {
    my ( $sth, $dbh, $rv ) = @_;
    if ($rv > 0) {
      $sf->dbID($self->last_insert_id('intron_supporting_evidence_id', undef, 'intron_supporting_evidence'));
    }
    return;
  });
  $sf->adaptor($self);
  
  if(!$sf->dbID()) {
    my $query = <<'SQL';
select intron_supporting_evidence_id 
from intron_supporting_evidence
where analysis_id =? 
and seq_region_id =? and seq_region_start =? and seq_region_end =?  and seq_region_strand =?
and hit_name =? 
SQL
    my $id = $self->dbc()->sql_helper()->execute_single_result(-SQL => $query, -PARAMS => $query_params);
    $sf->dbID($id);
  }
  
  return $sf->dbID();
}

=head2 store_transcript_linkage

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to link
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to link
  Arg[3]      : Integer an optional ID to give if the Transcript's own ID is possibly incorrect
  Example     : $isea->store_transcript_linkage($ise, $transcript);
                $isea->store_transcript_linkage($ise, $transcript, $tid);
  Description : Links a Transcript to a portion of Intron evidence
  Returntype  : None
  Exceptions  : Thrown if the given object is not a Transcript, if the 
                transcript is not stored, if the supporting evidence is not
                stored and for any DB exception. 

=cut

sub store_transcript_linkage {
  my ($self, $sf, $transcript, $transcript_id) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  
  throw "Cannot perform the link. The IntronSupportingEvidence must be persisted first" unless $sf->is_stored($self->db());
  
  my $insert_ignore = $self->insert_ignore_clause();
  my $sql = qq{
${insert_ignore} into transcript_intron_supporting_evidence
(transcript_id, intron_supporting_evidence_id, previous_exon_id, next_exon_id)
values (?,?,?,?)
  };

  my $intron = $sf->get_Intron($transcript);
  my ($previous_exon, $next_exon) = ($intron->prev_Exon(), $intron->next_Exon());
  $transcript_id ||= $transcript->dbID();
  
  my $params = [
    [$transcript_id, SQL_INTEGER],
    [$sf->dbID(), SQL_INTEGER],
    [$previous_exon->dbID(), SQL_INTEGER],
    [$next_exon->dbID(), SQL_INTEGER],
  ];
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params);
  
  return;
}

####### UPDATE

=head2 update

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to update
  Example     : $isea->update($ise);
  Description : Updates all attributes of an evidence object
  Returntype  : None
  Exceptions  : Thrown if the given object is not a IntronSupportingEvidence,
                if the object is not stored and for normal DB errors

=cut

sub update {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  if (! $sf->is_stored($self->db())) {
    throw "Cannot update the supporting evidence if it has not already been stored in this database";
  }
  
  my $params = [
    [$sf->analysis()->dbID(), SQL_INTEGER],
    [$sf->slice()->get_seq_region_id(), SQL_INTEGER],
    [$sf->start(), SQL_INTEGER],
    [$sf->end(), SQL_INTEGER],
    [$sf->strand(), SQL_INTEGER],
    [$sf->hit_name(), SQL_VARCHAR],
    [$sf->score(), SQL_FLOAT],
    [$sf->score_type(), SQL_VARCHAR],
    [$sf->is_splice_canonical() || 0, SQL_INTEGER],
    [$sf->dbID(), SQL_INTEGER],
  ];
  
  my $sql = <<'SQL';
UPDATE intron_supporting_evidence
SET analysis_id =?, seq_region_id =?, seq_region_start =?, 
seq_region_end =?, seq_region_strand =?, hit_name =?, score =?, score_type =?,
is_splice_canonical =?
WHERE intron_supporting_evidence_id =?
SQL
  
  $self->dbc()->sql_helper()->execute_update(-SQL => $sql, -PARAMS => $params);
  return;
}

####### DELETION

=head2 remove

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence
  Example			: $isea->remove($ise);
  Description	: Deletes the given IntronSupportingEvidence from the database. 
                This can only occur if the object has no linked transcripts
  Returntype 	: None
  Exceptions 	: Thrown if the IntronSupportingEvidence is not stored, if
                the object has linked transcripts and in the event of any
                database error

=cut

sub remove {
  my ($self, $sf) = @_;
  assert_ref($sf, 'Bio::EnsEMBL::IntronSupportingEvidence', 'intron_supporting_evidence');
  if (! $sf->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence if it has not already been stored in this database";
  }
  if($sf->has_linked_transcripts()) {
    throw sprintf('Cannot delete supporting evidence %d. It still has transcripts attached', $sf->dbID());
  }
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from intron_supporting_evidence where intron_supporting_evidence_id =?', 
    -PARAMS => [[$sf->dbID(), SQL_INTEGER]],
  );
  return;
}

=head2 remove_all_transcript_linkages

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence
  Example     : $isea->remove_all_transcript_linkages($ise);
  Description : Deletes the transcript links to the given IntronSupportingEvidence
  Returntype  : None
  Exceptions  : See remove_transcript_linkage

=cut

sub remove_all_transcript_linkages {
  my ($self, $sf) = @_;
  foreach my $transcript_id (@{$self->list_linked_transcript_ids($sf)}) {
    $self->_remove_transcript_linkage($sf, $transcript_id);
  }
  return;
}

=head2 remove_transcript_linkage

  Arg[1]      : Bio::EnsEMBL::IntronSupportingEvidence Evidence to unlink
  Arg[2]      : Bio::EnsEMBL::Transcript Transcript to unlink
  Example     : $isea->remove_transcript_linkages($ise, $transcript);
  Description : Deletes a transcript's link to the given IntronSupportingEvidence
  Returntype  : None
  Exceptions  : Thrown if the given object is not a Transcript, if the 
                transcript is not stored, if the supporting evidence is not
                stored and for any DB exception. 

=cut

sub remove_transcript_linkage {
  my ($self, $sf, $transcript) = @_;
  assert_ref($transcript, 'Bio::EnsEMBL::Transcript', 'transcript');
  if (! $transcript->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence to transcript linkage if the transcript has not already been stored in this database";
  }
  $self->_remove_transcript_linkage($sf, $transcript->dbID());
  return;
}

sub _remove_transcript_linkage {
  my ($self, $sf, $transcript_id) = @_;
  if (! $sf->is_stored($self->db())) {
    throw "Cannot delete the supporting evidence to transcript linkage if the evidence has not already been stored in this database";
  }
  $self->dbc()->sql_helper()->execute_update(
    -SQL => 'DELETE from transcript_intron_supporting_evidence where intron_supporting_evidence_id =? and transcript_id =?', 
    -PARAMS => [[$sf->dbID(), SQL_INTEGER], [$transcript_id, SQL_INTEGER]],
  );
  return;
}

1;
