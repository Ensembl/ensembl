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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor -
Performs database interaction related to PredictionTranscripts

=head1 SYNOPSIS

  # get a prediction transcript adaptor from the database
  $pta = $database_adaptor->get_PredictionTranscriptAdaptor();

  # get a slice on a region of chromosome 1
  $sa = $database_adaptor->get_SliceAdaptor();

  $slice = $sa->fetch_by_region( 'chromosome', 'x', 100000, 200000 );

  # get all the prediction transcripts from the slice region
  $prediction_transcripts = @{ $pta->fetch_all_by_Slice($slice) };

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor;

use vars qw( @ISA );
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Utils::Exception qw(deprecate throw warning);

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor );


# _tables
#
#  Arg [1]    : none
#  Example    : none
#  Description: Implements abstract superclass method to define the table used
#               to retrieve prediction transcripts from the database
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch

sub _tables {
  my $self = shift;

  return ['prediction_transcript', 'pt'];
}


# _columns

#  Arg [1]    : none
#  Example    : none
#  Description: Implements abstract superclass method to define the columns
#               retrieved in database queries used to create prediction 
#               transcripts.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : generic_fetch
#

sub _columns {
  my $self = shift;

  return qw( pt.prediction_transcript_id
             pt.seq_region_id
             pt.seq_region_start
             pt.seq_region_end
             pt.seq_region_strand
             pt.analysis_id
             pt.display_label);
}


=head2 fetch_by_stable_id

  Arg [1]    : string $stable_id
               The stable id of the transcript to retrieve
  Example    : $trans = $trans_adptr->fetch_by_stable_id('GENSCAN00000001234');
  Description: Retrieves a prediction transcript via its display_label.
               This method is called fetch_by_stable_id for polymorphism with
               the TranscriptAdaptor.  Prediction transcript display_labels are
               not necessarily stable in that the same identifier may be reused
               for a completely different prediction transcript in a subsequent
               database release.
  Returntype : Bio::EnsEMBL::PredictionTranscript
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my $self = shift;
  my $stable_id = shift;

  throw('Stable_id argument expected') if(!$stable_id);

  my $syn = $self->_tables()->[1];

  $self->bind_param_generic_fetch($stable_id,SQL_VARCHAR);
  my $pts = $self->generic_fetch("$syn.display_label = ?");

  return (@$pts) ? $pts->[0] : undef;
}



=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch transcripts on.
  Arg [3]    : (optional) boolean $load_exons
               if true, exons will be loaded immediately rather than
               lazy loaded later.
  Example    : $transcripts = $
  Description: Overrides superclass method to optionally load exons
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of transcripts whose
               exons are going to be used.
  Returntype : reference to list of transcripts
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Transcripts
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self  = shift;
  my $slice = shift;
  my $logic_name = shift;
  my $load_exons = shift;

  my $transcripts = $self->SUPER::fetch_all_by_Slice($slice,$logic_name);

  # if there are 0 transcripts still do lazy-loading
  if(!$load_exons || @$transcripts < 1) {
    return $transcripts;
  }

  # preload all of the exons now, instead of lazy loading later
  # faster than 1 query per transcript

  # get extent of region spanned by transcripts
  my ($min_start, $max_end);
  my $ext_slice;

  unless ($slice->is_circular()) {
    foreach my $t (@$transcripts) {
      if (!defined($min_start) || $t->seq_region_start() < $min_start) {
	$min_start = $t->seq_region_start();
      }
      if (!defined($max_end) || $t->seq_region_end() > $max_end) {
	$max_end = $t->seq_region_end();
      }
    }

    if ($min_start >= $slice->start() && $max_end <= $slice->end()) {
      $ext_slice = $slice;
    } else {
      my $sa = $self->db()->get_SliceAdaptor();
      $ext_slice = $sa->fetch_by_region($slice->coord_system->name(), $slice->seq_region_name(), $min_start, $max_end, $slice->strand(), $slice->coord_system->version());
    }

  } else {
    # feature might be crossing the origin of replication (i.e. seq_region_start > seq_region_end)
    # the computation of min_start|end based on seq_region_start|end is not safe
    # use feature start/end relative to the slice instead
    my ($min_start_feature, $max_end_feature);
    foreach my $t (@$transcripts) {
      if (!defined($min_start) || ($t->start() >= 0 && $t->start() < $min_start)) {
  	$min_start = $t->start();
  	$min_start_feature = $t;
      }
      if (!defined($max_end) || ($t->end() >= 0 && $t->end() > $max_end)) {
  	$max_end = $t->end();
  	$max_end_feature = $t;
      }
    }
    
    # now we can reassign min_start|end to seq_region_start|end of
    # the feature which spans the largest region
    $min_start = $min_start_feature->seq_region_start();
    $max_end = $max_end_feature->seq_region_end();

    my $sa = $self->db()->get_SliceAdaptor();
    $ext_slice = 
      $sa->fetch_by_region($slice->coord_system->name(), 
  			   $slice->seq_region_name(), 
  			   $min_start, 
  			   $max_end, 
  			   $slice->strand(), 
  			   $slice->coord_system->version());
  }


  # associate exon identifiers with transcripts

  my %tr_hash = map {$_->dbID => $_} @$transcripts;

  my $tr_id_str = '(' . join(',', keys %tr_hash) . ')';

  my $sth = $self->prepare
    ("SELECT prediction_transcript_id, prediction_exon_id, exon_rank " .
     "FROM   prediction_exon " .
     "WHERE  prediction_transcript_id IN $tr_id_str");

  $sth->execute();

  my ($ex_id, $tr_id, $rank);
  $sth->bind_columns(\$tr_id, \$ex_id, \$rank);

  my %ex_tr_hash;

  while($sth->fetch()) {
    $ex_tr_hash{$ex_id} ||= [];
    push @{$ex_tr_hash{$ex_id}}, [$tr_hash{$tr_id}, $rank];
  }

  $sth->finish();

  my $ea = $self->db()->get_PredictionExonAdaptor();
  my $exons = $ea->fetch_all_by_Slice($ext_slice);

  # move exons onto transcript slice, and add them to transcripts
  foreach my $ex (@$exons) {
    $ex = $ex->transfer($slice) if($slice != $ext_slice);

    if(!$ex) {
      throw("Unexpected. PredictionExon could not be transfered onto " .
            "PredictionTranscript slice.");
    }

    foreach my $row (@{$ex_tr_hash{$ex->dbID()}}) {
      my ($tr, $rank) = @$row;
      $tr->add_Exon($ex, $rank);
    }
  }

  return $transcripts;
}


=head2 fetch_by_prediction_exon_id

  Arg [1]    : Int $prediction_exon_id
               Unique database identifier for the prediction exon
               whose prediction transcript should be retrieved.
  Example    : $prediction_transcript = $prediction_transcript_adaptor->fetch_by_exon_id(1241);
  Description: Retrieves a prediction transcript from the database via the database identifier
               of one of its exons.
  Returntype : Bio::EnsEMBL::PredictionTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_prediction_exon_id {
  my ($self, $prediction_exon_id) = @_;

  # this is a cheap SQL call
  my $sth = $self->prepare(
  qq(
      SELECT pe.prediction_transcript_id
      FROM prediction_exon pe
      WHERE pe.prediction_exon_id = ?
  ));

  $sth->bind_param(1, $prediction_exon_id, SQL_INTEGER);
  $sth->execute();

  my ($prediction_transcript_id) = $sth->fetchrow_array();

  $sth->finish();

  return undef if (!defined $prediction_transcript_id);

  my $prediction_transcript = $self->fetch_by_dbID($prediction_transcript_id);
  return $prediction_transcript;
}





=head2 _objs_from_sth

  Arg [1]    : DBI:st $sth 
               An executed DBI statement handle
  Arg [2]    : (optional) Bio::EnsEMBL::Mapper $mapper 
               An mapper to be used to convert contig coordinates
               to assembly coordinates.
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice to map the prediction transcript to.   
  Example    : $p_transcripts = $self->_objs_from_sth($sth);
  Description: Creates a list of Prediction transcripts from an executed DBI
               statement handle.  The columns retrieved via the statement 
               handle must be in the same order as the columns defined by the
               _columns method.  If the slice argument is provided then the
               the prediction transcripts will be in returned in the coordinate
               system of the $slice argument.  Otherwise the prediction 
               transcripts will be returned in the RawContig coordinate system.
  Returntype : reference to a list of Bio::EnsEMBL::PredictionTranscripts
  Exceptions : none
  Caller     : superclass generic_fetch
  Status     : Stable

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();

  my @ptranscripts;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my (
    $prediction_transcript_id, $seq_region_id,     $seq_region_start,
    $seq_region_end,           $seq_region_strand, $analysis_id,
    $display_label);

  $sth->bind_columns(\(
                     $prediction_transcript_id, $seq_region_id,     $seq_region_start,
                     $seq_region_end,           $seq_region_strand, $analysis_id,
                     $display_label));

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

    # Finally, create the new PredictionTranscript.
    push( @ptranscripts,
          $self->_create_feature('Bio::EnsEMBL::PredictionTranscript', {
                                   '-start'    => $seq_region_start,
                                   '-end'      => $seq_region_end,
                                   '-strand'   => $seq_region_strand,
                                   '-adaptor'  => $self,
                                   '-slice'    => $slice,
                                   '-analysis' => $analysis,
                                   '-dbID' => $prediction_transcript_id,
                                   '-display_label' => $display_label
                                 } ) );

  }

  return \@ptranscripts;
}



=head2 store

  Arg [1]    : list of Bio::EnsEMBL::PredictionTranscript @pre_transcripts 
  Example    : $prediction_transcript_adaptor->store(@pre_transcripts);
  Description: Stores a list of given prediction transcripts in database. 
               Puts dbID and Adaptor into each object stored object.
  Returntype : none
  Exceptions : on wrong argument type 
  Caller     : general 
  Status     : Stable

=cut

sub store {
  my ( $self, @pre_transcripts ) = @_;

  my $ptstore_sth = $self->prepare
    (qq{INSERT INTO prediction_transcript (seq_region_id, seq_region_start,
                                           seq_region_end, seq_region_strand, 
                                           analysis_id, display_label)
        VALUES( ?, ?, ?, ?, ?, ?)});

  my $ptupdate_sth = $self->prepare
    (qq{UPDATE prediction_transcript SET display_label = ?
        WHERE  prediction_transcript_id = ?});

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();
  my $pexon_adaptor = $db->get_PredictionExonAdaptor();

  FEATURE: foreach my $pt (@pre_transcripts) {
    if(!ref($pt) || !$pt->isa('Bio::EnsEMBL::PredictionTranscript')) {
      throw('Expected PredictionTranscript argument not [' . ref($pt).']');
    }

    #skip prediction transcripts that have already been stored
    if($pt->is_stored($db)) {
      warning('Not storing already stored prediction transcript '. $pt->dbID);
      next FEATURE;
    }

    #get analysis and store it if it is not in the db
    my $analysis = $pt->analysis();
    if(!$analysis) {
      throw('Prediction transcript must have analysis to be stored.');
    }
    if(!$analysis->is_stored($db)) {
      $analysis_adaptor->store($analysis);
    }

    #ensure that the transcript coordinates are correct, they may not be,
    #if somebody has done some exon coordinate juggling and not recalculated
    #the transcript coords.
    $pt->recalculate_coordinates();

    my $original = $pt;
    my $seq_region_id;
    ($pt, $seq_region_id) = $self->_pre_store($pt);

    #store the prediction transcript
    $ptstore_sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $ptstore_sth->bind_param(2,$pt->start,SQL_INTEGER);
    $ptstore_sth->bind_param(3,$pt->end,SQL_INTEGER);
    $ptstore_sth->bind_param(4,$pt->strand,SQL_TINYINT);
    $ptstore_sth->bind_param(5,$analysis->dbID,SQL_INTEGER);
    $ptstore_sth->bind_param(6,$pt->display_label,SQL_VARCHAR);

    $ptstore_sth->execute();

    my $pt_id = $self->last_insert_id('prediction_transcript_id', undef, 'prediction_transcript');
    $original->dbID($pt_id);
    $original->adaptor($self);

    #store the exons
    my $rank = 1;
    foreach my $pexon (@{$original->get_all_Exons}) {
      $pexon_adaptor->store($pexon, $pt_id, $rank++);
    }

    # if a display label was not defined autogenerate one
    if(!defined($pt->display_label())) {
      my $zeros = '0' x (11 - length($pt_id));
      my $display_label = uc($analysis->logic_name()) . $zeros . $pt_id;
      $ptupdate_sth->bind_param(1,$display_label,SQL_VARCHAR);
      $ptupdate_sth->bind_param(2,$pt_id,SQL_INTEGER);
      $ptupdate_sth->execute();
      $original->display_label($display_label);
    }
  }
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::PredictionTranscript $pt 
  Example    : $prediction_transcript_adaptor->remove($pt);
  Description: removes given prediction transcript $pt from database. 
  Returntype : none
  Exceptions : throws if argument not a  Bio::EnsEMBL::PredictionTranscript
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $pre_trans = shift;

  if(!ref($pre_trans)||!$pre_trans->isa('Bio::EnsEMBL::PredictionTranscript')){
    throw('Expected PredictionTranscript argument.');
  }

  if(!$pre_trans->is_stored($self->db())) {
    warning('PredictionTranscript is not stored in this DB - not removing.');
    return;
  }

  #remove all associated prediction exons
  my $pexon_adaptor = $self->db()->get_PredictionExonAdaptor();
  foreach my $pexon (@{$pre_trans->get_all_Exons}) {
    $pexon_adaptor->remove($pexon);
  }

  #remove the prediction transcript
  my $sth = $self->prepare( "DELETE FROM prediction_transcript
                             WHERE prediction_transcript_id = ?" );
  $sth->bind_param(1,$pre_trans->dbID,SQL_INTEGER);
  $sth->execute();

  #unset the adaptor and internal id
  $pre_trans->dbID(undef);
  $pre_trans->adaptor(undef);
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$prediction_transcript_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all prediction transcript
               features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("prediction_transcript", undef, $ordered);
}

1;
